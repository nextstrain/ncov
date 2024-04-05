from itertools import product
from pathlib import Path
from typing import Any, List, Optional
import yaml

# Avoid writing aliases for identical objects.
# https://stackoverflow.com/a/30682604
yaml.Dumper.ignore_aliases = lambda *args: True


SUBSAMPLING_CONFIG_DIR = 'subsampling/'
Path(SUBSAMPLING_CONFIG_DIR).mkdir(exist_ok=True)


class Sample:
    group_by: Optional[List[str]]
    size: Optional[int]
    min_date: Optional[str]
    max_date: Optional[str]
    excludes: Optional[List[str]]
    disable_probabilistic_sampling: Optional[bool]
    priorities: Optional[Any]

    def __init__(self, name, **kwargs):
        self.name: str = name

        # Initialize instance attributes
        for option in self.__annotations__.keys():
            if option not in self.__annotations__:
                raise Exception(f'Option {option!r} not allowed.')
            # TODO: Check types
            value = kwargs[option] if option in kwargs else None
            self.__setattr__(option, value)

    def to_dict(self):
        options = dict()

        if self.group_by:
            options['group_by'] = self.group_by

        if self.size:
            options['max_sequences'] = self.size

        if self.min_date:
            options['min_date'] = self.min_date

        if self.max_date:
            options['max_date'] = self.max_date

        if self.disable_probabilistic_sampling:
            options['disable_probabilistic_sampling'] = True

        if self.excludes:
            options['exclude'] = self.excludes

        return options


class Config:
    samples: List[Sample]
    def __init__(self):
        self.samples = []

    def add(self, new_sample: Sample):
        if any(new_sample.name == sample.name for sample in self.samples):
            raise Exception(f"ERROR: A sample with the name {new_sample.name} already exists.")

        self.samples.append(new_sample)

    def to_dict(self):
        return {
            'samples': {
                sample.name: sample.to_dict() for sample in self.samples
            }
        }

    def to_file(self, path):
        print(f'Writing {path}.  n={sum(sample.size for sample in self.samples)}')
        with open(path, 'w') as f:
            yaml.dump(self.to_dict(), f, sort_keys=False)


def write_region_time_builds():
    TIMES = [
        '1M',
        '2M',
        '6M',
        'all-time',
    ]
    REGIONS = [
        'Africa',
        'Asia',
        'Europe',
        'Global',
        'North America',
        'Oceania',
        'South America',
    ]
    WEIGHT_EARLY, WEIGHT_RECENT = 1, 4
    WEIGHT_FOCAL, WEIGHT_CONTEXTUAL = 4, 1

    GROUP_BY_RECENT_TEMPORAL_RESOLUTION = {
        '1M': 'week',
        '2M': 'week',
        '6M': 'month',
        'all-time': 'month',
    }

    GROUP_BY_GEOGRAPHICAL_RESOLUTION = {
        'Africa': 'country',
        'Asia': 'country',
        'China': 'division',
        'India': 'division',
        'Europe': 'country',
        'South America': 'country',
        'North America': 'division',
        'Oceania': 'division',
    }

    for region, time in product(REGIONS, TIMES):

        # Generate by region followed by time.
        # For each region, 'all-time' gets special treatment because it does not need an early/recent split.
        build_name = f"{region.lower().replace(' ', '-')}_{time.lower()}"
        filename = Path(SUBSAMPLING_CONFIG_DIR, f"{build_name}.yaml")

        config = Config()

        # Global gets special treatment because it is not a region.
        if region == 'Global':
            target_size = 5150

            locations = [
                'Africa',
                'Asia',
                'China',
                'Europe',
                'India',
                'North America',
                'South America',
                'Oceania',
            ]

            # FIXME: add this to commit message instead
            # Deduced using https://www.mathwarehouse.com/calculators/ratio-simplifier-3-or-more-numbers.php
            weights = {
                'Africa': 30,
                'Asia': 40,
                'China': 35,
                'Europe': 25,
                'India': 35,
                'North America': 20,
                'South America': 18,
                'Oceania': 3,
            }

            sum_location_weights = sum(weights.values())
            assert sum_location_weights == 206

            excludes = {
                'Africa': ['region!=Africa'],
                'Asia': ['region!=Asia', 'country=China', 'country=India'],
                'China': ['country!=China'],
                'Europe': ['region!=Europe'],
                'India': ['country!=India'],
                'North America': ['region!=North America'],
                'South America': ['region!=South America'],
                'Oceania': ['region!=Oceania'],
            }

            if time == 'all-time':
                for location in locations:
                    config.add(Sample(
                        name=location.lower().replace(' ', '_'),
                        group_by=[
                            GROUP_BY_GEOGRAPHICAL_RESOLUTION[location],
                            'year',
                            'month',
                        ],
                        size=int(
                            target_size
                            * weights[location] / sum_location_weights
                        ),
                        excludes=excludes[location],
                    ))
            else:
                # Early sequences
                for location in locations:
                    config.add(Sample(
                        name=f"{location.lower().replace(' ', '_')}_early",
                        group_by=[
                            GROUP_BY_GEOGRAPHICAL_RESOLUTION[location],
                            'year',
                            'month',
                        ],
                        size=int(
                            target_size
                            * WEIGHT_EARLY / (WEIGHT_EARLY + WEIGHT_RECENT)
                            * weights[location] / sum_location_weights
                        ),
                        max_date=time,
                        excludes=excludes[location],
                    ))

                # Recent sequences
                for location in locations:
                    config.add(Sample(
                        name=f"{location.lower().replace(' ', '_')}_recent",
                        group_by=[
                            GROUP_BY_GEOGRAPHICAL_RESOLUTION[location],
                            GROUP_BY_RECENT_TEMPORAL_RESOLUTION[time],
                        ],
                        size=int(
                            target_size
                            * WEIGHT_RECENT / (WEIGHT_EARLY + WEIGHT_RECENT)
                            * weights[location] / sum_location_weights
                        ),
                        min_date=time,
                        excludes=excludes[location],
                    ))

        # Asia gets special treatment because two countries must be weighted differently.
        elif region == 'Asia':
            target_size = 4375

            locations = [
                'Asia',
                'China',
                'India',
            ]

            weights = {
                'Asia': 3,
                'China': 2,
                'India': 2,
            }

            sum_location_weights = sum(weights.values())
            assert sum_location_weights == 7

            excludes = {
                'Asia': ['region!=Asia', 'country=China', 'country=India'],
                'China': ['country!=China'],
                'India': ['country!=India'],
            }

            if time == 'all-time':
                # Focal sequences for region
                for location in locations:
                    config.add(Sample(
                        name=location.lower().replace(' ', '_'),
                        group_by=[
                            # TODO: use GROUP_BY_GEOGRAPHICAL_RESOLUTION?
                            'division',
                            'year',
                            'month',
                        ],
                        size=int(
                            target_size
                            * WEIGHT_FOCAL / (WEIGHT_FOCAL + WEIGHT_CONTEXTUAL)
                            * weights[location] / sum_location_weights
                        ),
                        excludes=excludes[location],
                    ))
                
                # Contextual sequences from the rest of the world
                config.add(Sample(
                    name='context',
                    group_by=[
                        'country',
                        'year',
                        'month',
                    ],
                    size=int(
                        target_size
                        * WEIGHT_CONTEXTUAL / (WEIGHT_FOCAL + WEIGHT_CONTEXTUAL)
                    ),
                    excludes=['region=Asia'],
                ))

            else:
                # Early focal sequences for region
                for location in locations:
                    config.add(Sample(
                        name=f"{location.lower().replace(' ', '_')}_early",
                        group_by=[
                            # TODO: use GROUP_BY_GEOGRAPHICAL_RESOLUTION?
                            'division',
                            'year',
                            'month',
                        ],
                        size=int(
                            target_size
                            * WEIGHT_EARLY / (WEIGHT_EARLY + WEIGHT_RECENT)
                            * WEIGHT_FOCAL / (WEIGHT_FOCAL + WEIGHT_CONTEXTUAL)
                            * weights[location] / sum_location_weights
                        ),
                        max_date=time,
                        excludes=excludes[location],
                    ))

                # Early contextual sequences from the rest of the world
                config.add(Sample(
                    name='context_early',
                    group_by=[
                        'country',
                        'year',
                        'month',
                    ],
                    size=int(
                        target_size
                        * WEIGHT_EARLY / (WEIGHT_EARLY + WEIGHT_RECENT)
                        * WEIGHT_CONTEXTUAL / (WEIGHT_FOCAL + WEIGHT_CONTEXTUAL)
                    ),
                    max_date=time,
                    excludes=['region=Asia'],
                ))

                # Recent focal sequences for region
                for location in locations:
                    config.add(Sample(
                        name=f"{location.lower().replace(' ', '_')}_recent",
                        group_by=[
                            # TODO: use GROUP_BY_GEOGRAPHICAL_RESOLUTION?
                            'division',
                            # TODO: use GROUP_BY_RECENT_TEMPORAL_RESOLUTION?
                            'year',
                            'month',
                        ],
                        size=int(
                            target_size
                            * WEIGHT_RECENT / (WEIGHT_EARLY + WEIGHT_RECENT)
                            * WEIGHT_FOCAL / (WEIGHT_FOCAL + WEIGHT_CONTEXTUAL)
                            * weights[location] / sum_location_weights
                        ),
                        min_date=time,
                        excludes=excludes[location],
                    ))

                # Recent contextual sequences from the rest of the world
                config.add(Sample(
                    name='context_recent',
                    group_by=[
                        'country',
                        # TODO: use GROUP_BY_RECENT_TEMPORAL_RESOLUTION?
                        'year',
                        'month',
                    ],
                    size=int(
                        target_size
                        * WEIGHT_RECENT / (WEIGHT_EARLY + WEIGHT_RECENT)
                        * WEIGHT_CONTEXTUAL / (WEIGHT_FOCAL + WEIGHT_CONTEXTUAL)
                    ),
                    min_date=time,
                    excludes=['region=Asia'],
                ))

        # Everything else is a "standard" region with dynamic geographical/temporal grouping.
        else:
            target_size = 4000

            if time == 'all-time':
                # Focal sequences for region
                config.add(Sample(
                    name='focal',
                    group_by=[
                        GROUP_BY_GEOGRAPHICAL_RESOLUTION[region],
                        'year',
                        'month',
                    ],
                    size=int(
                        target_size
                        * WEIGHT_FOCAL / (WEIGHT_FOCAL + WEIGHT_CONTEXTUAL)
                    ),
                    excludes=[f'region!={region}']
                ))

                # Contextual sequences from the rest of the world
                config.add(Sample(
                    name='context',
                    group_by=[
                        'country',
                        'year',
                        'month',
                    ],
                    size=int(
                        target_size
                        * WEIGHT_CONTEXTUAL / (WEIGHT_FOCAL + WEIGHT_CONTEXTUAL)
                    ),
                    excludes=[f'region={region}']
                ))
            else:
                # Early focal sequences for region
                config.add(Sample(
                    name='focal_early',
                    group_by=[
                        GROUP_BY_GEOGRAPHICAL_RESOLUTION[region],
                        'year',
                        'month',
                    ],
                    size=int(
                        target_size
                        * WEIGHT_EARLY / (WEIGHT_EARLY + WEIGHT_RECENT)
                        * WEIGHT_FOCAL / (WEIGHT_FOCAL + WEIGHT_CONTEXTUAL)
                    ),
                    max_date=time,
                    excludes=[f'region!={region}'],
                ))

                # Early contextual sequences from the rest of the world
                config.add(Sample(
                    name='context_early',
                    group_by=[
                        'country',
                        'year',
                        'month',
                    ],
                    size=int(
                        target_size
                        * WEIGHT_EARLY / (WEIGHT_EARLY + WEIGHT_RECENT)
                        * WEIGHT_CONTEXTUAL / (WEIGHT_FOCAL + WEIGHT_CONTEXTUAL)
                    ),
                    max_date=time,
                    excludes=[f'region={region}'],
                ))

                # Recent focal sequences for region
                config.add(Sample(
                    name='focal_recent',
                    group_by=[
                        GROUP_BY_GEOGRAPHICAL_RESOLUTION[region],
                        GROUP_BY_RECENT_TEMPORAL_RESOLUTION[time],
                    ],
                    size=int(
                        target_size
                        * WEIGHT_RECENT / (WEIGHT_EARLY + WEIGHT_RECENT)
                        * WEIGHT_FOCAL / (WEIGHT_FOCAL + WEIGHT_CONTEXTUAL)
                    ),
                    min_date=time,
                    excludes=[f'region!={region}'],
                ))

                # Recent contextual sequences from the rest of the world
                config.add(Sample(
                    name='context_recent',
                    group_by=[
                        'country',
                        GROUP_BY_RECENT_TEMPORAL_RESOLUTION[time],
                    ],
                    size=int(
                        target_size
                        * WEIGHT_RECENT / (WEIGHT_EARLY + WEIGHT_RECENT)
                        * WEIGHT_CONTEXTUAL / (WEIGHT_FOCAL + WEIGHT_CONTEXTUAL)
                    ),
                    min_date=time,
                    excludes=[f'region={region}'],
                ))

        # Double check the total sample size.
        total_size = sum(sample.size for sample in config.samples)
        assert target_size == total_size

        config.to_file(filename)


def write_reference_build():
    config = Config()
    config.add(Sample(
        name='clades',
        group_by=['Nextstrain_clade'],
        size=300,
    ))
    filename = Path(SUBSAMPLING_CONFIG_DIR, f"reference.yaml")
    config.to_file(filename)


def write_ci_build():
    config = Config()
    config.add(Sample(
        name='region',
        group_by=[
            'division',
            'year',
            'month',
        ],
        size=20,
        disable_probabilistic_sampling=True,
        excludes=['region!=Europe']
    ))
    config.add(Sample(
        name='global',
        group_by=[
            'year',
            'month',
        ],
        size=10,
        disable_probabilistic_sampling=True,
        excludes=['region=Europe'],
        # TODO: add Priority(type=proximity, focus=region)
    ))
    filename = Path(SUBSAMPLING_CONFIG_DIR, f"nextstrain_ci_sampling.yaml")
    config.to_file(filename)


def main():
    write_region_time_builds()
    write_reference_build()
    write_ci_build()
    # TODO: 100k builds
    # TODO: country builds


if __name__ == '__main__':
    main()
