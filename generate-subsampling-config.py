from itertools import product
from pathlib import Path
from typing import List, Optional
import yaml

# Avoid writing aliases for identical objects.
# https://stackoverflow.com/a/30682604
yaml.Dumper.ignore_aliases = lambda *args: True


SUBSAMPLING_CONFIG_DIR = 'subsampling/'
Path(SUBSAMPLING_CONFIG_DIR).mkdir(exist_ok=True)


class Filters:
    min_date: Optional[str]
    max_date: Optional[str]
    excludes: Optional[List[str]]

    def __init__(self, min_date=None, max_date=None, excludes=None):
        self.min_date = min_date
        self.max_date = max_date
        self.excludes = excludes
        # TODO: add more options

    def to_dict(self):
        options = dict()

        if self.min_date:
            options['min_date'] = self.min_date

        if self.max_date:
            options['max_date'] = self.max_date

        if self.excludes:
            options['exclude'] = self.excludes

        return options

File = str

class Sample:
    name: str
    filters: Optional[Filters]
    weight: Optional[int]
    weighted_sampling: Optional[File]
    uniform_sampling: Optional[List[str]]

    def __init__(self, name, weight=None, weighted_sampling=None, uniform_sampling=None, filters=None):
        self.name = name
        self.weight = weight
        self.weighted_sampling = weighted_sampling
        self.uniform_sampling = uniform_sampling
        self.filters = filters

    def to_dict(self):
        options = dict()

        if self.weight:
            options['weight'] = self.weight

        if self.filters:
            options['filters'] = self.filters.to_dict()

        if self.weighted_sampling:
            options['weighted_sampling'] = self.weighted_sampling

        if self.uniform_sampling:
            options['uniform_sampling'] = self.uniform_sampling

        return options


class Config:
    size: int
    samples: Optional[List[Sample]]

    def __init__(self, size, samples=None):
        if samples == None:
            samples = []

        self.size = size
        self.samples = samples

    def add(self, new_sample: Sample):
        if any(new_sample.name == sample.name for sample in self.samples):
            raise Exception(f"ERROR: A sample with the name {new_sample.name} already exists.")

        self.samples.append(new_sample)

    def to_dict(self):
        return {
            'size': self.size,
            'samples': {
                sample.name: sample.to_dict() for sample in self.samples
            }
        }

    def to_file(self, path):
        print(f'Writing {path}.')
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

        # Global gets special treatment because it is not a region.
        if region == 'Global':
            config = Config(size=5150)

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
                        uniform_sampling=[
                            GROUP_BY_GEOGRAPHICAL_RESOLUTION[location],
                            'year',
                            'month',
                        ],
                        weight=weights[location],
                        filters=Filters(excludes=excludes[location]),
                    ))
            else:
                # Early sequences
                for location in locations:
                    config.add(Sample(
                        name=f"{location.lower().replace(' ', '_')}_early",
                        uniform_sampling=[
                            GROUP_BY_GEOGRAPHICAL_RESOLUTION[location],
                            'year',
                            'month',
                        ],
                        weight=(WEIGHT_EARLY * weights[location]),
                        filters=Filters(
                            excludes=excludes[location],
                            max_date=time,
                        ),
                    ))

                # Recent sequences
                for location in locations:
                    config.add(Sample(
                        name=f"{location.lower().replace(' ', '_')}_recent",
                        uniform_sampling=[
                            GROUP_BY_GEOGRAPHICAL_RESOLUTION[location],
                            GROUP_BY_RECENT_TEMPORAL_RESOLUTION[time],
                        ],
                        weight=(WEIGHT_RECENT * weights[location]),
                        filters=Filters(
                            excludes=excludes[location],
                            min_date=time,
                        ),
                    ))

        # Asia gets special treatment because two countries must be weighted differently.
        elif region == 'Asia':
            config = Config(size=4375)

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
                        uniform_sampling=[
                            # TODO: use GROUP_BY_GEOGRAPHICAL_RESOLUTION?
                            'division',
                            'year',
                            'month',
                        ],
                        weight=(WEIGHT_FOCAL * weights[location]),
                        filters=Filters(
                            excludes=excludes[location],
                        ),
                    ))
                
                # Contextual sequences from the rest of the world
                config.add(Sample(
                    name='context',
                    uniform_sampling=[
                        'country',
                        'year',
                        'month',
                    ],
                    weight=WEIGHT_CONTEXTUAL,
                    filters=Filters(
                        excludes=['region=Asia'],
                    ),
                ))

            else:
                # Early focal sequences for region
                for location in locations:
                    config.add(Sample(
                        name=f"{location.lower().replace(' ', '_')}_early",
                        uniform_sampling=[
                            # TODO: use GROUP_BY_GEOGRAPHICAL_RESOLUTION?
                            'division',
                            'year',
                            'month',
                        ],
                        weight=(WEIGHT_EARLY * WEIGHT_FOCAL * weights[location]),
                        filters=Filters(
                            max_date=time,
                            excludes=excludes[location],
                        ),
                    ))

                # Early contextual sequences from the rest of the world
                config.add(Sample(
                    name='context_early',
                    uniform_sampling=[
                        'country',
                        'year',
                        'month',
                    ],
                    weight=(WEIGHT_EARLY * WEIGHT_CONTEXTUAL),
                    filters=Filters(
                        max_date=time,
                        excludes=['region=Asia'],
                    ),
                ))

                # Recent focal sequences for region
                for location in locations:
                    config.add(Sample(
                        name=f"{location.lower().replace(' ', '_')}_recent",
                        uniform_sampling=[
                            # TODO: use GROUP_BY_GEOGRAPHICAL_RESOLUTION?
                            'division',
                            # TODO: use GROUP_BY_RECENT_TEMPORAL_RESOLUTION?
                            'year',
                            'month',
                        ],
                        weight=(WEIGHT_RECENT * WEIGHT_FOCAL * weights[location]),
                        filters=Filters(
                            min_date=time,
                            excludes=excludes[location],
                        ),
                    ))

                # Recent contextual sequences from the rest of the world
                config.add(Sample(
                    name='context_recent',
                    uniform_sampling=[
                        'country',
                        # TODO: use GROUP_BY_RECENT_TEMPORAL_RESOLUTION?
                        'year',
                        'month',
                    ],
                    weight=(WEIGHT_RECENT * WEIGHT_CONTEXTUAL),
                    filters=Filters(
                        min_date=time,
                        excludes=['region=Asia'],
                    ),
                ))

        # Everything else is a "standard" region with dynamic geographical/temporal grouping.
        else:
            config = Config(size=4000)

            if time == 'all-time':
                # Focal sequences for region
                config.add(Sample(
                    name='focal',
                    uniform_sampling=[
                        GROUP_BY_GEOGRAPHICAL_RESOLUTION[region],
                        'year',
                        'month',
                    ],
                    weight=WEIGHT_FOCAL,
                    filters=Filters(
                        excludes=[f'region!={region}'],
                    ),
                ))

                # Contextual sequences from the rest of the world
                config.add(Sample(
                    name='context',
                    uniform_sampling=[
                        'country',
                        'year',
                        'month',
                    ],
                    weight=WEIGHT_CONTEXTUAL,
                    filters=Filters(
                        excludes=[f'region={region}'],
                    ),
                ))
            else:
                # Early focal sequences for region
                config.add(Sample(
                    name='focal_early',
                    uniform_sampling=[
                        GROUP_BY_GEOGRAPHICAL_RESOLUTION[region],
                        'year',
                        'month',
                    ],
                    weight=(WEIGHT_EARLY * WEIGHT_FOCAL),
                    filters=Filters(
                        max_date=time,
                        excludes=[f'region!={region}'],
                    ),
                ))

                # Early contextual sequences from the rest of the world
                config.add(Sample(
                    name='context_early',
                    uniform_sampling=[
                        'country',
                        'year',
                        'month',
                    ],
                    weight=(WEIGHT_EARLY * WEIGHT_CONTEXTUAL),
                    filters=Filters(
                        max_date=time,
                        excludes=[f'region={region}'],
                    ),
                ))

                # Recent focal sequences for region
                config.add(Sample(
                    name='focal_recent',
                    uniform_sampling=[
                        GROUP_BY_GEOGRAPHICAL_RESOLUTION[region],
                        GROUP_BY_RECENT_TEMPORAL_RESOLUTION[time],
                    ],
                    weight=(WEIGHT_RECENT * WEIGHT_FOCAL),
                    filters=Filters(
                        min_date=time,
                        excludes=[f'region!={region}'],
                    ),
                ))

                # Recent contextual sequences from the rest of the world
                config.add(Sample(
                    name='context_recent',
                    uniform_sampling=[
                        'country',
                        GROUP_BY_RECENT_TEMPORAL_RESOLUTION[time],
                    ],
                    weight=(WEIGHT_RECENT * WEIGHT_CONTEXTUAL),
                    filters=Filters(
                        min_date=time,
                        excludes=[f'region={region}'],
                    ),
                ))

        config.to_file(filename)


def write_reference_build():
    config = Config()
    config.add(Sample(
        name='clades',
        uniform_sampling=['Nextstrain_clade'],
        size=300,
    ))
    filename = Path(SUBSAMPLING_CONFIG_DIR, f"reference.yaml")
    config.to_file(filename)


def main():
    write_region_time_builds()
    # write_reference_build()
    # TODO: 100k builds
    # TODO: CI builds
    # TODO: country builds


if __name__ == '__main__':
    main()
