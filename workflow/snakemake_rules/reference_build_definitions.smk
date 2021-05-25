if config.get('us_state_builds'):
    for abbr, state in config.get('us_state_builds'):
        config['builds'][f'US-{abbr}'] = {
            "subsampling_scheme": "nextstrain_local",
            "region": "North America",
            "country": "USA",
            "division": state,
            "auspice_config": "nextstrain_profiles/nextstrain/north-america_auspice_config.json"
        }

if config.get("build_sizes"):
    from copy import deepcopy

    new_schemes = {}
    for scheme in config['subsampling']:
        for size, N in config["build_sizes"].items():
            tmp = deepcopy(config['subsampling'][scheme])
            for y in tmp:
                if 'max_sequences' in tmp[y] and tmp[y]['max_sequences']<1:
                    tmp[y]['max_sequences'] = int(tmp[y]['max_sequences']*N)

            new_schemes[scheme+'_'+size] = tmp

    config['subsampling'].update(new_schemes)

    new_builds = {}
    for build in config['builds']:
        for size, N in config["build_sizes"].items():
            tmp = deepcopy(config["builds"][build])
            tmp["subsampling_scheme"] += '_'+size
            if size=='standard':
                new_builds[build] = tmp
            else:
                new_builds[build+'-'+size] = tmp
    config['builds'] = new_builds
