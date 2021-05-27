if config.get('us_state_builds'):
    for abbr, state in config.get('us_state_builds'):
        config['builds'][f'US-{abbr}'] = {
            "subsampling_scheme": "nextstrain_local",
            "region": "North America",
            "country": "USA",
            "division": state,
            "auspice_config": "nextstrain_profiles/nextstrain/north-america_auspice_config.json"
        }

# If `build_sizes` are specified we want to create copies of
# the builds (one per build size) and create corresponding copies
# of the subsampling schemes and trait definitions.
if config.get("build_sizes"):
    from copy import deepcopy

    # Duplicate the subsampling schemes
    new_schemes = {}
    for scheme in config['subsampling']:
        for size, N in config["build_sizes"].items():
            tmp = deepcopy(config['subsampling'][scheme])
            for y in tmp:
                if 'max_sequences' in tmp[y] and tmp[y]['max_sequences']<1:
                    tmp[y]['max_sequences'] = int(tmp[y]['max_sequences']*N)

            new_schemes[scheme+'_'+size] = tmp

    config['subsampling'].update(new_schemes)

    # duplicate any trait definitions so that there's one per build
    for build_name in list(config.get("traits", {}).keys()):
        if build_name not in config['builds']:
            continue
        for size in config["build_sizes"]:
            if size=='standard':
                continue
            config['traits'][build_name+"-"+size] = config['traits'][build_name]
    
    # duplicate the builds themselves
    new_builds = {}
    for build in config['builds']:
        for size, N in config["build_sizes"].items():
            tmp = deepcopy(config["builds"][build])
            tmp["subsampling_scheme"] += '_'+size
            tmp["build_size"] = {"name": size, "n": N}
            if size=='standard':
                new_builds[build] = tmp
            else:
                new_builds[build+'-'+size] = tmp
    config['builds'] = new_builds
