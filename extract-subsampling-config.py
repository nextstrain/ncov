import argparse
from copy import deepcopy
import difflib
from pathlib import Path
import tempfile
# pip install "ruamel.yaml<0.18.0"
import ruamel.yaml

yaml=ruamel.yaml.YAML()


SUBSAMPLING_CONFIG_DIR = 'subsampling/'
Path(SUBSAMPLING_CONFIG_DIR).mkdir(exist_ok=True)


# https://stackoverflow.com/a/60099750
def recursive_delete_comments(d):
    if isinstance(d, dict):
        for k, v in d.items():
            recursive_delete_comments(k)
            recursive_delete_comments(v)
    elif isinstance(d, list):
        for elem in d:
            recursive_delete_comments(elem)
    try:
         # literal scalarstring might have comment associated with them
         attr = 'comment' if isinstance(d, ruamel.yaml.scalarstring.ScalarString) \
                  else ruamel.yaml.comments.Comment.attrib 
         delattr(d, attr)
    except AttributeError:
        pass


def recursive_replace(d, old, new):
    if isinstance(d, dict):
        for k, v in d.items():
            recursive_replace(k, old, new)
            if isinstance(v, str):
                d[k] = v.replace(old, new)
            else:
                recursive_replace(v, old, new)
    elif isinstance(d, list):
        for i, v in enumerate(d):
            if isinstance(v, str):
                d[i] = v.replace(old, new)
            else:
                recursive_replace(v, old, new)


def resolve_template(config, old, new):
    recursive_replace(config, old, new)


def write_subsampling_config(path, scheme):
    config = {
        'samples': scheme
    }

    with open(path, 'w') as f:
        yaml.dump(config, f)


def extract_from_workflow_config_builds(input_path, use_scheme_name_for_filename=False):
    print(f"Reading subsampling schemes from {input_path}. Configs that are already present and identical will be ignored.")

    with open(input_path) as f:
        workflow_config = yaml.load(f)
        recursive_delete_comments(workflow_config)

    # For each build entry, write the subsampling scheme as a file.
    for build_name, build_config in workflow_config['builds'].items():
        scheme_name = build_config['subsampling_scheme']
        if use_scheme_name_for_filename:
            output_path = Path(SUBSAMPLING_CONFIG_DIR, f"{scheme_name}.yaml")
        else:
            output_path = Path(SUBSAMPLING_CONFIG_DIR, f"{build_name}.yaml")

        # deepcopy for temporary inplace modifications
        scheme = deepcopy(workflow_config['subsampling'][scheme_name])

        if 'region' in build_config:
            resolve_template(scheme, '{region}', build_config['region'])
        if 'country' in build_config:
            resolve_template(scheme, '{country}', build_config['country'])
        # TODO: add other templates

        if output_path.exists():
            # Check that it is the same.
            new_config_path = tempfile.NamedTemporaryFile().name
            write_subsampling_config(new_config_path, scheme)
            with open(output_path) as existing_f, open(new_config_path) as new_f:
                diff = list(difflib.unified_diff(
                    existing_f.readlines(),
                    new_f.readlines(),
                ))
            if len(diff) != 0:
                print(f"ERROR: Subsampling config for {build_name} exists and differs.")
                for line in diff:
                    print(line, end="")
                exit(1)
        else:
            print(f"Writing new config to {output_path}...")
            write_subsampling_config(output_path, scheme)


def extract_from_workflow_config_subsampling(input_path):
    with open(input_path) as f:
        workflow_config = yaml.load(f)
        recursive_delete_comments(workflow_config)

    for name, scheme in workflow_config['subsampling'].items():
        output_path = Path(SUBSAMPLING_CONFIG_DIR, f"{name}.yaml")
        write_subsampling_config(output_path, scheme)


def main():
    # Extract one subsampling config per build in the following configfiles.
    extract_from_workflow_config_builds('nextstrain_profiles/100k/config-gisaid.yaml')
    extract_from_workflow_config_builds('nextstrain_profiles/100k/config-open.yaml')
    extract_from_workflow_config_builds('nextstrain_profiles/nextstrain-country/builds.yaml')
    extract_from_workflow_config_builds('nextstrain_profiles/nextstrain-open/builds.yaml')
    extract_from_workflow_config_builds('nextstrain_profiles/nextstrain-gisaid/builds.yaml')
    extract_from_workflow_config_builds('nextstrain_profiles/nextstrain-gisaid-21L/builds.yaml')

    # CI has a build definition but it's named "europe" which doesn't represent usage solely by CI.
    extract_from_workflow_config_builds('nextstrain_profiles/nextstrain-ci/builds.yaml',use_scheme_name_for_filename=True)

    # This file has no build definitions to extract from.
    extract_from_workflow_config_subsampling('defaults/parameters.yaml')



if __name__ == '__main__':
    main()
