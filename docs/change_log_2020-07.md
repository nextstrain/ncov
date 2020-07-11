# Change Log & Breaking Changes for July 2020 Repo Update

In order to make our repository more accessible for first-time users, and those who want to create their own customized build, we made a number of changes in July 2020.
Here, we outline those changes, and also list 'breaking changes' that may effect advanced users who were already running their own builds. 

You can check the lists below to see if there's any action you could take to make updating to the latest version of this repository is smooth!

## Summary of changes  

The goal of this release is to make it easier for new users to get their own `ncov` build up and running. This consists of a new, extensive tutorial hosted on github pages; a simplified repository structure; and more didactic file names. 

Summary of breaking changes:    
- All nextstrain-specific files now live under `nextstrain_profiles`  
- `config` is now named `defaults`  
- `envs`, `rules` and `schemas` now live under `workflow`  

### Additions and updates    
- README is updated  
- Adds in-depth tutorial to `docs/` and accompanying github pages setup files  
- `data/` now contains example data that should remain packaged with the repository  
- `narratives/` now contains only a `template_narrative.md` for new users (not to worry -- all finished Situation Reports are still hosted over in the `nextstrain/narratives` repo)  
- `my_profiles` now contains extensive `example` (formerly `king-county`) and `example_advanced_customization` (formerly `swiss`) profiles  

### File name changes  
- The main `config.yaml` file is now more precisely labeled as `parameters.yaml`  
- `reference.gb` > `reference_seq.gb`  
- `ordering.tsv` > `color_ordering.tsv`  
- `config` > `defaults`, with nextstrain-specific files (e.g., `description.md`) moved to `nextstrain_profiles/nextstrain/`  
- nextstrain-specific profiles labeled with `nextstrain_*`  


### Structural changes  
To improve repo clarity and approachability for new users, the following files have been reorganized.
- Workflow-related files now live in a `workflow` directory, rather than the top level  
  - `envs` > `workflow/envs`  
  - `schemas` > `workflow/schemas`  
  - `rules` > `workflow/snakemake_rules`  
- Nextstrain-specific profiles now live under `nextstrain_profiles`  

(Reiterated from last section because they're important :)
- `config` is now `defaults`, with nextstrain-specific files (e.g., `description.md`) moved to `nextstrain_profiles/nextstrain/`  
- `config.yaml` is now `parameters.yaml`

## Breaking Changes

If you've got custom builds running and regularly pull from the repository, you should check here to see if any of the changes you've made might make the merge harder.
Usually, just copying and renaming folders and/or files is enough to avoid too many merge conflicts!

First, you should check the file locations that have changed, above.
If you reference any of these files specifically in your pipelines or profile, be sure to update those links!

Be sure to also check any links to files that were previously in `profiles`. 
For example, if you were referencing the cluster information in `profiles/nextstrain-scicore/cluster.json` for example) - we've now got two profiles folders - you'd find `nextstrain-scicore` in the `nextstrain_profiles` folder now.

If you modified the `swiss`, `default`, or `king-county` profiles, this may cause a merge conflict, as these folders have been renamed.
To keep your builds as they are, we recommend copying these folders, and then after the merge, restoring the name you prefer.

We now have a `data` folder which comes with the repository and contains example build data.
If you had previously created a `data` folder, this may cause a merge conflict - usually this can be easily resolved.


## Developer Changes

For those who help run and maintain the Nextstrain.org builds, there's a few changes that might be different from what you're used to.

Just focusing on the files mostly accessed by the 'build shepherds,' here's what's good to know:
- The `config` folder is now `defaults`
- `orderings.tsv` is now `color_orderings.tsv`
- `rules/nextstrain_exports.smk` (the rules that do the 'finalizing' of the builds) isnow `workflow/snakemake_rules/export_for_nextstrain.smk`

The 'finishing builds' commands are also slightly different because the profiles have moved - they'll look familiar though!
You can already find these changed commands in the `export_for_nextstrain.smk` file, but here they are as well:

```
# To update ordering/lat_longs after AWS download:
#   snakemake --touch --forceall --profile nextstrain_profiles/nextstrain
#   snakemake --profile nextstrain_profiles/nextstrain clean_export_regions
#   snakemake --profile nextstrain_profiles/nextstrain export_all_regions
# When done adjusting lat-longs & orders, remember to run
#   snakemake --profile nextstrain_profiles/nextstrain all_regions
# to produce the final Auspice files!
```