Previously, we recommended using Snakemake profiles under a `my_profiles/` analysis directory. We now recommend using Snakemake config files directly via the `--configfile` parameter. You can still use existing profiles via `--configfile my_profiles/<profile_name>/builds.yaml`.

See [this guide](https://docs.nextstrain.org/projects/ncov/en/latest/tutorial/next-steps.html#create-analysis-directory) to create your own analysis directory.
