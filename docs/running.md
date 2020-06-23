# Running the analysis

>This section focuses on how to running the basic example build to give you a chance to practice and get a sense of how things work. The next section covers customizing and configuring your own build.

**To run our analyses, we need to:**  
1. Ensure our **sequence data and metadata is [properly formatted](data-prep.md)**  
2. **Specify which builds you want** to generate using a `builds.yaml` file  
3. **Execute the workflow**    
4. [Hopefully you don't have to] **troubleshoot**

## Step 1. Gather and format your data  

If you haven't done this step yet, check out our [data prep](data-prep.md) guide and come back when you're ready.  

## Step 2. Specify which builds to run    

In the orientation section, we learned that  
- [Nextstrain analyses are run using a workflow manager called Snakemake](orientation-workflow.md)  
- [A "build"](glossary.md#Build) is a bundle of input files, parameters, and commands   
- [Each build is configured by two primary files](orientation-files.md): `build.yaml` and `config.yaml`  

Let's start with defining a build in `./my_config/example/builds.yaml`.  
**We use the `builds.yaml` file to define what geographic areas of the world we want to focus on. Each block in this file will produce a separate output JSON for visualization**.  

The first block of the provided `build.yaml` file looks like this:  
```
# build focusing on King County (location) in Washington State (division) in the USA (country)

builds:
  usa_washington_king-county: # name of the build
    subsampling_scheme: location # what subsampling method to use (see default_config.yaml)
    geographic_scale: location # our focal area is a 'location' (e.g., county or city)
    region: North America
    country: USA
    division: Washington
    location: King County
```
Looking at this example, we can see that each build has a:  
- `build_name`, which is used for naming output files  
- `subsampling_scheme`, which specifies how sequences are selected. Default schemes exist for `region`, `country`, and `division`. Custom schemes [can be defined](###).
- `geographic_scale` specifies whether the focal area is a `location`, `division`, `country` or `region`.  
- `geographic_name` the name of the focal area; if not provided, defaults to the `build_name`  
- `region`, `country`, and `division`: specify the location of the sample, down to the specified `geographic_scale`.

The rest of the builds defined in this file serve as examples for division-, country- or region-focused analyses.

**To adapt this for your own analyses, copy `my_config/example` to `my_config/<my-new-name>` and open the `builds.yaml`
file in this directory.**

Go ahead and **swap out the values in this file with the geographic area of interest to you.** You can add, disable (comment out), or remove as many of these build definitions as you'd like.  

## Step 3: Run the workflow  

To actually execute the workflow, run:  

```bash  
sarscov2-tutorial$ snakemake --profile my_config/example -p  
```

`--profile` tells snakemake where to find your `builds.yaml` and `config.yaml` files.  
`-p` tells snakemake to print each command it runs to help you understand what it's doing.   

If you'd like to run a dryrun, try running with the `-np` flag, which will execute a dryrun. This prints out each command, but doesn't execute it.


## Step 4: Troubleshoot common issues

If you have a question which is not addressed here, please don't hestitate to [ask for help](index.md#Help)


#### My country / division does not show up on the map

This is most often a result of the country / division not being present in [the file defining the latitude & longitude of each deme](../config/lat_longs.tsv).
Adding it to that file (and rerunning the Snakemake rules downstream of this) should fix this.

#### My trait (e.g. division) is grey instead of colored

We generate the colors from the `colors` rule in the Snakefile, which uses the [ordering TSV](./default_config/ordering.tsv) to generate these. See ['customizing your analysis'](customizing-analysis.md) for more info.

_*A note about locations and colors:*_
Unless you want to specifically override the colors generated, it's usually easier to _add_ information to the default `ncov` files, so that you can benefit from all the information already in those files.

#### My genomes aren't included in the analysis

There are a few steps where sequences can be removed:

- During the `filter` step:
    - Samples which are included in [the exclude file](../default_config/exclude.tsv) are removed
    - Samples which fail the current filtering criteria, as defined in the `default_config.yaml` file, are removed. You can modify the snakefile as desired, but currently these are:
        - Minimum sequence length of 25kb
        - No ambiguity in (sample collection) date
    - Samples may be randomly removed during subsampling; see ['customizing your analysis'](customizing-analysis.md) for more info.
  - During the `refine` step, where samples that deviate more than 4 interquartile ranges from the root-to-tip vs time are removed

#### `Error: Where there's SAMPLING_TRAIT we should always have EXPOSURE_TRAIT`

This comes from an incomplete metadata file.
If you define (e.g.) `country` for a sample then you _must_ also define `country_exposure` for that sample.
If there is no (known) travel history, then you can set the same values for each.


#### Sequencing and alignment errors

Genome sequencing, bioinformatic processing of the raw data, and alignment of the sequences are all steps were errors can slip in.
Such errors can distort the phylogenetic analysis.
To avoid sequences with known problems to mess up the analysis, we keep a list of problematic sequences in `config/exclude.txt` and filter them out.
To facilitate spotting such problematic sequences, we added an additional quality control step that produces the files

 * `results/sequence-diagnostics.tsv`
 * `results/flagged-sequences.tsv`
 * `results/to-exclude.txt`

These files are the output of `scripts/diagnostics.py` and are produced by rule `diagnostic`.
The first file contains statistics for every sequence in the aligment, sorted by divergence worst highest to lowest.
The second file contains only those sequences with diagnostics exceeding thresholds each with their specific reason for flagging -- these are sorted by submission date (newest to oldest).
The third file contains only the names of the flagged sequences and mirrors the format of `config/exclude.txt`.
These names could be added to `config/exclude.txt` for permanent exclusion.
Note, however, that some sequences might look problematic due to alignment issues rather than intrinsic problems with the sequence.
The flagged sequences will be excluded from the current run.

To only run the sequence diagnostic, you can specify any of the three above files as target or run:
```bash
snakemake --profile my_config/<name> diagnostic
```

In addition, we provide rules to re-examine the sequences in `config/exclude.txt`.
By running
```bash
snakemake --profile my_config/<name> diagnose_excluded
```
the pipeline will produce

 * `results/excluded-sequence-diagnostics.tsv`
 * `results/excluded-flagged-sequences.tsv`
 * `results/check-exclusion.txt`

These files are meant to facilitate checking whether sequences in `config/exclude.txt` are excluded for valid reasons.

## [Previous Section: Orientation: which files should I touch?](orientation-files.md)
## [Next Section: Orientation: Customizing your analysis](customizing-analysis.md)
