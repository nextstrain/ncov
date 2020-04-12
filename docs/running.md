### Running a SARS-CoV-2 analysis

The pipeline described in this repo is designed primarily to run the phylogeographic analyses of SARS-CoV-2 data which are displayed on nexstrain.org, for instance [the global analysis](https://nextstrain.org/ncov/global), the [European subsampled build](https://nextstrain.org/ncov/europe), the [North American Build](https://nextstrain.org/ncov/north-america) etc.

It is also possible to run your own data through the analysis here.
Because Nextstrain is open-source, you may modify the analysis to suit your particular needs.

#### A few points before we dive in:

- If you haven't run an analysis using Nextstrain before, I **highly recommend** running the [zika tutorial](https://nextstrain.org/docs/tutorials/zika) first, which introduces these concepts in a gentler fashion.
- If you would like to use Nextstrain Groups, such as [this one](https://nextstrain.org/groups/blab/), to share your results through nextstrain.org then please get in touch! You will have control & ownership of the datasets and narratives, but they can be shared freely through nextstrain.org. These are also available in a private fashion for sensitive data.


#### This page consists of four parts:

- Curating your own input data
- Running the default build using Snakemake
- Understanding the parts of the analysis
- Debugging common issues


## Curating your own input data

Curating sequences and input data is a crucial part of running your own analysis.
Here we assume you have two data sources which you wish to analyse together:

1. Your own genomes & metadata (perhaps collected by your lab)
2. The global data, or a subset of this, so that you can see your data in context.

> As the global dataset grows, subsampling becomes important.
We will write guidance for subsampling in a future page, but you can investigate the `Snakefile_Regions` to see how we perform subsampling for our regional builds.

#### Global data

For the nextstrain.org analyses we use data obtained through [GISAID](https//gisaid.org).
Please see there for how you may obtain those genomic data for your own analysis as the terms of data sharing prevent us making the sequence data publicly available.
Included in this repository is [a curated list of metadata](../data/metadata.tsv) associated with those sequences.

#### Your own data

This should consist of
- a FASTA file with the (unaligned) genomes and the name must not contain characters such as spaces, or `()[]{}|#><` (except for the `>` character which starts each entry in a FASTA file).
- the metadata corresponding to each of those genomes.  Please see the [metadata documentation](./metadata.md) for details of the format of this metadata.


#### Combining the data

Let's assume you have now have four files:
1. `data/global_sequences.fasta` - genomes of worldwide data to provide phylogenetic context
2. `data/global_metadata.tsv` - Metadata of these (global) genomes. This could well be a copy of the `data/metadata.tsv` that's included with this repo.
3. `data/our_sequences.fasta` - Your own sequences, in FASTA format.
4. `data/our_metadata.tsv` - metadata of your genomes. Let's assume this follows the [same format](./metadata.md) as (2).

We can combine the two sets of genomes simply via
```bash
cat data/our_sequences.fasta data/global_sequences.fasta > data/sequences.fasta
```
And, as long as the metadata formats are the same, then we can add our metadata via:
```bash
cp ./data/global_metadata.tsv ./data/metadata.tsv
tail +2 ./data/our_metadata.tsv >> ./data/metadata.tsv
```
(Please double check the columns in this new metadata TSV match up. It's not a problem if there are more entries in the metadata than there are genomes.)



## Running the default build

If the data is in the correct formats (`./data/sequences.fasta` & `./data/metadata.tsv`) then we can generate the analysis by simply running
```bash
snakemake -p -s Snakefile --cores 2 auspice/ncov.json
```
Which will produce a `./auspice/ncov.json` file which you can visualise in Auspice via
```
auspice view --datasetDir auspice
```


## Understanding the parts of the analysis

The Snakemake analysis here consists of a number of rules, which are displayed below.
Not all of the rules here are essential, or may even be desirable for your analysis.
We maintain this snakefile primarily for our analyses, and thus your build may be able to be made a lot simpler!
The aim of this tutorial is to walk you through the rules in this basic analysis run by Nextstrain and give you the ability to change it to suit your needs.

> Note: this repo contains a few different Snakefiles, as we use them to automate a number of analyses, some of which are beyond the scope of this tutorial. 
This tutorial follows the main `Snakefile`.


Each snakemake rule produces output which may be used by other rules; these "dependencies" are indicated by grey arrows.
Additionally, there are input files which must be provided (e.g. the sequence data generated above, or other files which are part of this repo); these are indicated with red outlines.
Please inspect the `Snakefile` to see what each rule is doing in more detail and if there are any questions please [make an issue](https://github.com/nextstrain/ncov/issues/new) so that we can improve this page!

![snakemake_workflow](images/basic_snakemake_build.png)

## Debugging common issues

> If you have a question which is not addressed here, please consider [making an issue](https://github.com/nextstrain/ncov/issues/new) so that we can improve these documentation 👍

#### My country / division does not show up on the map

This is most often a result of the country / division not being present in [the file defining the latitude & longitdue of each deme](../config/lat_longs.tsv).
Adding it to that file (and rerunning the Snakemake rules downstream of this) should fix this.
You can rerun the appropriate parts of the build via:

```bash
snakemake -s Snakefile --cores 2 -p -f results/ncov_with_accessions.json
snakemake -s Snakefile --cores 2 -p -f auspice/ncov.tsv
```

#### My trait (e.g. division) is grey instead of colored

We generate the colors from the `colors` rule in the Snakefile, which uses the [ordering TSV](../ordering.tsv) to generate these.
Once you've modified this file, you can regenerate the appropriate parts of the analysis via:

```bash
snakemake -s Snakefile --cores 2 -p -f config/colors.tsv
snakemake -s Snakefile --cores 2 -p -f auspice/ncov.tsv
```


#### My genomes aren't included in the analysis

There are a few steps where sequences can be removed:

- During the `filter` step:
    - Samples which are included in [the exclude file](../config/exclude.tsv) are removed
    - Samples which fail the current filtering criteria, as defined in the Snakefile, are removed. You can modify the snakefile as desired, but currently these are:
        - Minimum sequence length of 25kb
        - No ambiguity in (sample collection) date
        - A random subsampling of 500 sequences per division per month
  - During the `refine` step, where samples that deviate more than 4 interquartile ranges from the root-to-tip vs time are removed

#### Error: Where there's SAMPLING_TRAIT we should always have EXPOSURE_TRAIT

This comes from an incomplete metadata file.
If you define (e.g.) `country` for a sample then you _must_ also define `country_exposure` for that sample.
If there is no (known) travel history, then you can set the same values for each.
