## Data Submitter's FAQ

We often recieve questions from data submittors about why their data is not visible on the [Nextstrain SARS-CoV-2 runs](https://nextstrain.org/ncov). 
This short FAQ highlights some of the main reasons why data may not be showing up on Nextstrain.

### Sequence Length & Number of N's

We currently only use full-genome sequences which are at least 27,000 bases in length. They also cannot have more than 3,000 bases that are 'N'.

### Subsampling

Nextstrain runs can be subsampled considerably. There are over >30,000 whole-genome sequences available on GISAID currently, but we typically include <5,000 in each of our runs. If the division your samples are from contains more than about 100 samples per month, they are likely to be downsampled. Be sure to check the appropriate regional build - these are sampled more heavily from the focal region, so there's a higher chance a sequence will be included in the run. We have regional builds for [North America](https://nextstrain.org/ncov/north-america), [South America](https://nextstrain.org/ncov/south-america), [Asia](https://nextstrain.org/ncov/asia), [Africa](https://nextstrain.org/ncov/africa/), [Europe](https://nextstrain.org/ncov/europe), and [Oceania](https://nextstrain.org/ncov/oceania).

### Missing Dates

We currently only include samples that have an **exact sampling date** (day, month, year). This is because we cannot accurately estimate the sample dates from the sequences at the moment, given the short duration of the pandemic so far, and the mutation rate.

If your sample has only year or only month and year as a sampling date, it will be automatically excluded from runs. If you have privacy/data sharing concerns, it's ok to slightly change the collection date randomly by +/- 1 or 2 days. Please do *not* use the sequencing or processing date, as these can negatively influence our runs. 

If you wish to add a corrected date to your samples, simply updating the sampling date in GISAID  will automatically update our system, and the sequence will be included in the next run!

### Many Samples with the Same Date

If we receive many samples that have identical dates as sample dates, we may exclude these manually. This is because this often indicates that the 'sample date' given is not actually the sample date, but the sequencing, processing, or uploading date. We try to email submitters when we do this to check whether the dates are truly the collection dates.

If you are genuinely submitting many sequences with identical dates, you can avoid us temporarily excluding them by emailing hello@nextstrain.org to let us know about the sequences and why they have the same date (ex: collected during investigation of a long-term care center).

### Missing USA State

We currently exclude samples from the USA which do not have a 'division' attribute (this is the USA state or territory where they were sampled). Adding a state/territory/division to your sample on GISAID will automatically update this on our system, and the sequence will appear in our next run.

### Divergence Issues

For quality control, we use a combination of automated and manual checks to ensure that sequences included seem to be free of sequencing and/or assembly error. If a sequenece is deemed to be far too divergent (has more mutations than we expect given the sampling date), or far too under-diverged (has far fewer mutations than we expect given the sampling date), it may be excluded. We cannot off direct help in these cases, but suggest you revisit the raw sequence files with the aid of someone with experience using your sequencing pipeline, in order to correct any sequencing and assembly errors.

