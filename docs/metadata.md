## Main Metadata File

Running the phylodynamics analysis described in this repo requires 2 main inputs: the sequence data, in FASTA format, and metadata in TSV format.
Here we describe the format of the metadata which the pipeline expects.
An example of this can be seen with the metadata [metadata.tsv](../data/metadata.tsv) for shared sequences which is what we use for the global analysis at [nextstrain.org/ncov/global](https://nextstrain.org/ncov/global). 

#### A few points before we dive in:
- The _order_ of the fields doesn't matter; but if you are going to join your metadata with the global collection then it's easiest to keep them in the same order!
- Not all fields are currently used, but this may change in the future.
- A tab character -- `\t` -- separates columns. Read more about the TSV format [here](https://en.wikipedia.org/wiki/Tab-separated_values).
- Data is case sensitive
- The "geographic" columns, such as "region" and "country" will be used to plot the samples on the map.
Adding a new value to these columns isn't a problem at all, but lat/long coordinates will also need to be defined for this in a seaparate file (documentation forthcoming).
- A lot of these fields will be available to interact with as "colorings" of the data in Auspice (the interactive visualisation used within Nextstrain). However which exact columns are used, and which colours are used for each value is completely customisable (documentation forthcoming).
- See also the [augur documentation](https://nextstrain-augur.readthedocs.io/en/stable/faq/metadata.html) relating to metadata formatting

### Missing data:

Missing data is to be expected for certain fields.
In general, missing data is represented by an empty string or a question mark character.
There is one important difference: if a discrete trait reconstruction (e.g. via `augur traits`) is to be run on this column, then a value of `?` will be inferred, whereas the empty string will be treated as missing data in the output. See below for how to represent uncertainty in sample collection date.

---


Here's an example of the first few columns of the metadata for a single strain, including the header row. _(Spacing between columns here is adjusted for clarity, and only the first 6 of 23 columns are shown)._
```
strain              virus  gisaid_epi_isl  genbank_accession   date        region   ...
NewZealand/01/2020  ncov   EPI_ISL_413490  ?                   2020-02-27  Oceania  ...
```

In total there are 23 columns of metadata for each genome.
Let's jump in and explore them, using the above example 👇

**Column 1: `strain`**

This needs to match the name of a sequence in the FASTA file exactly and must not contain characters such as spaces, or `()[]{}|#><`.
In our example we have a strain called "NewZealand/01/2020" so there should be a sequence in the FASTA file for ">NewZealand/01/2020" (sequence names in FASTA files always start with the `>` character, but this is not part of the name).

**Note that "strain" here carries no biological or functional significance** and should be thought of as synonymous with sample.

**Column 2: `virus`**

Currently unused.

**Column 3: `gisaid_epi_isl`**

If this genome is shared via [GISAID](https://www.gisaid.org/) then please include the EPI ISL here. In our example this is "EPI_ISL_413490".

**Column 4: `genbank_accession`**

If this genome is shared via [GenBank](https://www.ncbi.nlm.nih.gov/genbank/) then please include the accession number here. In our example this is "?" indicating that it hasn't (yet) been deposited in GenBank. (See above for more information on how to encode missing data.)

**Column 5: `date`** (really important!)

This describes the sample collection data (_not_ sequencing date!) and must be formated according as `YYYY-MM-DD`. 
Our example was collected on Feb 27, 2020 and is therefore represented as "2020-02-27".

You can specify unknown dates or month by replacing the respected values by `XX` (ex: `2013-01-XX` or `2011-XX-XX`) and completely unknown dates can be shown with `20XX-XX-XX` (which does not restrict the sequence to being in the 21st century - they could be earlier).
Please be aware that our current pipeline will filter out any genomes with an unknown date, however you can change this for your pipeline!

**Column 6: `region`**

The region the sample was collected in -- for our example this is "Oceania".
Please use one of "Africa", "Asia", "Europe", "North America", "Oceania" or "South America".
If you sequence a genome from Antartica please get in touch!


**Column 7: `country`**

The country the sample was collected in. Our example, "NewZealand/01/2020", was collected in ....... New Zealand.
You can run `tail +2 data/metadata.tsv | cut -f 7 | sort | uniq` to see all the countries currently present in the metadata.
As of April 10 there were 64! 🌎

**Column 8: `division`**

Division currently doesn't have a precise definition and we use it differently for different regions.
For instance, for samples in the USA we are detailing the state in which the sample was collected here. For other countries, it might be a county, region, or other administrative sub-division.
To see the divisions which are currently set for your country you can run the following command (replace "New Zealand" with your country):
```bash
tail +2 data/metadata.tsv | cut -f 7,8 | grep "^New Zealand" | cut -f 2 | sort | uniq
```

**Column 9: `location`**

Similarly to `division`, but for a smaller geographic resolution. This data is often unavailable, and missing data here is typically represented by an empty field or the same value as `division` is used.
In our example the division is "Auckland", which conveniently or confusingly is both a province of New Zealand and a city.

**Column 10: `region_exposure`**

If the sample has a known travel history and infection is thought to have occured in this location, then represent this here.
In our example, which represents New Zealand's first known case, the patient had recently arrived from Iran, thus the value here is "Asia".
Specifying these travel histories helps inform the model we use to reconstruct the geographical movements of the virus.

If there is no travel history then set this to be the same value as `region`.


**Column 11: `country_exposure`**

Analogous to `region_exposure` but for `country`.
In our example, given the patient's travel history, this is set to "Iran".

**Column 12: `division_exposure`**

Analogous to `region_exposure` but for `division`. If we don't know the exposure division, we may specify the value for `country_exposure` here as well.

**Column 13: `segment`**

Unused. Typically the value "genome" is set here.

**Column 14: `length`**

Genome length (numeric value).

**Column 15: `host`**

Host from which the sample was collected.
Currently we have multiple values in the dataset, including "Human", "Canine", "Manis javanica" and "Rhinolophus affinis".

**Column 16: `age`**

Numeric age of the patient from whom the sample was collected.
We round to an integer value.
This will show up in auspice when clicking on the tip in the tree which brings up an info box and 

**Column 17: `sex`**

Sex of the patient from whom the sample was collected.
This will show up in auspice when clicking on the tip in the tree which brings up an info box.

**Column 18: `originating_lab`**

Please see [GISAID](https://www.gisaid.org/help/publish-with-gisaid-references/) for more information.

**Column 19: `submitting_lab`**

Please see [GISAID](https://www.gisaid.org/help/publish-with-gisaid-references/) for more information.

**Column 20: `authors`**

Author of the genome sequence, or the paper which announced this genome.
Typically written as "LastName et al".
In our example, this is "Storey et al".
This will show up in auspice when clicking on the tip in the tree which brings up an info box.


**Column 21: `url`**

The URL, if available, pointing to the genome data.
For most SARS-CoV-2 data this is [https://www.gisaid.org](https://www.gisaid.org).

**Column 22: `title`**

The URL, if available, of the publication announcing these genomes. 

**Column 23: `date_submitted`**

Date the genome was submitted to a public database (most often GISAID).
In `YYYY-MM-DD` format (see `date` for more information on this formatting).

