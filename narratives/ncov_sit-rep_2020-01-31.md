---
title: Genomic analysis of nCoV spread. Situation report 2020-01-30.
authors: "Trevor Bedford, Richard Neher, James Hadfield, Emma Hodcroft, Misja Ilcisin, Nicola Müller"
authorLinks: "https://nextstrain.org"
affiliations: "Fred Hutch, Seattle, USA and Biozentrum, Basel, Switzerland"
date: "2020 January 30"
dataset: "https://nextstrain.org/ncov/2020-01-30?d=map"
abstract: "This report uses publicly shared novel coronavirus (nCoV) genomic data from GISAID and Genbank to estimate rates and patterns of viral epidemic spread. We plan to issue updated situation reports as new data is produced and shared. This website is optimized for display on desktop browsers."
---

# [Executive summary](https://nextstrain.org/ncov/2020-01-30/2020-01-30)

```auspiceMainDisplayMarkdown
## Executive summary

Using 42</tag> publicly shared novel coronavirus (nCoV) genomes, we examined genetic diversity to infer date of common ancestor and rate of spread.
We find:
* the 42</tag> sampled genomes are very similar, differing from the consensus by 0-7 mutations
* This lack of genetic diversity has a parsimonious explanation that the outbreak descends either from a single introduction into the human population or a small number of animal-to-human transmissions of very similar viruses.
* This event most likely occurred in November or early December 2019.
* There has been ongoing human-to-human spread since this point resulting in observed cases.
* Using estimates of total case count from Imperial College London of several thousand cases, we infer a reproductive number between 1.8 and 3.5 indicating rapid growth in the November 2019-Janurary 2020 period.
```

# [Coronaviruses](https://nextstrain.org/ncov/2020-01-30)

### Further Reading:

* General information on coronaviruses on [Wikipedia](https://en.wikipedia.org/wiki/Coronavirus) _2020-01-30_
* Summary of the nCov outbreak on [Wikipedia](https://en.wikipedia.org/wiki/2019%E2%80%9320_Wuhan_coronavirus_outbreak) _2020-01-30_
* Material provided by the [US CDC](https://www.cdc.gov/coronavirus/index.html) _2020-01-29_
* Organization and genome on [ViralZone](https://viralzone.expasy.org/764?outline=all_by_species) _2020-01-23_
* Interactive risk analysis by [MOBS-lab](https://datastudio.google.com/reporting/3ffd36c3-0272-4510-a140-39e288a9f15c/page/U5lCB) _2010-01-29_
* Interactive risk analysis by [ROCS-lab](http://rocs.hu-berlin.de/corona/) _2010-01-29_

```auspiceMainDisplayMarkdown

## Different human coronaviruses

Coronaviruses (CoV) are members of a diverse species of positive-sense single-stranded RNA ((+)ssRNA) viruses which have a history of causing respiratory infections in humans.
Some variants of coronaviruses are associated with outbreaks, others are continuously circulating and cause mostly mild respiratory infections (e.g. the common cold).

#### SARS-CoV & MERS-CoV
The most well known of these coronaviruses is [SARS-CoV](https://en.wikipedia.org/wiki/Severe_acute_respiratory_syndrome) ("severe acute respiratory syndrome"), which in a November 2002 to July 2003 outbreak spread around the world and resulted in [over 8000 cases and 774 deaths](https://www.theguardian.com/world/2017/dec/10/sars-virus-bats-china-severe-acute-respiratory-syndrome), with a case fatality rate of around 9–11%.

In 2012, a novel coronavirus, [MERS-CoV](https://en.wikipedia.org/wiki/Middle_East_respiratory_syndrome) ("Middle East respiratory syndrome"), causing severe respiratory symptoms was identified. MERS has resulted in fatalities comparable to SARS, however the transmission route of MERS is very different. Whereas SARS was efficiently spread from one human to another, human MERS infections were generally a result of independent zoonoses (animal to human transmissions) from camels (see [Dudas _et al._](https://elifesciences.org/articles/31257) for more information). This has lead to a self-limiting outbreak largely restricted to the Arabian Peninsula.


#### Seasonal CoV
However, not all coronaviruses are as deadly as SARS-CoV and MERS-CoV.
There are four "seasonal" coronaviruses that commonly infect humans each year.
Compared with SARS, these seasonal coronavirus strains are ["much more prevalent, much less severe, and common causes of influenza‐like illness (ILI)"](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5820427/).
In fact, [5](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2879166/)–[12](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5820427/)% of all ILI cases test positive for coronaviruses, so they are rather common, resulting in millions of infections every year with low severity.
These seasonal coronaviruses are the results of separate spillovers from the bat animal reservoir into humans in the past ~100 years, in which after spillover, each seasonal virus established itself and spread widely in the human population.


#### Animal reservoirs
Coronaviruses infect a wide range of animals, and the human outbreaks described above are a result of one or more "jumps" from these animal reservoirs into the human population.
SARS is believed to have arrived in the human population from [horseshoe bats via a masked palm civet intermediary](https://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1006698).


#### Human-to-human transmission
The ability for different lineages to be transmitted between humans is extremely important to understand the potential development of an outbreak.
Due to the ability of SARS to spread between humans and the high case fatality rate, SARS (or a SARS-like virus) is considered a [global public health threat](https://www.who.int/whr/2007/overview/en/index1.html) by the WHO.

```

# [Novel coronavirus (nCoV) 2019-2020](https://nextstrain.org/ncov/2020-01-30)

### Further Reading:

* New China virus: Five questions scientists are asking  [Nature news](https://www.nature.com/articles/d41586-020-00166-6) _2020-01-22_
* China virus latest: first US case confirmed  [Nature news](https://www.nature.com/articles/d41586-020-00154-w) _2020-01-21_
* New virus surging in Asia rattles scientists  [Nature news](https://www.nature.com/articles/d41586-020-00129-x) _2020-01-20_
* New virus identified as likely cause of mystery illness in China [Nature news](https://www.nature.com/articles/d41586-020-00020-9) _2020-01-08_

```auspiceMainDisplayMarkdown

## Recent outbreak of a novel coronavirus
In December 2019, a new illness was first detected in Wuhan, China.
We now know this to be another outbreak of coronavirus in humans (the 7th), and it is provisionally being called nCoV (novel coronavirus).

As of January 30th, 2020 over 7,914 cases and 170 deaths [have been reported](https://en.wikipedia.org/wiki/2019%E2%80%9320_outbreak_of_novel_coronavirus_(2019-nCoV)).
It's still too early to know the case fatality rate, but early indications are that it is significantly less than SARS-CoV.
The case counts are dramatically rising in part due to increased surveillance and testing.

While the outbreak seems to be centered in Wuhan, which is now [under quarantine](https://twitter.com/PDChina/status/1220060879112282117), the virus has spread throughout China and abroad, including Hong Kong, Singapore, Japan, and Thailand, as well as Europe, North America, South Asia, the Middle East, and Australia. Limited local transmission outside of China has been reported.

The origin of the virus is still unclear, however [genomic analyis](https://virological.org/t/ncovs-relationship-to-bat-coronaviruses-recombination-signals-no-snakes/331) suggests nCoV is most closely related to viruses previously identified in bats.
It is plausible that there were other intermediate animal transmissions before the introduction into humans.
There is no evidence of snakes as an intermediary.

#### Nextstrain narratives

The following pages contain analysis performed using [Nextstrain](https://nextstrain.org).
Scrolling through the left hand sidebar will reveal paragraphs of text with a corresponding visualization of the genomic data on the right hand side.

To have full genomes of a novel and large RNA virus this quickly is a remarkable achievement.
These analyses have been made possible by the rapid and open sharing of genomic data and interpretations by scientists all around the world (see the final slide for a visualization of sequencing authorship).
```

# [How to interpret the phylogenetic trees](https://nextstrain.org/ncov/2020-01-30)

### Further Reading:

* [Exploring interactive phylogenies with Auspice](https://neherlab.org/201901_krisp_auspice.html) _2019-01-24_

```auspiceMainDisplayMarkdown
## Transmission trees vs phylogenetic trees

Pathogens spread through rapid replication in one host followed by transmission to another host.
An epidemic can only take off when one infection results in more than one subsequent infections.

As the pathogen replicates and spreads, its genome needs to be replicated many times and random mutations (copying mistakes)  will accumulate in the genome.
Such random mutations can help to track the spread of the pathogen and learn about its transmission routes and dynamics.

<div>
  <img alt="cartoon showing how transmission tree and phylogenetic tree relate" width="500" src="https://neherlab.org/talk_images/infection_tree_combined.png"/>
</div>

The illustration above shows a sketch of a transmission tree with a subset of cases that were sampled (blue).
In practice, the transmission tree is unknown and typically only rough estimates of case counts are available.
Genome sequences allow us to infer parts of the transmission tree.
In this example, three mutations (little diamonds) are indicated on the tree.
Sequences that have the same mutations are more closely related, so these mutations allow us to group samples into clusters of closely related viruses that belong to the same transmission chains.

### Reading a Phylogenetic Tree

Below, we see an illustration with a phylogenetic tree on the left, where mutations are shown as colored circles. On the right are the corresponding sequences, also with mutations shown as colored circles.
We can see that sequences that share the same mutations group together.
When sequences appear linked by a flat vertical line, like A and B, this means there are no differences between them – their sequences are identical.

When a sequence sits on a long line on its own, like C or E, this means it has unique mutations not found in other sequences. The longer the line, the more mutations.
A and B also have unique mutations (the green circle) not shared by the other sequences, but they are identical to each other.

<div>
  <img alt="cartoon of phylogenetic tree and corresponding alignment, with samples labelled A-E" width="500" src="https://data.nextstrain.org/toy_alignment_tree.png"/>
</div>

At the moment, the novel coronavirus (nCoV) phylogeny may not look much like a 'tree'.
Many of the sequences are identical – they sit together on vertical lines like A and B (some are on the left-most part of the tree).
Others have unique or shared mutations and so sit on lines, or 'branches', going to the right.
You can see how many mutations a branch has by hovering your mouse over it.
```

# [Phylogenetic analysis](https://nextstrain.org/ncov/2020-01-30?m=div&d=tree)

Here we present a phylogeny of 42</tag> strains of nCoV that have been publicly shared.
Information on how the analysis was performed is available [in this GitHub repository](github.com/nextstrain/ncov).

<br>

The colours represent the within-country region or US-state of isolation, with the x-axis representing nucleotide divergence.

<br>

Divergence is measured as the number of changes (mutations) in the genome.
Several sequences have zero mutations -- meaning they are all identical to the root (center) of the tree.
Other viruses have between one and seven mutations.

<br>

Sequencing the genome of a large novel RNA virus in an evolving outbreak situation is challenging.
Some of the differences observed in these sequences may be sequencing errors rather than actual mutations.
Insertions, deletions, and differences at the ends of the genome are more likely to be errors and so we masked these for the purposes of this analysis.

# [Phylogenetic Interpretation](https://nextstrain.org/ncov/2020-01-30?m=div&d=tree)

We currently see little genetic diversity across the nCoV sequences, with 11</tag> out of 42</tag> sequences having no unique mutations.

<br>

Low genetic diversity across these sequences suggests that the most recent common ancestor of all nCoV sequences was fairly recent, since mutations accumulate slowly compared to other RNA viruses at a rate of around 1-2 mutations per month for coronaviruses.
Generally, repeated introductions from an animal reservoir will show significant diversity (this has been true for Lassa, Ebola, MERS-CoV and avian flu).
The observation of such strong clustering of human infections can be explained by an outbreak that descends from a single zoonotic introduction event into the human population followed by human-to-human epidemic spread.

<br>

We are starting to see groups of sequences that share mutations.
One cluster contains sequences from Guangdong and four isolates from the US.
Other clusters contain two to four isolates.
Sequences in these clusters tend to be from more recent samples, suggesting that the virus has started to accumulate mutations as it spread in Wuhan and subsequently to other cities.
There is currently no evidence that these mutations change how the virus behaves -- it is expected that RNA viruses mutate.

# [Within-family transmission 1](https://nextstrain.org/ncov/2020-01-30?m=div&d=tree&f_location=Zhuhai)

There are three genetically-identical isolates from Zhuhai (Southeastern China, Guangdong Province) which form a cluster, sharing one unique mutation seen in no other isolate (you can hover your mouse over the branches to see which mutations are present).

<br>

Two of these cases (ending 028 and 040) are [known to come from a single family](https://twitter.com/JingLu_LuJing/status/1220143773532880896), again indicating human-to-human transmission.
We don't have information about the third case.


# [Within-family transmission 2](https://nextstrain.org/ncov/2020-01-30?m=div&d=tree&f_location=Shenzhen)

Of the six isolates from Guangdong Province (which includes the city of Shenzhen) we see four isolates which are genetically identical.
These sequences differ by 3 mutations from the root of the tree.

<br>

Three of the sequences from Guangdong (ending F025, F013, and F012) are [known to come from a single family](https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(20)30154-9/fulltext), and almost certainly represent human-to-human transmission.


<br>

# [Within-family transmission 2 - shared mutations](https://nextstrain.org/ncov/2020-01-30?m=div&d=tree&f_location=Shenzhen,Los%20Angeles,Orange%20County,Seattle,Chicago,Phoenix)

The three mutations found in this cluster are also present in the Arizona, USA isolate, and two of the mutations are found in three other USA isolates.


<br>

# [Within-family transmission 3](https://nextstrain.org/ncov/2020-01-30?m=div&d=tree&f_location=Paris)

Finally, the two sequences from France are identical, sharing one unique mutation, and one mutation also found in one of the USA isolates and the Taiwanese isolate.

<br>

The two French sequences are [known to be from the same family](https://www.thelocal.fr/20200129/coronavirus-in-france-what-you-need-to-know) - a Chinese couple from Wuhan.


# [Cases outside China](https://nextstrain.org/ncov/2020-01-30?c=country&d=tree&m=div)

There are reported diagnostically confirmed nCoV cases in many East and South-East Asian countries, USA, Australia, the Middle East, and Europe.
Vietnam, Japan, and Germany have reported transmission within the country, albeit always with a known link to Wuhan, China.

<br>

The only currently available sequence data for cases outside of China are two cases from Thailand, five from the USA, two from France, and one from Taiwan.
The Thai samples are genetically identical to nine Chinese sequences, including seven isolated in Wuhan.
Four sequences from the USA share two mutations with the cluster of sequences from Shenzhen.
The remaining sequence from the USA shares a mutation with the sequence from Taiwan and the two from France.

<br>

The most parsimonious explanation for the observed pattern of mutation sharing between the US and Shenzhen sequences is a virus variant with the two shared mutations was circulating in Wuhan and was independently exported to Shenzhen and multiple times to the USA.
There is no evidence for a link between US sequences other than a link to Wuhan.

# [Dating the time of the most recent common ancestor](https://nextstrain.org/ncov/2020-01-30?d=tree)
The high similarity of the genomes suggests they share a recent common ancestor (i.e. that they have descended from the same ancestral virus recently).
Otherwise, we would expect a higher number of differences between the samples.

<br>

Previous research on related coronavirus suggests that these viruses accumulate between 1 and 3 changes in their genome per month (rates of 3 &times; 10<sup>-4</sup> to 2 &times; 10<sup>-3</sup> per site per year).

<br>

On the right, we explore how different assumptions about the rate of change, and the observed genetic diversity, give us estimates about when all sequenced cases last shared a common ancestor.

```auspiceMainDisplayMarkdown
## Date of the common ancestor of outbreak viruses
The time of the most recent common ancestor (or tMRCA) of a set of sequenced cases denotes when these sequenced cases last shared a common ancestor.
This time can be as early as the time when a virus first entered the human population, but can also be substantially later as show in the figure below.

<div>
  <img alt="Example phylogeny where the time of the initial zoonosis is different from the most recent common ancestor of several sequenced cases" width="500" src="https://data.nextstrain.org/zoonosis.png"/>
</div>

Using different approaches, the tMRCA can be estimated either jointly with the rate of evolution or one can fix this rate.
Below, we estimate this tMRCA using different approaches.

1) With the additional sequences shared during the past week, the tree now shows several distinct clusters such that our analysis from 2020-01-25 assuming a star-like topology is no longer appropriate.


We reproduce here our analysis based on data available up to 2020-01-25
assuming a star-like phylogeny structure along with a Poisson distribution of mutations through time to estimate the time of the most recent common ancestor ('TMRCA') of sequenced viruses.
**We found that the common ancestor most likely existed between mid-November and the beginning of December 2019. The biggest source of uncertainty is the substitution rate.**

<div>
  <img alt="graph of TMRCA estimates based on different mutation rates" width="500" src="https://data.nextstrain.org/ncov_poisson-tmrca.png"/>
</div>

Using the entire data set, the nextstrain analysis pipeline estimates that the common ancestor most likely existed between late-Nov and the beginning of December 2019.

2) With more and more sequences being made available publicly, we can start to estimate the rate of evolution and the tMRCA from sequence data directly.
There, however, remains a great deal of uncertainty in these estimates and the results have to therefore be treated as very uncertain.
[Analyzing](https://nicfel.github.io/nCov-Nicola/ExponentialCoalescent_20200130.html) this data assuming once a constant population of infected individuals and once and exponentially growing population estimates the tMRCA to be anywhere between the end of November and the middle of December.

<div>
  <img alt="TMRCA estimates inferred from genetic sequence data and sampling times in BEAST" width="500" src="https://data.nextstrain.org/beast_coal-tmrca_20200130.png"/>
</div>


More detailed modeling of the onset of the outbreak are ongoing.

```

# [Estimating the growth rate](https://nextstrain.org/ncov/2020-01-30?d=tree)

An important quantity in the spread of a pathogen is the average number of secondary cases each infection produces.

<br>

This number is known as R0 ("R-zero" or "R-nought").
One the right, we present simple estimates of R0.

```auspiceMainDisplayMarkdown
## Estimates of epidemic growth rate
Scientists at Imperial College London have used the number of cases observed outside of China to estimate the [total number of cases](https://www.imperial.ac.uk/mrc-global-infectious-disease-analysis/news--wuhan-coronavirus/) and suggested that there have been at least several thousand cases by 2020-01-22.
With the additional exported cases since and the continued growth of confirmed cases in China, we currently have to expect at least 50000 cases to date.
Together with our previous estimates of the age of the outbreak and information on the infectious period, we can estimate plausible ranges of R0 using a branching process model.

**We find plausible estimates of R0 between 1.8 and 3.5.**

If we assume the outbreak started at the beginning of November 2019 (12 weeks ago), we find that R0 should range between 1.8 and 2.5, depending on how large ('n') the outbreak is now.
<div>
  <img alt="graph of R0 estimates with epidemic start 12 weeks ago" width="500" src="https://data.nextstrain.org/ncov_branching-R0-early_2020-01-29.png"/>
</div>

If we assume a more recent start, at the beginning of December 2019 (8 weeks ago), the estimates for R0 range between 2.2 and 3.5:
<div>
  <img alt="graph of R0 estimates with epidemic start 8 weeks ago" width="500" src="https://data.nextstrain.org/ncov_branching-R0-recent_2020-01-29.png"/>
</div>
This estimates are broadly consistent with those by other scientists which mostly fall between R0=2-3, see for example <a href="https://www.biorxiv.org/content/10.1101/2020.01.25.919787v1">this preprint</a>.
Importantly, R0 is a quantity that depends strongly on the socio-economic context and infection control measures.
```

# [Scientific credit](https://nextstrain.org/ncov/2020-01-30?d=map&c=author)

We would like to acknowledge the amazing and timely work done by all scientists involved in this outbreak, but particularly those working in China.
Only through the rapid sharing of genomic data and metadata are analyses such as these possible.

<br>

The nCoV genomes were generously shared by scientists at the:

 * Shanghai Public Health Clinical Center & School of Public Health, Fudan University, Shanghai, China
 * National Institute for Viral Disease Control and Prevention, China CDC, Beijing, China
 * Institute of Pathogen Biology, Chinese Academy of Medical Sciences & Peking Union Medical College, Beijing, China
 * Wuhan Institute of Virology, Chinese Academy of Sciences, Wuhan, China
 * Department of Microbiology, Zhejiang Provincial Center for Disease Control and Prevention, Hangzhou, China
 * Guangdong Provincial Center for Diseases Control and Prevention
 * Department of Medical Sciences, National Institute of Health, Nonthaburi, Thailand
 * Division of Viral Diseases, Centers for Disease Control and Prevention, USA
 * Centers for Disease Control, R.O.C., Taipei, Taiwan
 * Institut Pasteur, Paris, France

# [Detailed scientific credit](https://nextstrain.org/ncov/2020-01-30?d=map&c=author)

These data were shared via [GISAID](https://gisaid.org).
We gratefully acknowledge their contributions.

<br>

To the right we give specific sequences shared by each lab.

```auspiceMainDisplayMarkdown

The nCoV genomes were generously shared by scientists at the

 * Shanghai Public Health Clinical Center & School of Public Health, Fudan University, Shanghai, China
   - Wuhan-Hu-1/2019
 * National Institute for Viral Disease Control and Prevention, China CDC, Beijing, China
   - Wuhan/IVDC-HB-01/2019
   - Wuhan/IVDC-HB-04/2020
   - Wuhan/IVDC-HB-05/2019)
 * Institute of Pathogen Biology, Chinese Academy of Medical Sciences & Peking Union Medical College, Beijing, China
   - Wuhan/IPBCAMS-WH-01/2019
   - Wuhan/IPBCAMS-WH-02/2019
   - Wuhan/IPBCAMS-WH-03/2019
   - Wuhan/IPBCAMS-WH-04/2019
 * Wuhan Institute of Virology, Chinese Academy of Sciences, Wuhan, China
   - Wuhan/WIV02/2019
   - Wuhan/WIV04/2019
   - Wuhan/WIV05/2019
   - Wuhan/WIV06/2019
   - Wuhan/WIV07/2019
 * Department of Microbiology, Zhejiang Provincial Center for Disease Control and Prevention, Hangzhou, China
   - Zhejiang/WZ-01/2020
   - Zhejiang/WZ-02/2020
 * Guangdong Provincial Center for Diseases Control and Prevention
   - Guangdong/20SF001/2020
   - Guangdong/20SF012/2020
   - Guangdong/20SF013/2020
   - Guangdong/20SF014/2020
   - Guangdong/20SF025/2020
   - Guangdong/20SF028/2020
   - Guangdong/20SF040/2020
   - Guangdong/20SF174/2020
   - Guangdong/20SF206/2020
   - Guangdong/20SF207/2020
   - Foshan/20SF207/2020
   - Foshan/20SF210/2020
   - Foshan/20SF211/2020

 * Department of Medical Sciences, National Institute of Health, Nonthaburi, Thailand
   - Nonthaburi/61/2020
   - Nonthaburi/74/2020
 * Division of Viral Diseases, Centers for Disease Control and Prevention, USA
   - USA-WA1/2020
   - USA/AZ1/2020
   - USA/IL1/2020
   - USA/CA1/2020
   - USA/CA2/2020
 * Centers for Disease Control, R.O.C., Taipei, Taiwan
   - Taiwan/2/2020
 * Institut Pasteur, Paris, France
   - France/IDF0372/2020
   - France/IDF0373/2020

```
