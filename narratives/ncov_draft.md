---
title: DRAFT nCoV narrative
authors: "Bedford, Neher, Hadfield, Mueller et al (TODO)"
authorLinks: "https://nextstrain.org"
affiliations: "Fred Hutch, Seattle, USA and Biozentrum Basel, Switzerland"
date: "23rd January 2020"
dataset: "https://nextstrain.org/ncov?d=map"
abstract: "This is a DRAFT narrative containing our commentary of the ongoing novel (wuhan) Coronavirus."
---


# [Coronaviruses](https://nextstrain.org/ncov)

### Further Reading:

* General information on Coronaviruses on [wikipedia](https://en.wikipedia.org/wiki/Coronavirus) _2020-01-23_
* Material provided by the [US CDC](https://www.cdc.gov/coronavirus/index.html) _2020-01-23_
* Organization and genome on [viralzone](https://viralzone.expasy.org/764?outline=all_by_species) _2020-01-23_

```auspiceMainDisplayMarkdown

## Different human coronaviruses

Coronaviruses are members of a diverse species of (+)ssRNA viruses which have a history of causing respiratory infections in humans.
Some variants of coronaviruses are associated with outbreaks, others are continuously circulating and cause mostly mild respiratory infections (e.g. the common cold).

#### SARS-CoV
The most well known of these is the SARS (SARS-CoV) outbreak in 2002-2003 which spread around the world and resulted in [over 8000 cases and 774 deaths](https://www.theguardian.com/world/2017/dec/10/sars-virus-bats-china-severe-acute-respiratory-syndrome), with a fatality rate of around 9-11%.

#### Seasonal CoV
However it's not fair to categorize all coronaviruses as deadly as SARS-CoV and MERS-CoV.
There are now [seven different groups of coronaviruses](https://www.cdc.gov/coronavirus/types.html) which can cause infections in humans, with four of them considered "seasonal" and commonly infecting humans each year.
Compared with SARS, these seasonal coronavirus strains are ["much more prevalent, much less severe, and common causes of influenza‐like illness (ILI)"](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5820427/).
In fact, [5](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2879166/) - [12](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5820427/)% of all ILI cases test positive for coronaviruses, so they are rather common.

#### MERS-CoV

In 2012, a novel coronavirus (MERS-CoV) causing severe respiratory symptoms was identified. MERS has resulted in fatalities comparable to SARS, however the transmission route of MERS is very different. Whereas SARS was efficiently spread from one human to another, human MERS infections were generally a result of independent zoonoses (animal to human transmissions) from camels (see [Dudas _et al._](https://elifesciences.org/articles/31257) for more information).

#### Animal reservoirs
Coronaviruses infect a wide range of animals, and the human outbreaks described above are a result of one or more "jumps" from these animal reservoirs into the human population.
SARS was originally believed to have arrived in the human population from Masked palm civets, but is now considered to have arrived [from horseshoe bats](https://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1006698).


#### Human to Human transmission
The ability for different lineages to be transmitted between humans is extremely important to understand the potential ramifications of an outbreak.
Due to the ability of SARS to spread between humans and the high case fatality rate, SARS (or a SARS-like virus) is considered a [global public health threat](https://www.who.int/whr/2007/overview/en/index1.html) by the WHO.

```


# [novel coronavirus (nCoV) 2019-2020](https://nextstrain.org/ncov)

### Further Reading:

* New China virus: Five questions scientists are asking  [Nature news](https://www.nature.com/articles/d41586-020-00166-6) _2020-01-22_
* China virus latest: first US case confirmed  [Nature news](https://www.nature.com/articles/d41586-020-00154-w) _2020-01-21_
* New virus surging in Asia rattles scientists  [Nature news](https://www.nature.com/articles/d41586-020-00129-x) _2020-01-20_
* New virus identified as likely cause of mystery illness in China [Nature news](https://www.nature.com/articles/d41586-020-00020-9) _2020-01-08_

```auspiceMainDisplayMarkdown

## Recent outbreak of a novel coronavirus
In December 2019, a new illness was first detected in Wuhan, China.
We now know this to another outbreak of coronavirus in humans (the 7th), and is being called nCoV (novel coronavirus).

As of January 23rd, over 624 cases and 17 deaths [have been reported](https://en.wikipedia.org/wiki/2019%E2%80%9320_outbreak_of_novel_coronavirus_(2019-nCoV)).
It's still too early to know the case fatality rate, but early indications are that it is less than SARS-CoV.
The case counts are dramatically rising in part due to increased surveillance and testing, with 131 new cases reported yesterday (Jan 22nd).

While the outbreak seems to be centered in Wuhan, which is now [under quarantine](https://twitter.com/PDChina/status/1220060879112282117), the virus has spread throughout China and abroad, including Hong Kong, Macau, Thailand, USA, Japan and South Korea.


The origin of the virus is still unclear, however [genomic analyis](http://virological.org/t/ncovs-relationship-to-bat-coronaviruses-recombination-signals-no-snakes/331) suggests nCoV is most closely related to viruses previously identified in bats.
It is plausible that there were other intermediate animal transmissions before the introduction into humans.

#### Nextstrain narratives

The following pages contain analysis performed using [nextstrain](https://nextstrain.org).
Scrolling through the left hand sidebar will reveal paragraphs of text with a corresponding visualization of the genomic data on the right hand side.

To have full genomes of a novel & large RNA virus this quickly is a remarkable achievement.
These analyses have been made possible by the rapid and open sharing of genomic data and interpretations by scientists all around the world (see the final slide for a visualization of sequencing authorship).


```







# [How to interpret the phylogenetic trees](https://nextstrain.org/ncov)

```auspiceMainDisplayMarkdown
## Transmission trees vs phylogenetic trees

Pathogens spread by through rapid replication in one host followed by transmission to another host and an outbreak can only take off when one infection results in more than one subsequent infections.
As the pathogen replicates and spreads, its genome needs to be replicated many times and random mutations -- copying mistakes -- will accumulate in the genome.
Such random mutations can help to track the spread of the pathogen and learn about its transmission routes and dynamics.

<div>
  <img alt="tmrca" width="500" src="https://neherlab.org/talk_images/infection_tree_combined.png"/>
</div>

The illustration above shows a sketch of transmission tree with a subset of cases that were sampled (blue).
In practice, the transmission tree is unknown and typically only rough estimates of case counts are available.
Genome sequences allow us to infer parts of the transmission tree.
In this example, three mutations (little diamonds) are indicated on the tree and these mutations allow us to group samples into clusters of closely related viruses that belong to the same transmission chains.

```

# [Phylogenetic analysis](https://nextstrain.org/ncov?d=tree)

Here we present a maximum likelihood phylogeny of 24 strains of nCoV that have been publicly shared.
Information on how the analysis was performed is available [in this GitHub repository](github.com/nextstrain/ncov).
Sequencing the genome of a large novel RNA virus in an evolving outbreak situation is challenging and not all differences observed in these sequences might represent actual mutations but could instead be due to sequencing errors.
Insertions and deletions as well as differences at the ends of the genome are more likely to be errors and we therefore masked those for the purposes of this analysis.

The colours represent the city of isolation, with the x-axis representing nucleotide divergence (you can hover over the branches to see which mutations are present on a given branch.)
Divergence is measured as the number of changes per site.
Since the genome is 29'000 bases long, one mutation corresponds to a distance 1/29'000 = 0.0000335.

* We currently see little genetic diversity across the nCov sequences with 8 out of 24 sequences being genetically identical.

* Low genetic diversity across these sequences suggest that the most recent common ancestor of all nCov sequences was fairly recent.

* At the moment, most mutations that can be observed are singletons -- that is they are unique to individual genomes. Only sequences form two clusters from Guangdong share mutations -- we will explore in the following slides.

Previous estimates of the rate at which mutations accumulate in Coronaviruses suggest a rate of about 1-3 changes per month.
The sequence data available to date therefore points towards an common ancestor of the sampled cases between early November to early December 2019.
While the sequences provide allow to make inferences about the age of the outbreak, they don't allow us to differentiate between multiple transmission from a homogeneous animal reservoir and a single animal-human transmission followed by spread among humans.


# [Potential within-family transmission I](https://nextstrain.org/ncov?d=tree&f_city=Shenzhen)

Of the four isolates from Shenzhen (Southeastern China) we see three isolates which are genetically identical and share three mutations unique to those three samples.
These three samples are [known to come from a single family](https://twitter.com/JingLu_LuJing/status/1220143773532880896) and likely represent human to human transmission.

The fourth sample does not seem to be related to the other three or to any other of the available sequences.
Its genome has one mutation not seen in any other genome.

# [Potential within-family transmission II](https://nextstrain.org/ncov?d=tree&f_city=Zhuhai)

Similarly, there are two genomically-identical isolates from Zhuhai (Southeastern China) which form a cluster, sharing one unique mutation seen in no other isolate.
These two cases are also [known to come from a single family](https://twitter.com/JingLu_LuJing/status/1220143773532880896) again indicating likely human to human transmission.

# [Cases outside china](https://nextstrain.org/ncov?c=country&d=tree,map)

There are reported nCoV cases in Thailand, USA, Japan and South Korea.
These cases are all linked to Wuhan and there we are not aware of evidence for transmission of nCov in these countries.


The only currently available sequence data for cases outside of China is for the two cases from Thailand, which are coloured here in Red.
These samples are genetically identical to sequences isolated in Wuhan, as well as 6 other genomes sequences from Chinese cases.


# [Dating the zoonosis](https://nextstrain.org/ncov?d=tree)
The high similarity if the genomes available to date suggest they have a recent common ancestor -- otherwise we would expect that a higher number of differences.
Previous research on related coronavirus suggest that these viruses accumulate between 1 and 3 changes in their genome per month (rates of 0.003 - 0.001/site and year).
On the right, we explore the implications of the observed genetic diversity and these rate estimates for the timing of the outbreak.

```auspiceMainDisplayMarkdown

Here, we use this star-like structure along with a Poisson distribution of mutations through time to estimate the time of the most recent common ancestor of sequenced viruses:

<div>
  <img alt="tmrca" width="500" src="http://data.nextstrain.org/ncov_poisson-tmrca.png"/>
</div>

```

# [Scientific credit](https://nextstrain.org/ncov?d=map&c=author)

We would like to acknowledge the amazing and timely work done by all scientists involved in this outbreak, but particularly those working in China.
Only through the rapid sharing of genomic data and metadata are analyses such as these possible.

.

The nCoV genomes were generously shared by scientists at the Shanghai Public Health Clinical Center & School of Public Health, Fudan University, Shanghai, China (Wuhan-Hu-1/2019), at the National Institute for Viral Disease Control and Prevention, China CDC, Beijing, China (Wuhan/IVDC-HB-01/2019, Wuhan/IVDC-HB-04/2020, Wuhan/IVDC-HB-05/2019) at the Institute of Pathogen Biology, Chinese Academy of Medical Sciences & Peking Union Medical College, Beijing, China (Wuhan/IPBCAMS-WH-01/2019, Wuhan/IPBCAMS-WH-02/2019, Wuhan/IPBCAMS-WH-03/2019, Wuhan/IPBCAMS-WH-04/2019), at the Wuhan Institute of Virology, Chinese Academy of Sciences, Wuhan, China (Wuhan/WIV02/2019, Wuhan/WIV04/2019, Wuhan/WIV05/2019, Wuhan/WIV06/2019, Wuhan/WIV07/2019), at the Department of Microbiology, Zhejiang Provincial Center for Disease Control and Prevention, Hangzhou, China (Zhejiang/WZ-01/2020, Zhejiang/WZ-02/2020), at the Guangdong Provincial Center for Diseases Control and Prevention (Guangdong/20SF012/2020, Guangdong/20SF013/2020, Guangdong/20SF014/2020, Guangdong/20SF025/2020, Guangdong/20SF028/2020, Guangdong/20SF040/2020) and at the Department of Medical Sciences, National Institute of Health, Nonthaburi, Thailand (Nonthaburi/61/2020, Nonthaburi/74/2020) via GISAID. We gratefully acknowledgement their contributions.
