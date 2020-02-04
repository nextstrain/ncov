---
title: Genomic analysis of nCoV spread. Situation report 2020-02-03.
authors: "Trevor Bedford, Richard Neher, James Hadfield, Emma Hodcroft, Misja Ilcisin, Nicola Müller"
authorLinks: "https://nextstrain.org"
affiliations: "Fred Hutch, Seattle, USA and Biozentrum, Basel, Switzerland"
date: "2020 February 03"
dataset: "https://nextstrain.org/ncov?d=map"
abstract: "This report uses publicly shared novel coronavirus (nCoV) genomic data from GISAID and Genbank to estimate rates and patterns of viral epidemic spread. We plan to issue updated situation reports as new data is produced and shared. This website is optimized for display on desktop browsers."
---
<!-- Translators: Only text after : in the above ^ needs to be translated -->
<!-- Comment tags like these do not need to be translated, they are only to help you! -->

<!-- This is left-side text -->
# [Executive summary](https://nextstrain.org/ncov)

<!-- It seems the ending '.' is needed or the whole slide 'breaks'... -->
Since the last report we have updated:
* Discussion of the latest [German sequence](http://localhost:4000/narratives/ncov/sit-rep/2020-03-03/newformat?n=7).
* Estimates of [TMRCA](http://localhost:4000/narratives/ncov/sit-rep/2020-02-03/newformat?n=9).
* Estimates of [R0](http://localhost:4000/narratives/ncov/sit-rep/2020-02-03/newformat?n=10).

_Please see the right for background links_

<!-- This is right-side text -->
```auspiceMainDisplayMarkdown
## Executive summary

Using 52</tag> publicly shared novel coronavirus (nCoV) genomes, we examined genetic diversity to infer date of common ancestor and rate of spread.
We find:
* the 52</tag> sampled genomes are very similar, differing from the consensus by 0-7 mutations
* This lack of genetic diversity has a parsimonious explanation that the outbreak descends either from a single introduction into the human population or a small number of animal-to-human transmissions of very similar viruses.
* All the sequenced cases included in this analysis likely share a common ancestor sometime between mid November and mid December 2019.
* There has been ongoing human-to-human spread since this point resulting in observed cases.
* Using estimates of total case count from Imperial College London of several thousand cases, we infer a reproductive number between 1.8 and 3.5 indicating rapid growth in the November 2019-Janurary 2020 period.


<div>
  <a href="nextstrain.org"><img alt="microscopy image of coronaviruses" width="100" src="https://nextstrain.org/static/ncov_narrative-76cfd610d11ef708d213a3170de9519f.png"/> Background on Coronaviruses </a>

  <a href="nextstrain.org"><img alt="illustration of a coronavirus" width="100" src="http://emmahodcroft.com/images/nCoV-CDC.jpg"/> Recent nCoV Outbreak </a>

  <a href="nextstrain.org"><img alt="cartoon of a phylogenetic tree" width="100" src="http://emmahodcroft.com/images/toy_alignment_mini.png"/> How to Read Phylogenies</a>


</div>

```



<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [Novel coronavirus (nCoV) 2019-2020](https://nextstrain.org/ncov)

### Further Reading:

* New China virus: Five questions scientists are asking  [Nature news](https://www.nature.com/articles/d41586-020-00166-6) _2020-01-22_
* China virus latest: first US case confirmed  [Nature news](https://www.nature.com/articles/d41586-020-00154-w) _2020-01-21_
* New virus surging in Asia rattles scientists  [Nature news](https://www.nature.com/articles/d41586-020-00129-x) _2020-01-20_
* New virus identified as likely cause of mystery illness in China [Nature news](https://www.nature.com/articles/d41586-020-00020-9) _2020-01-08_

<!-- This is right-side text -->
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

## Conspiracy Theories

A number of false ideas have been circulated about the origins of the novel coronavirus.
During outbreaks like this one, the spread of information that's known to be incorrect can lead to more panic, and cause people not to trust scientists and governments, meaning they are likely to follow advisories and take appropriate precautions.

In an effort to try and explain why these views are incorrect, we have 'debunked' these theories at the pages below:

<div>

  <a href="nextstrain.org"><img alt="picture of a snake" width="100" src="http://emmahodcroft.com/images/snake-freeToUse.jpg"/> 'Snake' Origins of nCoV </a>
  <a href="nextstrain.org"><img alt="illustration of HIV" width="100" src="http://emmahodcroft.com/images/HIV-wiki.jpg"/> 'HIV Engineering' Idea</a>


</div>

```

<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [Phylogenetic analysis](https://nextstrain.org/ncov?m=div&d=tree)

Here we present a phylogeny of 52</tag> strains of nCoV that have been publicly shared.
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
<!-- There is NO right-side text -->


<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [Phylogenetic Interpretation](https://nextstrain.org/ncov?m=div&d=tree)

We currently see little genetic diversity across the nCoV sequences, with 11</tag> out of 52</tag> sequences having no unique mutations. (Where is the 11 calculated?, seems like it should be more from the divergence tree)

<br>

Low genetic diversity across these sequences suggests that the most recent common ancestor of all nCoV sequences was fairly recent and was from many independent introductions from an animal reservoir.
<!-- since mutations accumulate slowly compared to other RNA viruses at a rate of around 1-2 mutations per month for coronaviruses.
Generally, repeated introductions from an animal reservoir will show significant diversity (this has been true for Lassa, Ebola, MERS-CoV and avian flu).
The observation of such strong clustering of human infections can be explained by an outbreak that descends from a single zoonotic introduction event into the human population followed by human-to-human epidemic spread. -->

<br>

Clusters seem to contain more recent samples, and sequences from outside of China often fall into clusters, though not always. This suggests the virus began accumulating changes as it spread in Wuhan and China, and these have subsequently been exported to other cities and countries.

<br>

<!-- We are starting to see groups of sequences that share mutations.
One cluster contains sequences from Guangdong and four isolates from the US.
Other clusters contain two to four isolates.
Sequences in these clusters tend to be from more recent samples, suggesting that the virus has started to accumulate mutations as it spread in Wuhan and subsequently to other cities. -->
There is currently no evidence that these mutations change how the virus behaves -- it is expected that RNA viruses mutate.
<!-- There is NO right-side text -->

<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [Large Cluster](https://nextstrain.org/ncov/2020-02-03?m=div&d=tree)

Perhaps this slide would zoom in on the largest cluster, maybe with map view.

<br>

We could perhaps discuss what it means that both Chinese and exported samples are in this cluster, and share different combinations of mutations.

<br>

Three of the sequences from Guangdong (ending F025, F013, and F012) are [known to come from a single family](https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(20)30154-9/fulltext), and almost certainly represent human-to-human transmission.
<!-- There is NO right-side text -->

<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [Cluster comprised entirely from sequences outside of China](https://nextstrain.org/ncov/2020-02-03?m=div&d=tree)

On the left, a cluster of 5 sequences consisting of isolates from the USA, France, Taiwan, and Australian is shown.

<br>

Sequences that are sampled at the same point of an outbreak are more likely to share a recent common ancestor than sequences sampled at different points in an epidemic.
With a lot of newer sequences originating outside of China, it is expected that clusters of sequences that are entirely comprised of sequences from outside of China are observed.
The most parsimonious explanation to why these sequences would cluster together, is the absence of sequences from China, rather than direct transmission between them outside of China.
In fact, the two French sequences in this cluster are [known to be from the same family](https://www.thelocal.fr/20200129/coronavirus-in-france-what-you-need-to-know) - a Chinese couple from Wuhan.

<br>

<!-- There is NO right-side text -->

<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [New German Sequence](https://nextstrain.org/ncov/2020-02-03?m=div&d=tree&s=Germany/BavPat1/2020)

Here we could say something about the latest German sequence

<br>

This came from one of the 8 employees of a small company in Bavaria who were infected by a visitor from Wuhan who was not symptomatic during the visit.

None of those infected have had severe symptoms, but this does provide some suggestion that the infection can be transmitted before and/or without symptoms appear.

<!-- There is NO right-side text -->

<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [Cases outside China](https://nextstrain.org/ncov/2020-02-03?c=country&d=tree&m=div)

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
<!-- There is NO right-side text -->

<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [Dating the time of the most recent common ancestor](https://nextstrain.org/ncov/2020-02-03?d=tree)
The high similarity of the genomes suggests they share a recent common ancestor (i.e. that they have descended from the same ancestral virus recently).
Otherwise, we would expect a higher number of differences between the samples.

<br>

Previous research on related coronavirus suggests that these viruses accumulate between 1 and 4 changes in their genome per month (rates of 3 &times; 10<sup>-4</sup> to 2 &times; 10<sup>-3</sup> per site per year).

<br>

On the right, we explore how different assumptions about the rate of change, and the observed genetic diversity, give us estimates about when all sequenced cases last shared a common ancestor.

<!-- This is right-side text -->
```auspiceMainDisplayMarkdown
## Date of the common ancestor of outbreak viruses
The time of the most recent common ancestor (or tMRCA) of a set of sequenced cases denotes when these sequenced cases last shared a common ancestor.
This time can be as early as the time when a virus first entered the human population, but can also be substantially later as show in the figure below.

<div>
  <img alt="Example phylogeny where the time of the initial zoonosis is different from the most recent common ancestor of several sequenced cases" width="500" src="https://raw.githubusercontent.com/nicfel/nCov-Nicola/master/figures/zoonosis.png"/>
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
[Analyzing](https://nicfel.github.io/nCov-Nicola/ExponentialCoalescent_2020203.html) this data assuming once a constant population of infected individuals and once and exponentially growing population estimates the tMRCA to be anywhere between the mid-November and the mid-December.

<div>
  <img alt="TMRCA estimates inferred from genetic sequence data and sampling times in BEAST" width="500" src="https://raw.githubusercontent.com/nicfel/nCov-Nicola/master/figures/beast_coal-tmrca_2020203.png"/>
</div>


More detailed modeling of the onset of the outbreak are ongoing and with more and more sequence data being made publicly available, the quality of these estimates will improve.

```

<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [Estimating the growth rate](https://nextstrain.org/ncov/2020-02-03?d=tree)

An important quantity in the spread of a pathogen is the average number of secondary cases each infection produces.

<br>

This number is known as R0 ("R-zero" or "R-nought").
One the right, we present simple estimates of R0.

<!-- This is right-side text -->
```auspiceMainDisplayMarkdown
## Estimates of epidemic growth rate
Scientists at Imperial College London have used the number of cases observed outside of China to estimate the [total number of cases](https://www.imperial.ac.uk/mrc-global-infectious-disease-analysis/news--wuhan-coronavirus/) and suggested that there have been at least several thousand cases by 2020-01-22.
With the additional exported cases since and the continued growth of confirmed cases in China, we currently have to expect at least 50000 cases to date.
Together with our previous estimates of the age of the outbreak and information on the infectious period, we can estimate plausible ranges of R0 a [birth-death model](https://en.wikipedia.org/wiki/Birth%E2%80%93death_process).
These estimates can be sensitive to, amongst other factors, the date of the onset of the outbreak, as well as to the initial number of infected and should therefore be treated with some caution.

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

<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [Scientific credit](https://nextstrain.org/ncov/2020-02-03?d=map&c=author)

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
 <!-- There is NO right-side text -->

<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [Detailed scientific credit](https://nextstrain.org/ncov/2020-02-03?d=map&c=author)

These data were shared via [GISAID](https://gisaid.org).
We gratefully acknowledge their contributions.

<br>

To the right we give specific sequences shared by each lab.

<!-- This is right-side text -->
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
