---
title: Genomic analysis of COVID-19 spread. Situation report 2020-03-02.
authors: "Trevor Bedford, Richard Neher, James Hadfield, Emma Hodcroft, Misja Ilcisin, Nicola Müller"
authorLinks: "https://nextstrain.org"
affiliations: "Fred Hutch, Seattle, USA and Biozentrum, Basel, Switzerland"
date: "2020 March 03"
dataset: "https://nextstrain.org/ncov/2020-01-30"
abstract: "This report uses publicly shared COVID-19 genomic data from GISAID and Genbank to estimate rates and patterns of viral epidemic spread. We plan to issue updated situation reports as new data is produced and shared. This website is optimized for display on desktop browsers."
---

# [Executive summary](https://nextstrain.org/ncov/2020-01-30)

```auspiceMainDisplayMarkdown
## Executive summary

Using 146</tag> publicly shared COVID-19 genomes, we examined genetic diversity to characterize the spread of COVID-19 in different areas and infer date of common ancestor.
We find:
* the genetic sequences support the hypothesis of sustained spread of COVID-19 in the greater Seattle area since a few weeks.
* COVID-19 was introduced into Italy at least twice with subsequent community spread.
* All the sequenced cases included in this analysis likely share a common ancestor sometime between mid November and mid December 2019.
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

As of March 2, 2020 there are 90,912 cases and 3,117 deaths [have been reported](https://gisanddata.maps.arcgis.com/apps/opsdashboard/index.html#/bda7594740fd40299423467b48e9ecf6).
It's still too early to know the case fatality rate, but early indications are that it is significantly less than SARS-CoV.
The case counts are dramatically rising in part due to increased surveillance and testing.

While the outbreak seems to be centered in Wuhan, which is now [under quarantine](https://twitter.com/PDChina/status/1220060879112282117), the virus has spread throughout China and abroad, including Hong Kong, Singapore, Japan, and Thailand, as well as Europe, North America, South Asia, the Middle East, and Australia. Limited local transmission outside of China has been reported.

The Nextstrain team has been fielding a lot of questions about the origins of this virus and what steps individuals and families should take to prepare themselves. See the end of this narrative for a list of our most commonly asked questions and answers from sources we trust.

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

### Reading a typed Phylogenetic Tree

Phylogenetic trees often contain additional information, such as where geographically individual sequences were isolated from.
Additionally, possible locations of internal nodes can be inferred using mathematical models as well.
Interpreting these should, however, be done with caution, as the sampling and sequencing or lack thereof can significantly influence the interpretation,

In the following example, we first show fully sampled phylogenetic tree, with samples from two different locations denoted by orange and blue.
In the fully sampled case on the right, our interpretation of what happened, was that there were three different introductions from orange to blue.
When removing the one orange sequence in the middle, our interpretation is now that there was one introduction into blue that happened much earlier.
In the last example, we have only one sequence from orange, which could lead us to think that there was one introduction from orange into blue.


<div>
  <img alt="Example phylogeny where all or only a subset of cases are included in the final phylogeny" width="500" src="https://raw.githubusercontent.com/nicfel/nCov-Nicola/master/figures/introductions.png"/>
</div>

Overall, the inferred locations of where a lineage has been in the past, should be considered as highly uncertain.

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

We currently have sequenced cases from five different continents.
While the early cases where all directly linked to cases in Wuhan, associated with the seafood market outbreak, we now observe various different cases that show evidence for community spread or were imported from sources outside China.

<br>
Generally, repeated introductions from an animal reservoir will show significant diversity (this has been true for Lassa, Ebola, MERS-CoV and avian flu).
The observation of such strong clustering of human infections can be explained by an outbreak that descends from a single zoonotic introduction event into the human population followed by human-to-human epidemic spread.

<br>

# [At least two introductions with community spread into italy](https://nextstrain.org/ncov/2020-02-03?m=div&d=tree)

We currently have 3 sequences from Italy, two of which from the Rome area and one from Lombardy in northern Italy.
These 3 sequences share a common ancestor early in the epidemic, which strongly suggests that there were at least two introductions with community spread into Italy.

<br>

The two sequences from the Rome are cluster together and therefore indicate local spread.
The sequence from Lombardy (Italy/CDG1/2020) clusters together with sequences with known travel history to Italy that were most likely infected in Italy.

<br>

# [Export of COVID-19 from Italy to other countries in Europe and Brazil](https://nextstrain.org/ncov/2020-02-03?m=div&d=tree)

There are several cases of COVID-19 with known travel history to Lombardy.


<br>

We can observe these patterns in the phylogenetic tree.
Brazil/SPBR-01/2020 and Finland/FIN-25/2020 for example have a known travel history to Italy and cluster together with a sequence from Lombardy, Italy/CDG1/2020.

<br>


<!-- This is left-side text -->
# [Likely spread of SARS-Cov 2 in the Seattle area](https://nextstrain.org/ncov/2020-02-03?m=div&d=tree)

There are now several cases of COVID-19 reported in the greater Seattle area and the US in general.
The newly isolated and sequenced case from Snohomish County (north of Seattle), is genetically closely related to a case isolated at the end of January in the same area.

<br>

There are two explanations to why this is the case.
Either, the virus was introduced at least twice into the greater Seattle area from a common source in China.
The far more likely explanation, however, is that the virus was circulating undetected in the area for a while.

<br>



# [Dating the time of the most recent common ancestor](https://nextstrain.org/ncov/2020-01-30?d=tree)

The time of the most recent common ancestor (or tMRCA) of a set of sequenced cases denotes when these sequenced cases last shared a common ancestor.
This time can be as early as the time when a virus first entered the human population, but can also be substantially later as show in the figure below.

<div>
  <img alt="Example phylogeny where the time of the initial zoonosis is different from the most recent common ancestor of several sequenced cases" width="500" src="https://raw.githubusercontent.com/nicfel/nCov-Nicola/master/figures/zoonosis.png"/>
</div>



```auspiceMainDisplayMarkdown
## Date of the common ancestor of outbreak viruses

Several research groups and people have estimated the time of the most recent common ancestor.
http://virological.org/t/phylodynamic-analysis-129-genomes-24-feb-2020/356
http://virological.org/t/evolutionary-epidemiological-analysis-of-93-genomes/405/2

The common ancestor of all sequences is most likely between mid November and mid December.
This would be consistent with all currently sequenced cases descending from the [initial cluster of cases at the Wuhan seafood market](http://virological.org/t/phylodynamic-analysis-129-genomes-24-feb-2020/356).

<div>
  <img alt="graph of TMRCA estimates based on different mutation rates" width="500" src="https://data.nextstrain.org/ncov_poisson-tmrca.png"/>
</div>

Using the entire data set, the nextstrain analysis pipeline estimates that the common ancestor most likely existed between late-Nov and the beginning of December 2019.

There is a [confirmed case in Wuhan with onset date of December 1, 2019](https://twitter.com/trvrb/status/1220749265380593664), which would put an upper bound on the date of most recent common ancestor.
The common ancestor of viruses sequenced to date might be later than this date though.

More detailed modeling of the onset of the outbreak are ongoing.
Despite considerable uncertainty, our best guess is remains late November/early December 2019.

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

# [The Nextstrain team's take on COVID-19 frequent questions](https://nextstrain.org/ncov/2020-01-30/2020-01-30)

As a team, we’ve been fielding a lot of questions about COVID-19 and how individuals and families should be preparing. We have compiled a list of the most common questions we receive and answered them with information and recommendations from resources that we trust. 

With all of this, it is important to note that as COVID-19 is a novel virus, we are still learning about disease characteristics including who is most likely to become infected or experience serious illness as a result of it. Current clinical knowledge is based primarily on data from the China outbreak. We believe everyone should be taking steps to protect themselves and their families. We hope the resources below will help you to do that. 

```auspiceMainDisplayMarkdown

## COVID-19 FAQs

#### How quickly is COVID-19 spreading?
Based on modeling estimates and phylodynamic analysis, the current epidemic doubling time is about 7 days. There is uncertainty with these estimates and this rate will vary based on community contact patterns. Nonpharmaceutical interventions can be effective in reducing this rate, as we have [seen in China](https://twitter.com/jodigraphics15/status/1232484344872071168?s=20). Sudden jumps in case counts are likely due to increased testing capacity catching cases that were not previously missed by surveillance systems.

Current case counts, deaths, and recoveries are being tracked in a number of places. These are updated as regularly as possible, but may not always be immediately up to date. One dashboards that we particularly like is from the JohnsCenter for Systems Science and Engineering:
[Johns Hopkins CSSE Dashboard](https://gisanddata.maps.arcgis.com/apps/opsdashboard/index.html#/bda7594740fd40299423467b48e9ecf6).

#### How is COVID-19 spread?
COVID-19 is a respiratory virus, which spread in 3 main ways:
1. Between people who are in close contact with each other (within about 6 feet, or 2 metres), through small aerosolized particles that can be produced when a person talks, coughs, or sneezes. 
2. In (relatively) larger droplets that are produced someone coughs or sneezes (these can travel over similar distances)
3. Contact with surfaces or objects infected with the virus (also known as fomites). If a person touches an infected surface and then touches their eyes, nose, or mouth without washing their hands they may become sick by introducing the virus into their body.

For COVID-19 there is evidence of transmission in clinically mild and asymptomatic cases, which means people who are infected with COVID-19 and able to transmit the virus to others, may not appear, or even feel, sick.

#### Who is at the highest risk for COVID-19?
Studies of patients in China have shown that most COVID-19 cases are adults, just 2.1% of 44,672 cases in China were in individuals below 20 years of age and just 1% of those are in children under 10. Approximately 80% of people experience mild illness, 14% experience severe disease, and 5% were critically ill. There is evidence of asymptomatic cases. 

Based on case reports, individuals over 60 and those who have weakened immune systems or comorbidities like chronic lung disease, diabetes, and cardiovascular disease are at a higher risk for severe disease. Children do not seem to be at a high risk for either contracting COVID-19 or having severe disease as a result of an infection. A great graphic from [Jessie Bloom and his team](https://research.fhcrc.org/bloom/en.html that illustrates this point can be viewed [here] (https://raw.githubusercontent.com/jbloom/CoV_vs_flu_CFR/master/CFR-stats.jpg).

More detailed information on clinical characteristics can be found in the original journal article [Characteristics of and Important Lessons From the Coronavirus Disease 2019 (COVID-19) Outbreak in China by Dr. Zunyou Wu and Dr. Jennifer M. McGoogan](https://jamanetwork.com/journals/jama/fullarticle/2762130).


#### What are the symptoms of COVID-19?
The most commonly reported symptoms are fever, dry cough, and shortness of breath.[^1] Again, these symptoms will be mild for about 80% of cases. COVID-19 is a lower respiratory infection (so most symptoms will appear “below” the voicebox) but has flu-like symptoms and how its symptoms present varies from person to person. If you are concerned about symptoms you or someone else may be having, we recommend calling your doctor to discuss next steps.

[^1] Again, more detailed information on clinical characteristics can be found in the original journal article [Characteristics of and Important Lessons From the Coronavirus Disease 2019 (COVID-19) Outbreak in China by Zunyou Wu, MD, PhD and Jennifer M. McGoogan, PhD](https://jamanetwork.com/journals/jama/fullarticle/2762130).

The full report of the WHO-China Joint Mission on Coronavirus Disease 2019 (COVID-19O), can be found [here](https://www.who.int/docs/default-source/coronaviruse/who-china-joint-mission-on-covid-19-final-report.pdf). 

#### What can be done to prevent or reduce the spread of COVID-19?
COVID-19 is a respiratory virus and nonpharmaceutical interventions are the most effective tool we currently have to prevent and slow the spread of COVID-19 at this time. There is not currently a vaccine or specific treatment for COVID-19. Some ways to implement nonpharmaceutical interventions include:

- Practice [social distancing](https://en.wikipedia.org/wiki/Social_distancing), such as limiting attendance at events with large groups of people
- Stay home, especially if you are feeling ill.
- Implement good hand washing practices - it is extremely important to wash hands regularly  
- Cover coughs and sneezes in your elbow or tissue
- Avoid touching your eyes, nose, and mouth with unwashed hands
- Disinfect frequently touched surfaces, such as doorknobs
- Begint to take your temperature every day, if you develop a fever, stay home (if possible, separated from other people you live with) and call your doctor
- Start preparations in anticipation of social distancing and potential supply chain shortages. This includes ensuring you have sufficient supplies of prescription medicine and have about a 2 week supply of food and other necessary household goods. If possible, discuss a plan to work from home with your employer or how classes could be attended remotely. 

With these steps  in mind, it is important _to not panic and panic buy_. Panic buying unnecessarily increases strain on supply chains and can make it difficult to ensure that everyone is able to get supplies that they need. For an excellent explanation on why we should all be having a proactive, but measured response to the situation, we recommend reading [this twitter thread](https://twitter.com/firefoxx66/status/1233666678841597952?s=20) by Nextstrain team member Emma Hodcroft, PhD. 

Some other resources and articles we have found extremely useful during our own preparations are: 

[Past Time to Tell the Public: “It Will Probably Go Pandemic, and We Should All Prepare Now" by Jody Lanard and Peter M. Sandman](https://virologydownunder.com/past-time-to-tell-the-public-it-will-probably-go-pandemic-and-we-should-all-prepare-now/).

[Ready - Pandemic Preparedness Page](https://www.ready.gov/pandemic) - a helpful public service campaign on disaster preparedness from the US Department of Homeland Security.  

[So you think you’re about to be in a pandemic? by Ian McKay, PhD and Katherine E Arden, PhD](https://virologydownunder.com/so-you-think-youve-about-to-be-in-a-pandemic/).

[US CDC Guidelines specifically on prevention of COVID-19 spread](https://www.cdc.gov/coronavirus/2019-ncov/community/index.html). 

[US CDC Guidelines on Nonpharmaceutical Interventions](https://www.cdc.gov/nonpharmaceutical-interventions/index.html).

#### What should you do if you think you might be sick with COVID-19?
If you are experiencing influenza like illness, call your doctor before going to the clinic or emergency room. If your healthcare providers are concerned that you might have COVID-19, they will have specific steps for you to follow when you arrive at the clinic that will reduce the risk of getting other sicks or contracting an infection from someone already at the clinic. If you show up to a clinic or hospital without giving advance warning, you will expose yourself to a significant number of sick people and may expose others to your illness, potentially increasing transmission. 

We also recommend that people limit unnecessary trips to the emergency room. As the number of cases increases, there is strain put on the healthcare system. Reducing the number of people arriving at clinics and hospitals will decrease this strain and health keep you away from people who may be infected with COVID-19 or other respiratory viruses. If you are experiencing a health issue that you are concerned about, we recommend contacting your doctor’s office to determine where and when the best place for you to seek care is. 

#### What is the origin of the COVID-19 virus?
The origin of the virus is still unclear, however [genomic analyis](https://virological.org/t/ncovs-relationship-to-bat-coronaviruses-recombination-signals-no-snakes/331) suggests nCoV is most closely related to viruses previously identified in bats. It is plausible that there were other intermediate animal transmissions before the introduction into humans. There is no evidence of snakes as an intermediary.

There is not evidence that the COVID-19 outbreak is the result of a lab escape. 
The data we have on the COVID-19 outbreak is consistent with a zoonotic origin and inconsistent with a lab escape scenario. For a full analysis the evidencing supporting a zoonotic origin over a lab escape origin, please read [this twitter thread](https://twitter.com/trvrb/status/1230634136102064128?s=20) by Nextstrain team member Trevor Bedford, PhD, which discusses this issue in-depth.

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
