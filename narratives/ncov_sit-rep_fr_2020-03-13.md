---
title: Analyse génomique de la propagation du SARS-Cov-2. Rapport de situation 2020-03-13.
authors: "Emma Hodcroft, Nicola Müller, Cassia Wagner, Misja Ilcisin, James Hadfield, Sidney M. Bell, Richard Neher, Trevor Bedford. (Translation: Maxime Morin, Meriem El Karoui, Etienne Simon-Loriere)"
authorLinks: "https://nextstrain.org"
affiliations: "Fred Hutch, Seattle, USA; Biozentrum, Basel, Switzerland; CZI, CA, USA"
date: "13 Mars 2020"
dataset: "https://nextstrain.org/ncov/2020-03-13?d=map&legend=closed"
abstract: "Ce rapport utilise des données génomiques partagées publiquement pour suivre la propagation du SARS-Cov-2. Nos rapports sont mis à jour chaque semaine."
---
<!-- Translators: Only text after : in the above ^ needs to be translated -->
<!-- Comment tags like these do not need to be translated, they are only to help you! -->
<!-- Ensure that links always end in a 'letter' (. counts) If some kind of text doesn't follow them, it breaks the slide. -->
<!-- numbers can be tagged ilke this: 161</tag> - this is just for us to help find them to update! Just leave in the </tag> bit. -->

<!-- This is left-side text -->
# [Table des matières](https://nextstrain.org/ncov/2020-03-13?d=tree,map&p=grid)

* [Ressources sur COVID-19](https://nextstrain.org/narratives/ncov/sit-rep/fr/2020-03-13?n=2).
* [Une note sur les questions d'échantillonnage](https://nextstrain.org/narratives/ncov/sit-rep/fr/2020-03-13?n=3).
* [Circulation en Europe](https://nextstrain.org/narratives/ncov/sit-rep/fr/2020-03-13?n=4).
* [Transmission locale dans les îles britanniques et en Irlande](https://nextstrain.org/narratives/ncov/sit-rep/fr/2020-03-13?n=5).
* [Propagation du SRAS-CoV-2 en provenance d'Iran](https://nextstrain.org/narratives/ncov/sit-rep/fr/2020-03-13?n=6).
* [Introductions aux Etats-Unis](https://nextstrain.org/narratives/ncov/sit-rep/fr/2020-03-13?n=7).
* [Propagation du SRAS-CoV-2 dans l'État de Washington (USA)](https://nextstrain.org/narratives/ncov/sit-rep/fr/2020-03-13?n=8).
* [Propagation du SRAS-CoV-2 en Californie](https://nextstrain.org/narratives/ncov/sit-rep/fr/2020-03-13?n=9).
* [Ce que vous pouvez faire](https://nextstrain.org/narratives/ncov/sit-rep/fr/2020-03-13?n=10).
* [FAQs & Idées reçues](https://nextstrain.org/narratives/ncov/sit-rep/fr/2020-03-13?n=11).
* [Crédit scientifique](https://nextstrain.org/narratives/ncov/sit-rep/fr/2020-03-13?n=12).

<!-- This is right-side text -->
```auspiceMainDisplayMarkdown
# Résumé exécutif

Ici, nous avons analysé 410 </tag> génomes du SARS-Cov-2 partagés publiquement. En comparant ces génomes viraux entre eux, nous pouvons caractériser comment Sars-Cov-2 évolue et se déplace dans le monde.

Pour un aperçu actuel du nombre de cas de coronavirus dans le monde, voir [Our World In Data (en anglais)](https://ourworldindata.org/coronavirus).

Dans ce rapport, nous montrons que le virus circule largement à travers le monde, avec des preuves de transmission locale sur plusieurs continents. En ce moment, nous demandons instamment aux autorités de se  concentrer sur les efforts pour ralentir la propagation au sein des communautés; les interdictions de voyager sont moins susceptibles d'être efficaces.

Dans les mises à jour de cette semaine, nous rapportons:  

* COVID-19 circule largement à travers l'Europe, avec des mouvements importants entre les pays.

* Nous identifions au moins 4 introductions au Royaume-Uni, certaines avec transmission communautaire continue.

* Il y a eu un certain nombre de cas liés aux voyages reliant l'Iran à d'autres parties du monde.

* Il y a eu à ce jour de nombreuses introductions aux États-Unis, ce qui a résulté en des chaînes de transmission locales dans plusieurs États.

* L'épidémie continue de croître dans l'État de Washington; certains cas sont étroitement liés à ceux du navire de croisière Grand Princess.

* Il y a une circulation locale de COVID-19 en Californie.

* Des mesures de distanciation sociale doivent être adoptées rapidement pour alléger la charge pesant sur les systèmes de santé et protéger les personnes vulnérables.
```


<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [Ressources sur COVID-19](https://nextstrain.org/ncov/2020-03-05)

Nous avons préparé quelques ressources qui valent la peine d'être lues pour vous familiariser avec COVID-19 et le virus qui le provoque, SARS-CoV-2.

Ces informations faciliteront l'interprétation des données que nous présentons dans ce document. Si vous n'êtes pas familier avec les arbres phylogénétiques, nous vous encourageons à consulter le ["Comment lire les phylogénie" en anglais](https://nextstrain.org/narratives/trees-background) et à revenir quand vous serez prêt.

<!-- This is right-side text -->
```auspiceMainDisplayMarkdown

## Contexte

<div>
  <a href="https://nextstrain.org/help/coronavirus/human-CoV"><img alt="image de microscopie des coronavirus" width="100" src="https://nextstrain.org/static/ncov_narrative-76cfd610d11ef708d213a3170de9519f.png"/> Informations sur les coronavirus </a>

  <a href="https://nextstrain.org/help/coronavirus/SARS-CoV-2"><img alt="illustration d'un coronavirus" width="100" src="http://data.nextstrain.org/img_nCoV-CDC.jpg"/> Contexte de l'émergence de COVID-19 </a>

  <a href="https://nextstrain.org/help/general/how-to-read-a-tree"><img alt="illustration d'un arbre phylogénétique" width="100" src="http://data.nextstrain.org/img_toy_alignment_mini.png"/> Comment lire les phylogénies </a>

</div>

## Lectures additionnelles:

* Informations générales sur les coronavirus sur [Wikipedia](https://fr.wikipedia.org/wiki/Coronavirus)
* Documents fournis par l'[US CDC](https://www.cdc.gov/coronavirus/index.html)

## Récits Nextstrain

Les pages suivantes contiennent des analyses effectuées à l'aide de [Nextstrain](https://nextstrain.org).
En faisant défiler la barre latérale gauche, vous découvrirez des paragraphes de texte avec sur le côté droit une visualisation correspondant aux données génomiques.

Obtenir des génomes complets d'un nouveau et grand virus à ARN si rapidement est une réalisation remarquable. Ces analyses ont été rendues possibles grâce au partage rapide et ouvert de données génomiques et aux interprétations de scientifiques du monde entier (voir la diapositive finale pour une visualisation des crédits des séquençages).
```

<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [Une note sur les questions d'échantillonnage](https://nextstrain.org/ncov/2020-03-13?c=country&r=country&d=map&p=grid&legend=closed)
Nous avons actuellement des séquences d'échantillons prélevés dans 30 pays sur 5 continents. C'est un exploit incroyable: séquencer un virus à ARN inconnu de grande taille au milieu d'une pandémie est difficile et n'est possible que grâce au travail incroyable et au partage en temps opportun de données par des scientifiques et des médecins du monde entier.
<br><br>
Bien que ces données nous permettent d'inférer de nombreuses caractéristiques utiles de l'épidémie et de suivre sa propagation en temps réel, il est important de souligner que nos conclusions sont limitées par les données disponibles.
<br><br>
Par exemple, la carte montre très peu de séquences des pays du Sud. Ce n'est PAS parce que COVID-19 ne circule pas dans ces zones, ou que ces cas ne sont pas aussi cruciaux à comprendre; c'est simplement que nous n'avons pas beaucoup de données disponibles dans ces régions. La taille de chaque cercle sur la carte indique la quantité de données actuellement disponibles dans cette région, et non pas la taille réelle de l'épidémie.

<!-- There is NO right-side text -->


<!-- ############ SLIDE BREAK ############# -->
<!-- This is left-side text -->
# [Circulation en Europe](https://nextstrain.org/ncov/2020-03-13?c=country&legend=closed&f_country=Belgium,France,Germany,Ireland,Italy,Netherlands,Portugal,Spain,Sweden,Switzerland,United%20Kingdom,Ireland&label=clade:A2&m=div&d=map,tree&p=grid)
Ici, nous observons un grand clade de séquences d'Europe. Il est important de noter que des séquences de nombreux pays différents s'intercalent, ce qui indique COVID-19 circule déjà assez largement à travers l'Europe.
<br><br>
En zoomant sur la carte, on observe qu'il existe de nombreux liens entre l'Italie et d'autres zones ; il est toutefois important de garder à l'esprit que la direction de ces liens ne peut pas toujours être déduite avec certitude. D'autres hypothèses peuvent également expliquer ces données (par exemple, si un cas non échantillonné a infecté à la fois un cas secondaire séquencé en Italie et un cas secondaire séquencé ailleurs).

<!-- There is no right side text -->

# [Transmission locale dans les îles britanniques et en Irlande](https://nextstrain.org/ncov/2020-03-13?c=country&legend=closed&d=tree&f_country=United%20Kingdom,Ireland&p=full)
En prenant l'exemple des îles britanniques et de l'Irlande, on peut observer plusieurs cas où des virus étroitement liés à des échantillons provenant d'autres pays apparaissent dans les îles britanniques et en Irlande.
<br><br>
Cela correspond à 4 introductions ou plus en provenance d'autres endroits.
<br><br>
Nous voyons également des situations où, après une introduction, il y a plusieurs cas étroitement liés au même endroit. Ceci est cohérent avec de la transmission intracommunautaire à partir de plusieurs de ces introductions.
<!-- There is no right side text -->

<!-- ############ SLIDE BREAK ############# -->
<!-- This is left-side text -->
# [Propagation du SRAS-CoV-2 en provenance d'Iran](https://nextstrain.org/ncov/2020-03-13?d=tree,map&label=clade:A3&p=grid&legend=closed&m=div)
Un certain nombre de génomes ont été séquencés à partir de patients ayant des antécédents de voyage en Iran. Ces génomes sont tous extrêmement similaires et indiquent que l'épidémie en Iran est peut-être le résultat d'une seule introduction qui a ensuite été transmise à de nombreux autres endroits.
<br><br>
Notez qu'il n'y a pas de génomes complets disponibles provenant de patients en Iran.
<!-- There is NO right-side text -->

<!-- ############ SLIDE BREAK ############# -->
<!-- This is left-side text -->
# [Introductions aux Etats-Unis](https://nextstrain.org/ncov/2020-03-13?d=tree,map&f_country=USA&m=div&p=full&legend=closed)
Ici, nous pouvons observer que le virus a été introduit aux États-Unis à de multiples occasions indépendantes.
<br><br>
Ici, nous pouvons observer que le virus a été introduit aux États-Unis à plusieurs reprises de manière indépendante. La plupart de ces introductions ne sont associées à aucun autre cas échantillonné aux États-Unis, nous ne savons donc pas si ces introductions ont conduit à des épidémies locales. Cependant, étant donné que la capacité de test diagnostic n'a pas encore augmenté suffisamment dans la plupart des régions, nous nous attendons à ce qu'il y ait de nombreux cas non signalés.
<br><br>
Pour l'état de Washington et la Californie, cependant, nous observons des groupes de cas qui sont étroitement liés. Cela suggère une transmission continue et une propagation locale dans ces deux États.
<!-- There is no right side text -->

<!-- ############ SLIDE BREAK ############# -->
<!-- This is left-side text -->
# [Propagation du SRAS-CoV-2 dans l'État de Washington (USA)](https://nextstrain.org/ncov/2020-03-13?c=division&r=division&d=tree,map&f_country=USA&label=clade:B1&m=div&p=grid&legend=closed)
Ici, nous observons un grand regroupement de cas dans l'État de Washington qui sont tous étroitement liés. Nous concluons donc qu'il existe une large propagation locale dans cet État.
<br><br>
Il est intéressant de noter que les échantillons de l'Etat de Washington sont intercalés avec des échantillons du navire de croisière Grand Princess. Nous ne savons pas encore si le virus s'est propagé du navire de croisière à Washington ou l'inverse; à mesure que nous obtiendrons plus de données, nous mettrons à jour notre analyse.
<!-- There is NO right-side text -->

<!-- This is left-side text -->
# [Propagation du SRAS-CoV-2 en Californie](https://nextstrain.org/ncov/2020-03-13?c=country&r=division&d=tree,map&f_division=California&m=div&p=grid&legend=closed)
En regardant les échantillons de Californie, nous observons la preuve d'introductions multiples. Plus importants encore, nous observons au moins un groupe de cas étroitement liés, tous échantillonnés en Californie sur une courte période (cliquez sur ['Explore the Data'](https://nextstrain.org/ncov) et recherchez 'CA9' pour voir l'exemple).
<br><br>
Cela suggère fortement qu'il y a une transmission locale en cours en Californie.
<!-- There is NO right-side text -->

<!-- ############ SLIDE BREAK ############# -->
<!-- This is left-side text -->
# [A retenir](https://nextstrain.org/ncov/2020-03-13?c=country&d=map&p=full)
- Le virus a été introduit à plusieurs reprises dans de nombreuses régions du globe. Toutes les introductions n'entraînent pas une transmission locale.
<br><br>
- Nous observons des preuves de transmission locale à travers l'Europe, certaines parties des États-Unis, la Chine et l'Asie du Sud-Est.
<br><br>
- Le contrôle des épidémies locales par la distanciation sociale est crucial pour protéger les personnes vulnérables.

<!-- This is the right-side text -->

```auspiceMainDisplayMarkdown
# Ce que vous pouvez faire

La distanciation sociale - c'est-à-dire la diminution du nombre de personnes que vous rencontrez chaque jour - cela peut être difficile, mais cette mesure est extrêmement bénéfique pour le bien public.  

Si chacun diminuait ses contacts quotidiens de 25 %, on pourrait s'attendre à une diminution de 50 % du nombre cumulé de cas au cours du mois prochain ([Klein et al., 2020-03-13](https://institutefordiseasemodeling.github.io/COVID-public/reports/Working%20paper%20%E2%80%93%20model-based%20estimates%20of%20COVID-19%20burden%20in%20King%20and%20Snohomish%20counties%20through%20April%207.pdf)). Vous n'êtes pas sûr de ce que signifie la distanciation sociale ? [Consultez ce guide pratique en anglais](https://www.theatlantic.com/family/archive/2020/03/coronavirus-what-does-social-distancing-mean/607927/).
<div>
  <img src="https://github.com/nextstrain/ncov/raw/master/figures/social-distancing-efficacy.png" width="70%">
</div>

## Mesures que vous pouvez prendre en tant que particulier
* Réduisez le nombre de personnes avec lesquelles vous êtes en contact chaque jour, surtout si vous faites partie d'un groupe vulnérable (par exemple, les personnes âgées et les personnes souffrant de maladies préexistantes).
* N'oubliez pas que même si vous ne faites pas partie de ces personnes vulnérables, de nombreuses personnes autour de vous le sont ; suivez ces pratiques pour protéger les autres.
* Lavez-vous les mains "comme si vous veniez de couper un piment et que vous deviez changer une lentille de contact".
* Restez à la maison si vous êtes malade ; soyez prêt à vous procurer quelques provisions supplémentaires au cas où vous auriez besoin de vous mettre en quarantaine.
* Si vous êtes un employeur, encouragez vos employés à rester chez eux lorsqu'ils sont malades (et soutenez-les financièrement).

## Mesures que vous pouvez prendre en tant que vous pouvez prendre en tant qu'autorité responsable
* Rendre les tests gratuits et largement disponibles.  
* Mettre en place des mesures de distanciation sociale.
* Soutenir financièrement les personnes touchées par les mesures de distanciation sociale (par exemple, les travailleurs payés a l'heure, les personnes ayant des responsabilités de soins aux personnes âgées ou aux enfants, les petites entreprises, etc.)
```

<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [FAQs & Idées reçues](https://nextstrain.org/ncov/2020-03-05)

#### Nous savons que beaucoup de gens ont des questions sur COVID-19.

#### [Nous avons mis en place un guide pour tenter de répondre aux questions les plus fréquemment posées](https://nextstrain.org/help/coronavirus/FAQ).

#### La Fédération des scientifiques américains (Federation of American Scientists) maintient également [une excellente FAQ](https://covid19.fas.org/l/en).

<!-- This is right-side text -->
```auspiceMainDisplayMarkdown
# Lectures additionnelles (en anglais):

* "Don't believe the conspiracy theories you hear about coronavirus & HIV" [article](https://massivesci.com/notes/wuhan-coronavirus-ncov-sars-mers-hiv-human-immunodeficiency-virus/) _2020-01-31_

* "Baseless Conspiracy Theories Claim New Coronavirus Was Bioengineered" [article](https://www.factcheck.org/2020/02/baseless-conspiracy-theories-claim-new-coronavirus-was-bioengineered/) _2020-02-07_

* "No, The Wuhan Coronavirus Was Not Genetically Engineered To Put Pieces Of HIV In It" [article](https://www.forbes.com/sites/victoriaforster/2020/02/02/no-coronavirus-was-not-bioengineered-to-put-pieces-of-hiv-in-it/#5d339e8e56cb) _2020-02-02_

* "Busting coronavirus myths" [AFP Fact Check](https://factcheck.afp.com/busting-coronavirus-myths) _2020-02-19_


# Idées reçues

Un certain nombre d'idées fausses ont été diffusées sur les origines du nouveau coronavirus.
Lors d'épidémies comme celle-ci, la diffusion d'informations connues pour être incorrectes peut conduire à davantage de panique et amener les gens à ne pas faire confiance aux scientifiques et aux gouvernements, ce qui signifie qu'ils sont moins susceptibles de suivre les avis et de prendre les précautions appropriées.

Afin d'essayer d'expliquer pourquoi ces vues sont incorrectes, les scientifiques ont abordé ces théories dans les pages ci-dessous (en anglais):

<div>

  <a href="http://virological.org/t/ncovs-relationship-to-bat-coronaviruses-recombination-signals-no-snakes-no-evidence-the-2019-ncov-lineage-is-recombinant/331"><img alt="photo d'un serpent" width="100" src="http://data.nextstrain.org/img_snake-freeToUse.jpg"/> 'Snake' Origins of SARS-CoV-2 (Technical) </a>
  <a href="https://twitter.com/trvrb/status/1223666856923291648"><img alt="illustration du VIH" width="100" src="http://data.nextstrain.org/img_HIV-wiki.jpg"/> 'HIV Engineering' Idea (Twitter thread)</a>


</div>


```


<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [Crédit scientifique](https://nextstrain.org/ncov/2020-03-05?d=map&c=author)

Nous tenons à souligner le superbe travail accompli si rapidement par tous les scientifiques impliqués dans cette épidémie, et en particulier ceux qui travaillent en Chine. Ce n'est que par le partage rapide des données génomiques et des métadonnées que de telles analyses sont possibles.

<br>

Nous remercions également [GISAID](https://gisaid.org) d'avoir fourni la plate-forme à travers laquelle ces données peuvent être téléchargées et partagées.

<!-- Do not need to translate institutions names -->
<!-- This is right-side text -->
```auspiceMainDisplayMarkdown

Nous sommes reconnaissants pour les données recueillies par ces laboratoires d'origine:

* Arizona Department of Health Services
* Auckland Hospital
* BCCDC Public Health Laboratory
* Bamrasnaradura Hospital
* Bundeswehr Institute of Microbiology
* CNR Virus des Infections Respiratoires - France SUD
* CR&WISCO GENERAL HOSPITAL
* California Department of Health
* California Department of Public Health
* Center of Medical Microbiology, Virology, and Hospital Hygiene
* Center of Medical Microbiology, Virology, and Hospital Hygiene, University of Duesseldorf
* Centers for Disease Control, R.O.C. (Taiwan)
* Centre for Human and Zoonotic Virology (CHAZVY), College of Medicine University of Lagos/Lagos University Teaching Hospital (LUTH), part of the Laboratory Network of the Nigeria Centre for Disease Control (NCDC)
* Centre for Infectious Diseases and Microbiology - Public Health
* Centre for Infectious Diseases and Microbiology Laboratory Services
* Centro Hospital do Porto, E.P.E. - H. Geral de Santo Antonio
* Centro Hospitalar e Universitario de Sao Joao, Porto
* Charite Universitatsmedizin Berlin, Institute of Virology; Institut fur Mikrobiologie der Bundeswehr, Munich
* Department of Infectious Diseases, Istituto Superiore di Sanita, Roma , Italy
* Department of Infectious and Tropical Diseases, Bichat Claude Bernard Hospital, Paris
* Department of Internal Medicine, Triemli Hospital
* Department of Laboratory Medicine, National Taiwan University Hospital
* Department of Microbiology, Institute for Viral Diseases, College of Medicine, Korea University
* Department of Pathology, Toshima Hospital
* Department of Virology III, National Institute of Infectious Diseases
* Department of Virology and Immunology, University of Helsinki and Helsinki University Hospital, Huslab Finland
* Department of microbiology laboratory,Anhui Provincial Center for Disease Control and Prevention
* Dept. of Pathology, National Institute of Infectious Diseases
* Dept. of Virology III, National Institute of Infectious Diseases
* Dienst Gezondheid & Jeugd Zuid-Holland Zuid
* Division of Infectious Diseases, Department of Internal Medicine, Korea University College of Medicine
* Division of Infectious Diseases, University Hospital Zurich
* Division of Viral Diseases, Center for Laboratory Control of Infectious Diseases, Korea Centers for Diseases Control and Prevention
* Dutch COVID-19 response team
* ErasmusMC
* Foundation Elisabeth-Tweesteden Ziekenhuis
* Foundation Pamm
* Fujian Center for Disease Control and Prevention
* General Hospital of Central Theater Command of People's Liberation Army of China
* Guangdong Provincial Center for Diseases Control and Prevention; Guangdong Provincial Public Health
* Guangdong Provincial Center for Diseases Control and Prevention; Guangdong Provinical Public Health
* Guangdong Provincial Center for Diseases Control and Prevention;Guangdong Provincial Institute of Public Health
* Guangdong Provincial Institution of Public Health, Guangdong Provinical Center for Disease Control and Prevention
* HUS Diagnostiikkakeskus, Hallinto
* Hangzhou Center for Disease Control and Prevention
* Hangzhou Center for Disease and Control Microbiology Lab
* Harborview Medical Center
* Hong Kong Department of Health
* Hospital Israelita Albert Einstein
* IL Department of Public Health Chicago Laboratory
* INMI Lazzaro Spallanzani IRCCS
* Indian Council of Medical Research - National Institute of Virology
* Indian Council of Medical Research-National Institute of Virology
* Institute of Pathogen Biology, Chinese Academy of Medical Sciences & Peking Union Medical College
* Institute of Viral Disease Control and Prevention, China CDC
* Instituto Nacional de Enfermedades Respiratorias
* KU Leuven, Clinical and Epidemiological Virology
* Klinik Hirslanden Zurich
* Korea Centers for Disease Control & Prevention (KCDC) Center for Laboratory Control of Infectious Diseases Division of Viral Diseases
* Laboratoire National de Sante
* Laboratoire de Virologie, HUG
* Laboratorio di Microbiologia e Virologia, Universita Vita-Salute San Raffaele, Milano
* Laboratory Medicine
* Lapland Central Hospital
* MHC Brabant Zuidoost
* MHC Drente
* MHC Flevoland
* MHC Gooi & Vechtstreek
* MHC Haaglanden
* MHC Kennemerland
* MHC Rotterdam-Rijnmond
* MHC Utrecht
* MHC West-Brabant
* MSHS Clinical Microbiology Laboratories
* Massachusetts Department of Public Health
* Monash Medical Centre
* NHC Key laboratory of Enteric Pathogenic Microbiology, Institute of Pathogenic Microbiology
* National Centre for Infectious Diseases
* National Influenza Center - National Institute of Hygiene and Epidemiology (NIHE)
* National Influenza Centre, National Public Health Laboratory, Kathmandu, Nepal
* National Institute for Viral Disease Control and Prevention, China CDC
* National Public Health Laboratory
* National Public Health Laboratory, National Centre for Infectious Diseases
* Pathology Queensland
* Providence Regional Medical Center
* Public Health Ontario Laboratory
* RIVM
* Respiratory Virus Unit, Microbiology Services Colindale, Public Health England
* Seattle Flu Study
* Serology, Virology and OTDS Laboratories (SAViD), NSW Health Pathology Randwick
* Servicio Microbiologia. Hospital Clinico Universitario. Valencia.
* Shenzhen Key Laboratory of Pathogen and Immunity, National Clinical Research Center for Infectious Disease, Shenzhen Third People's Hospital
* Singapore General Hospital
* Sorbonne Universite, Inserm et Assistance Publique-Hopitaux de Paris (Pitie Salpetriere)
* State Health Office Baden-Wuerttemberg
* Taiwan Centers for Disease Control
* Texas Department of State Health Services
* The Central Hospital Of Wuhan
* The National Institute of Public Health Center for Epidemiology and Microbiology
* The University of Hong Kong - Shenzhen Hospital
* Tianmen Center for Disease Control and Prevention
* UCD National Virus Reference Laboratory
* University of Washington Virology Lab
* Union Hospital of Tongji Medical College, Huazhong University of Science and Technology
* Valley Medical Center
* Virology Department, Sheffield Teaching Hospitals NHS Foundation Trust
* Virology Unit, Institut Pasteur du Cambodge.
* Wales Specialist Virology Centre
* Washington State Department of Health
* Washington State Public Health Lab
* Weifang Center for Disease Control and Prevention
* West of Scotland Specialist Virology Centre, NHSGGC
* Wisconsin Department of Health Services
* Wuhan Fourth Hospital
* Wuhan Jinyintan Hospital
* Wuhan Lung Hospital
* Yongchuan District Center for Disease Control and Prevention
* Zhejiang Provincial Center for Disease Control and Prevention
* Zhongxian Center for Disease Control and Prevention

```


<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [Crédit scientifique détaillé](https://nextstrain.org/ncov/2020-03-05?d=map&c=author)

Ces données ont été partagées par [GISAID](https://gisaid.org). Nous remercions chaleureusement leurs contributions.

<br>

Sur la droite nous indiquons les séquences partagées par chaque laboratoire.

<!-- This is right-side text -->
```auspiceMainDisplayMarkdown

Les génomes du SRAS-CoV-2 ont été généreusement partagés par les scientifiques de ces laboratoires:
* Arizona Department of Health Services
	* USA/AZ1/2020

* Auckland Hospital
	* NewZealand/01/2020

* BCCDC Public Health Laboratory
	* Canada/BC_37_0-2/2020

* Bamrasnaradura Hospital
	* Nonthaburi/61/2020
	* Nonthaburi/74/2020

* Beijing Institute of Microbiology and Epidemiology
	* pangolin/Guangdong/P2S/2019
	* pangolin/Guangxi/P1E/2017
	* pangolin/Guangxi/P2V/2017
	* pangolin/Guangxi/P3B/2017
	* pangolin/Guangxi/P4L/2017
	* pangolin/Guangxi/P5E/2017
	* pangolin/Guangxi/P5L/2017

* Bundeswehr Institute of Microbiology
	* Germany/BavPat2/2020
	* Germany/BavPat3/2020

* CNR Virus des Infections Respiratoires - France SUD
	* France/RA739/2020

* CR&WISCO GENERAL HOSPITAL
	* Wuhan/HBCDC-HB-05/2020

* California Department of Health
	* USA/CA3/2020
	* USA/CA4/2020
	* USA/CA5/2020

* California Department of Public Health
	* USA/CA-CDPH-UC1/2020
	* USA/CA-CDPH-UC2/2020
	* USA/CA-CDPH-UC3/2020
	* USA/CA-CDPH-UC4/2020
	* USA/CA-CDPH-UC5/2020
	* USA/CA-CDPH-UC6/2020
	* USA/CA-CDPH-UC7/2020
	* USA/CA-CDPH-UC8/2020
	* USA/CA-CDPH-UC9/2020
	* USA/CA1/2020
	* USA/CA2/2020
	* USA/CA6/2020
	* USA/CA7/2020
	* USA/CA8/2020
	* USA/CA9/2020
	* USA/UC-CDPH-UC11/2020

* Center of Medical Microbiology, Virology, and Hospital Hygiene
	* Germany/NRW-01/2020
	* Germany/NRW-02-1/2020
	* Germany/NRW-03/2020
	* Germany/NRW-04/2020

* Center of Medical Microbiology, Virology, and Hospital Hygiene, University of Duesseldorf
	* Germany/NRW-011/2020
	* Germany/NRW-05/2020
	* Germany/NRW-06/2020
	* Germany/NRW-07/2020
	* Germany/NRW-08/2020
	* Germany/NRW-09/2020
	* Germany/NRW-10/2020

* Centers for Disease Control, R.O.C. (Taiwan)
	* Taiwan/2/2020

* Centre for Human and Zoonotic Virology (CHAZVY), College of Medicine University of Lagos/Lagos University Teaching Hospital (LUTH), part of the Laboratory Network of the Nigeria Centre for Disease Control (NCDC)
	* Nigeria/Lagos01/2020

* Centre for Infectious Diseases and Microbiology - Public Health
	* Australia/NSW10/2020
	* Australia/NSW12/2020
	* Australia/NSW13/2020
	* Australia/NSW14/2020

* Centre for Infectious Diseases and Microbiology Laboratory Services
	* Australia/NSW01/2020
	* Australia/NSW05/2020
	* Australia/NSW06/2020
	* Australia/NSW07/2020
	* Australia/NSW08/2020
	* Australia/NSW09/2020
	* Sydney/2/2020

* Centre for Infectious Diseases and Microbiology- Public Health
	* Australia/NSW11/2020

* Centro Hospital do Porto, E.P.E. - H. Geral de Santo Antonio
	* Portugal/CV62/2020

* Centro Hospitalar e Universitario de Sao Joao, Porto
	* Portugal/CV63/2020

* Charite Universitatsmedizin Berlin, Institute of Virology; Institut fur Mikrobiologie der Bundeswehr, Munich
	* Germany/BavPat1/2020

* Department of Infectious Diseases, Istituto Superiore di Sanita, Roma , Italy
	* Italy/CDG1/2020

* Department of Infectious Diseases, Istituto Superiore di Sanita, Rome, Italy
	* Italy/SPL1/2020

* Department of Infectious and Tropical Diseases, Bichat Claude Bernard Hospital, Paris
	* France/IDF0372-isl/2020
	* France/IDF0372/2020
	* France/IDF0373/2020
	* France/IDF0386-islP1/2020
	* France/IDF0386-islP3/2020
	* France/IDF0515-isl/2020
	* France/IDF0515/2020
	* France/IDF0571/2020

* Department of Internal Medicine, Triemli Hospital
	* Switzerland/1000477102/2020
	* Switzerland/1000477377/2020

* Department of Laboratory Medicine, National Taiwan University Hospital
	* Taiwan/NTU01/2020
	* Taiwan/NTU02/2020
	* Taiwan/NTU03/2020

* Department of Microbiology, Institute for Viral Diseases, College of Medicine, Korea University
	* SouthKorea/KUMC01/2020
	* SouthKorea/KUMC02/2020
	* SouthKorea/KUMC04/2020
	* SouthKorea/KUMC06/2020

* Department of Pathology, Toshima Hospital
	* Japan/TK/20-31-3/2020

* Department of Virology III, National Institute of Infectious Diseases
	* Japan/AI/I-004/2020

* Department of Virology and Immunology, University of Helsinki and Helsinki University Hospital, Huslab Finland
	* Finland/FIN01032020/2020
	* Finland/FIN03032020A/2020
	* Finland/FIN03032020B/2020
	* Finland/FIN03032020C/2020

* Department of microbiology laboratory,Anhui Provincial Center for Disease Control and Prevention
	* Anhui/SZ005/2020

* Dept. of Pathology, National Institute of Infectious Diseases
	* Japan/NA-20-05-1/2020
	* Japan/OS-20-07-1/2020

* Dept. of Virology III, National Institute of Infectious Diseases
	* Japan/KY-V-029/2020
	* Japan/TY-WK-012/2020
	* Japan/TY-WK-501/2020
	* Japan/TY-WK-521/2020

* Dienst Gezondheid & Jeugd Zuid-Holland Zuid
	* Netherlands/Hardinxveld_Giessendam_1364806/2020

* Division of Infectious Diseases, Department of Internal Medicine, Korea University College of Medicine
	* SouthKorea/KUMC03/2020
	* SouthKorea/KUMC05/2020

* Division of Infectious Diseases, University Hospital Zurich
	* Switzerland/1000477796/2020
	* Switzerland/1000477797/2020
	* Switzerland/1000477806/2020

* Division of Viral Diseases, Center for Laboratory Control of Infectious Diseases, Korea Centers for Diseases Control and Prevention
	* SouthKorea/KCDC05/2020
	* SouthKorea/KCDC06/2020
	* SouthKorea/KCDC07/2020
	* SouthKorea/KCDC12/2020
	* SouthKorea/KCDC24/2020

* Dutch COVID-19 response team
	* Netherlands/Gelderland_1/2020
	* Netherlands/Limburg_2/2020
	* Netherlands/Limburg_3/2020
	* Netherlands/Limburg_4/2020
	* Netherlands/Limburg_5/2020
	* Netherlands/Limburg_6/2020
	* Netherlands/NoordBrabant_1/2020
	* Netherlands/NoordBrabant_10/2020
	* Netherlands/NoordBrabant_11/2020
	* Netherlands/NoordBrabant_12/2020
	* Netherlands/NoordBrabant_13/2020
	* Netherlands/NoordBrabant_14/2020
	* Netherlands/NoordBrabant_15/2020
	* Netherlands/NoordBrabant_16/2020
	* Netherlands/NoordBrabant_17/2020
	* Netherlands/NoordBrabant_18/2020
	* Netherlands/NoordBrabant_19/2020
	* Netherlands/NoordBrabant_2/2020
	* Netherlands/NoordBrabant_20/2020
	* Netherlands/NoordBrabant_21/2020
	* Netherlands/NoordBrabant_22/2020
	* Netherlands/NoordBrabant_23/2020
	* Netherlands/NoordBrabant_24/2020
	* Netherlands/NoordBrabant_25/2020
	* Netherlands/NoordBrabant_26/2020
	* Netherlands/NoordBrabant_27/2020
	* Netherlands/NoordBrabant_28/2020
	* Netherlands/NoordBrabant_29/2020
	* Netherlands/NoordBrabant_3/2020
	* Netherlands/NoordBrabant_30/2020
	* Netherlands/NoordBrabant_31/2020
	* Netherlands/NoordBrabant_32/2020
	* Netherlands/NoordBrabant_33/2020
	* Netherlands/NoordBrabant_34/2020
	* Netherlands/NoordBrabant_35/2020
	* Netherlands/NoordBrabant_36/2020
	* Netherlands/NoordBrabant_37/2020
	* Netherlands/NoordBrabant_38/2020
	* Netherlands/NoordBrabant_39/2020
	* Netherlands/NoordBrabant_4/2020
	* Netherlands/NoordBrabant_5/2020
	* Netherlands/NoordBrabant_6/2020
	* Netherlands/NoordHolland_1/2020
	* Netherlands/NoordHolland_2/2020
	* Netherlands/Overijssel_1/2020
	* Netherlands/Overijssel_2/2020
	* Netherlands/Utrecht_1/2020
	* Netherlands/Utrecht_10/2020
	* Netherlands/Utrecht_11/2020
	* Netherlands/Utrecht_12/2020
	* Netherlands/Utrecht_13/2020
	* Netherlands/Utrecht_14/2020
	* Netherlands/Utrecht_15/2020
	* Netherlands/Utrecht_16/2020
	* Netherlands/Utrecht_2/2020
	* Netherlands/Utrecht_3/2020
	* Netherlands/Utrecht_4/2020
	* Netherlands/Utrecht_5/2020
	* Netherlands/Utrecht_6/2020
	* Netherlands/Utrecht_7/2020
	* Netherlands/Utrecht_8/2020
	* Netherlands/ZuidHolland_1/2020
	* Netherlands/ZuidHolland_10/2020
	* Netherlands/ZuidHolland_11/2020
	* Netherlands/ZuidHolland_13/2020
	* Netherlands/ZuidHolland_14/2020
	* Netherlands/ZuidHolland_15/2020
	* Netherlands/ZuidHolland_16/2020
	* Netherlands/ZuidHolland_17/2020
	* Netherlands/ZuidHolland_18/2020
	* Netherlands/ZuidHolland_19/2020
	* Netherlands/ZuidHolland_2/2020
	* Netherlands/ZuidHolland_20/2020
	* Netherlands/ZuidHolland_21/2020
	* Netherlands/ZuidHolland_22/2020
	* Netherlands/ZuidHolland_23/2020
	* Netherlands/ZuidHolland_24/2020
	* Netherlands/ZuidHolland_5/2020
	* Netherlands/ZuidHolland_6/2020
	* Netherlands/ZuidHolland_7/2020
	* Netherlands/ZuidHolland_8/2020
	* Netherlands/ZuidHolland_9/2020

* ErasmusMC
	* Netherlands/Nieuwendijk_1363582/2020
	* Netherlands/Rotterdam_1363790/2020

* Foundation Elisabeth-Tweesteden Ziekenhuis
	* Netherlands/Tilburg_1363354/2020
	* Netherlands/Tilburg_1364286/2020

* Foundation Pamm
	* Netherlands/Berlicum_1363564/2020

* Fujian Center for Disease Control and Prevention
	* Fujian/13/2020
	* Fujian/8/2020

* General Hospital of Central Theater Command of People's Liberation Army of China
	* Wuhan/WH01/2019
	* Wuhan/WH02/2019
	* Wuhan/WH03/2020
	* Wuhan/WH04/2020

* Guangdong Provincial Center for Diseases Control and Prevention; Guangdong Provincial Public Health
	* Foshan/20SF207/2020
	* Foshan/20SF210/2020
	* Foshan/20SF211/2020
	* Guangdong/20SF012/2020
	* Guangdong/20SF013/2020
	* Guangdong/20SF014/2020
	* Guangdong/20SF025/2020
	* Guangdong/20SF028/2020
	* Guangdong/20SF040/2020

* Guangdong Provincial Center for Diseases Control and Prevention; Guangdong Provinical Public Health
	* Guangdong/20SF174/2020
	* Guangzhou/20SF206/2020

* Guangdong Provincial Center for Diseases Control and Prevention;Guangdong Provincial Institute of Public Health
	* Guangdong/20SF201/2020

* Guangdong Provincial Institution of Public Health, Guangdong Provinical Center for Disease Control and Prevention
	* Guangdong/2020XN4239-P0034/2020
	* Guangdong/2020XN4243-P0035/2020
	* Guangdong/2020XN4273-P0036/2020
	* Guangdong/2020XN4276-P0037/2020
	* Guangdong/2020XN4291-P0038/2020
	* Guangdong/2020XN4373-P0039/2020
	* Guangdong/2020XN4433-P0040/2020
	* Guangdong/2020XN4448-P0002/2020
	* Guangdong/2020XN4459-P0041/2020
	* Guangdong/2020XN4475-P0042/2020
	* Guangdong/DG-S2-P0054/2020
	* Guangdong/DG-S41-P0056/2020
	* Guangdong/DG-S6-P0055/2020
	* Guangdong/DG-S9-P0045/2020
	* Guangdong/FS-S29-P0051/2020
	* Guangdong/FS-S30-P0052/2020
	* Guangdong/FS-S34-P0015/2020
	* Guangdong/FS-S42-P0046/2020
	* Guangdong/FS-S48-P0047/2020
	* Guangdong/FS-S50-P0053/2020
	* Guangdong/GD2020012-P0022/2020
	* Guangdong/GD2020016-P0011/2020
	* Guangdong/GD2020080-P0010/2020
	* Guangdong/GD2020085-P0043/2020
	* Guangdong/GD2020086-P0021/2020
	* Guangdong/GD2020087-P0008/2020
	* Guangdong/GD2020115-P0009/2020
	* Guangdong/GD2020134-P0031/2020
	* Guangdong/GD2020139-P0007/2020
	* Guangdong/GD2020227-P0029/2020
	* Guangdong/GD2020233-P0027/2020
	* Guangdong/GD2020234-P0023/2020
	* Guangdong/GD2020241-P0013/2020
	* Guangdong/GD2020246-P0028/2020
	* Guangdong/GD2020258-P0018/2020
	* Guangdong/GDFS2020052-P0025/2020
	* Guangdong/GDFS2020054-P0005/2020
	* Guangdong/GDFS2020056-P0044/2020
	* Guangdong/GDFS2020127-P0026/2020
	* Guangdong/GDSZ202004-P0004/2020
	* Guangdong/GDSZ202008-P0020/2020
	* Guangdong/GDSZ202009-P0032/2020
	* Guangdong/GDSZ202013-P0014/2020
	* Guangdong/GDSZ202015-P0019/2020
	* Guangdong/GZ-S6-P0050/2020
	* Guangdong/JM-S1-P0062/2020
	* Guangdong/MM-S1-P0048/2020
	* Guangdong/SZ-N128-P0057/2020
	* Guangdong/SZ-N59-P0049/2020
	* Guangdong/ZH-N22-P0059/2020
	* Guangdong/ZH-S33-P0058/2020
	* Guangdong/ZQ-S2-P0061/2020
	* Guangdong/ZS-S6-P0060/2020

* HUS Diagnostiikkakeskus, Hallinto
	* Finland/FIN-25/2020

* Hangzhou Center for Disease Control and Prevention
	* Hangzhou/HZCDC0001/2020

* Hangzhou Center for Disease and Control Microbiology Lab
	* Hangzhou/HZ-1/2020

* Harborview Medical Center
	* USA/WA3-UW1/2020
	* USA/WA9-UW6/2020

* Hong Kong Department of Health
	* HongKong/VB20024950/2020
	* HongKong/VB20026565/2020
	* HongKong/VM20001061/2020
	* HongKong/case42_VM20002493/2020
	* HongKong/case48_VM20002507/2020
	* HongKong/case52_VM20002582/2020
	* HongKong/case78_VM20002849/2020
	* HongKong/case85_VM20002868/2020
	* HongKong/case90_VM20002907/2020
	* canine/HongKong/20-02756/2020

* Hospital Israelita Albert Einstein
	* Brazil/SPBR-01/2020
	* Brazil/SPBR-02/2020
	* Brazil/SPBR-03/2020

* Hospital Sao Joaquim Beneficencia Portuguesa
	* Brazil/SPBR-04/2020
	* Brazil/SPBR-05/2020
	* Brazil/SPBR-06/2020

* IL Department of Public Health Chicago Laboratory
	* USA/IL1/2020
	* USA/IL2/2020

* INMI Lazzaro Spallanzani IRCCS
	* Italy/INMI1-cs/2020
	* Italy/INMI1-isl/2020

* Indian Council of Medical Research - National Institute of Virology
	* India/1-27/2020

* Indian Council of Medical Research-National Institute of Virology
	* India/1-31/2020

* Institute of Pathogen Biology, Chinese Academy of Medical Sciences & Peking Union Medical College
	* Wuhan/IPBCAMS-WH-01/2019
	* Wuhan/IPBCAMS-WH-02/2019
	* Wuhan/IPBCAMS-WH-03/2019
	* Wuhan/IPBCAMS-WH-04/2019
	* Wuhan/IPBCAMS-WH-05/2020

* Institute of Viral Disease Control and Prevention, China CDC
	* Wuhan/IVDC-HB-envF13-20/2020
	* Wuhan/IVDC-HB-envF13-21/2020
	* Wuhan/IVDC-HB-envF13/2020
	* Wuhan/IVDC-HB-envF54/2020

* Instituto Nacional de Enfermedades Respiratorias
	* Mexico/CDMX/InDRE_01/2020

* Jingzhou Center for Disease Control and Prevention
	* Jingzhou/HBCDC-HB-01/2020

* KU Leuven, Clinical and Epidemiological Virology
	* Belgium/GHB-03021/2020

* Klinik Hirslanden Zurich
	* Switzerland/1000477757/2020

* Korea Centers for Disease Control & Prevention (KCDC) Center for Laboratory Control of Infectious Diseases Division of Viral Diseases
	* SouthKorea/KCDC03/2020

* Laboratoire National de Sante
	* Luxembourg/Lux1/2020

* Laboratoire de Virologie, HUG
	* Switzerland/AG0361/2020
	* Switzerland/BL0902/2020
	* Switzerland/GE3121/2020
	* Switzerland/GE3895/2020
	* Switzerland/GE5373/2020
	* Switzerland/GE9586/2020
	* Switzerland/TI9486/2020
	* Switzerland/VD5615/2020

* Laboratorio di Microbiologia e Virologia, Universita Vita-Salute San Raffaele, Milano
	* Italy/UniSR1/2020

* Laboratory Medicine
	* Taiwan/CGMH-CGU-01/2020

* Lapland Central Hospital
	* Finland/1/2020

* MHC Brabant Zuidoost
	* Netherlands/Eindhoven_1363782/2020

* MHC Drente
	* Netherlands/Dalen_1363624/2020

* MHC Flevoland
	* Netherlands/Zeewolde_1365080/2020

* MHC Gooi & Vechtstreek
	* Netherlands/Blaricum_1364780/2020
	* Netherlands/Naarden_1364774/2020

* MHC Haaglanden
	* Netherlands/Nootdorp_1364222/2020

* MHC Hart voor Brabant
	* Netherlands/Oisterwijk_1364072/2020

* MHC Kennemerland
	* Netherlands/Haarlem_1363688/2020

* MHC Rotterdam-Rijnmond
	* Netherlands/Rotterdam_1364040/2020

* MHC Utrecht
	* Netherlands/Utrecht_1363564/2020
	* Netherlands/Utrecht_1363628/2020
	* Netherlands/Utrecht_1364066/2020

* MHC West-Brabant
	* Netherlands/Andel_1365066/2020
	* Netherlands/Helmond_1363548/2020

* MSHS Clinical Microbiology Laboratories
	* USA/NY1-PV08001/2020

* Massachusetts Department of Public Health
	* USA/MA1/2020

* Monash Medical Centre
	* Australia/VIC01/2020

* NHC Key laboratory of Enteric Pathogenic Microbiology, Institute of Pathogenic Microbiology
	* Jiangsu/JS01/2020
	* Jiangsu/JS02/2020
	* Jiangsu/JS03/2020

* National Centre for Infectious Diseases
	* Singapore/12/2020
	* Singapore/13/2020
	* Singapore/14/2020
	* Singapore/3/2020
	* Singapore/4/2020

* National Influenza Center - National Institute of Hygiene and Epidemiology (NIHE)
	* Vietnam/VR03-38142/2020

* National Influenza Centre, National Public Health Laboratory, Kathmandu, Nepal
	* Nepal/61/2020

* National Institute for Viral Disease Control and Prevention, China CDC
	* Beijing/IVDC-BJ-005/2020
	* Chongqing/IVDC-CQ-001/2020
	* Henan/IVDC-HeN-002/2020
	* Jiangsu/IVDC-JS-001/2020
	* Jiangxi/IVDC-JX-002/2020
	* Shandong/IVDC-SD-001/2020
	* Shanghai/IVDC-SH-001/2020
	* Sichuan/IVDC-SC-001/2020
	* Wuhan/IVDC-HB-01/2019
	* Wuhan/IVDC-HB-04/2020
	* Wuhan/IVDC-HB-05/2019
	* Yunnan/IVDC-YN-003/2020

* National Public Health Laboratory
	* Singapore/11/2020

* National Public Health Laboratory, National Centre for Infectious Diseases
	* Singapore/10/2020
	* Singapore/7/2020
	* Singapore/8/2020
	* Singapore/9/2020

* Pathology Queensland
	* Australia/QLD01/2020
	* Australia/QLD02/2020
	* Australia/QLD03/2020
	* Australia/QLD04/2020
	* Australia/QLD09/2020

* Providence Regional Medical Center
	* USA/WA1/2020

* Public Health Ontario Laboratory
	* Canada/ON-PHL2445/2020
	* Canada/ON-VIDO-01/2020

* RIVM
	* Netherlands/Delft_1363424/2020
	* Netherlands/Diemen_1363454/2020
	* Netherlands/Loon_op_zand_1363512/2020
	* Netherlands/Oss_1363500/2020
	* NetherlandsL/Houten_1363498/2020

* Respiratory Virus Unit, Microbiology Services Colindale, Public Health England
	* England/01/2020
	* England/02/2020
	* England/09c/2020
	* England/200641094/2020
	* England/200690245/2020
	* England/200690300/2020
	* England/200690306/2020
	* England/200690756/2020
	* England/200940527/2020
	* England/200960041/2020
	* England/200960515/2020
	* England/200981386/2020
	* England/200990002/2020
	* England/200990006/2020
	* England/200990660/2020
	* England/200990723/2020
	* England/200990724/2020
	* England/200990725/2020
	* England/200991076/2020
	* England/201000003/2020
	* England/201040081/2020
	* England/201040141/2020

* Seattle Flu Study
	* USA/WA-S2/2020
	* USA/WA-S3/2020

* Second Hospital of Anhui Medical University
	* Hefei/2/2020

* Serology, Virology and OTDS Laboratories (SAViD), NSW Health Pathology Randwick
	* Sydney/3/2020

* Servicio Microbiologia. Hospital Clinico Universitario. Valencia.
	* Spain/Valencia1/2020
	* Spain/Valencia2/2020

* Shenzhen Key Laboratory of Pathogen and Immunity, National Clinical Research Center for Infectious Disease, Shenzhen Third People's Hospital
	* Shenzhen/SZTH-002/2020
	* Shenzhen/SZTH-003/2020
	* Shenzhen/SZTH-004/2020

* Shenzhen Third People's Hospital
	* Shenzhen/SZTH-001/2020

* Singapore General Hospital
	* Singapore/1/2020
	* Singapore/2/2020

* Singapore General Hospital, Molecular Laboratory, Division of Pathology
	* Singapore/5/2020
	* Singapore/6/2020

* Sorbonne Universite, Inserm et Assistance Publique-Hopitaux de Paris (Pitie Salpetriere)
	* France/IDF0626/2020

* South China Agricultural University
	* pangolin/Guandong/1/2019

* State Health Office Baden-Wuerttemberg
	* Germany/Baden-Wuerttemberg-1/2020

* Taiwan Centers for Disease Control
	* Taiwan/3/2020
	* Taiwan/4/2020

* Texas Department of State Health Services
	* USA/TX1/2020

* The Central Hospital Of Wuhan
	* Wuhan/HBCDC-HB-02/2020

* The National Institute of Public Health Center for Epidemiology and Microbiology
	* CzechRepublic/951/2020

* The University of Hong Kong - Shenzhen Hospital
	* Shenzhen/HKU-SZ-002/2020
	* Shenzhen/HKU-SZ-005/2020

* Tianmen Center for Disease Control and Prevention
	* Tianmen/HBCDC-HB-07/2020

* UCD National Virus Reference Laboratory
	* Ireland/COR-20134/2020

* UW Virology Lab
	* USA/WA-UW15/2020
	* USA/WA-UW16/2020
	* USA/WA-UW17/2020
	* USA/WA-UW18/2020
	* USA/WA-UW19/2020
	* USA/WA-UW20/2020
	* USA/WA-UW21/2020
	* USA/WA11-UW7/2020
	* USA/WA12-UW8/2020
	* USA/WA13-UW9/2020
	* USA/WA14-UW10/2020
	* USA/WA15-UW11/2020
	* USA/WA16-UW12/2020
	* USA/WA17-UW13/2020
	* USA/WA18-UW14/2020

* Union Hospital of Tongji Medical College, Huazhong University of Science and Technology
	* Wuhan/HBCDC-HB-03/2020
	* Wuhan/HBCDC-HB-04/2020

* Unknown
	* Netherlands/Coevorden_1363618/2020

* Valley Medical Center
	* USA/WA8-UW5/2020

* Virology Department, Sheffield Teaching Hospitals NHS Foundation Trust
	* England/Sheff01/2020
	* England/Sheff02/2020

* Virology Unit, Institut Pasteur du Cambodge.
	* Cambodia/0012/2020

* WA State Department of Health
	* USA/WA1-A12/2020

* Wales Specialist Virology Centre
	* Wales/PHW03/2020
	* Wales/PHW05/2020
	* Wales/PHW1/2020
	* Wales/PHW2/2020

* Washington State Department of Health
	* USA/WA1-F6/2020
	* USA/WA2/2020

* Washington State Public Health Lab
	* USA/WA4-UW2/2020
	* USA/WA6-UW3/2020
	* USA/WA7-UW4/2020

* Weifang Center for Disease Control and Prevention
	* China/WF0001/2020
	* China/WF0002/2020
	* China/WF0003/2020
	* China/WF0004/2020
	* China/WF0006/2020
	* China/WF0009/2020
	* China/WF0012/2020
	* China/WF0014/2020
	* China/WF0015/2020
	* China/WF0016/2020
	* China/WF0017/2020
	* China/WF0018/2020
	* China/WF0019/2020
	* China/WF0020/2020
	* China/WF0021/2020
	* China/WF0023/2020
	* China/WF0024/2020
	* China/WF0026/2020
	* China/WF0028/2020
	* China/WF0029/2020

* West of Scotland Specialist Virology Centre, NHSGGC
	* Scotland/CVR01/2020
	* Scotland/CVR02/2020
	* Scotland/CVR03/2020
	* Scotland/CVR04/2020
	* Scotland/CVR05/2020

* Wisconsin Department of Health Services
	* USA/WI1/2020

* Wuhan Fourth Hospital
	* Wuhan/WH05/2020

* Wuhan Institute of Virology, Chinese Academy of Sciences
	* bat/Yunnan/RaTG13/2013

* Wuhan Jinyintan Hospital
	* Wuhan/HBCDC-HB-01/2019
	* Wuhan/HBCDC-HB-02/2019
	* Wuhan/HBCDC-HB-03/2019
	* Wuhan/HBCDC-HB-04/2019
	* Wuhan/WIV02/2019
	* Wuhan/WIV04/2019
	* Wuhan/WIV05/2019
	* Wuhan/WIV06/2019
	* Wuhan/WIV07/2019

* Wuhan Lung Hospital
	* Wuhan/HBCDC-HB-06/2020

* Yongchuan District Center for Disease Control and Prevention
	* Chongqing/YC01/2020

* Zhejiang Provincial Center for Disease Control and Prevention
	* Zhejiang/WZ-01/2020
	* Zhejiang/WZ-02/2020

* Zhongxian Center for Disease Control and Prevention
	* Chongqing/ZX01/2020

```
