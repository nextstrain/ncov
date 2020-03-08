---
title: Analyse génomique de la propagation du nCoV. Rapport de situation au 2020-01-30.
authors: "Trevor Bedford, Richard Neher, James Hadfield, Emma Hodcroft, Misja Ilcisin, Nicola Müller, Etienne Simon-Lorière, Pierre Barrat-Charlaix"
authorLinks: "https://nextstrain.org"
affiliations: "Fred Hutch, Seattle, États-Unis; Biozentrum, Bâle, Suisse; et Institut Pasteur, Paris, France"
date: "2020 Jan 30"
dataset: "https://nextstrain.org/ncov/2020-01-30?d=map"
abstract: "Ce rapport utilise des données génomiques du nouveau coronavirus (nCoV) rendues publiques sur GISAID et Genbank afin d'estimer la vitesse et les modalités de propagation de l'épidémie virale. Nous prévoyons d'émettre des rapports de situation à mesure que de nouvelles données sont produites et partagées. Ce site web est optimisé pour être affiché sur un ordinateur de bureau."
---

# [Résumé exécutif](https://nextstrain.org/ncov/2020-01-30)

```auspiceMainDisplayMarkdown
## Résumé exécutif

Nous avons examiné la diversité génétique de 42</tag> génomes du nouveau coronavirus (nCoV) rendus publics pour estimer la date de l'ancêtre commun et le taux de propagation.
Principaux résultats :
* Les 42</tag> génomes sont très similaires, différant du consensus de 0 à 7 mutations.
* L'explication parcimonieuse de cette faible diversité est que l'épidémie est issue d'une unique introduction du virus dans la population humaine, ou d'un faible nombre de transmissions animal-homme de virus très similaires. 
* Cette transmission a probablement eu lieu en novembre ou début décembre 2019.  
* La transmission d'homme à homme a eu lieu en continu depuis cette date, aboutissant aux cas observés. 
* À l'aide des estimations du nombre total de cas de l'Imperial College de Londres, nous avons inféré un taux de reproduction entre 1,8 et 3,5, indiquant une croissance rapide sur la période de novembre à janvier. 
```

# [Coronavirus](https://nextstrain.org/ncov/2020-01-30)

### Lectures additionnelles:

* Informations générales sur le coronavirus sur [Wikipedia](https://fr.wikipedia.org/wiki/Coronavirus) _2020-01-30_
* Résumé de l'épidémie du nCoV sur [Wikipedia](https://fr.wikipedia.org/wiki/%C3%89pid%C3%A9mie_de_coronavirus_de_2019-2020) _2020-01-30_
* Documents fournis par l'[US CDC](https://www.cdc.gov/coronavirus/index.html) _2020-01-29_
* Caractéristiques du génome sur [ViralZone](https://viralzone.expasy.org/764?outline=all_by_species) _2020-01-23_
* Analyse interactive du risque par [MOBS-lab](https://datastudio.google.com/reporting/3ffd36c3-0272-4510-a140-39e288a9f15c/page/U5lCB) _2010-01-29_
* Analyse interactive du risque par [ROCS-lab](http://rocs.hu-berlin.de/corona/) _2010-01-29_

```auspiceMainDisplayMarkdown

## Différents coronavirus humains

Les coronavirus (CoV) regroupent diverses espèces de virus à ARN de sens positif (ARNsb (+)) causant des infections respiratoires chez l'Homme.
Certains variants de coronavirus sont associées à des épidémies alors que d'autres circulent en permanence et causent principalement des infections respiratoires légères (par exemple le rhume).

#### SARS-CoV et MERS-CoV
Le plus connu de ces coronavirus est le [SRAS-CoV](https://fr.wikipedia.org/wiki/Syndrome_respiratoire_aigu_s%C3%A9v%C3%A8re_li%C3%A9_au_coronavirus) ("syndrome respiratoire aigu sévère"), responsable d'une épidémie mondiale entre novembre 2002 et juillet 2003 ayant causé [plus 8000 cas et 774 décès](https://www.theguardian.com/world/2017/dec/10/sars-virus-bats-china-severe-acute-respiratory-syndrome), soit un taux de mortalité de 9–11%.

En 2012, un nouveau coronavirus causant de graves symptômes respiratoires, [MERS-CoV](https://fr.wikipedia.org/wiki/Coronavirus_du_syndrome_respiratoire_du_Moyen-Orient) ("Middle East respiratory syndrome"), a été identifié. Le MERS a causé un nombre de décès comparable au SRAS, bien que son mode de transmission soit très différent. Alors que le SRAS se transmettait efficacement d'homme à homme, les infections par le MERS étaient généralement l'effet de zoonoses indépendantes (transmissions de l'animal à l'Homme) depuis le chameau (voir [Dudas _et al._](https://elifesciences.org/articles/31257) pour plus d'informations). L'épidémie était donc intrinsèquement limitée et essentiellement restreinte à la péninsule arabique. 


#### CoV saisonniers

Cependant, tous les coronavirus ne sont pas aussi mortels que le SRAS-CoV et le MERS-CoV.
Quatre coronavirus "saisonniers" infectent habituellement l'homme chaque année.
Par rapport au SRAS, ces coronavirus saisonniers sont ["bien plus prévalents, beaucoup moins graves, et sont une cause fréquente de syndromes grippaux"](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5820427/).
En effet, [5](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2879166/)–[12](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5820427/)% de tous les syndromes grippaux sont testés positivement pour un coronavirus. Ils sont donc plutôt communs et causent des millions d'infections légères chaque année.
Les coronavirus saisonniers sont le résultat de plusieurs zoonoses depuis le réservoir animal des chauves-souris vers les humains dans la dernière centaine d'années. Après ces contagions inter-espèces, chaque virus saisonnier s'est répandu et établi dans la population humaine.


#### Réservoirs animaux
Les coronavirus infectent une grande diversité d'animaux, et les épidémies humaines décrites plus haut sont le résultat d'un ou plusieurs "sauts" de ces réservoirs animaux vers la population humaine.
On pense que le SRAS est arrivé dans la population humaine depuis [la chauves-souris rhinolophe par l'intermédiaire de la civette masquée](https://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1006698).


#### Transmission de personne à personne
La capacité de différents lignages à se transmettre d'homme à homme est extrêmement importante pour comprendre le développement potentiel d'une épidémie. 
À cause de sa capacité à se propager chez l'humain et de son fort taux de mortalité, le SRAS (ou un virus semblable au SRAS) est considéré comme une [menace de santé publique globale](https://www.who.int/whr/2007/overview/fr/index1.html) par l'OMS.

```

# [Nouveau coronavirus (nCoV) 2019-2020](https://nextstrain.org/ncov/2020-01-30)

### Lectures additionnelles (en anglais):

* New China virus: Five questions scientists are asking  [Nature news](https://www.nature.com/articles/d41586-020-00166-6) _2020-01-22_
* China virus latest: first US case confirmed  [Nature news](https://www.nature.com/articles/d41586-020-00154-w) _2020-01-21_
* New virus surging in Asia rattles scientists  [Nature news](https://www.nature.com/articles/d41586-020-00129-x) _2020-01-20_
* New virus identified as likely cause of mystery illness in China [Nature news](https://www.nature.com/articles/d41586-020-00020-9) _2020-01-08_

```auspiceMainDisplayMarkdown

## Épidémie récente d'un nouveau coronavirus
En décembre 2019, une nouvelle maladie a été détectée pour la première fois à Wuhan, en Chine.
Nous savons maintenant qu'il s'agit d'une nouvelle épidémie de coronavirus chez l'homme (la 7ème), et ce virus est provisoirement appelé nCoV (nouveau coronavirus).

Au 30 janvier, plus de 7,914 cas et 170 decès [ont été signalés](https://en.wikipedia.org/wiki/2019%E2%80%9320_outbreak_of_novel_coronavirus_(2019-nCoV)).
Il est encore trop tôt pour connaître le taux de mortalié, mais les premiers chiffres indiquent qu'il est nettement inférieur à celui du SRAS-CoV.
Le nombre de cas augmente considérablement en partie en raison de la surveillance et des tests accrus.

Alors que l'épidémie semble être concentrée à Wuhan, qui est maintenant [en quarantaine](https://twitter.com/PDChina/status/1220060879112282117), le virus s'est propagé dans toute la Chine et à l'étranger, notamment à Hong Kong, Macao, Thaïlande, Japon, Corée du Sud, États-Unis et France. Une transmission locale limitée en dehors de la Chine a été signalée.

L'origine du virus n'est pas encore claire, cependant [l'analyse génomique](https://virological.org/t/ncovs-relationship-to-bat-coronaviruses-recombination-signals-no-snakes/331) suggère que le nCoV est proche de virus précédemment identifiés chez des chauves-souris.
Il est plausible qu'il y ait eu d'autres transmissions animales intermédiaires avant l'introduction chez l'homme.
Il n'y a aucune preuve que l'intermédiaire soit des serpents.

#### récits Nextstrain

Les pages suivantes contiennent des analyses effectuées à l'aide de [Nextstrain](https://nextstrain.org).
En faisant défiler la barre latérale gauche, vous découvrirez des paragraphes de texte avec sur le côté droit une visualisation correspondant aux données génomiques.

Avoir des génomes complets d'un nouveau et  grand virus à ARN si rapidement est une réalisation remarquable.
Ces analyses ont été rendues possibles grâce au partage rapide et ouvert de données génomiques et aux interprétations de scientifiques du monde entier (voir la diapositive finale pour une visualisation des crédits des séquençages).
```

# [How to interpret the phylogenetic trees](https://nextstrain.org/ncov/2020-01-30)

### Further Reading:

* [Exploring interactive phylogenies with Auspice](https://neherlab.org/201901_krisp_auspice.html) _2019-01-24_

```auspiceMainDisplayMarkdown
## Arbres de transmission vs arbres phylogénétiques

Les agents pathogènes se propagent par une réplication rapide dans un hôte, suivie d'une transmission à un autre hôte.
Une épidémie ne peut décoller que lorsqu'une infection entraîne ensuite plusieurs nouvelles infections.

À mesure que l'agent pathogène se réplique et se propage, son génome doit se répliquer de nombreuses fois et des mutations aléatoires (erreurs de copie) vont s'accumuler dans le génome.
Ces mutations aléatoires peuvent aider à suivre la propagation de l'agent pathogène et à en apprendre davantage sur ses voies de transmission et sa dynamique.

<div>
  <img alt="schéma montrant la relation entre l'arbre de transmission et l'arbre phylogénétique" width="500" src="https://neherlab.org/talk_images/infection_tree_combined.png"/>
</div>

L'illustration ci-dessus montre un croquis d'un arbre de transmission avec un sous-ensemble de cas qui ont été échantillonnés (bleu).
En pratique, l'arbre de transmission est inconnu et généralement, seules des estimations approximatives du nombre de cas sont disponibles.
Les séquences du génome nous permettent de déduire des parties de l'arbre de transmission.
Dans cet exemple, trois mutations (petits losanges) sont indiquées sur l'arbre.
Les séquences qui ont les mêmes mutations sont plus étroitement liées, de sorte que ces mutations nous permettent de regrouper les échantillons en groupes de virus étroitement liés qui appartiennent aux mêmes chaînes de transmission.

### Lecture d'un arbre phylogénétique

Ci-dessous, nous voyons une illustration avec un arbre phylogénétique sur la gauche, où les mutations sont représentées par des cercles colorés. À droite se trouvent les séquences correspondantes, également avec des mutations représentées par des cercles colorés.
Nous pouvons voir que les séquences qui partagent les mêmes mutations se regroupent.
Lorsque des séquences apparaissent liées par une ligne verticale plate, comme A et B, cela signifie qu'il n'y a aucune différence entre elles - leurs séquences sont identiques.

Lorsqu'une séquence se trouve seule sur une longue ligne, comme C ou E, cela signifie qu'elle a des mutations uniques qui ne se trouvent pas dans d'autres séquences. Plus la ligne est longue, plus il y a de mutations.
A et B ont également des mutations uniques (le cercle vert) non partagées par les autres séquences, mais elles sont identiques l'une à l'autre.

<div>
  <img alt="schéma de l'arbre phylogénétique et de l'alignement correspondant, avec des échantillons étiquetés A-E" width="500" src="https://data.nextstrain.org/toy_alignment_tree.png"/>
</div>

Pour le moment, la nouvelle phylogénie du coronavirus (nCoV) peut ne pas ressembler beaucoup à un «arbre».
Beaucoup de séquences sont identiques - elles sont positionnées ensemble sur des lignes verticales comme A et B (certaines sont sur la partie la plus à gauche de l'arbre).
D'autres ont des mutations uniques ou partagées et se trouvent donc sur des lignes, ou des «branches», allant vers la droite.
Vous pouvez voir le nombre de mutations d'une branche en passant votre souris dessus.
```

# [Analyse phylogénétique](https://nextstrain.org/ncov/2020-01-30?m=div&d=tree)

Nous présentons ici une phylogénie de 42 </tag> souches de nCoV qui ont été partagées publiquement.
Des informations sur la façon dont l'analyse a été effectuée sont disponibles [dans ce recueil GitHub](github.com/nextstrain/ncov).

<br>

Les couleurs représentent la région d'isolement à l'intérieur d’un pays ou d’un état des États-Unis, avec l'axe des x représentant la divergence nucléotidique.

<br>

La divergence est mesurée comme le nombre de changements (mutations) dans le génome.
Plusieurs séquences n'ont aucune mutation, ce qui signifie qu'elles sont toutes identiques à la racine (centre) de l'arbre.
D'autres virus ont entre une et sept mutations.

<br>

Le séquençage du génome d'un nouveau et grand virus à ARN dans une situation d'épidémie est difficile.
Certaines des différences observées dans ces séquences peuvent être des erreurs de séquençage plutôt que des mutations réelles.
Les insertions, les suppressions et les différences aux extrémités du génome sont plus susceptibles d'être des erreurs et nous les avons donc masquées pour cette analyse.

# [Interprétation phylogénétique](https://nextstrain.org/ncov/2020-01-30?m=div&d=tree)

Nous constatons actuellement peu de diversité génétique parmi les séquences nCoV, avec 11 </tag> sur 42 </tag> séquences n'ayant pas de mutations uniques.

<br>

La faible diversité génétique à travers ces séquences suggère que l'ancêtre commun le plus récent de toutes les séquences nCoV était assez récent, car les mutations s'accumulent lentement par rapport aux autres virus à ARN, à un taux d'environ 1-2 mutations par mois pour les coronavirus.
En général, les introductions répétées à partir d'un réservoir animal montrent une diversité importante (cela a été vrai pour Lassa, Ebola, MERS-CoV et la grippe aviaire).
L'observation d'une telle concentration d'infections humaines groupées peut s'expliquer par une épidémie qui découle d'un seul événement d'introduction zoonotique dans la population humaine suivie d'une propagation d'épidémie interhumaine.

<br>

Nous commençons à voir des groupes de séquences qui partagent des mutations.
Un cluster contient des séquences de Guangdong et quatre isolats des États-Unis.
D'autres groupes contiennent de deux à quatre isolats.
Les séquences de ces groupes ont tendance à provenir d'échantillons plus récents, ce qui suggère que le virus a commencé à accumuler des mutations lors de sa propagation à Wuhan puis dans d'autres villes.
Il n'y a actuellement aucune preuve que ces mutations modifient le comportement du virus - il est attendu que les virus à ARN mutent.

# [Transmission intra-familiale 1](https://nextstrain.org/ncov/2020-01-30?m=div&d=tree&f_location=Zhuhai)

Il existe trois isolats génétiquement identiques de Zhuhai (sud-est de la Chine, province du Guangdong) qui forment un groupe partageant une mutation unique, observée dans aucun autre isolat (vous pouvez passer votre souris sur les branches pour voir quelles mutations sont présentes).

<br>

Deux de ces cas (se terminant par 028 et 040) sont [connus pour provenir d'une même famille](https://twitter.com/JingLu_LuJing/status/1220143773532880896), indiquant à nouveau une transmission interhumaine.
Nous n'avons pas d'informations sur le troisième cas.


# [Transmission intra-familiale 2](https://nextstrain.org/ncov/2020-01-30?m=div&d=tree&f_location=Shenzhen)

Des six isolats de la province du Guangdong (qui comprend la ville de Shenzhen), nous voyons quatre isolats génétiquement identiques.
Ces séquences diffèrent par 3 mutations de la racine de l'arbre.

<br>

Trois des séquences de Guangdong (se terminant par F025, F013 et F012) sont [connues pour provenir d'une même famille](https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(20)30154-9/fulltext), et représentent presque certainement une transmission interhumaine.

<br>

# [Transmission intra-familiale 2 - mutations partagées](https://nextstrain.org/ncov/2020-01-30?m=div&d=tree&f_location=Shenzhen,Los%20Angeles,Orange%20County,Seattle,Chicago,Phoenix)

Les trois mutations trouvées dans ce groupe sont également présentes dans l'isolat de l'Arizona, aux États-Unis, et deux des mutations se retrouvent dans trois autres isolats des États-Unis.

<br>

# [Transmission intra-familiale 3](https://nextstrain.org/ncov/2020-01-30?m=div&d=tree&f_location=Paris)

Enfin, les deux séquences de France sont identiques, partageant une mutation unique, et une mutation également trouvée dans l'un des isolats américains et l'isolat taiwanais.

<br>

Les deux séquences françaises sont [connues pour être de la même famille](https://www.thelocal.fr/20200129/coronavirus-in-france-what-you-need-to-know) - un couple chinois de Wuhan.

<br>

# [Cas hors de Chine](https://nextstrain.org/ncov/2020-01-30?c=country&d=tree&m=div)

Des cas de nCoV confirmés par diagnostic ont été signalés dans de nombreux pays d'Asie de l'Est et d’Asie du Sud-Est, aux États-Unis, en Australie, au Moyen-Orient et en Europe.
Le Vietnam, le Japon et l'Allemagne ont signalé une transmission dans le pays, bien que toujours avec un lien connu avec Wuhan, en Chine.

<br>

Les seules données de séquence actuellement disponibles pour les cas hors de Chine sont les deux cas de Thaïlande, cinq des États-Unis, deux de France et un de Taïwan.
Les échantillons thaïlandais sont génétiquement identiques à neuf séquences chinoises, dont sept isolées à Wuhan.
Quatre séquences des États-Unis partagent deux mutations avec le groupe de séquences de Shenzhen.
La séquence restante des États-Unis partage une mutation avec la séquence de Taïwan et les deux de France.

<br>

L'explication la plus parcimonieuse du schéma observé de mutations retrouvées parmi les séquences américaines et de Shenzhen est qu’un variant du virus avec les deux mutations partagées circulait à Wuhan et a été exporté indépendamment à Shenzhen et plusieurs fois aux États-Unis.
Il n'y a aucune preuve d'un lien entre des séquences américaines autre qu'un lien avec Wuhan.

# [Datation de l'ancêtre commun le plus récent](https://nextstrain.org/ncov/2020-01-30?d=tree)
La grande similitude des génomes suggère qu'ils partagent un ancêtre commun récent (c'est-à-dire qu'ils sont descendus du même virus ancestral récemment). Sinon, nous nous attendrions à un plus grand nombre de différences entre les échantillons.

<br>

Des recherches antérieures sur des coronavirus apparentés suggèrent que ces virus accumulent entre 1 et 3 changements dans leur génome par mois (taux de 3 &times; 10<sup>-4</sup> à 1 &times; 10<sup>-3</sup> par site et par an).

<br>

sur la droite, nous explorons comment différentes hypothèses sur le taux de changement et la diversité génétique observée nous donnent des estimations du timing de l'épidémie.

```auspiceMainDisplayMarkdown
## Date de l'ancêtre commun des virus épidémiques
Avec les séquences supplémentaires partagées au cours de la semaine dernière, l'arbre montre maintenant plusieurs groupes distincts, de sorte que notre analyse du 2020-01-25 qui supposait une topologie en étoile n'est plus appropriée.

Nous reproduisons ici notre analyse sur la base des données disponibles au 2020-01-25
en supposant une structure de phylogénie en forme d'étoile ainsi qu'une distribution de Poisson des mutations dans le temps pour estimer le timing de l'ancêtre commun le plus récent («TMRCA») des virus séquencés.
** Nous trouvons que l'ancêtre commun existait probablement entre mi-novembre et début décembre 2019. La plus grande source d'incertitude est le taux de substitution. **

<div>
  <img alt="graphique des estimations du TMRCA basées sur différents taux de mutation" width="500" src="https://data.nextstrain.org/ncov_poisson-tmrca.png"/>
</div>

En utilisant l'ensemble des données, le pipeline d'analyse de la souche suivante estime que l'ancêtre commun existait probablement entre fin novembre et début décembre 2019.

Il y a un [cas confirmé à Wuhan dont la date de survenue est le 1er décembre 2019](https://twitter.com/trvrb/status/1220749265380593664), qui mettrait une limite supérieure à la date du dernier ancêtre commun.
L'ancêtre commun des virus séquencés à ce jour pourrait cependant être postérieur à cette date.

Une modélisation plus détaillée du début de l'épidémie est en cours.
Malgré une incertitude considérable, notre meilleure estimation reste fin novembre / début décembre.

```

# [Estimation du taux de croissance](https://nextstrain.org/ncov/2020-01-30?d=tree)

Une mesure importante dans la propagation d'un pathogène est le nombre moyen de cas secondaires produits par chaque infection.

<br>

Ce nombre est appelé R0 ("R-zéro").
À droite, nous présentons des estimations simples de R0.

```auspiceMainDisplayMarkdown
## Estimates of epidemic growth rate
## Estimations du taux de croissance épidémique
Les scientifiques de l'Imperial College de Londres ont utilisé le nombre de cas observés hors de Chine pour estimer le [nombre total de cas](https://www.imperial.ac.uk/mrc-global-infectious-disease-analysis/news--wuhan-coronavirus/) et ont suggéré qu'il y avait eu au moins plusieurs milliers de cas au 2020-01-22.
Avec les nouveaux cas exportés depuis et la croissance continue des cas confirmés en Chine, nous devons actuellement nous attendre à au moins 50 000 cas à ce jour.
Avec nos estimations précédentes de l'âge de l'épidémie et des informations sur la période infectieuse, nous pouvons estimer les plages plausibles de R0 en utilisant un modèle de ramifications.

** Nous trouvons des estimations plausibles de R0 entre 1,8 et 3,5. **

Si nous supposons que l'épidémie a commencé au début de novembre 2019 (il y a 12 semaines), nous trouvons que R0 devrait se situer entre 1,8 et 2,5, selon la taille ('n') de l'épidémie.
<div>
  <img alt="graphique des estimations de R0 avec début d'épidémie il y a 12 semaines" width="500" src="https://data.nextstrain.org/ncov_branching-R0-early_2020-01-29.png"/>
</div>

Si nous supposons un démarrage plus récent, début décembre 2019 (il y a 8 semaines), les estimations de R0 se situent entre 2,2 et 3,5:
<div>
  <img alt="graphique des estimations de R0 avec début d'épidémie il y a 8 semaines" width="500" src="https://data.nextstrain.org/ncov_branching-R0-recent_2020-01-29.png"/>
</div>
Ces estimations sont globalement cohérentes avec celles d'autres scientifiques, qui se situent principalement entre R0 = 2-3, voir par exemple <a href="https://www.biorxiv.org/content/10.1101/2020.01.25.919787v1">cette prépublication</a>.
De manière importante, R0 est une variable qui dépend fortement du contexte socio-économique et des mesures de contrôle des infections.
```

# [Crédit scientifique](https://nextstrain.org/ncov/2020-01-30?d=map&c=author)

Nous tenons à souligner le superbe travail accompli si rapidement par tous les scientifiques impliqués dans cette épidémie, et en particulier ceux qui travaillent en Chine.
Ce n'est que par le partage rapide des données génomiques et des métadonnées que de telles analyses sont possibles.

<br>

Les génomes nCoV ont été généreusement partagés par les scientifiques de:

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

# [Crédit scientifique détaillé](https://nextstrain.org/ncov/2020-01-30?d=map&c=author)

Ces données ont été partagées par [GISAID](https://gisaid.org).
Nous remercions chaleureusement leurs contributions.

<br>

Sur la droite nous indiquons les séquences partagées par chaque laboratoire.

```auspiceMainDisplayMarkdown

Les génomes nCoV ont été généreusement partagés par les scientifiques de:

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
