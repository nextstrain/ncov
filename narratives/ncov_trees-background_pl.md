---
title: Jak interpretować drzewa filogenetyczne
authors:
  - Emma Hodcroft
  - Nicola Müller
  - James Hadfield
  - Sidney M. Bell
  - Richard Neher
  - Trevor Bedford
authorLinks:
  - https://neherlab.org/emma-hodcroft.html
  - https://bedford.io/team/nicola-mueller/
  - https://bedford.io/team/james-hadfield/
  - https://twitter.com/sidneymbell
  - https://neherlab.org/richard-neher.html
  - https://bedford.io/team/trevor-bedford/
affiliations: "Fred Hutch, Seattle, USA; Biozentrum, Basel, Switzerland; CZI, CA, USA"
translators:
  - Anna Fijarczyk
  - Marta Niedzicka
translatorLinks:
  - https://twitter.com/afijarczyk
  - https://www.researchgate.net/profile/Marta_Niedzicka
date: "13 marca 2020"
dataset: "https://nextstrain.org/ncov/2020-03-11?d=tree&legend=open&c=country"
abstract: "Ta prezentacja pokazuje jak czytać i interpretować drzewa filogenetyczne w epidemiologii genomowej. Prezentacja jest dostosowana do wyświetlania w przeglądarkach komputerowych."
---
<!-- Translators: Only text after : in the above ^ needs to be translated -->
<!-- Comment tags like these do not need to be translated, they are only to help you! -->
<!-- Ensure that links always end in a 'letter' (. counts) If some kind of text doesn't follow them, it breaks the slide. -->
<!-- numbers can be tagged ilke this: 161</tag> - this is just for us to help find them to update! Just leave in the </tag> bit. -->

<!-- This is left-side text -->
<!--abstract: "This narrative explains how to read and interpret the phylogenetic trees that inform genomic epidemiology. This website is optimized for display on desktop browsers."-->
# [Spis treści](https://nextstrain.org/ncov/2020-03-11?d=tree&legend=open&c=country)

* [Jak sieci transmisji są powiązane z drzewami filogenetycznymi](https://nextstrain.org/narratives/trees-background/pl?n=2)?  
* [Jak czytać drzewo](https://nextstrain.org/narratives/trees-background/pl?n=3)?  
* [Jak panel "zmienności" odnosi się do drzewa](https://nextstrain.org/narratives/trees-background/pl?n=4)?   
* [Mierzenie różnic przy pomocy dywergencji genetycznej](https://nextstrain.org/narratives/trees-background/pl?n=5).  
* [Mierzenie różnic w czasie](https://nextstrain.org/narratives/trees-background/pl?n=6).  
* [Datowanie początku wybuchu epidemii](https://nextstrain.org/narratives/trees-background/pl?n=7).  
* [Jak należy interpretować cechy (kolory) na drzewie](https://nextstrain.org/narratives/trees-background/pl?n=8)?  
* [Jak mapa odnosi się do drzewa](https://nextstrain.org/narratives/trees-background/pl?n=9)?  
* [Lektura zaawansowana: miary niepewności drzewa](https://nextstrain.org/narratives/trees-background/pl?n=10).  
* [O zestawie danych](https://nextstrain.org/narratives/trees-background/pl?n=11).  

<!-- No right-side text -->


<!-- Table of Contents(https://nextstrain.org/ncov/2020-03-11?d=tree&legend=open&c=country)
* [How are transmission networks related to phylogenetic trees](https://nextstrain.org/narratives/trees-background?n=1)?  
* [How do I read a tree](https://nextstrain.org/narratives/trees-background?n=2)?  
* [How does the "diversity" panel relate to the tree](https://nextstrain.org/narratives/trees-background?n=3)?   
* [Measuring differences with genetic divergence](https://nextstrain.org/narratives/trees-background?n=4).  
* [Measuring differences over time](https://nextstrain.org/narratives/trees-background?n=5).  
* [Dating the start of an outbreak](https://nextstrain.org/narratives/trees-background?n=6)?  
* [How should I interpret traits (colors) on the tree](https://nextstrain.org/narratives/trees-background?n=7)?  
* [How does the map relate to the tree](https://nextstrain.org/narratives/trees-background?n=8)?  
* [Advanced reading: uncertainty in trees](https://nextstrain.org/narratives/trees-background?n=9).  
* [About the dataset](https://nextstrain.org/narratives/trees-background?n=10).  -->



<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [Jak sieci transmisji są powiązane z drzewami filogenetycznymi?](https://nextstrain.org/ncov/2020-03-11?d=tree&p=full)
Patogeny rozprzestrzeniają się poprzez szybką replikację w jednym gospodarzu, a nastepnie transmisję do kolejnego gospodarza. Epidemia może wybuchnąć jednynie wtedy gdy pojedyncza infekcja prowadzi do więcej niż jednej kolejnej infekcji. 
<br><br>
Podczas gdy patogen namnaża i rozprzestrzenia się, jego genom podlega wielokrotnej replikacji co prowadzi do akumulacji losowych mutacji (błędów kopiowania) w genomie; jest to naturalny proces. Takie losowe mutacje pomagają w śledzeniu rozprzestrzeniania się patogenu i dowiadywaniu się o drogach i historii transmisji.

<!-- How are transmission networks related to phylogenetic trees?
Pathogens spread through rapid replication in one host followed by transmission to another host. An epidemic can only take off when one infection results in more than one subsequent infections.
As the pathogen replicates and spreads, its genome needs to be replicated many times and random mutations (copying mistakes) will accumulate in the genome; this is normal. Such random mutations can help to track the spread of the pathogen and learn about its transmission routes and dynamics.-->

<!-- This is right-side text -->
```auspiceMainDisplayMarkdown
# Przykład
<div width="50%" margin="auto">
<p>
<img width="500px" alt="cartoon showing how transmission tree and phylogenetic tree relate" src="https://github.com/nextstrain/nextstrain.org/raw/master/static-site/content/help/01-general/figures/infection_tree_combined.png"/>
</p>
<p>
Powyższa ilustracja pokazuje schemat drzewa transmisji. Każde kółko przedstawia jeden przypadek (zainfekowaną osobę), gdzie linie poziome wskazują czas trwania infekcji. Pionowo połączone kółka przedstawiają transmisję z jednej osoby na drugą. 
<br> <br>
W tym przypadku, widzimy pełny obraz drzewa transmisji. W praktyce, mamy dostęp jedynie do podzbioru przypadków (niebieskie); drzewo transmisji jest nieznane i zazwyczaj znamy tylko przybliżone oszacowania liczby przypadków. Sekwencje genomów pozwalają nam rekonstruować niektóre części drzewa transmisji. W tym przypadku widzimy trzy mutacje (romby) na drzewie. Sekwencje, które współdzielą te same mutacje, są bliżej ze sobą spokrewnione, w związku z tym te mutacje pozwalają nam pogrupować próbki w grupy blisko spokrewnionych wirusów pochodzących z tego samego ciągu transmisji.
</p>
</div>
```

<!-- The illustration above shows a sketch of a transmission tree. Each circle represents a case (infected person), with horizontal lines indicating the duration of their infection. Connected cases represent transmissions from one person to the next.
<br> <br>
Here, we see the full picture of the transmission tree. In practice, however, only a subset of cases are sampled (blue); the transmission tree is unknown and typically only rough estimates of case counts are available. Genome sequences allow us to infer parts of the transmission tree. In this example, three mutations (little diamonds) are indicated on the tree. Sequences that have the same mutations are more closely related, so these mutations allow us to group samples into clusters of closely related viruses that belong to the same transmission chains. -->


<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [Jak czytać drzewo?](https://nextstrain.org/ncov/2020-03-11)

Oś x drzewa przedstawia stopień zróżnicowania (w skali czasu lub dywergencji genetycznej -- dojdziemy do tego później). Oś y jedynie pomaga romieścić wszystkie gałęzie tak aby były widoczne; nie ma żadnych jednostek pomiaru.
<br><br>
Na krańcach gałęzi znajdują się próby (tzn. niebieskie przypadki z poprzedniego slajdu). Węzły wewnątrz drzewa przedstawiają przypadki, dla których nie mamy prób, ale które uważamy, że były przodkami wszystkich przypadków, które od nich odchodzą (tzn. czerwone węzły z poprzedniego slajdu). Te relacje są wnioskowane na podstawie analizy wzorów mutacji obserwowanych wśród próbkowanych przypadków.

<!-- The x axis of a tree represents the degree of difference (in time or genetic divergence -- we'll get to that next). The y axis just helps spread things out so we can see everything; it doesn't have any units of measurement.
The tips of the tree represent samples (i.e., blue cases from the last slide). The internal nodes represent cases that weren't sampled, but that we think were the source of all the cases descendant from them (i.e., the red nodes from the last slide). These relationships are inferred by analyzing the pattern of mutations observed in the sampled cases. -->

<!-- This is right-side text -->
```auspiceMainDisplayMarkdown
## Przykład
<div width="50%" margin="auto">
<p>
<img width="700px" alt="Example phylogeny where all or only a subset of cases are included in the final phylogeny" src="https://github.com/nextstrain/nextstrain.org/raw/master/static-site/content/help/01-general/figures/toy_alignment_tree.png"/>
</p>
<p>
Powyżej po lewej stronie, widzimy ilustrację drzewa filogenetycznego, na którym mutacje są oznaczone kółkami. Po prawej znajdują się sekwencje odpowiadające poszczególnym gałęziom drzewa (próbkom), również z mutacjami zaznaczonymi kółkami. Widzimy, że sekwencje, które dzielą te same mutacje grupują się razem. Gdy pokrewne przypadki zdają się być połączone pionową linią na drzewie, jak np. A i B, to oznacza, że nie ma różnic między nimi - ich sekwencje są identyczne.
<br><br>
Kiedy sekwencja znajduje się na długiej gałęzi, np. C czy E, to oznacza, że posiada unikatowe mutacje nie występujące w żadnej innej sekwencji. Im dłuższa linia, tym więcej mutacji. A i B również mają unikatowe mutacje (zielone kółko) nie współdzielone z żadną inną sekwencją, ale względem siebie są identyczne.
<br><br>
Na podstawie tego drzewa, możemy wyciągnąć wniosek, że A i B są blisko ze sobą spokrewnione, oraz D i E są blisko ze sobą spokrewnione. A i B są bliżej spokrewnione z C niż do D i E <!-- tu byl chyba blad, powinno byc Di E zamiast D i C? -->
</p>

<!-- Above, we see an illustration with a phylogenetic tree on the left, where mutations are shown as colored circles. On the right are the corresponding sequences, also with mutations shown as colored circles. We can see that sequences that share the same mutations group together. When sequences appear linked by a flat vertical line, like A and B, this means there are no differences between them – their sequences are identical.
<br><br>
When a sequence sits on a long line on its own, like C or E, this means it has unique mutations not found in other sequences. The longer the line, the more mutations.
A and B also have unique mutations (the green circle) not shared by the other sequences, but they are identical to each other.
<br><br>
Based on this tree, we would conclude that A & B closely related to each other, and D & E are closely related to each other. A & B are more closely related to C than they are to D & C. -->


### Materiały dodatkowe
* [How to read a tree: instruktaż z Arctic Network](https://artic.network/how-to-read-a-tree.html).  
* [How to read a tree: wideo z Khan academy](https://www.khanacademy.org/science/high-school-biology/hs-evolution/hs-phylogeny/a/phylogenetic-trees).  

</div>

```


<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [Jak panel "zmienności" odnosi się do drzewa?](https://nextstrain.org/ncov/2020-03-11?d=tree,entropy&c=gt-ORF1b_314&legend=open)

Spójrzmy na pierwsze 169</tag> sekwencji wirusa SARS-CoV-2 (wirus, który wywołuje COVID-19), które zostały publicznie udostępnione. Tak jak na poprzedniej stronie, przyrównaliśmy te sekwencje wirusa do siebie (możesz sprawdzić jak wszystkie analizy zostały wykonane [na GitHub](https://github.com/nextstrain/ncov)).
<br><br>
U góry znajduje się drzewo filogenetyczne, a poniżej wykres słupkowy przedstawiający zmienność (tzn. mutacje) w genomie.
Bez tych mutacji, nie bylibyśmy w stanie zbudować drzewa, wobec tego te dwa elementy są ze sobą ściśle powiązane.
<br><br>
W tym panelu "zmienności", oś pozioma przedstawia kolejne pozycje w genomie wirusa (w sumie około trzydzieści tysięcy!). 
Oś pionowa pokazuje jak wiele zmienności jest w każdej pozycji.
<br><br>
Drzewo pokolorowaliśmy na podstawie jednej z tych mutacji -- w tym przypadku mutacji w kodonie 314 w genie "ORF1b".
Nie ma powodu myśleć zawczasu, że ta mutacja jest funkcjonalna (tzn. wywołuje jakąś biologiczną zmianę).
Dokładnie takie mutacje wykorzystujemy do definiowania związków pomiędzy sekwencjami i konstruowania drzewa. 

<!-- There is NO right-side text -->

<!-- Let's take a look at the first 169</tag> strains of SARS-CoV-2 (the virus that causes COVID-19) that have been publicly shared. Just as on the last page, we built an alignment of these viral sequences (you can see how all of the analyses mentioned here were done [on GitHub](https://github.com/nextstrain/ncov)).
<br><br>
Here we're displaying the phylogenetic tree above a bar-chart showing the variation (i.e. mutations) in the genome.
Without these mutations we couldn't build the tree, so the two are intimitely connected.
<br><br>
In this "divisity" panel, the horizontal axis is each site in the viral genome (all thirty thousand or so of them!).
The vertical axis indicates how much variability there is at each site.
<br><br>
We've coloured the tree according to one of these mutations -- in this case codon 314 in the gene "ORF1b".
There's no a priori reason to think this mutation is a functional mutation (i.e. conferring any biological change).
It is precisely mutations such as this which we use to define the relationships between sequences and construct the tree. -->


<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [Mierzenie różnic przy pomocy dywergencji genetycznej](https://nextstrain.org/ncov/2020-03-11?c=num_date&d=tree&m=div)
To jest drzewo filogenetycze zbudowane na podstawie pierwszych 169</tag> sekwencji wirusa SARS-CoV-2 (który wywołuje COVID-19), które zostały publicznie udostępnione.
<br><br>
Oś pozioma odpowiada dywergencji, czyli liczbie zmian (mutacji) w genomie, względem korzenia drzewa (tzn. początku wybuchu epidemii). 
Niektóre sekwencje mają 0 mutacji -- co oznacza że są identyczne z korzeniem (środkiem) drzewa.
Inne wirusy mają między 1 a 11 mutacji.
<br><br>
Obecnie nasze drzewo mało przypomina prawdziwe 'drzewo'. Wiele sekwencji jest identycznych - znajdują się obok siebie na pionowych liniach, jak A i B (niektóre z nich położone są po lewej stronie drzewa).
Inne mają unikatowe lub współdzielone mutacje i są położone na poziomym liniach, lub 'gałęziach', skierowanych na prawo.
Możesz sprawdzić ile mutacji znajduje się na gałęzi przesuwając nad nią myszkę. 

<!-- There is NO right-side text -->

<!-- This is a phylogeny of the first 169</tag> strains of SARS-CoV-2 (the virus that causes COVID-19) that have been publicly shared.
<br><br>
Here, the horizontal axis indicates divergence, which is the number of changes (mutations) in the genome, relative to the root of the tree (i.e., the start of the outbreak).
Some sequences may have zero mutations -- meaning they are all identical to the root (center) of the tree.
Other viruses have between one and eleven mutations.
<br><br>
At the moment, this may not look much like a ‘tree’. Many of the sequences are identical – they sit together on vertical lines like A and B (some are on the left-most part of the tree).
Others have unique or shared mutations and so sit on lines, or ‘branches’, going to the right.
You can see how many mutations a branch has by hovering your mouse over it. -->


<!-- ############ SLIDE BREAK ############# -->
# [Mierzenie różnic w czasie](https://nextstrain.org/ncov/2020-03-11?c=num_date&d=tree&legend=open)
Możemy również pokazać jak wirus rozprzestrzenia się w czasie używając daty pobrania próby jako osi x.
W tym przypadku, oś x przedstawia czas uzyskania próby. Krańce gałęzi odpowiadają dacie pobrania danej próby. Daty dla węzłów wewnątrz drzewa -- "brakująch przypadków" -- są wnioskowane na podstawie tego, kiedy ich potomkowie zostali próbkowani oraz tempa z jakim wirus mutuje. 
<br><br>
Zwróć uwagę, jak wiele sekwencji, które wcześniej były położone w jednej pionowej linii (wskazującej na identyczne genomy) są teraz rozłożone w czasie. 
Tak dzieje się, kiedy tempo, z jakim wirus mutuje, jest niewiele wolniejsze od tempa, z jakim się rozprzestrzenia. 
Przewiń do poprzedniego slajdu i porównaj z obecnym, aby zobaczyć jak drzewo się zmienia.

<!-- There is NO right-side text -->

<!-- We can also visualize how the virus has spread over time by using the sampling date as the x axis.
Here, the x axis represents the sampling date of each virus. The tips' positions reflect the date those samples were taken. The dates of internal nodes -- the "missing cases" -- are inferred based on when their descendants were sampled and the rate at which the virus mutates.
<br><br>
Notice how many sequences that previously sat in a line (indicating identical genomes) are now spread apart in time.
This happens when the rate at which the virus mutates is slightly slower than the rate at which is spreads.
You can scroll up and down between the previous slide and this one, to see how the tree changes. -->


<!-- ############ SLIDE BREAK ############# -->
# [Datowanie początku wybuchu epidemii](https://nextstrain.org/ncov/2020-03-11?c=num_date&d=tree&legend=open)

Genomikę możemy również wykorzystać w celu oszacowania daty wybuchu epidemii, nawet jeśli to było zanim zdaliśmy sobie z niej sprawę. 
Ponieważ możemy przypisać datę do każdej próbki i węzła na drzewie, możemy to wykorzystać do oszacowania daty pojawienia się "korzenia" drzewa. Korzeń przedstawia "ostatniego wspólnego przodka" wszystkich sekwencji SARS-CoV-2, które dotychczas uzyskaliśmy. Np. twoi dziadkowie są "ostatnimi wspólnymi przodkami" dla ciebie i wszystkich twoich kuzynów pierwszego stopnia.
<br><br>
Przesuń myszkę na najbardziej wysuniętą na lewo linię, a zobaczysz, że szacowana data tej epidemii jest pomiędzy połową listopada a połową grudnia 2019.

<!-- There is NO right-side text -->

<!-- We can also use genomics to assign a date to when an outbreak started, even if this was before we realized it was happening.
Because we can assign dates to each sample and node in the tree, we can use this to infer the date of the 'root' of the tree. This represents the "most recent common ancestor" of all the SARS-CoV-2 sequences we have so far. E.g., your grandparents are the "most recent common ancestors" of you and all of your first cousins.
<br><br>
If you mouse over the leftmost vertical line, you can see that the inferred start date is between mid-November and mid-December of 2019 for this particular outbreak. -->


<!-- ############ SLIDE BREAK ############# -->
<!-- This is left-side text -->
# [Jak należy interpretować cechy (kolory) na drzewie?](https://nextstrain.org/ncov/2020-03-11)
Drzewa filogenetyczne często zawierają dodatkowe informacje, takie jak lokalizacja kolekcji prób. Na tej podstawie możemy wnioskować o lokalizacji wewnętrznych węzłów drzewa (hipotetyczne pośrednie, niepróbkowane przypadki) używając modeli matematycznych. To nam może pomóc w zrozumieniu jak wirus się przemieszcza z jednego miejsca na inne.
<br><br>
Takie wnioski należy jednak interpretować z ostrożnością, ponieważ próbkowanie i sekwencjonowanie lub ich brak może istotnie wpłynąć na interpretację wyników.


<!-- Phylogenetic trees often contain additional information, such as the location of each sample collection. From this, we can infer the locations of internal nodes (hypothesized intermediate, unsampled cases) using mathematical models. This can help us understand how the virus is moving from one location to the next.
<br><br>
Interpreting these should, however, be done with caution, as the sampling and sequencing or lack thereof can significantly influence the interpretation.-->

<!-- This is right-side text -->
```auspiceMainDisplayMarkdown
# Przykład
<div width="50%" margin="auto">
<p>
<img width="700px" alt="Illustration showing how sampling effects interpretation of viral spread" src="https://github.com/nextstrain/nextstrain.org/raw/master/static-site/content/help/01-general/figures/introductions.png"/>
</p>
<p>
Po lewej stronie, pokazujemy hipotetyczne drzewo filogenetyczne z wszystkimi próbami, gdzie próby z dwóch lokalizacji zaznaczone są na pomarańczowo i niebiesko. Gdy podążamy w dół drzewa, widzimy trzy momenty kiedy kolor (lokalizacja) zmienia się z pomarańczowego na niebieski. Z tego możemy wnioskować, że nastąpiły trzy różne introdukcje z pomarańczowej lokalizacji do niebieskiej. 
<br><br>
Taka interpretacja silnie zależy od skali pobranych prób: w środkowym drzewie usunęliśmy jedną pomarańczową próbkę. Teraz widzimy tylko jedno przejście z pomarańczowego na niebieski, co sugeruje, że nastąpiła tylko jedna introdukcja do niebieskiej lokalizacji, i że nastąpiła ona wcześniej niż na pierwszym drzewie.
<br><br>
Na ostatnim przykładzie, mamy tylko jedną pomarańczową próbkę, co może świadczyć o jednej introdukcji z pomarańczowej do niebieskiej lokalizacji.
<br><br>
Podsumowując, takie wnioskowanie może być bardzo cenne, aczkolwiek musi być interpretowane z ostrożnością.
</p>
```
<!-- On the left, we show a fully sampled phylogenetic tree, with samples from two different locations denoted by orange and blue. As we walk down the tree, we observe three instances where the color (location) switches from orange to blue. From this, we would conclude that there were three different introductions from the orange location to the blue location.
<br><br>
But, this interpretation relies on sampling: in the middle tree, we've removed one orange sample. We now observe only one switch from orange to blue, suggesting that there was only one introduction into blue that happened much earlier.
<br><br>
In the last example, we have only one sequence from orange, which could lead us to think that there was one introduction from orange into blue.
<br><br>
Thus, while these inferences can be invaluable, they also must be interpreted with caution.-->

<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
<!-- ############ SLIDE BREAK ############# -->
<!-- This is left-side text -->
# [Jak mapa odnosi się do drzewa?](https://nextstrain.org/ncov/2020-03-11?d=tree,map&legend=closed)

Na tym slajdzie pokazujemy drzewo pokolorowane według lokalizacji pobranych próbek (i szacowanych lokalizacji w wewnętrznych węzłach).
Jeśli kilkniesz na ['Explore the data'](https://nextstrain.org/ncov), możesz obejrzeć animację pokazującą jak, według oszacowań, wirus rozprzestrzeniał się od momentu wybuchu epidemii.

<!-- Here, we show the tree colored by the location of each sample (and inferred location for each internal node).
If you click ['Explore the data'](https://nextstrain.org/ncov), you can play an animation of how the inferred spread of the virus over the course of the outbreak. -->

<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [Lektura zaawansowana: miary niepewności drzewa](https://nextstrain.org/ncov/2020-03-11)
Wcześniej pisaliśmy, jak wewnętrzne węzły przedstawiają _hipotetyczne_ niepróbkowane przypadki. W rzeczywistości, wszystkie drzewa są _hipotezami_ opisującymi jak patogen ewoluuje i przemieszcza się w czasie. Drzewa, które prezentujemy na Nextstrain są oszacowaniami punktowymi<!-- point estimates? --> -- tzn. wersjami historii, dla których prawdopodobieństwo obserwowania danych, jakimi dysponujemy, jest największe. <!-- nie jestem pewna czy 'zmaksymalizowane' nie jest zbyt techniczne-->
<br><br>
Niemniej jednak, takie oszacowania zawsze obarczone są niepewnością. Mówiąc ogólnie, części drzewa, które są gęściej próbkowane, są też bardziej wiarygodne; obszary słabo próbkowane są mniej pewne.

<!-- Earlier, we talked about how internal nodes represent _hypothesized_ unsampled cases. In fact, all trees represent _hypotheses_ about how a pathogen has evolved and moved over time. The trees we present on Nextstrain are point estimates -- that is, the version of this history that maximizes the probability of observing the data that we do.
<br><br>
However, there is always uncertainty in these estimates. Generally speaking, parts of the tree that are densely sampled are more certain; areas that are sparsely sampled are less certain.-->


```auspiceMainDisplayMarkdown
# Przykład
<div width="50%" margin="auto">
<p>
<img width="700px" alt="Illustration of the uncertainty inherent in tree reconstruction" src="https://github.com/nextstrain/nextstrain.org/raw/c69bfd0750c284ff12f33682f8d82848e13d9e15/static-site/content/help/01-general/figures/hcov_densitree.png"/>
</p>
</div>
```

<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [Podziękowania](https://nextstrain.org/ncov/2020-03-05?d=map&c=author)

Wyrażamy podziękowania dla wszystkich naukowców zaangażowanych w ten wybuch pandemii za wspaniałą i śpieszną pracę. W szczególności dziękujemy tym, którzy pracują w Chinach. Analizy takie jak ta nie byłyby możliwe gdyby nie natychmiastowe publikowanie zsekwencjonowanych genomów oraz metadanych z nimi związanych.

<br><br>

Ponadto szczególnie dziękujemy [GISAID](https://gisaid.org) za udostępnienie platformy do wymiany tego typu danych.

<!-- Do not need to translate insitutions names -->
<!-- This is right-side text -->
```auspiceMainDisplayMarkdown

Jesteśmy wdzięczni za dane zebrane w tych ośrodkach badawczych:

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
# [Szczegółowe podziękowania](https://nextstrain.org/ncov/2020-03-05?d=map&c=author)

Te dane były udostępnione poprzez platformę [GISAID](https://gisaid.org).
Szczególnie dziękujemy za wkład ich drużyny. 

<br><br>

Po prawej stronje znajduje się szczegółowy spis sekwencji opublikowanych przez poszczególne instytucje badawcze. 

<!-- This is right-side text -->
```auspiceMainDisplayMarkdown

Genomy wirusa SARS-CoV-2 zostały udostępnione przez naukowców pracujących w niżej wymienionych instytucjach badawczych:

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
