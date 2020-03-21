---
title: Come interpretare gli alberi filogenetici
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
affiliations: "Fred Hutch, Seattle, USA; Biozentrum, Basel, Switzerland; Chan Zuckerberg Initiative, CA, USA"
translators:
  - Alice Ledda
  - Alessia Lepore
  - Sara Pilia
date: "13 marzo 2020"
dataset: "https://nextstrain.org/ncov/2020-03-11?d=tree&legend=open&c=country"
abstract: "Questo testo spiega come leggere e interpretare gli alberi filogenetici che caratterizzano l'epidemiologia genomica. Questo sito è ottimizzato per visualizzazione da browser desktop."
---
<!-- Translators: Only text after : in the above ^ needs to be translated -->
<!-- Comment tags like these do not need to be translated, they are only to help you! -->
<!-- Ensure that links always end in a 'letter' (. counts) If some kind of text doesn't follow them, it breaks the slide. -->
<!-- numbers can be tagged ilke this: 161</tag> - this is just for us to help find them to update! Just leave in the </tag> bit. -->

<!-- This is left-side text -->
<!-- # [Table of Contents](https://nextstrain.org/ncov/2020-03-11?d=tree&legend=open&c=country) -->
# [Indice](https://nextstrain.org/ncov/2020-03-11?d=tree&legend=open&c=country)

<!-- * [How are transmission networks related to phylogenetic trees](https://nextstrain.org/narratives/trees-background?n=1)? -->
* [Qual è la relazione tra le reti di trasmissione e gli alberi filogenetici](https://nextstrain.org/narratives/trees-background/it?n=2)?
<!-- * [How do I read a tree](https://nextstrain.org/narratives/trees-background?n=2)?  -->
* [Come faccio a leggere un albero](https://nextstrain.org/narratives/trees-background/it?n=3)?
<!-- * [How does the "diversity" panel relate to the tree](https://nextstrain.org/narratives/trees-background?n=3)?  -->
* [Qual è la relazione tra il pannello "diversità" e l'albero](https://nextstrain.org/narratives/trees-background/it?n=4)?
<!-- * [Measuring differences with genetic divergence](https://nextstrain.org/narratives/trees-background?n=4).  -->
* [Misurare le differenze usando le distanze genetiche](https://nextstrain.org/narratives/trees-background/it?n=5).
<!-- * [Measuring differences over time](https://nextstrain.org/narratives/trees-background?n=5).  -->
* [Misurare le differenze nel tempo](https://nextstrain.org/narratives/trees-background/it?n=6).
<!-- * [Dating the start of an outbreak](https://nextstrain.org/narratives/trees-background?n=6) -->
* [Datare l'inizio dell'epidemia](https://nextstrain.org/narratives/trees-background/it?n=7).
<!-- * [How should I interpret traits (colors) on the tree](https://nextstrain.org/narratives/trees-background?n=7)?  -->
* [Come devo interpretare i tratti (colori) sull'albero](https://nextstrain.org/narratives/trees-background/it?n=8)?
<!-- * [How does the map relate to the tree](https://nextstrain.org/narratives/trees-background?n=8)?  -->
* [Che relazione c'è tra la cartina e l'albero](https://nextstrain.org/narratives/trees-background/it?n=9)?
<!-- * [Advanced reading: uncertainty in trees](https://nextstrain.org/narratives/trees-background?n=9).  -->
* [Approfondimenti: incertezze sull'albero](https://nextstrain.org/narratives/trees-background/it?n=10).
<!-- * [About the dataset](https://nextstrain.org/narratives/trees-background?n=10). -->
* [Dettagli sul dataset](https://nextstrain.org/narratives/trees-background/it?n=11).


<!-- No right-side text -->


<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
<!-- * [How are transmission networks related to phylogenetic trees](https://nextstrain.org/narratives/trees-background?n=1)? -->
# [Qual è la relazione tra le reti di trasmissione e gli alberi filogenetici](https://nextstrain.org/ncov/2020-03-11?d=tree&p=full)

<!-- Pathogens spread through rapid replication in one host followed by transmission to another host. An epidemic can only take off when one infection results in more than one subsequent infections. -->
I patogeni si diffondono attraverso una replicazione rapida in un organismo ospite, seguita dalla trasmissione ad un altro organismo ospite. <!--Un'epidemia può verificarsi solo quando un'infezione risulta in più di una infezione successiva. -->
Parliamo di un'epidemia solo quando un'infezione in organismo ospite dà succesivamente luogo a più di un'infezione in altri organismi.
<br><br>
<!-- As the pathogen replicates and spreads, its genome needs to be replicated many times and random mutations (copying mistakes) will accumulate in the genome; this is normal. Such random mutations can help to track the spread of the pathogen and learn about its transmission routes and dynamics. -->
Quando il patogeno si replica e si diffonde, il suo genoma ha bisogno di essere replicato molte volte e mutazioni casuali (errori di replicazione) si accumulano sul genoma; questo è un fatto normale. Tali mutazioni casuali possono aiutare a tracciare la diffusione del patogeno e imparare cose sulle sue rotte di trasmissione e la sua dinamica.

<!-- This is right-side text -->
```auspiceMainDisplayMarkdown
<!-- # An example -->
# Un esempio
<div width="50%" margin="auto">
<p>
<!-- <img width="500px" alt="cartoon showing how transmission tree and phylogenetic tree relate" src="https://github.com/nextstrain/nextstrain.org/raw/master/static-site/content/help/01-general/figures/infection_tree_combined.png"/> -->
<img width="500px" alt="disegno che mostra la relazione tra l'albero di trasmissione e l'albero filogenetico" src="https://github.com/nextstrain/nextstrain.org/raw/master/static-site/content/help/01-general/figures/infection_tree_combined.png"/>
</p>
<p>
<!-- The illustration above shows a sketch of a transmission tree. Each circle represents a case (infected person), with horizontal lines indicating the duration of their infection. Connected cases represent transmissions from one person to the next. -->
L'immagine qui sopra mostra un'illustrazione schematica di un albero di trasmissione. Ciascun cerchio rappresenta un caso (una persona infetta), e le linee orizzontali indicano la durata della loro infezione. Casi connessi indicano la trasmissione da una persona alla successiva.
<br> <br>
<!-- Here, we see the full picture of the transmission tree. In practice, however, only a subset of cases are sampled (blue); the transmission tree is unknown and typically only rough estimates of case counts are available. Genome sequences allow us to infer parts of the transmission tree. In this example, three mutations (little diamonds) are indicated on the tree. Sequences that have the same mutations are more closely related, so these mutations allow us to group samples into clusters of closely related viruses that belong to the same transmission chains. -->
Qui vediamo una immagine completa dell'albero di trasmissione. In pratica, comunque, solo un sottoinsieme dei casi vengono censiti (identificati in blu); l'albero di trasmissione non è noto e sono disponibili solo stime rozze sul numero dei casi. Le sequenze genomiche ci permettono di inferire statisticamente alcune parti dell'albero di trasmissione. In questo esempio abbiamo evidenziato tre mutazioni sull'albero (piccoli rombi). Siccome le sequenze che hanno le stesse mutazioni sono relazionate più strettamente, queste mutazioni ci permettono di raggruppare i campioni in gruppi di virus simili che appartengono alle stesse catene di trasmissione.
</p>
</div>
```


<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
<!-- # [How do I read a tree?](https://nextstrain.org/ncov/2020-03-11) -->
# [Come faccio a leggere un albero?](https://nextstrain.org/ncov/2020-03-11)

<!-- The x axis of a tree represents the degree of difference (in time or genetic divergence -- we'll get to that next). The y axis just helps spread things out so we can see everything; it doesn't have any units of measurement. -->
In un albero, sull'asse delle x viene riportato il grado di differenza (in unità di tempo o di distanza genetica - ci arriveremo). L'asse y serve solo a allargare il disegno in modo che si possano vedere meglio i dettagli; non ha unità di misura.
<!-- AL: (in unità di tempo o di distanza genetica - lo spiegeremo più avanti) -->
<br><br>
<!-- The tips of the tree represent samples (i.e., blue cases from the last slide). The internal nodes represent cases that weren't sampled, but that we think were the source of all the cases descendant from them (i.e., the red nodes from the last slide). These relationships are inferred by analyzing the pattern of mutations observed in the sampled cases. -->
Le foglie dell'albero rappresentano i campioni raccolti (per esempio i casi blu della slide precedente). I nodi interni rappresentano casi che non sono stati osservati, ma che si pensa siano l'antenato di tutti i casi che discendono da loro (per esempio i nodi rossi nella slide precedente). Queste relazioni vengono inferite statisticamente usando la struttura delle mutazioni nei casi osservati.

<!-- This is right-side text -->
```auspiceMainDisplayMarkdown
<!-- ## An example -->
## Un esempio
<div width="50%" margin="auto">
<p>
<!-- <img width="700px" alt="Example phylogeny where all or only a subset of cases are included in the final phylogeny" src="https://github.com/nextstrain/nextstrain.org/raw/master/static-site/content/help/01-general/figures/toy_alignment_tree.png"/> -->
<img width="700px" alt="Esempio di un albero filogenetico in cui vengono inclusi tutti o solo una parte dei casi" src="https://github.com/nextstrain/nextstrain.org/raw/master/static-site/content/help/01-general/figures/toy_alignment_tree.png"/>
</p>
<p>
<!-- Above, we see an illustration with a phylogenetic tree on the left, where mutations are shown as colored circles. On the right are the corresponding sequences, also with mutations shown as colored circles. We can see that sequences that share the same mutations group together. When sequences appear linked by a flat vertical line, like A and B, this means there are no differences between them – their sequences are identical. -->
In alto, sulla sinistra, si vede un'illustrazione con un albero filogenetico in cui le mutazioni sono indicate come cerchi colorati. Sulla destra vengono riportate le sequenze corrispondenti, con le mutazioni evidenziate nello stesso colore dei cerchi. Si vede chiaramente che le sequenze che condividono le stesse mutazioni stanno nello stesso gruppo. Quando le sequenze appaiono legate de una riga verticale, come A e B, significa che non ci sono differenze tra loro - queste sequenze sono identiche.
<br><br>
<!--When a sequence sits on a long line on its own, like C or E, this means it has unique mutations not found in other sequences. The longer the line, the more mutations. -->
Quando una sequenza si trova da sola in fondo ad una lunga linea, come C o E, significa che ha mutazioni uniche che non sono state trovate in altre sequenze. Quanto più è lunga la linea, più è grande il numero di mutazioni.
<!-- AL: significa che ha mutazioni uniche -->
<!-- A and B also have unique mutations (the green circle) not shared by the other sequences, but they are identical to each other.-->
Anche A e B hanno mutazioni uniche (il cerchio verde) che non sono condivise da altre sequenze, ma sono identiche tra di loro.
<br><br>
<!-- Based on this tree, we would conclude that A & B closely related to each other, and D & E are closely related to each other. A & B are more closely related to C than they are to D & C. -->
Basandoci su questo albero, concluderemmo che A e B sono strettamente relazionate tra di loro, così come D ed E. A e B sono più strettamente relazionate con C di quanto non lo siano D ed E.
</p>

<!--### Further reading-->
### Approfondimenti

<!-- * [How to read a tree: tutorial from Arctic Network](https://artic.network/how-to-read-a-tree.html). -->
* [Come leggere un albero filogenetico: guida di Arctic Network](https://artic.network/how-to-read-a-tree.html).

<!-- * [How to read a tree: video from Khan academy](https://www.khanacademy.org/science/high-school-biology/hs-evolution/hs-phylogeny/a/phylogenetic-trees). -->
* [Come leggere un albero filogenetico: video di Khan academy](https://www.khanacademy.org/science/high-school-biology/hs-evolution/hs-phylogeny/a/phylogenetic-trees).


</div>

```


<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
<!-- # [How does the "diversity" panel relate to the tree?](https://nextstrain.org/ncov/2020-03-11?d=tree,entropy&c=gt-ORF1b_314&legend=open) -->
# [Qual è la relazione tra il pannello "diversità" e l'albero?](https://nextstrain.org/ncov/2020-03-11?d=tree,entropy&c=gt-ORF1b_314&legend=open)

<!-- Let's take a look at the first 169</tag> strains of SARS-CoV-2 (the virus that causes COVID-19) that have been publicly shared. Just as on the last page, we built an alignment of these viral sequences (you can see how all of the analyses mentioned here were done [on GitHub](https://github.com/nextstrain/ncov)). -->
Diamo un'occhiata ai primi 169</tag> ceppi virali del SARS-CoV-2 (il virus che causa COVID-19) che sono stati pubblicamente condivisi. Come abbiamo fatto anche nella pagina precedente, abbiamo costruito un allineamento di queste sequenze virali (potete vedere tutte le analisi che sono state fatte [su GitHub](https://github.com/nextstrain/ncov)).
<br><br>
<!-- Here we're displaying the phylogenetic tree above a bar-chart showing the variation (i.e. mutations) in the genome. -->
Qui mostriamo l'albero filogenetico sopra un diagramma che mostra le variazioni (mutazioni) lungo il genoma.
<!-- Without these mutations we couldn't build the tree, so the two are intimitely connected. -->
Senza queste mutazioni non potremmo costruire l'albero, quindi i due sono strettamente connessi.
<br><br>
<!--In this "divisity" panel, the horizontal axis is each site in the viral genome (all thirty thousand or so of them!). -->
In questo pannello di diversità, l'asse orizzontale rappresenta i trentamila siti lungo il genoma virale.
<!-- The vertical axis indicates how much variability there is at each site.-->
L'asse verticale indica quanta variabilità è stata riscontrata in ciascun sito.
<!-- We've coloured the tree according to one of these mutations -- in this case codon 314 in the gene "ORF1b". -->
Abbiamo colorato l'albero seguendo una di queste mutazioni -- in questo caso il codone 314 nel gene  "ORF1b".
<!-- There's no a priori reason to think this mutation is a functional mutation (i.e. conferring any biological change). -->
Non c'è nessuna ragione per pensare a priori che questa sia una mutazione funzionale (cioè che conferisca un cambiamento biologico).
<!-- It is precisely mutations such as this which we use to define the relationships between sequences and construct the tree. -->
Mutazioni come questa sono precisamente quelle che usiamo per definire una relazione tra le sequenze e costruire l'albero.

<!-- There is NO right-side text -->


<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [Misurare le differenze usando le distanze genetiche](https://nextstrain.org/ncov/2020-03-11?c=num_date&d=tree&m=div)
<!--This is a phylogeny of the first 169</tag> strains of SARS-CoV-2 (the virus that causes COVID-19) that have been publicly shared.-->
Questo è un albero filogenetico dei primi 169</tag> ceppi virali di SARS-CoV-2 (il virus che causa COVID-19) che sono stati condivisi pubblicamente.
<br><br>
<!-- AL: Altri virus hanno tra zero e undici mutazioni.(ossia è per altri virus generici e non degli altri virus in particolare) -->
<!-- Here, the horizontal axis indicates divergence, which is the number of changes (mutations) in the genome, relative to the root of the tree (i.e., the start of the outbreak). -->
Qui, l'asse orizzontale indica la divergenza, che è data dal numero di cambiamenti (mutazioni) nel genoma, relativamente alla radice dell'albero (che indica l'inizio dell'epidemia).
<!-- Some sequences may have zero mutations -- meaning they are all identical to the root (center) of the tree. -->
Alcune sequenze possono avere zero mutazioni -- significa che sono uguali alla radice al centro dell'albero.
<!-- Other viruses have between one and eleven mutations. -->
Gli altri virus hanno tra zero e undici mutazioni.
<br><br>
<!-- At the moment, this may not look much like a ‘tree’. Many of the sequences are identical – they sit together on vertical lines like A and B (some are on the left-most part of the tree). -->
Al momento, questo potrebbe non assomigliare molto ad un 'albero'. Molte di queste sequenze sono identiche - stanno insieme su linee verticali, come A e B (alcune sono nella parte più a sinistra dell'albero).
<!-- Others have unique or shared mutations and so sit on lines, or ‘branches’, going to the right. -->
Altre hanno mutazioni uniche o condivise, per questo le troviamo su linee orizzontali, o 'rami' dirette verso destra.
<!-- You can see how many mutations a branch has by hovering your mouse over it.-->
Si puo' vedere quante mutazioni sono contenute in un ramo passandoci sopra con il mouse.

<!-- There is NO right-side text -->

<!-- ############ SLIDE BREAK ############# -->
<!-- # [Measuring differences over time](https://nextstrain.org/ncov/2020-03-11?c=num_date&d=tree&legend=open) -->
# [Misurare le differenze nel tempo](https://nextstrain.org/ncov/2020-03-11?c=num_date&d=tree&legend=open)

<!--We can also visualize how the virus has spread over time by using the sampling date as the x axis. -->
Possiamo visualizzare come il virus si è diffuso nel tempo, usando la data di raccolta dei campioni sull'asse delle x.
<!-- Here, the x axis represents the sampling date of each virus. The tips' positions reflect the date those samples were taken. The dates of internal nodes -- the "missing cases" -- are inferred based on when their descendants were sampled and the rate at which the virus mutates.-->
Qui l'asse x rappresenta la data in cui ciascun virus è stato isolato. La posizione delle foglie identifica la data in cui i campioni sono stati raccolti. Le date dei nodi interni -- i 'casi mancanti' -- sono state inferite statisticamente basandosi sulle date in cui i loro discendenti sono stati identificati e sulla velocità a cui il virus muta.
<!-- AL: ..-sono state dedotte basandosi sulle date...e sulla velocità a cui il virus muta
<br><br>
<!-- Notice how many sequences that previously sat in a line (indicating identical genomes) are now spread apart in time. -->
È importante notare quante sequenze che prima erano sulla stessa linea (che indicava genomi identici) siano ora distanziate nel tempo.
<!-- This happens when the rate at which the virus mutates is slightly slower than the rate at which is spreads. -->
Questo succede quando la velocità a cui il virus muta è lievemente minore di quella a cui si diffonde.
<!-- You can scroll up and down between the previous slide and this one, to see how the tree changes. -->
Potete scorrere tra questa slide e la precedente per vedere meglio come cambia l'albero.
<!-- There is NO right-side text -->


<!-- ############ SLIDE BREAK ############# -->
<!-- # [Dating the start of an outbreak](https://nextstrain.org/ncov/2020-03-11?c=num_date&d=tree&legend=open) -->
# [Datare l'inizio dell'epidemia](https://nextstrain.org/ncov/2020-03-11?c=num_date&d=tree&legend=open)

<!-- We can also use genomics to assign a date to when an outbreak started, even if this was before we realized it was happening. -->
Possiamo anche usare la genetica per assegnare una data all'inizio dell'epidemia, anche se questo è successo ancora prima che ci accorgessimo dell'epidemia.
<!-- Because we can assign dates to each sample and node in the tree, we can use this to infer the date of the 'root' of the tree. This represents the "most recent common ancestor" of all the SARS-CoV-2 sequences we have so far. E.g., your grandparents are the "most recent common ancestors" of you and all of your first cousins.-->
Dato che possiamo assegnare delle date ad ogni campione e nodo nell'albero, possiamo usare questo fatto per inferire statisticamente la data della 'radice' dell'albero. Questo rappresenta 'l'antenato comune più recente' di tutte le sequenze del virus SARS-CoV-2 che abbiamo finora. Per esempio: i vostri nonni sono 'l'antenato comune più recente' vostro e dei vostri cugini di primo grado.
<!-- AL: del virus SARS-CoV-2 che abbiamo finora.
<!-- AL: cugini di primo grado.
<br><br>
<!-- If you mouse over the leftmost vertical line, you can see that the inferred start date is between mid-November and mid-December of 2019 for this particular outbreak. -->
Posizionando il mouse sulla linea verticale più a sinistra, si può vedere che la data d'inizio inferita per questa specifica epidemia è tra metà novembre e metà dicembre del 2019.

<!-- There is NO right-side text -->


<!-- ############ SLIDE BREAK ############# -->
<!-- This is left-side text -->
<!-- # [How should I interpret traits (colors) on the tree?](https://nextstrain.org/ncov/2020-03-11) -->
# [Come devo interpretare i tratti (colori) sull'albero?](https://nextstrain.org/ncov/2020-03-11)

<!-- Phylogenetic trees often contain additional information, such as the location of each sample collection. From this, we can infer the locations of internal nodes (hypothesized intermediate, unsampled cases) using mathematical models. This can help us understand how the virus is moving from one location to the next. -->
Gli alberi filogenetici contengono spesso informazioni aggiuntive, quali per esempio la provenienza geografica del campione. Da questo si può inferire statisticamente la localizzazione geografica dei nodi interni (ipotetici casi intermedi e non identificati) usando dei modelli matematici. Questo può aiutare a capire meglio come il virus si stia muovendo da un luogo ad un altro.
<!-- AL: ...dedurre la localizzazione geografica...
<br><br>
<!--Interpreting these should, however, be done with caution, as the sampling and sequencing or lack thereof can significantly influence the interpretation.-->
L'interpretazione di queste inferenze comunque deve essere fatta con cautela dato che i campioni e le sequenze, o la loro mancanza, possono fortemente influenzarne l'interpretazione.
<!-- AL: L'interpretazione di queste inferenze deve comunque essere fatta con cautela dato che i campioni e le sequenze, o la loro mancanza, può fortemente influenzarne l'interpretazione.

<!-- This is right-side text -->
```auspiceMainDisplayMarkdown
<!-- # An example -->
# Un esempio
<div width="50%" margin="auto">
<p>
<!-- <img width="700px" alt="Illustration showing how sampling effects interpretation of viral spread" src="https://github.com/nextstrain/nextstrain.org/raw/master/static-site/content/help/01-general/figures/introductions.png"/> -->
<img width="700px" alt="Illustrazione che mostra come il modo in cui vengono raccolti i campioni influenzi l'interpretazione della diffusione del virus" src="https://github.com/nextstrain/nextstrain.org/raw/master/static-site/content/help/01-general/figures/introductions.png"/> -->
</p>
<p>
<!-- On the left, we show a fully sampled phylogenetic tree, with samples from two different locations denoted by orange and blue. As we walk down the tree, we observe three instances where the color (location) switches from orange to blue. From this, we would conclude that there were three different introductions from the orange location to the blue location. -->
A sinistra, mostriamo un albero filogenetico completo, con campioni provenienti da due diverse località, indicate in blu e arancione. Scorrendo l'albero vediamo tre punti in cui il colore (località) cambia da arancione a blu. Perciò possiamo concludere che ci siano state tre diverse introduzioni dalla località arancione a quella blu.
<br><br>
<!-- But, this interpretation relies on sampling: in the middle tree, we've removed one orange sample. We now observe only one switch from orange to blue, suggesting that there was only one introduction into blue that happened much earlier. -->
Ma questa interpretazione si basa sulla disponibilità dei dati: abbiamo rimosso un campione arancione in mezzo all'albero. Ora osserviamo solo un cambiamento di colore, suggerendo una sola introduzione nella regione blu, avvenuta molto prima.
<br><br>
<!-- In the last example, we have only one sequence from orange, which could lead us to think that there was one introduction from orange into blue. -->
Nell'ultimo esempio abbiamo solo una sequenza dalla regione arancione, che ci potrebbe portare a pensare che ci sia stata una introduzione dall'arancione al blu.
<br><br>
<!-- Thus, while these inferences can be invaluable, they also must be interpreted with caution. -->
Quindi, per quanto queste inferenze statistiche siano importantissime, devono essere interpretate con cauzione.
</p>
```
<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
<!-- ############ SLIDE BREAK ############# -->
<!-- This is left-side text -->
<!-- # [How does the map relate to the tree?](https://nextstrain.org/ncov/2020-03-11?d=tree,map&legend=closed) -->
# [Che relazione c'è tra la mappa e l'albero?](https://nextstrain.org/ncov/2020-03-11?d=tree,map&legend=closed)

<!--Here, we show the tree colored by the location of each sample (and inferred location for each internal node). -->
Qui mostriamo l'albero colorato secondo la provenienza di ogni campione (e provenienza inferita statisticamente, per i nodi interni).
<!-- AL: ...(e provenienza dedotta... -->
<!-- AL: ... diffusione dedotta.. -->
<!-- If you click ['Explore the data'](https://nextstrain.org/ncov), you can play an animation of how the inferred spread of the virus over the course of the outbreak. -->
Cliccando su ['Esplora i dati'](https://nextstrain.org/ncov) si può vedere un'animazione della diffusione inferita del virus nel corso dell'epidemia.


<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
<!-- # [Advanced reading: uncertainty in trees](https://nextstrain.org/ncov/2020-03-11) -->
# [Approfondimenti: incertezze negli alberi](https://nextstrain.org/ncov/2020-03-11)

<!-- Earlier, we talked about how internal nodes represent _hypothesized_ unsampled cases. In fact, all trees represent _hypotheses_ about how a pathogen has evolved and moved over time. The trees we present on Nextstrain are point estimates -- that is, the version of this history that maximizes the probability of observing the data that we do. -->
Abbiamo detto prima che i nodi interni rappresentano _ipotetici_ casi non identificati. Di fatto, tutti i casi rappresentano delle _ipotesi_ su come il patogeno si sia evoluto e spostato nel tempo. Gli alberi che presentiamo su Nextstrain sono stime puntuali -- cioè le versioni di questa storia che massimizzano la probabilità di osservare quello che osserviamo in realtà.
<br><br>
<!-- However, there is always uncertainty in these estimates. Generally speaking, parts of the tree that are densely sampled are more certain; areas that are sparsely sampled are less certain. -->
Comunque c'è sempre dell'incertezza in queste stime. In generale, le parti dell'albero in cui abbiamo più campioni sono più affidabili; le parti che contengono meno campioni sono meno affidabili.

```auspiceMainDisplayMarkdown
<!-- # An illustration -->
Un'illustrazione
<div width="50%" margin="auto">
<p>
<!-- <img width="700px" alt="Illustration of the uncertainty inherent in tree reconstruction" src="https://github.com/nextstrain/nextstrain.org/raw/c69bfd0750c284ff12f33682f8d82848e13d9e15/static-site/content/help/01-general/figures/hcov_densitree.png"/> -->
<img width="700px" alt="Illustrazione dell'incertezza inerente alla ricostruzione dell'albero" src="https://github.com/nextstrain/nextstrain.org/raw/c69bfd0750c284ff12f33682f8d82848e13d9e15/static-site/content/help/01-general/figures/hcov_densitree.png"/>
</p>
</div>
```

<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
<!-- # [Scientific credit](https://nextstrain.org/ncov/2020-03-05?d=map&c=author) -->
# [Riconoscimenti scientifici](https://nextstrain.org/ncov/2020-03-05?d=map&c=author)

<!-- We would like to acknowledge the amazing and timely work done by all scientists involved in this outbreak, but particularly those working in China. -->
Ringraziamo tutti gli scienziati al lavoro su questa epidemia, e in particolare quelli che lavorano in Cina, per lo straordinario e tempestivo lavoro fatto.
<!-- Only through the rapid sharing of genomic data and metadata are analyses such as these possible. -->
È solo attraverso la rapida condivisione di dati e metadati genomici che analisi come questa sono possibili.
<br><br>

<!-- We also gratefully acknowledge [GISAID](https://gisaid.org) for providing the platform through which these data can be uploaded and shared. -->
Siamo anche grati a [GISAID](https://gisaid.org) per aver fornito la piattaforma su cui questi dati sono stati caricati e condivisi.

<!-- Do not need to translate insitutions names -->
<!-- This is right-side text -->
```auspiceMainDisplayMarkdown

<!--We are grateful for the data gathered by these originating labs:-->
Siamo grati per i dati raccolti da questi laboratori di origine:

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
<!-- # [Detailed scientific credit](https://nextstrain.org/ncov/2020-03-05?d=map&c=author) -->
# [Riconoscimenti scientifici dettagliati](https://nextstrain.org/ncov/2020-03-05?d=map&c=author)

<!-- These data were shared via [GISAID](https://gisaid.org). -->
Questi dati sono condivisi via [GISAID](https://gisaid.org).
<!-- We gratefully acknowledge their contributions. -->
Riconosciamo con gratitudine i loro contributi.
<br><br>

<!--  To the right we give specific sequences shared by each lab. -->
Sulla destra le sequenze specifiche condivise da ogni laboratorio.

<!-- This is right-side text -->
```auspiceMainDisplayMarkdown

<!-- The SARS-CoV-2 genomes were generously shared by scientists at these submitting labs: -->
I genomi del SARS-CoV-2 sono stati generosamente condivisi dagli scienziati presso questi laboratori:

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
