---
title: Analisi genomica della diffusione del COVID-19. Aggiornamento al 13 marzo 2020.
authors:
  - Emma Hodcroft
  - Nicola Müller
  - Cassia Wagner
  - Misja Ilcisin
  - James Hadfield
  - Sidney M. Bell
  - Richard Neher
  - Trevor Bedford
authorLinks:
  - https://neherlab.org/emma-hodcroft.html
  - https://bedford.io/team/nicola-mueller/
  - https://bedford.io/team/cassia-wagner/
  - https://bedford.io/team/misja-ilcisin/
  - https://bedford.io/team/james-hadfield/
  - https://twitter.com/sidneymbell
  - https://neherlab.org/richard-neher.html
  - https://bedford.io/team/trevor-bedford/
affiliations: "Fred Hutch, Seattle, USA; Biozentrum, Basel, Switzerland; CZI, CA, USA"
translators:
  - Andrea Claudi
  - Alice Ledda
  - Alessia Lepore
  - Sara Pilia
date: "2020 March 13"
dataset: "https://nextstrain.org/ncov/2020-03-13?d=map&legend=closed"
abstract: "Questo report usa dati genomici pubblicamente condivisi per ricostruire la diffusione del COVID-19. Questi report sono aggiornati settimanalmente"
---
<!-- Translators: Only text after : in the above ^ needs to be translated -->
<!-- Comment tags like these do not need to be translated, they are only to help you! -->
<!-- Ensure that links always end in a 'letter' (. counts) If some kind of text doesn't follow them, it breaks the slide. -->
<!-- numbers can be tagged ilke this: 161</tag> - this is just for us to help find them to update! Just leave in the </tag> bit. -->

<!-- This is left-side text -->
# [Indice dei contenuti](https://nextstrain.org/ncov/2020-03-13?d=tree,map&p=grid)
<!-- # [Table of Contents](https://nextstrain.org/ncov/2020-03-13?d=tree,map&p=grid)-->

* [Risorse di base](https://nextstrain.org/narratives/ncov/sit-rep/it/2020-03-13?n=2).
<!-- * [Background resources](https://nextstrain.org/narratives/ncov/sit-rep/2020-03-13?n=2). -->
* [Una nota sulla raccolta dei campioni](https://nextstrain.org/narratives/ncov/sit-rep/it/2020-03-13?n=3).
<!-- * [A note on sampling](https://nextstrain.org/narratives/ncov/sit-rep/2020-03-13?n=3). -->
* [Circolazione in Europa](https://nextstrain.org/narratives/ncov/sit-rep/it/2020-03-13?n=4).
<!-- * [Circulation in Europe](https://nextstrain.org/narratives/ncov/sit-rep/2020-03-13?n=4).  -->
* [Trasmissione locale in Regno Unito](https://nextstrain.org/narratives/ncov/sit-rep/it/2020-03-13?n=5).
<!-- * [Local transmission in the U.K.](https://nextstrain.org/narratives/ncov/sit-rep/2020-03-13?n=5). -->
* [Diffusione del virus SARS-CoV-2 dall'Iran](https://nextstrain.org/narratives/ncov/sit-rep/it/2020-03-13?n=6).
<!-- * [Spread of SARS-CoV-2 from Iran](https://nextstrain.org/narratives/ncov/sit-rep/2020-03-13?n=6). -->
* [Introduzioni negli Stati Uniti](https://nextstrain.org/narratives/ncov/sit-rep/it/2020-03-13?n=7).
<!-- *[Introductions to the U.S.](https://nextstrain.org/narratives/ncov/sit-rep/2020-03-13?n=7). -->
* [Diffusione del virus SARS-CoV-2 nello Stato di Washington](https://nextstrain.org/narratives/ncov/sit-rep/it/2020-03-13?n=8).
<!-- *[Spread of SARS-CoV-2 in Washington state](https://nextstrain.org/narratives/ncov/sit-rep/2020-03-13?n=8). -->
* [Diffusione del virus SARS-CoV-2 in California](https://nextstrain.org/narratives/ncov/sit-rep/it/2020-03-13?n=9).
<!-- * [Spread of SARS-CoV-2 in California](https://nextstrain.org/narratives/ncov/sit-rep/2020-03-13?n=9). -->
* [Cosa puoi fare](https://nextstrain.org/narratives/ncov/sit-rep/it/2020-03-13?n=10).
<!-- * [What you can do](https://nextstrain.org/narratives/ncov/sit-rep/2020-03-13?n=10). -->
* [FAQ & idee errate](https://nextstrain.org/narratives/ncov/sit-rep/it/2020-03-13?n=11).
<!-- * [FAQ & common misconceptions](https://nextstrain.org/narratives/ncov/sit-rep/2020-03-13?n=11). -->
* [Ringraziamenti scientifici](https://nextstrain.org/narratives/ncov/sit-rep/it/2020-03-13?n=12).
<!-- * [Scientific credit](https://nextstrain.org/narratives/ncov/sit-rep/2020-03-13?n=12). -->


<!-- This is right-side text -->
```auspiceMainDisplayMarkdown
<!-- # Executive summary -->
# Sommario

<!-- Here, we analyzed 410</tag> publicly shared COVID-19 genomes. By comparing these viral genomes to each other, we can characterize how COVID-19 is evolving and moving around the world. -->
Di seguito, abbiamo analizzato 410</tag> genomi del COVID-19 pubblicamente condivisi. Confrontando tra loro questi genomi virali possiamo determinare come il COVID-19 si sta evolvendo e spostando nel mondo.

<!-- For a current snapshot of the number of coronavirus cases around the world, see [Our World In Data](https://ourworldindata.org/coronavirus). -->
Per seguire in tempo reale il numero di casi di coronavirus nel mondo, vai a [Our World In Data](https://ourworldindata.org/coronavirus).

<!-- In this report, we show that the virus is widely circulating across the globe, with evidence of local transmission on multiple continents. -->
In questo report, osserviamo come il virus stia circolando in tutto il mondo con prove di trasmissione locale in più continenti.
<!-- At this time, we urge focus on efforts to slow the spread within communities; travel bans are less likely to be effective. -->
In questo momento, riteniamo necessario che gli sforzi siano concentrati sul rallentare la diffusione all'interno delle comunità; i divieti di viaggiare hanno meno probabilità di essere efficaci.

<!-- In this week's updates, we report: -->
Nell'aggiornamento di questa settimana:

<!-- * COVID-19 is circulating widely across Europe, with significant movement between countries. -->
* Il COVID-19 sta circolando in tutta Europa, spostandosi in modo significativo tra i diversi Stati.

<!-- * We identify at least 4 introductions to the UK, some with onward community transmission. -->
* Abbiamo identificato almeno 4 diverse introduzioni del virus nel Regno Unito, alcune delle quali seguite dalla trasmissione all'interno della comunità.

<!-- * There have been a number of travel-related cases linking Iran with other parts of the world. -->
* Ci sono stati un certo numero di casi collegati a viaggi effettuati tra l'Iran e altre parti del mondo.

<!-- * There have been many introductions into the U.S. to date, resulting in local transmission chains in multiple states. -->
* Ad oggi, ci sono state molte introduzioni del virus negli Stati Uniti che hanno causato focolai in più Stati.

<!-- * The outbreak continues to grow in Washington state; some cases are closely related to those from the Grand Princess cruise ship. -->
* L'epidemia continua a crescere nello Stato di Washington; alcuni di questi casi sono strettamenti connessi a quelli della nave da crociera Grand Princess.

<!-- * There is local circulation of COVID-19 in California. -->
* C'è una diffusione locale del COVID-19 in California.

<!-- * Social distancing measures should be enacted swiftly to ease the burden on healthcare systems and protect the vulnerable. -->
* Misure di limitazione delle interazioni sociali dovrebbero essere messe in atto velocemente per alleviare la pressione sui diversi sistemi sanitari e proteggere i soggetti più vulnerabili.

```


<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
<!-- # [COVID-19 Resources](https://nextstrain.org/ncov/2020-03-05) -->
# [Risorse sul COVID-19](https://nextstrain.org/ncov/2020-03-05)
<!-- We've prepared some resources that are worth reading to familiarize yourself with COVID-19 and the virus that causes it, SARS-CoV-2. -->
Di seguito elenchiamo alcune letture utili per familiarizzare con i concetti generali relativi al COVID-19 e al virus che ne è la causa, il SARS-CoV-2.

<!-- This information will make interpreting the data we present in this narrative easier; if you aren't familiar with phylogenetic trees, we encourage you to check out the ['How to Read Phylogenies' narrative](https://nextstrain.org/narratives/trees-background) and come back when you're ready. -->
Queste informazioni renderanno più semplice interpretare i dati che presenteremo più avanti; se non conosci gli alberi filogenetici, ti consigliamo di leggere ['Come leggere la filogenetica'](https://nextstrain.org/narratives/trees-background/it) e di continuare a leggere il resto più tardi.

<!-- This is right-side text -->
```auspiceMainDisplayMarkdown

<!-- ## Background -->
## Informazioni di base

<div>
 <!-- <a href="https://nextstrain.org/help/coronavirus/human-CoV"><img alt="microscopy image of coronaviruses" width="100" src="https://nextstrain.org/static/ncov_narrative-76cfd610d11ef708d213a3170de9519f.png"/> Background on Coronaviruses </a> -->
  <a href="https://nextstrain.org/help/coronavirus/human-CoV"><img alt="Immagine di un coronavirus al microscopio" width="100" src="https://nextstrain.org/static/ncov_narrative-76cfd610d11ef708d213a3170de9519f.png"/> Background on Coronaviruses </a>

 <!-- <a href="https://nextstrain.org/help/coronavirus/SARS-CoV-2"><img alt="illustration of a coronavirus" width="100" src="http://data.nextstrain.org/img_nCoV-CDC.jpg"/> Recent COVID-19 Outbreak Background </a> -->
  <a href="https://nextstrain.org/help/coronavirus/SARS-CoV-2"><img alt="illustrazione di un coronavirus" width="100" src="http://data.nextstrain.org/img_nCoV-CDC.jpg"/> Recent COVID-19 Outbreak Background </a>

 <!-- <a href="https://nextstrain.org/narratives/trees-background"><img alt="cartoon of a phylogenetic tree" width="100" src="http://data.nextstrain.org/img_toy_alignment_mini.png"/> How to Read Phylogenies</a> -->
  <a href="https://nextstrain.org/narratives/trees-background/it"><img alt="disegno di un albero filogenetico" width="100" src="http://data.nextstrain.org/img_toy_alignment_mini.png"/> Come leggere la filogenetica </a>
</div>

<!-- ## Further Reading -->
## Ulteriori letture

<!-- * Summary of the SARS-CoV-2 outbreak on [Wikipedia](https://en.wikipedia.org/wiki/2019%E2%80%9320_Wuhan_coronavirus_outbreak). -->
* Informazioni generali sull'epidemia di SARS-CoV-2 su [Wikipedia](https://it.wikipedia.org/wiki/Pandemia_di_COVID-19_del_2019-2020).

<!-- * Material provided by the [US CDC](https://www.cdc.gov/coronavirus/index.html). -->
* Materiale fornito dal Centro per il Controllo e la Prevenzione delle Malattie degli Stati Uniti [US CDC](https://www.cdc.gov/coronavirus/index.html) (in inglese).


<!-- ## Nextstrain narratives -->
## Analisi Nextstrain

<!-- The following pages contain analysis performed using [Nextstrain](https://nextstrain.org). -->
Le pagine che seguono contengono alcune analisi eseguite utilizzando [Nextstrain](https://nextstrain.org).
<!-- Scrolling through will reveal paragraphs of text with a corresponding visualization of the genomic data. -->
Scorrendo lungo la barra laterale sinistra sarà possibile leggere il testo relativo alla visualizzazione dei dati genomici.

<!-- To have full genomes of a novel and large RNA virus this quickly is a remarkable achievement. -->
Avere a disposizione così rapidamente i genomi completi per un virus a RNA nuovo e così grande è un risultato notevole.
<!-- These analyses have been made possible by the rapid and open sharing of genomic data and interpretations by scientists all around the world (see the final slide for a visualization of sequencing authorship). -->
Queste analisi sono state rese possibili grazie alla rapida e libera diffusione dei dati genomici e alla loro interpretazione da parte degli scienziati di tutto il mondo (nella pagina finale è possibile leggere i nomi degli autori del sequenziamento).
```

<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
<!-- # [A note about sampling](https://nextstrain.org/ncov/2020-03-13?c=country&r=country&d=map&p=grid&legend=closed) -->
# [Una nota su come sono state raccolte le sequenze](https://nextstrain.org/ncov/2020-03-13?c=country&r=country&d=map&p=grid&legend=closed)
<!-- We currently have sequences from samples taken in 30 countries across 5 continents. This is an incredible feat -- sequencing an unknown, large RNA virus in the midst of a pandemic is difficult, and is only possible through the incredible work and timely sharing of data by scientists and physicians around the world. -->
Al momento abbiamo sequenze da campioni raccolti in 30 Paesi e 5 continenti. Questa è un'impresa incredibile -- sequenziare un virus a RNA così grande e sconosciuto nel bel mezzo di una pandemia è difficile ed è stato possibile solamente grazie all'incredibile sforzo e la tempestiva condivisione dei dati da parte di scienziati e medici in tutto il mondo.
<br><br>
<!-- While this data enables us to infer many useful characteristics of the outbreak and track its spread in real time, it's important to emphasize that our conclusions are limited by the available data. -->
Anche se questi dati ci permettono di dedurre molte utili caratteristiche dell'epidemia e di tracciare la sua evoluzione in tempo reale, è però importante enfatizzare che le nostre conclusioni sono limitate dalla disponibilità di dati.
<br><br>
<!-- For example, the map shows very few sequences from the global south. This is NOT because COVID-19 isn't circulating in these areas, or that these cases are not as crucial to understand; rather, we just don't have much data available from these areas. The size of each circle on the map indicates how much data is currently available from that area, rather than the true size of the outbreak. -->
Per esempio, la cartina mostra molte poche sequenze dal sud del mondo. Questo NON è dovuto al fatto che il COVID-19 non stia circolando in queste aree, o che questi casi non siano cruciali per la comprensione; in realtà, ciò è dovuto al fatto che non abbiamo abbastanza dati disponibili da queste aree. La dimensione di ogni cerchio sulla cartina indica quanti dati sono disponibili da una certa area, e non la reale dimensione del focolaio.
<!-- There is NO right-side text -->


<!-- ############ SLIDE BREAK ############# -->
<!-- This is left-side text -->
<!-- # [Circulation across Europe](https://nextstrain.org/ncov/2020-03-13?c=country&legend=closed&f_country=Belgium,France,Germany,Ireland,Italy,Netherlands,Portugal,Spain,Sweden,Switzerland,United%20Kingdom,Ireland&label=clade:A2&m=div&d=map,tree&p=grid)-->
# [Circolazione del virus in Europa](https://nextstrain.org/ncov/2020-03-13?c=country&legend=closed&f_country=Belgium,France,Germany,Ireland,Italy,Netherlands,Portugal,Spain,Sweden,Switzerland,United%20Kingdom,Ireland&label=clade:A2&m=div&d=map,tree&p=grid)

<!-- Here, we see a large clade of sequences from Europe. -->
Qui osserviamo un grande ramo che contiene sequenze provenienti dall'Europa.
<!-- Notably, sequences from many different countries intercalate, indicating that COVID-19 is already circulating quite widely across Europe. -->
È interessante notare che le sequenze di molti Paesi diversi risultano mescolate, e ciò indica che il COVID-19 sta già circolando piuttosto ampiamente in tutta Europa.
<br><br>

<!-- Zooming in on the map, we see that there are many links between Italy and other areas; however, it is important to keep in mind that the directionality of these links can't always be confidently inferred. Other hypotheses can also explain these data (e.g., if an unsampled case infected both a secondary case sequenced in Italy and a secondary case sequenced elsewhere). -->
Ingrandendo la cartina, si può osservare che ci sono molte connessioni tra l'Italia e altre aree; tuttavia, è importante ricordare che la direzionalità di queste connessioni non può essere sempre dedotta in modo sicuro. Ci sono anche altre ipotesi che possono spiegare questi dati (ad esempio, se un caso non campionato ha infettato sia un caso secondario sequenziato in Italia che un caso secondario sequenziato altrove).
<!-- There is no right side text -->


<!-- # [Local transmission in the British Isles & Ireland](https://nextstrain.org/ncov/2020-03-13?c=country&legend=closed&d=tree&f_country=United%20Kingdom,Ireland&p=full) -->
# [Trasmissione locale nelle Isole Britanniche e in Irlanda](https://nextstrain.org/ncov/2020-03-13?c=country&legend=closed&d=tree&f_country=United%20Kingdom,Ireland&p=full)

<!-- Looking at the British Isles & Ireland as an example, we can see several instances where viruses that are closely related to samples from other countries appear in the British Isles & Ireland.-->
Prendendo le Isole Britanniche e l'Irlanda come esempio, osserviamo diversi casi in cui in questi luoghi compaiono virus molto simili a quelli presenti in altri Paesi.
<br><br>

<!-- This is consistent with 4 or more introductions from other locations. -->
Questo è coerente con l'ipotesi che ci siano state 4 o più introduzioni del virus da altre località.
<br><br>

<!-- We also see instances where after an introduction, there are several closely-related cases from the same location. This is consistent with local community transmission from more than one of these introductions. -->
Abbiamo osservato anche casi in cui, dopo un ingresso del virus, ci sono diversi casi strettamente connessi tra loro nella stessa località. Questo è coerente con l'ipotesi di una trasmissione del virus a livello locale derivata da uno o più di questi ingressi.

<!-- There is no right side text -->

<!-- ############ SLIDE BREAK ############# -->
<!-- This is left-side text -->
<!-- # [Spread of SARS-CoV-2 from Iran](https://nextstrain.org/ncov/2020-03-13?d=tree,map&label=clade:A3&p=grid&legend=closed&m=div) -->
# [Diffusione del virus SARS-CoV-2 dall'Iran](https://nextstrain.org/ncov/2020-03-13?d=tree,map&label=clade:A3&p=grid&legend=closed&m=div)

<!-- A number of genomes have been sequenced from patients reporting travel history to Iran. These genomes are all extremely similar, and indicate that the outbreak in Iran may be the result of a single transmission which has subsequently been transmitted to many other places. -->
Numerosi genomi sono stati sequenziati da pazienti che hanno dichiarato di aver viaggiato in Iran. Questi genomi sono tutti estremamente simili tra loro ed indicano che l'epidemia in Iran potrebbe essere il risultato di una singola trasmissione che è stata successivamente diffusa a molte altre località.
<br><br>

<!-- Note that there are no full genomes available from patients in Iran. -->
Si noti che non ci sono genomi completi disponibili da pazienti in Iran.

<!-- There is NO right-side text -->

<!-- ############ SLIDE BREAK ############# -->
<!-- This is left-side text -->
<!--  # [Introductions to the U.S.](https://nextstrain.org/ncov/2020-03-13?d=tree,map&f_country=USA&m=div&p=full&legend=closed) -->
# [Introduzione negli Stati Uniti](https://nextstrain.org/ncov/2020-03-13?d=tree,map&f_country=USA&m=div&p=full&legend=closed)

<!-- Here, we can see that the virus has been introduced to the U.S. on multiple independent occasions. -->
Qui possiamo osservare come il virus sia stato introdotto negli Stati Uniti in più occasioni indipendenti tra loro.
<br><br>

<!-- Most of these introductions aren't associated with any other sampled cases from the U.S., so we're not sure if these introductions led to local outbreaks. -->
La maggior parte di queste introduzioni non sono associate con nessuno dei campioni raccolti negli Stati Uniti, perciò non siamo sicuri se queste introduzioni abbiano causato focolai locali.
<!-- However, given that testing capacity is not yet ramped up in most areas, we expect there are many unreported cases. -->
Tuttavia, dato che il numero di test effettuato in molte aree non è ancora alto, prevediamo l'esistenza di molti casi che non sono stati segnalati.
<br><br>

<!-- For Washington and California, though, we do see clusters of cases that are closely related. -->
Nello Stato di Washington e in California, tuttavia, vediamo gruppi di casi geneticamente molto simili.
<!-- This suggests ongoing transmission and local spread within these two states. -->
Ciò suggerisce che ci siano dei focolai attivi e stia avvenendo diffusione locale del virus all’interno di questi due Stati.
<!-- There is no right side text -->

<!-- ############ SLIDE BREAK ############# -->
<!-- This is left-side text -->
<!-- # [Spread of SARS-CoV-2 in Washington state](https://nextstrain.org/ncov/2020-03-13?c=division&r=division&d=tree,map&f_country=USA&label=clade:B1&m=div&p=grid&legend=closed) -->
# [Diffusione del virus SARS-CoV-2 nello Stato di Washington](https://nextstrain.org/ncov/2020-03-13?c=division&r=division&d=tree,map&f_country=USA&label=clade:B1&m=div&p=grid&legend=closed)

<!-- Here, we see a large cluster of cases from Washington that are all closely related. -->
Qui, osserviamo un grande gruppo di casi da Washington che sono tutti molto simili geneticamente.
<!-- From this, we conclude that there is extensive local spread within Washington state. -->
Da ciò, concludiamo che esiste un grosso focolaio locale del virus nello Stato di Washington.
<br><br>

<!-- Interestingly, the Washington samples intercalate with samples from the Grand Princess cruise ship. -->
È interessante notare che i campioni di Washington si mescolano nell'albero con quelli della nave da crociera Grand Princess.
<!-- We aren't sure yet whether the virus spread from the cruise ship to Washington or the other way around; as we get more data, we'll update our analysis. -->
Non è ancora chiaro se il virus si sia diffuso dalla nave da crociera a Washington o viceversa; man mano che otteremo più dati, aggiorneremo la nostra analisi.

<!-- There is NO right-side text -->

<!-- This is left-side text -->
<!-- # [Spread of SARS-CoV-2 in California](https://nextstrain.org/ncov/2020-03-13?c=country&r=division&d=tree,map&f_division=California&m=div&p=grid&legend=closed) -->
# [Diffusione del virus SARS-CoV-2 in California](https://nextstrain.org/ncov/2020-03-13?c=country&r=division&d=tree,map&f_division=California&m=div&p=grid&legend=closed)

<!-- Looking at samples from California, we see evidence for multiple introductions. -->
Osservando i campioni provenienti dalla California, troviamo prove di introduzioni multiple.
<!-- More importantly, we see at least one cluster of closely related cases, all sampled in California over a short time period (click on ['Explore the Data'](https://nextstrain.org/ncov) and search for 'CA9' to see on example). -->
Ancora più importante, vediamo almeno un gruppo di casi molto simili geneticamente in campioni rilevati in California durante un breve periodo di tempo (clicca su ['Esplora i dati'](https://nextstrain.org/ncov) e cerca 'CA9' per vedere un esempio).
<br><br>

<!-- This strongly suggests that there is ongoing local transmission within California. -->
Questo suggerisce fortemente che in California sia attivo un focolaio locale del virus.
<!-- There is NO right-side text -->

<!-- ############ SLIDE BREAK ############# -->
<!-- This is left-side text -->
<!-- # [Takeaways](https://nextstrain.org/ncov/2020-03-13?c=country&d=map&p=full) -->
# [Cose utili da sapere](https://nextstrain.org/ncov/2020-03-13?c=country&d=map&p=full)
<!-- - The virus has been introduced to many parts of the globe multiple times. Not all introductions result in local transmission. -->
- Il virus è stato introdotto in molte parti del mondo diverse volte. Non tutte le volte le introduzioni del virus hanno poi causato focolai locali.
<br><br>

<!-- - We see evidence of local transmission across Europe, parts of the United States, China, and Southeast Asia. -->
- Abbiamo prove di trasmissione locale del virus in Europa, parte degli Stati Uniti, Cina e del Sud-Est dell'Asia.
<br><br>

<!-- - Controlling local outbreaks through social distancing is crucial to protect the vulnerable. -->
- Controllare le epidemie locali attraverso la limitazione delle interazioni sociali è fondamentale per proteggere le categorie più vulnerabili.


<!-- This is the right-side text -->

```auspiceMainDisplayMarkdown
<!-- # What you can do -->
# Cosa puoi fare

<!-- Social distancing -- that is, decreasing the number of people you encounter each day -- can be challenging, but is hugely beneficial to the public good.  -->
Il distanziamento sociale -- cioè diminuire il numero di persone che si incontrano ogni giorno -- può essere difficile, ma è estremamente vantaggioso per il bene colletivo.

<!-- If everyone decreased their daily contacts by 25%, we would expect to see a 50% decrease in the cumulative number of cases over the next month ([Klein et al., 2020-03-13](https://institutefordiseasemodeling.github.io/COVID-public/reports/Working%20paper%20%E2%80%93%20model-based%20estimates%20of%20COVID-19%20burden%20in%20King%20and%20Snohomish%20counties%20through%20April%207.pdf)). Not sure what social distancing means? [Check out this helpful guide](https://www.theatlantic.com/family/archive/2020/03/coronavirus-what-does-social-distancing-mean/607927/). -->
Se tutti diminuissero del 25% i loro contatti giornalieri con altre persone, ci potremmo aspettare una diminuizione del 50% nel numero dei casi totali nel mese successivo ([Klein et al., 2020-03-13](https://institutefordiseasemodeling.github.io/COVID-public/reports/Working%20paper%20%E2%80%93%20model-based%20estimates%20of%20COVID-19%20burden%20in%20King%20and%20Snohomish%20counties%20through%20April%207.pdf)). Non sei sicuro di cosa significhi distanziamento sociale? [Leggi questa utile guida (in inglese)](https://www.theatlantic.com/family/archive/2020/03/coronavirus-what-does-social-distancing-mean/607927/).
<div>
  <img src="https://github.com/nextstrain/ncov/raw/master/figures/social-distancing-efficacy.png" width="70%">
</div>


<!-- ## Steps individuals can take -->
## Cosa si puo fare individualmente

<!-- * Reduce the number of people you are in contact with each day, especially if you are in a vulnerable group (e.g., seniors and those with pre-existing conditions). -->
* Ridurre il numero di persone con cui si è in contatto ogni giorno, in particolare se fai parte di una categoria vunerabile (ad esempio: anziani e persone con patologie mediche pre-esistenti).

<!-- * Remember that even if you are not super vulnerable, many people around you are; follow these practices to protect others. -->
* Ricorda che anche se non fai parte delle categorie vulnerabili, molte persone attorno a te lo sono; segui queste regole per proteggere gli altri.

<!-- * Wash your hands "like you just chopped a jalapeno and have to change a contact lens".  -->
* Lava le mani come se avessi appena tagliato un peperoncino e dovessi cambiarti le lenti a contatto.

<!-- * Stay home if you are sick; be prepared with a few extra supplies in case you need to self-quarantine.  -->
* Se stai male, stai a casa; prepara qualche provvista extra nel caso avessi bisogno di metterti in auto-quarantena.

<!-- * If you are an employer, encourage your employees to stay home when sick (and financially support them to do so).  -->
* Se sei un datore di lavoro, incoraggia i tuoi impiegati a stare a casa quando sono malati (e fornisci loro il supporto economico per farlo).

<!-- ## Steps officials can take  -->
## Quello che i governi possono fare

<!--  * Make testing free and broadly available. -->
* Far sì che i test siano gratuiti e ampiamente disponibili.

<!-- * Put social distancing measures in place. -->
* Implementare misure di distanziamento sociale.

<!-- * Financially support those impacted by social distancing measures (e.g., hourly workers, those with elder or childcare responsibilities, small businesses, etc.). -->
* Garantire supporto finanziario a chi è affetto dalla misure di distanziamento sociale (ad esempio: lavoratori ad ore, persone con anziani o bambini a carico, piccole attività commerciali, etc. )

```



<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
<!-- # [FAQs & Misconceptions](https://nextstrain.org/ncov/2020-03-05) -->
# [Domande frequenti e convinzioni errate](https://nextstrain.org/ncov/2020-03-05)

<!-- #### We know that a lot of people have questions about COVID-19. -->
#### Sappiamo che molte persone hanno domande sul COVID-19.

<!-- #### [We've set up a guide to try and answer the most frequently asked questions](https://nextstrain.org/help/coronavirus/FAQ). -->
#### [Per questo abbiamo creato una guida per provare a rispondere ad alcune delle domande più frequenti](https://nextstrain.org/help/coronavirus/FAQ).

<!-- #### The Federation of American Scientists also maintains [a great resource for FAQs](https://covid19.fas.org/l/en). -->
#### Anche la Federazione degli Scienziati Americani mantiene [un'ottima risorsa per delle FAQ](https://covid19.fas.org/l/en).

<!-- This is right-side text -->
```auspiceMainDisplayMarkdown
# Ulteriori letture (in inglese)
<!-- # Further reading -->

* "Don't believe the conspiracy theories you hear about coronavirus & HIV" [article](https://massivesci.com/notes/wuhan-coronavirus-ncov-sars-mers-hiv-human-immunodeficiency-virus/) _2020-01-31_

* "Baseless Conspiracy Theories Claim New Coronavirus Was Bioengineered" [article](https://www.factcheck.org/2020/02/baseless-conspiracy-theories-claim-new-coronavirus-was-bioengineered/) _2020-02-07_

* "No, The Wuhan Coronavirus Was Not Genetically Engineered To Put Pieces Of HIV In It" [article](https://www.forbes.com/sites/victoriaforster/2020/02/02/no-coronavirus-was-not-bioengineered-to-put-pieces-of-hiv-in-it/#5d339e8e56cb) _2020-02-02_

* "Busting coronavirus myths" [AFP Fact Check](https://factcheck.afp.com/busting-coronavirus-myths) _2020-02-19_

<!-- # Misconceptions -->
# Convinzioni errate

<!-- A number of misconceptions have been circulated about the origins of the novel coronavirus. -->
Stanno circolando diverse idee errate sulle origini del nuovo coronavirus.
<!-- During outbreaks like this one, the spread of information that's known to be incorrect can lead to more panic, and cause people not to trust scientists and governments, meaning they are less likely to follow advisories and take appropriate precautions. -->
Durante epidemie come questa, la diffusione incontrollata di informazioni note per essere errate può generare panico e diminuire la fiducia della popolazione negli scienziati e nel governo, e per questo motivo ad essere meno disposti a seguire i consigli e a prendere le precauzioni appropriate per evitare il contagio.

<!-- In an effort to try and explain why these views are incorrect, scientists have addressed these theories at the pages below: -->
Nel tentativo di spiegare perché queste convinzioni non sono corrette, gli scienziati hanno affrontato queste teorie nelle pagine che seguono (in inglese):

<div>

  <a href="http://virological.org/t/ncovs-relationship-to-bat-coronaviruses-recombination-signals-no-snakes-no-evidence-the-2019-ncov-lineage-is-recombinant/331"><img alt="picture of a snake" width="100" src="http://data.nextstrain.org/img_snake-freeToUse.jpg"/> 'Snake' Origins of SARS-CoV-2 (Technical) </a>
  <a href="https://twitter.com/trvrb/status/1223666856923291648"><img alt="illustration of HIV" width="100" src="http://data.nextstrain.org/img_HIV-wiki.jpg"/> 'HIV Engineering' Idea (Twitter thread)</a>


</div>


```


<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
<!-- # [Scientific credit](https://nextstrain.org/ncov/2020-03-05?d=map&c=author) -->
# [Ringraziamenti scientifici](https://nextstrain.org/ncov/2020-03-05?d=map&c=author)

<!-- We would like to acknowledge the amazing and timely work done by all scientists involved in this outbreak, and particularly those working in China. -->
Ringraziamo tutti gli scienziati al lavoro su questa epidemia, e in particolare quelli che lavorano in Cina, per lo straordinario e tempestivo lavoro fatto.
<!-- Only through the rapid sharing of genomic data and metadata are analyses such as these possible. -->
È solo attraverso la rapida condivisione di dati e metadati genomici che analisi come questa sono possibili.

<br>

<!-- We also gratefully acknowledge [GISAID](https://gisaid.org) for providing the platform through which these data can be uploaded and shared. -->
Siamo anche grati a [GISAID](https://gisaid.org) per aver fornito la piattaforma su cui questi dati sono stati caricati e condivisi.

<!-- Do not need to translate institutions names -->
<!-- This is right-side text -->
```auspiceMainDisplayMarkdown

<!-- We are grateful for the data gathered by these originating labs: -->
Ringraziamo i seguenti laboratori per i dati raccolti:

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
# [Ringraziamenti scientifici dettagliati](https://nextstrain.org/ncov/2020-03-05?d=map&c=author)


<!-- These data were shared via [GISAID](https://gisaid.org). -->
Quest dati sono stati condivisi via [GISAID](https://gisaid.org).
<!-- We gratefully acknowledge their contributions. -->
Riconosciamo con gratitudine i loro contributi.

<br>

<!--  To the right we give specific sequences shared by each lab. -->
Sulla destra si possono vedere le sequenze specifiche condivise da ogni laboratorio.

<!-- This is right-side text -->
```auspiceMainDisplayMarkdown

<!-- The SARS-CoV-2 genomes were generously shared by scientists at these submitting labs: -->
I genomi del SARS-CoV-2 sono stati generosamente condivisi dagli scienziati di questi laboratori:

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
