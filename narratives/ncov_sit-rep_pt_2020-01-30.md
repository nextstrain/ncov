---
title: Análise genômica do espalhamento do coronavírus nCoV. Relatório da situação até 30/01/2020.
authors: "Trevor Bedford, Richard Neher, James Hadfield, Emma Hodcroft, Misja Ilcisin, Nicola Müller"
authorLinks: "https://nextstrain.org"
affiliations: "Fred Hutch, Seattle, USA and Biozentrum, Basel, Switzerland"
date: "2020 Jan 30"
dataset: "https://nextstrain.org/ncov/2020-01-30?d=map"
abstract: "Esse relatório usa dados genômicos publicamente compartilhados sobre o novo coronavírus (nCov) da GISAID e Genbank, para estimar taxas e padrões de espalhamento de epidemia viral. Nós planejamos liberar relatórios atualizados assim que novos dados sejam produzidos e compartilhados. Esse site é otimizado para visualização em navegadores de computadores. Tradução por: Glaucio Santos & Anderson Brito"
---

# [Resumo](https://nextstrain.org/ncov/2020-01-30)

```auspiceMainDisplayMarkdown
## Resumo

Examinamos a diversidade genética dos 42</tag> genomas do novo coronavírus (nCoV) compartilhados publicamente, a fim de inferir a data do ancestral comum e taxas de espalhamento.
Concluímos que:
* Os 42</tag> genomas amostrados são muito similares, diferindo em apenas 0-5 mutações
* Essa falta de diversidade genética tem uma explicação cautelosa de que o surto ou descende de uma única introdução do vírus na população humana, ou descende de um pequeno número de transmissões de animais para humanos.
* Esse evento provavelmente ocorreu em Novembro ou no início de Dezembro de 2019.
* Desde então a transmissão de humanos para humanos tem ocorrido, resultando em casos observáveis da doença .
* Usando o total estimado pelo Imperial College London — que aponta milhares de casos —, inferimos um número reprodutivo básico entre 1,8 e 3,5, indicando rápido crescimento do surto viral no período entre Novembro e Janeiro.
```

# [Coronavírus](https://nextstrain.org/ncov/2020-01-30)

### Leituras adicionais:

* Informações gerais sobre os coronavírus na [Wikipédia](https://pt.wikipedia.org/wiki/Coronavirus) 30/01/2020
* Resumo do surto de nCoV na [Wikipédia](https://pt.wikipedia.org/wiki/Epidemia_de_pneumonia_por_novo_coronav%C3%ADrus_de_2019%E2%80%932020) 25/01/2020
* Material disponibilizado pelo [Centro de Controle e Prevenção de Doenças dos Estados Unidos (CDC)](https://www.cdc.gov/coronavirus/index.html) 29/01/2020
* Características do genoma em [ViralZone](https://viralzone.expasy.org/764?outline=all_by_species) 23/01/2020
* Análise interativa de risco do [MOBS-lab](https://datastudio.google.com/reporting/3ffd36c3-0272-4510-a140-39e288a9f15c/page/U5lCB) 29/01/2020
* Análise interativa de risco do [ROCS-lab](http://rocs.hu-berlin.de/corona/) 29/01/2020

```auspiceMainDisplayMarkdown

## Diferentes coronavírus humanos

Os coronavírus (CoV) são membros de um conjunto numeroso de especíes de vírus de RNA de sentido positivo ((+) ssRNA), com histórico de causar infecções respiratórias em humanos.
Alguns tipos de coronavírus estão associados com surtos, já outros estão continuamente circulando e causando principalmente infecções respiratórias leves (como por exemplo, o resfriado comum).



#### SARS-CoV e MERS-CoV
O mais conhecido dos coronavírus é o [SARS-CoV](https://pt.wikipedia.org/wiki/S%C3%ADndrome_respirat%C3%B3ria_aguda_grave) (Síndrome respiratória aguda grave), que foi responsável por um surto global entre Novembro de 2002 e Julho de 2003, que resultou em [mais de 8000 casos e 774 mortes](https://www.theguardian.com/world/2017/dec/10/sars-virus-bats-china-severe-acute-respiratory-syndrome), com uma taxa de letalidade entre 9 e 11%.

Em 2012, um novo coronavírus foi identificado, o [MERS-CoV](https://pt.wikipedia.org/wiki/Coronav%C3%ADrus_da_s%C3%ADndrome_respirat%C3%B3ria_do_Oriente_M%C3%A9dio) ("Síndrome respiratória do Oriente Médio"), o qual causa graves sintomas respiratórios. O MERS resultou em um número de mortes comparável ao SARS, porém sua via de transmissão é muito diferente. Enquanto o SARS é transmitido eficientemente de um humano para outro, infecções humanas com MERS são geralmente resultante de contato com animas, isto é, por camelos (veja [Dudas _et al._](https://elifesciences.org/articles/31257) para mais informações). Esse fator fez com que o surto de MERS se auto-restringisse a Península Arábica.




#### CoVs Sazonais
Contudo, nem todos os coronavírus são tão mortais quanto o SARS-CoV e MERS-CoV.
Há quatro coronavírus "sazonais" que  normalmente infectam humanos todo ano.
Comparado com o SARS, as cepas sazonais são ["muito mais prevalentes, muito menos graves, e geralmente causam doenças parecidas com a gripe (ILI, segundo a sigla em inglês)"](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5820427/).
Na verdade, [5](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2879166/)–[12%](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5820427/) de todos os casos de ILI testam positivo para coronavírus, portanto os sazonais são muito comuns, o que resulta em milhões de infecções leves todo ano.
Estes coronavírus sazonais são resultados de transmissões independentes de vírus provenientes de morcegos (reservatórios naturais de vírus) para humanos ao longo de aproximadamente 100 anos. Após tais transmissões entre espécies, cada vírus se estabeleceu e se espalhou amplamente na população humana.




#### Reservatórios animais
Coronavírus infectam uma ampla gama de animas, e os surtos virais em humanos descritos anteriormente são resultado de um ou mais "saltos" dos vírus desses reservatórios animais para populações humanas.
Acredita-se que o SARS chegou até a população humana por meio de [morcegos-de-ferradura via uma civetas (Paguma larvata) como intermediário](https://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1006698).


#### Transmissão de Pessoa para Pessoa
A capacidade de transmitir diferentes linhagens virais entre humanos é extremamente importante para entender o crescimento potencial de um surto.
Devido a sua habilidade de se transmitir entre humanos, acrescido de sua alta taxa de letalidade, o SARS (ou um vírus similar a ele) é considerado pela OMS uma [ameaça global de saúde pública](https://www.who.int/whr/2007/overview/en/index1.html).

```




# [Novo coronavírus (nCoV) 2019-2020](https://nextstrain.org/ncov/2020-01-30)

### Leituras Adicionais:

* Novo vírus na China: 5 Perguntas que cientistas estão fazendo (em inglês) [Nature news](https://www.nature.com/articles/d41586-020-00166-6) _2020-01-22_
* Última atualização do vírus da China: primeiro caso confirmado nos EUA (em inglês) [Nature news](https://www.nature.com/articles/d41586-020-00154-w) 21/01/2020
* Novo vírus surgindo na Ásia preocupa cientistas (em inglês) [Nature news](https://www.nature.com/articles/d41586-020-00129-x) 20/01/2020
* Novo vírus identificado como provável causa de doença misteriosa na China (em inglês) [Nature news](https://www.nature.com/articles/d41586-020-00020-9) 08/01/2020

```auspiceMainDisplayMarkdown




## Surto recente de um novo coronavírus
Em Dezembro de 2019, uma nova doença foi detectada pela primeira vez em Wuhan, China. Com sua identificação, ficou estabelecido que trata-se de mais um surto de coronavírus em humanos (o 7º), e o vírus passou a ser provisoriamente chamado de nCoV (novo coronavírus).

Até 30 de janeiro, mais de 7,914 casos e 170 mortes [foram reportadas](https://pt.wikipedia.org/wiki/Epidemia_de_pneumonia_por_novo_coronav%C3%ADrus_de_2019%E2%80%932020). Ainda é muito cedo para saber a taxa de letalidade, mas segundo indicações preliminares é menor que a do SARS-CoV. O número de casos estão aumentando dramaticamente, em partes graças ao aumento de vigilância epidemiológica e de testes.

Enquanto o surto parece ter seu epicentro em Wuhan, que atualmente está [sob quarentena](https://twitter.com/PDChina/status/1220060879112282117), o vírus tem se espalhado por toda a China, e tem sido detectado no exterior, como em Hong Kong, Singapura, Japão, e Tailândia, assim como na Europa, América do Norte, Ásia Meridional, Oriente Médio, e Austrália. Fora da China, a transmissão local (de pessoa para pessoa) foi reportada, mas é limitada até o momento.

A origem do vírus ainda não está clara, entretanto [análises genômicas](https://virological.org/t/ncovs-relationship-to-bat-coronaviruses-recombination-signals-no-snakes/331) sugerem que o nCoV é provavelmente relacionado a vírus anteriormente identificados em morcegos.
É plausível que a transmissão ocorreu tendo outros animais como intermediários antes de ser transmitido aos humanos. Porém, não há evidências de cobras como intermediário.




#### Narrativas da Nextstrain

As páginas seguintes contém análises realizadas usando o [Nextstrain](https://nextstrain.org). Rolando a página no menu a esquerda, textos correspondentes aos dados genômicos visualizados na direita surgirão na tela.

Ter acesso a longos genomas completos, de um novo vírus de RNA, e tão rapidamente, é uma conquista notável. Essas análises têm sido possíveis devido ao rápido compartilhamento público de dados genômicos e interpretações por cientistas de todo o mundo (veja a última página para conhecer os autores do sequenciamento).
```




# [Como interpretar as árvores filogenéticas?](https://nextstrain.org/ncov/2020-01-30)

### Leituras adicionais:

* [Explorando filogenias interativas com o Auspice](https://neherlab.org/201901_krisp_auspice.html) 24/01/2019

```auspiceMainDisplayMarkdown
## Árvores de transmissão versus árvores filogenéticas

Os patógenos se propagam através da replicação rápida em um hospedeiro, seguido da sua transmissão para outro hospedeiro. Uma epidemia só começa quando uma infecção resulta em mais de uma infecção subsequente.

Conforme o patógeno replica e se espalha, seu genoma precisa ser copiado muitas vezes e isto causa mutações aleatórias ("erros de cópia") que se acumularão no genoma. Tais mutações aleatórias podem ajudar a rastrear o espalhamento do patógeno, e aprender sobre suas vias e dinâmicas de transmissão.

<div>
  <img alt="cartoon showing how transmission tree and phylogenetic tree relate" width="500" src="https://neherlab.org/talk_images/infection_tree_combined.png"/>
</div>

A ilustração acima mostra um esboço de uma árvore de transmissão com um subconjunto de casos que foram amostrados (em azul). Na prática, a árvore de transmissão é desconhecida, e assim sendo, só temos disponíveis as estimativas do número total de casos.
Sequências genômicas nos permitem inferir partes da árvore de transmissão.
Nesse exemplo, três mutações (losangos pequenos) estão indicadas na árvore. As sequências que têm as mesmas mutações estão intimamente relacionadas, logo, essas mutações nos permitem agrupar amostras em grupos de vírus próximos, os quais pertencem às mesmas cadeias de transmissão.

### Lendo uma árvore filogenética

Abaixo, vemos uma ilustração com uma árvore filogenética a esquerda, onde mutações são mostradas como círculos coloridos. A direita estão as sequências correspondentes, também com mutações mostradas como círculos coloridos. Podemos ver que sequências que compartilham as mesmas mutações se agrupam juntas.
Quando sequências aparecem ligadas por uma linha vertical, como A e B, isso significa que não há diferenças entre eles: suas sequências são idênticas.

Quando uma sequência está sozinha em uma linha longa, como C ou E, isso indica a presença de mutações únicas, que não se encontram em outras sequências. Quanto mais longa a linha, mais mutações. A e B também tem mutações únicas (círculo verde) não compartilhadas por outras sequências, mas ambos vírus são idênticos.

<div>
  <img alt="cartoon of phylogenetic tree and corresponding alignment, with samples labelled A-E" width="500" src="https://nextstrain-data.s3.amazonaws.com/toy_alignment_tree.png"/>
</div>

No momento, a filogenia do novo coronavírus (nCov) pode não parecer muito como uma "árvore". Muitas das sequências são idênticas: aparecem juntas em linhas verticais, como A e B. Outras têm mutações únicas ou compartilhadas e, portanto, aparecem nas linhas ("ramificações") que vão para a direita. Você pode ver quantas  mutações uma ramificação tem passando seu mouse sobre ela (verifique na próxima página).
```




# [Análise filogenética](https://nextstrain.org/ncov/2020-01-30?m=div&d=tree)

Aqui apresentamos uma filogenia de 42</tag> cepas de nCoV que foram compartilhadas publicamente. Informações sobre como a análise foi realizada está disponível em [nosso repositório no GitHub](github.com/nextstrain/ncov).

<br>

As cores representam a região ou estado onde foram coletadas as amostras, e o eixo X representa a distância genética (medida em nucleotídeos).

<br>

A divergência é medida conforme o número de mudanças (mutações) no genoma.
Muitas sequências têm 0 mutações, o que significa que são todas idênticas até a raíz (base) da árvore. Outros vírus têm entre uma e sete mutações.

<br>

Sequenciar o genoma de um novo vírus de RNA longo durante um surto é desafiante. Algumas das diferenças observadas nessas sequências podem ser resultado de erros de sequenciamento, e não mutações reais. Inserções, deleções, e diferenças nos extremos do genoma têm grandes chances de serem erros, por isso as ignoramos visando cumprir os propósitos dessa análise.




# [Interpretação filogenética](https://nextstrain.org/ncov/2020-01-30?m=div&d=tree)

Atualmente vemos uma diversidade genética mínima nas sequências do nCoV: 11</tag> das 42</tag> sequências não têm mutações exclusivas.

<br>

A baixa diversidade genética nessas sequências sugerem que o ancestral comum de todas as sequências de nCoV é bastante recente, já que as mutações se acumulam lentamente comparado a outros vírus de RNA (com taxas em torno de 1-2 mutações por mês para coronavírus). Geralmente, quando as introduções por reservatório são múltiplas, é possível observar diversidade genética significativa (como tem sido com a febre de Lassa, o Ebola, a MERS-CoV e a gripe aviária).
O fato de podermos observar uma agrupamento tão forte em infecções humanas pode indicar que o surto descende de uma única introdução zoonótica (vinda de um animal) para a população humana, seguida do espalhamento epidêmico de pessoa para pessoa.

<br>

Nesse momento começamos a ver grupos de sequências que compartilham mutações. 
Um grupo inclui sequências de Guangdong, e quatro isolados dos Estados Unidos. 
Outro grupo contém de dois a quatro isolados. 
Sequências nestes grupos tendem a ser de amostras mais recentes, sugerindo que o vírus começou a acumular mutações a medida que ele se espalha in Wuhan e para outras cidades. 
Atualmente não há qualquer evidência de que tais mutações mudaram a forma como o vírus se comporta: é normal e esperado que vírus de RNA sofram mutações.




# [Transmissão dentro da família 1](https://nextstrain.org/2020-01-30/ncov?m=div&d=tree&f_location=Zhuhai)

Há três isolados de Zhuhai (sudeste da China, na província de Guangdong) que são geneticamente idênticos, e que formam um grupo, compartilhando uma mutação única que não é vista em outros isolamentos (passe o mouse sobre as ramificações para ver quais mutações estão presentes).

<br>

Dois desses casos (identificados por 028 and 040) [vêm de uma mesma família](https://twitter.com/JingLu_LuJing/status/1220143773532880896), indicando novamente que ocorreu uma transmissão de humano para humano. 
Não temos informações a cerca do terceiro caso.




# [Transmissão dentro da família 2](https://nextstrain.org/ncov/2020-01-30?m=div&d=tree&f_location=Shenzhen)

Dos seis (6) isolados da província de Guangdong (onde se localiza a cidade de Shenzhen), quatro (4) são geneticamente idênticos. 
Essas sequências compartilham as mesmas 3 mutações únicas. 

<br>

Três das sequências de Guangdong (identificadas por F025, F013, e F012) são [provenientes de uma mesma família](https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(20)30154-9/fulltext), e muito provavelmente representam novamente transmissão de humano para humano.



# [Transmissão dentro da família 2 - mutações compartilhadas](https://nextstrain.org/ncov/2020-01-30?m=div&d=tree&f_location=Shenzhen,Los%20Angeles,Orange%20County,Seattle,Chicago,Phoenix)

As três mutações encontradas nesse grupo estão também presentes no Arizona, isolado dos EUA, e duas das mutações são encontradas em outros três isolados dos EUA.




# [Transmissão dentro da família 3](https://nextstrain.org/ncov/2020-01-30?m=div&d=tree&f_location=Paris)

Por fim, as duas sequências da França são idênticas, compartilhando uma única mutação, e uma mutação também encontrada em um dos isolados dos EUA e de Taiwan.

<br>

As duas sequências francesas são [originários da mesma família](https://www.thelocal.fr/20200129/coronavirus-in-france-what-you-need-to-know) - um casal chinês de Wuhan.




# [Casos fora da China](https://nextstrain.org/ncov/2020-01-30?c=country&d=tree&m=div)

Existem casos diagnosticados e confirmados de nCoV em muitos países do Sudeste Asiático e Ásia Oriental, UEUA, Austrália, Oriente Médio, e Europa.
Vietnã, Japão, e Alemanha registraram transmissões locais dentro desses países, muito embora sempre com uma conexão com Wuhan, China.

<br>

As únicas sequências genômicas atualmente disponível de casos fora da China são de: dois casos da Tailândia; cinco dos EUA, dois da França, e um de Taiwan.
As sequências tailandesas são geneticamente idênticas a nove sequências chinesas, incluindo sete de Wuhan.

Quatro sequências dos EUA compartilham duas mutações com grupos de sequências de Shenzhen.
As sequências restantes dos EUA compartilham uma mutação com a sequência de Taiwan, e duas da França.

<br>

A explicação mais cautelosa sobre o padrão observado de mutações compartilhadas entre sequências dos EUA e de Shenzhen indica que o variante viral com duas mutações compartilhadas estava circulando em Wuhan, e foi independentemente exportado para Shenzhen, e introduzido múltiplas vezes nos EUA.
Não há evidências de outras conexões de casos nos EUA que não sejam ligadas a Wuhan.




# [Calculando a provável data do ancestral comum mais recente](https://nextstrain.org/ncov/2020-01-30?d=tree)
A grande similaridade dos genomas sugere que os vírus compartilham um ancestral comum recente (ou seja, eles descenderam do mesmo vírus ancestral recentemente). Caso contrário, esperaríamos um número maior de diferenças genéticas entre as amostras.

<br>

Pesquisas anteriores com coronavírus similares sugerem que esses vírus acumulam entre 1 e 3 mutações em seus genomas por mês (taxas de 3 &times; 10<sup>-4</sup> e 1 &times; 10<sup>-3</sup> substituições por sítio por ano).

<br>

A direita nós exploramos como diferentes proposições sobre as taxas de substituição (“mutações”) e a diversidade genética observada nos dá estimativas temporais sobre o surto viral.




```auspiceMainDisplayMarkdown
## Provável data do ancestral comum dos vírus causadores do surto

Com as sequência adicionais compartilhadas ao longo da última semana, a filogenia agora mostra distintos grupos, de forma que nossa análise prévia supondo uma topologia em forma de estrela não é mais apropriada.

De toda forma, mostramos abaixo o resultado de nossa análise baseada nos dados disponíveis até dia 25/01/2020, assumindo uma filogenia em forma de estrela, com uma distribuição Poisson de mutações ao longo do tempo, para estimar o tempo até o ancestral comum mais recente (conhecido em inglês pela sigla “tMRCA”) dos vírus sequenciados.

**Usando o conjunto de dados completo, a análise feita pelo nextstrain estima que o ancestral comum de nCoV provavelmente existiu entre o fim de Novembro e princípios de Dezembro de 2019. A taxa de substituição ainda se impõem como a principal fonte de incertezas nessa análise.**

<div>
  <img alt="graph of TMRCA estimates based on different mutation rates" width="500" src="https://nextstrain-data.s3.amazonaws.com/ncov_poisson-tmrca.png"/>
</div>

Existe um [caso confirmado em Wuhan com início de sintomas em 01/12/2019](https://twitter.com/trvrb/status/1220749265380593664), o qual colocaria um limite máximo para a data de ancestralidade comum.
O ancestral comum dos vírus sequenciados até hoje deve, no entanto, deve ter existido antes da data citada acima.

Modelagem mais detalhadas sobre o início do surto estão em andamento.
Apesar das incertezas consideráveis, nova melhor suposição continua sendo que a epidemia começou entre fins de Novembro e princípios de Dezembro.

```




# [Estimando a taxa de crescimento](https://nextstrain.org/ncov/2020-01-30?d=tree)

Uma importante medida durante o espalhamento de um patógeno é o número médio de novos casos secundários que cada pessoa infectada gera.

<br>

Esse número é conhecido como R0 ("R-zero").
A direita, apresentamos estimativas simples do R0.

```auspiceMainDisplayMarkdown
## Estimativas sobre taxas de crescimento epidêmico

Cientistas do Imperial College London utilizaram o número de casos observados fora da China para estimar o [número total de casos](https://www.imperial.ac.uk/mrc-global-infectious-disease-analysis/news--wuhan-coronavirus/), e sugeriram a ocorrência de pelo menos alguns milhares de casos até o dia 22/01/2020.

Com os adicionais casos importados, e o contínuado crescimento em número de casos na China, não atualmente estimamos pelo menos 5000 casos até hora. 

Juntamente com nossas estimativas prévias sobre a data de surgimento do surto viral, e informações sobre a duração do período de infecção, podemos estimar possíveis faixas de R0 usando um modelo matémático de processos de ramificação.

**Nossas estimativas apontam um R0 provável entre 1,8 e 3,5.**

Se assumimos que o surto começou no início de novembro de 2019 (12 semanas atrás), cremos que o R0 deveria oscilar entre 1,8 e 2,5, dependendo de quão grande ('n') o surto esteja agora.
<div>
  <img alt="Gráfico de estimativas de R0 com o início da epidemia 12 semanas atrás" width="500" src="https://nextstrain-data.s3.amazonaws.com/ncov_branching-R0-early_2020-01-29.png"/>
</div>

Se consideramos um início mais recente, no começo de dezembro de 2019 (8 semanas atrás), as estimativas são de que o R0 oscile entre 2,2 e 3,5:
<div>
  <img alt="Gráfico de estimativas de R0 com o início da epidemia 8 semanas atrás" width="500" src="https://nextstrain-data.s3.amazonaws.com/ncov_branching-R0-recent_2020-01-29.png"/>
</div>

Estas estimativas são de maneira geral consistentes com aquelas feitas por outros cientistas, os quais em grande parte mostram um R0 entre 2 e 3, veja por exemplo <a href="https://www.biorxiv.org/content/10.1101/2020.01.25.919787v1">esta pré-publicação</a>. 
Importante salientar: o R0 é uma medida que depende fortemente do contexto sócio-econômico, e de medidas de controle de infecções.
```




# [Crédito científico](https://nextstrain.org/ncov/2020-01-30?d=map&c=author)

Gostaríamos de reconhecer o trabalho incrível e oportuno realizado por todos os cientistas envolvidos nesse surto, sobretudo os que estão trabalhando na China.
Somente com o compartilhamento rápido de dados genômicos e de metadados que tais análises foram possíveis.

<br>

Os genomas do nCoV foram generosamente compartilhados por cientistas de várias instituições. São elas:

 * Shanghai Public Health Clinical Center & School of Public Health, Fudan University, Xangai, China
 * National Institute for Viral Disease Control and Prevention, China CDC, Pequim, China
 * Institute of Pathogen Biology, Chinese Academy of Medical Sciences & Peking Union Medical College, Pequim, China
 * Wuhan Institute of Virology, Chinese Academy of Sciences, Wuhan, China
 * Department of Microbiology, Zhejiang Provincial Center for Disease Control and Prevention, Hancheu, China
 * Guangdong Provincial Center for Diseases Control and Prevention
 * Department of Medical Sciences, National Institute of Health, Nonthaburi, Tailândia
 * Division of Viral Diseases, Centers for Disease Control and Prevention, Estados Unidos
 * Centers for Disease Control, R.O.C., Taipé, Taiwan
 * Institut Pasteur, Paris, França



# [Crédito cientifíco detalhado](https://nextstrain.org/ncov/2020-01-30?d=map&c=author)

Esses dados foram compartilhados através do [GISAID](https://gisaid.org).
Agradecemos as suas contribuições com profunda gratidão.

<br>

A direita, especificamos quais sequências foram compartilhadas por quais laboratórios.

```auspiceMainDisplayMarkdown

Os genomas do nCoV foram generosamente cedidos por cientistas de várias instituições. São elas:

 * Shanghai Public Health Clinical Center & School of Public Health, Fudan University, Xangai, China
   - Wuhan-Hu-1/2019
 * National Institute for Viral Disease Control and Prevention, China CDC, Pequim, China
   - Wuhan/IVDC-HB-01/2019
   - Wuhan/IVDC-HB-04/2020
   - Wuhan/IVDC-HB-05/2019)
 * Institute of Pathogen Biology, Chinese Academy of Medical Sciences & Peking Union Medical College, Pequim, China
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
 * Department of Microbiology, Zhejiang Provincial Center for Disease Control and Prevention, Hancheu, China
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
 * Department of Medical Sciences, National Institute of Health, Nonthaburi, Tailândia
   - Nonthaburi/61/2020
   - Nonthaburi/74/2020
 * Division of Viral Diseases, Centers for Disease Control and Prevention, Estados Unidos
   - USA-WA1/2020
   - USA/AZ1/2020
   - USA/IL1/2020
   - USA/CA1/2020
   - USA/CA2/2020
 * Centers for Disease Control, R.O.C., Taipé, Taiwan
   - Taiwan/2/2020
 * Institut Pasteur, Paris, França
   - France/IDF0372/2020
   - France/IDF0373/2020

```
