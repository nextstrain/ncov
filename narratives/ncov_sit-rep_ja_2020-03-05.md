---
title: 新型コロナウイルス拡散の遺伝的解析. 状況報告 2020-01-23.
authors: "Trevor Bedford, Richard Neher, James Hadfield, Emma Hodcroft, Misja Ilcisin, Nicola Müller, Takeshi Sato, Fengjun Zhang"
authorLinks: "https://nextstrain.org"
affiliations: "Fred Hutch, Seattle, USA and Biozentrum, Basel, Switzerland"
date: "2020 March 05"
dataset: "https://nextstrain.org/ncov/2020-03-05"
abstract: "このレポートは、GISAIDとGenbankが公的に共有している新型コロナウイルス(nCOV)の遺伝情報を用い、ウイルスの疫学的拡散の速度とパターンを評価しています。データが新たに取得、共有され次第状況報告を更新する予定です。このウェブサイトはデスクトップブラウザ上での表示に最適化されています."
---
<!-- Translators: Only text after : in the above ^ needs to be translated -->
<!-- Comment tags like these do not need to be translated, they are only to help you! -->
<!-- Ensure that links always end in a 'letter' (. counts) If some kind of text doesn't follow them, it breaks the slide. -->
<!-- numbers can be tagged ilke this: 161</tag> - this is just for us to help find them to update! Just leave in the </tag> bit. -->

<!-- This is left-side text -->
# [エグゼクティブサマリ](https://nextstrain.org/ncov/2020-03-05)

<!-- This is right-side text -->
```auspiceMainDisplayMarkdown
## エグゼクティブサマリ

公的にシェアされたCOVID-19の 169</tag> のゲノム情報を用いて 、私たちは様々な地域でのCOVID-19の拡散を特徴付け、共通祖先の出現日を推測するために、遺伝的多様性を調査しました.



私たちが発見したこととしては、
* COVID-19はイタリアに少なくとも2度持ち込まれ、それに続きコミュニティで拡散した([ページを見る](https://nextstrain.org/narratives/ncov/sit-rep/ja/2020-03-05?n=7)).
* イタリアから持ち込まれた症例が出現したとされる6つの異なる国由来のシーケンスのクラスターがこれに含まれる([これを見る](https://nextstrain.org/narratives/ncov/sit-rep/ja/2020-03-05?n=9))
* 遺伝子シーケンスのデータは1月中頃からgreater SeattleエリアでCOVID-19の未検出でかつ持続的に拡散しているという仮設を支持する([これを見る](https://nextstrain.org/narratives/ncov/sit-rep/ja/2020-03-05?n=10))
* この解析で用いられたものも含み、シーケンスが行われたすべての症例は、2019年11月中頃から12月中頃までの期間に共通祖先を持つ可能性がある.([これを見る](https://nextstrain.org/narratives/ncov/sit-rep/ja/2020-03-05?n=11))

```


<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [コロナウイルス](https://nextstrain.org/ncov/2020-03-05)

### 参考文献:

* [Wikipedia](https://en.wikipedia.org/wiki/2019%E2%80%9320_Wuhan_coronavirus_outbreak) _2020-01-30_ のコロナウイルス感染症流行の要旨
* [米CDC](https://www.cdc.gov/coronavirus/index.html) _2020-01-29_ が提供する資料

<!-- This is right-side text -->
```auspiceMainDisplayMarkdown

## COVID-19 資料

下記に、COVID-19とその症状を引き起こすウイルスであるSARS-CoV-2の詳細を知るために一読する値する資料を用意しました.
この情報によって、このNarrativeの中で私たちが提示するデータを解釈するのがより容易になります.

<div>
  <a href="https://nextstrain.org/help/coronavirus/human-CoV"><img alt="microscopy image of coronaviruses" width="100" src="https://nextstrain.org/static/ncov_narrative-76cfd610d11ef708d213a3170de9519f.png"/> コロナウイルスの背景 </a>

  <a href="https://nextstrain.org/help/coronavirus/SARS-CoV-2"><img alt="illustration of a coronavirus" width="100" src="http://data.nextstrain.org/img_nCoV-CDC.jpg"/> 最近のCOVID-19流行の背景 </a>

  <a href="https://nextstrain.org/help/general/how-to-read-a-tree"><img alt="cartoon of a phylogenetic tree" width="100" src="http://data.nextstrain.org/img_toy_alignment_mini.png"/> 系統発生の読み方 </a>

</div>

## Nextstrain Narrative

このあとのページが[Nextstrain](https://nextstrain.org) を使って行われた解析です
左側のサイドバーをスクロールしていくと、右側に対応するゲノムデータの可視化とともに段落としてテキストが表示されます.

新型でかつ大型のRNAウイルスを早い段階でこのような形にするのは優れた業績です.
このような解析は素早くかつオープンなゲノムデータの共有と世界中の科学者による解釈によって可能になりました(シーケンシングのオーサシップについては最後のスライドを参照してください).


```


<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [よくある質問と誤解](https://nextstrain.org/ncov/2020-03-05)

### 参考文献:

* "Baseless Conspiracy Theories Claim New Coronavirus Was Bioengineered" [記事](https://www.factcheck.org/2020/02/baseless-conspiracy-theories-claim-new-coronavirus-was-bioengineered/) _2020-02-07_

* "No, The Wuhan Coronavirus Was Not Genetically Engineered To Put Pieces Of HIV In It" [記事](https://www.forbes.com/sites/victoriaforster/2020/02/02/no-coronavirus-was-not-bioengineered-to-put-pieces-of-hiv-in-it/#5d339e8e56cb) _2020-02-02_

* "Busting coronavirus myths" [AFP ファクトチェック](https://factcheck.afp.com/busting-coronavirus-myths) _2020-02-19_

<!-- This is right-side text -->
```auspiceMainDisplayMarkdown

## FAQと誤解

### FAQ

多くの方々がCOVID-19について疑問を持っているということは承知しています。
よくある質問について回答を試みるために[こちら](https://nextstrain.org/help/coronavirus/FAQ)にガイドを用意しました

<div>

  <a href="https://nextstrain.org/help/coronavirus/FAQ"><img alt="picture of a question mark" width="100" src="http://data.nextstrain.org/img_question-mark.jpg"/> COVID-19 FAQ </a>

</div>


### 誤解

新型コロナウイルスの由来については多くの誤解が流布されているます。
このような大きな感染症流行の間は、不適切であるとわかる情報の拡散はより大きなパニックを引き起こし、人々に科学者と政府を不信を抱かせることになり、これは勧告に従ったり、適切な警戒をしなくなるようになることを意味します.

次のような見方がなぜ適切でないかを説明しようと、科学者は下記のページで理論について言及しました：

<div>

  <a href="http://virological.org/t/ncovs-relationship-to-bat-coronaviruses-recombination-signals-no-snakes-no-evidence-the-2019-ncov-lineage-is-recombinant/331"><img alt="picture of a snake" width="100" src="http://data.nextstrain.org/img_snake-freeToUse.jpg"/> SARS-CoV-2の起源 "蛇"  (専門的) </a>
  <a href="https://twitter.com/trvrb/status/1223666856923291648"><img alt="illustration of HIV" width="100" src="http://data.nextstrain.org/img_HIV-wiki.jpg"/> "HIVエンジニアリング"という考え(Twitter)</a>


</div>


```


<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [系統発生解析](https://nextstrain.org/ncov/2020-03-05?d=tree)


ここに示すのが、SARS-CoV-2(COVID-19を引き起こすウイルス)の169</tag>系統の発生解析です。
解析方法の情報にについては[このGitHubリポジトリ](https://github.com/nextstrain/ncov)で閲覧可能です.

<br>

色は単離された国内での地域あるいはアメリカの州を表し、x軸は取得された日付を表します。
y軸はどのようにシーケンスがつながるかを表し、計測における単位はありません。

<br>

サンプル取得日は有益ですが、2つのシーケンスが遺伝的に関連敷いているかを確かに示しているとは限りません - 2つの遺伝シーケンスはサンプル取得日が違ったとしても同一ですし、この表示方法ではかなり異なって見えてしまっています。

<!-- There is NO right-side text -->


<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [系統解析における "多様性"](https://nextstrain.org/ncov/2020-03-05?d=tree&m=div)

サンプル取得日の代わりにx軸に"多様性"を表すように表示を変えました

<br>

さきほどの表示とは違った見え方をするシーケンスの数の多さに注意してください、この画面では垂直方向に一列に並んでいます。
前のスライドとこのスライドの間をスクロールで行き来して、系統樹の変化を見ることができます。

<br>

多様性はゲノム中の変化(変異)の数で計測しています。
変異がゼロのものもあります -- 系統樹の根(中心)と全く同一であることを意味します。
その他のウイルスは、変異が1-11の間です.

<br>

流行が今まさに起きている状況での、大型かつ新タガタのRNAウイルスのゲノムのシーケンスは難しいものです.
これらのシーケンス結果中の観測された違いには、実際の変異というよりは、シーケンス作業そのものエラーであるものも含まれるかもしれません.
ゲノム末端での挿入、欠損、変異はエラーである可能性が高く、この解析の目的に従い、ゲノム末端での違いを除外しています.

<br>

系統樹を"時間" ビューで表示する場合もありますし、 "多様性" ビューで表示する場合もあります。私たちが何を強調しようとしているかによって変わります。

<!-- There is NO right-side text -->


<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [系統解析の解釈](https://nextstrain.org/ncov/2020-03-05?d=tree)

私たちは現在5つの異なる大陸で採取されたサンプル由来のシーケンスを所有しています.
早期の症例は、海産食品市場での感染流行と関連した、武漢での症例と直接的に結びついていましたが、現在では、コミュニティでの拡散を示すものや、中国国外の感染源から持ち込まれた、様々な症例が見られます.

<br>

一般的に、動物を宿主とし何度も持ち込まれるウイルスは大きな多様性を示します(ラサ、エボラ、MERSや鳥インフルエンザで当てはまります).
このようにヒトでの感染でクラスタ化が強いことが観測されるということは、動物からヒトへののウイルスの持ち込みが一度きりの感染流行で、ヒトからヒトへの疫学的拡散が続いたと説明することができます.

<!-- There is NO right-side text -->


<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [イタリアへの少なくとも2度の持ち込み, コミュニティでの拡散](https://nextstrain.org/ncov/2020-03-05?d=tree&f_country=Italy)

現時点でイタリアからのシーケンスが3つあり、そのうち2つはローマ(Rome)から1つは北イタリアのロンバルディア(Lombardy)由来です.

<br>

これらの3つのシーケンスは感染流行の初期に共通の祖先を持ちます(左側、系統樹の元の近く),
このことはイタリアには少なくとも2度ウイルスが持ち込まれたことを強く示唆します。

<br>

Nuno Faria 博士らは [ここで](http://virological.org/t/first-cases-of-coronavirus-disease-covid-19-in-brazil-south-america-2-genomes-3rd-march-2020/409)ブラジルや世界の他のシーケンスが "北イタリアの感染流行が一つの感染源由来ではなく、複数回当該の地域に持ち込まれた結果である可能性が高い" という優れた分析を提示しています.



<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [イタリア内での検出されていない伝染の可能性](https://nextstrain.org/ncov/2020-03-05?d=tree&label=clade:A1a&m=div)

ローマRome(2020/01/29)由来の2つのシーケンスは互いに直接的に結びついていて、どちらの患者の方も中国への渡航歴がありました.

<br>

"多様性" ビューに切り替えて、イタリアの2つのシーケンスが同一であることを示しました。
一方で、(イングランド、ブラジル、スイス、アメリカ、中国由来の)近いシーケンスは実際そのイタリアのシーケンスとは2-4つの変異の差があります。

<br>

しかしながら、ブラジルのシーケンス(SPBR-02)はミラノ(Milan)のロンバルディア(Lombardy)への渡航歴があり、このクラスタ内でのスイスのシーケンスもイタリアへの渡航歴があったと考えられています.アメリカのサンプルの渡航歴については不明です。England/09cの1件は中国からの直接的な持ち込みです.

<br>
(中国からの)イングランドのサンプルは早期のローマからのイタリアのサンプルとイタリアへの渡航が基地のサンプル(スイス、ブラジル)またはその可能性がある(アメリカ)のサンプルの間にあります.これはこれらの早期のイタリアのシーケンスが最近のイタリアに関係するシーケンスと直接結びついていると想定すべきではないということを意味します.最近のサンプルはイタリアに独立して別に持ち込まれた持ち込まれた可能性があります。

<!-- There is NO right-side text -->



<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [イタリア, ロンバルディア(Lombardy)から地球規模への拡散](https://nextstrain.org/ncov/2020-03-05?d=tree&label=clade:A2)

ロンバルディア(Lombardy) (Italy/CDG1/2020)由来のシーケンスはイタリアへの渡航歴があり、そこで感染した可能性が高い既知のシーケンスとクラスタを形成しています: メキシコ、ドイツ、ブラジルそしてフィンランドからのシーケンスです.

<br>

ドイツの "BavPat1" は感染流行のかなり早期に中国から持ち込まれたものの一部です。もう一方のシーケンスに対するその類似度(変異はたった1つしか変わりません)は早期のドイツのクラスタに由来する、ヨーロッパでの未検出("cryptic"な)伝染を示唆しえます.

<br>

ヨーロッパに2度独立して持ち込まれた結果であるということもありえます - まだ採取されていないどこかのシーケンスが、"BavPat1 と他方のクラスタの間に入る、このとき、私たちはどちらのシナリオが確かに正しいかを言うことができなくなります.

<!-- There is NO right-side text -->


<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [シアトル(Seattle)でのSARS-CoV-2拡散の可能性](https://nextstrain.org/ncov/2020-03-05?label=clade:B1%20&d=tree)

COVID-19のいくつかの症例がシアトル大都市圏(greater Seattle area)とアメリカ全域で報告されています.
新たに分離されシーケンスが行われた症例は、1月半ばに同じ地域で単離されたケースと遺伝的に強く関連します.

<br>

これには2つの説明の可能性があります.
ウイルスは少なくとも2度シアトル大都市圏(greater Seattle area)に中国の共通の感染源から持ち込まれた可能性があります.
他方の説明としては、一定期間の間、ウイルスは未検出の状態で伝搬していたというものです.

<br>

Trevor Bedford(Nextstrainの共同設立者)はこれらの可能性については素晴らしい投稿をしています.[ここで](https://bedford.io/blog/ncov-cryptic-transmission/)閲覧することができます.

<br>

残ったワシントン(Washington)のシーケンスからは別なことがわかります: シアトル大都市圏(greater Seattle area)由来のシーケンスとともにクラスタを形成します.
このことはコミュニティでの拡散があること、SARS-CoV-2 ウイルスが当該のエリアですでにある一定の期間伝搬していることを強く示唆します.

<!-- There is NO right-side text -->



<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [最近共通祖先時間の推定](https://nextstrain.org/ncov/2020-03-05?label=clade:B1%20&d=tree)

シーケンスが行われた症例の集合の、最近共通祖先時間(つまりtMRCA)は、これらのシーケンスされた症例が最後に共通の祖先を持つのがいつかを示します.
下の図に示すように、この時間はウイルスがヒト社会に初めて侵入した時間と同じ程度であることがありますが、それよりもかなり遅いということもありえます.

<div>
  <img alt="Example phylogeny where the time of the initial zoonosis is different from the most recent common ancestor of several sequenced cases" width="500" src="https://raw.githubusercontent.com/nicfel/nCov-Nicola/master/figures/zoonosis.png"/>
</div>


<!-- This is right-side text -->
```auspiceMainDisplayMarkdown

## 感染が流行したウイルスの共通祖先の出現日

いくつかの研究グループや個人が最も近い共通祖先時間を評価しています。[Rambautによるこの投稿](http://virological.org/t/phylodynamic-analysis-of-sars-cov-2-update-2020-03-06/420)か[T Stadlerによるこの投稿](http://virological.org/t/evolutionary-epidemiological-analysis-of-93-genomes)を参照してください.

全シーケンスの共通の祖先は11月中頃と12月中頃の間に現れた可能性が最も高いです.
ことことは、現在のところ[武漢海産食品市場での初期のクラスタ](http://virological.org/t/phylodynamic-analysis-of-sars-cov-2-update-2020-03-06/420)以来、シーケンスされた全症例で一貫しているでしょう.


<div>
  <img alt="estimate of the tMRCA using Bayesian phylogenetics" width="500" src="https://raw.githubusercontent.com/nicfel/nCov-Nicola/master/figures/beast_coal-tmrca_2020303.png"/>
</div>

```





<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [クレジット](https://nextstrain.org/ncov/2020-03-05?d=map&c=author)

私たちはこの感染症流行の中、関係したすべての科学者によるすばらしくかつタイムリーな研究に、とくに中国内でのものに謝意を表します.
このようなゲノム情報とメタ情報の迅速な共有を通してのみ、解析は実現可能です.

<br>

私たちはデータをアップロード・共有するプラットフォームを提供している [GISAID](https://gisaid.org) にも謝意を表します.

<!-- Do not need to translate insitutions names -->
<!-- This is right-side text -->
```auspiceMainDisplayMarkdown

私たちはこれらの研修室によって集積されたデータに感謝しています:

* Centre for Infectious Diseases and Microbiology Laboratory Services
* Pathology Queensland
* Monash Medical Centre
* National Institute for Viral Disease Control and Prevention, China CDC
* KU Leuven, Clinical and Epidemiological Virology
* Hospital Israelita Albert Einstein
* Virology Unit, Institut Pasteur du Cambodge.
* BCCDC Public Health Laboratory
* Yongchuan District Center for Disease Control and Prevention
* Zhongxian Center for Disease Control and Prevention
* Respiratory Virus Unit, Microbiology Services Colindale, Public Health England
* Lapland Central Hospital
* HUS Diagnostiikkakeskus, Hallinto
* Guangdong Provincial Center for Diseases Control and Prevention; Guangdong Provincial Public Health
* Department of Infectious and Tropical Diseases, Bichat Claude Bernard Hospital, Paris
* Sorbonne Universite, Inserm et Assistance Publique-Hopitaux de Paris (Pitie Salpetriere)
* CNR Virus des Infections Respiratoires - France SUD
* Fujian Center for Disease Control and Prevention
* State Health Office Baden-Wuerttemberg
* Charite Universitatsmedizin Berlin, Institute of Virology; Institut fur Mikrobiologie der Bundeswehr, Munich
* Guangdong Provincial Center for Diseases Control and Prevention; Guangdong Provinical Public Health
* Guangdong Provincial Center for Diseases Control and Prevention;Guangdong Provincial Institute of Public Health
* Hangzhou Center for Disease and Control Microbiology Lab
* Hangzhou Center for Disease Control and Prevention
* Second Hospital of Anhui Medical University
* Hong Kong Department of Health
* Department of Infectious Diseases, Istituto Superiore di Sanita, Roma , Italy
* INMI Lazzaro Spallanzani IRCCS
* Department of Infectious Diseases, Istituto Superiore di Sanita, Rome, Italy
* Department of Virology III, National Institute of Infectious Diseases
* Dept. of Virology III, National Institute of Infectious Diseases
* Dept. of Pathology, National Institute of Infectious Diseases
* NHC Key laboratory of Enteric Pathogenic Microbiology, Institute of Pathogenic Microbiology
* Jingzhou Center for Disease Control and Prevention
* Division of Viral Diseases, Center for Laboratory Control of Infectious Diseases, Korea Centers for Diseases Control and Prevention
* Instituto Nacional de Enfermedades Respiratorias
* National Influenza Centre, National Public Health Laboratory, Kathmandu, Nepal
* Bamrasnaradura Hospital
* The University of Hong Kong - Shenzhen Hospital
* Shenzhen Third People's Hospital
* Shenzhen Key Laboratory of Pathogen and Immunity, National Clinical Research Center for Infectious Disease, Shenzhen Third People's Hospital
* Singapore General Hospital
* National Public Health Laboratory, National Centre for Infectious Diseases
* National Public Health Laboratory
* National Centre for Infectious Diseases
* Singapore General Hospital, Molecular Laboratory, Division of Pathology
* Korea Centers for Disease Control & Prevention (KCDC) Center for Laboratory Control of Infectious Diseases Division of Viral Diseases
* Serology, Virology and OTDS Laboratories (SAViD), NSW Health Pathology Randwick
* Centers for Disease Control, R.O.C. (Taiwan)
* Taiwan Centers for Disease Control
* Laboratory Medicine
* Department of Laboratory Medicine, National Taiwan University Hospital
* Tianmen Center for Disease Control and Prevention
* Arizona Department of Health Services
* California Department of Public Health
* California Department of Health
* IL Department of Public Health Chicago Laboratory
* Massachusetts Department of Public Health
* Texas Department of State Health Services
* WA State Department of Health
* Washington State Department of Health
* Providence Regional Medical Center
* Wisconsin Department of Health Services
* National Influenza Center - National Institute of Hygiene and Epidemiology (NIHE)
* Wuhan Jinyintan Hospital
* The Central Hospital Of Wuhan
* Union Hospital of Tongji Medical College, Huazhong University of Science and Technology
* CR & Wisco General Hospital
* Wuhan Lung Hospital
* Institute of Pathogen Biology, Chinese Academy of Medical Sciences & Peking Union Medical College
* Institute of Viral Disease Control and Prevention, China CDC
* General Hospital of Central Theater Command of People's Liberation Army of China
* Wuhan Fourth Hospital
* Zhejiang Provincial Center for Disease Control and Prevention
* Wuhan Institute of Virology, Chinese Academy of Sciences
* Shandong First Medical University & Shandong Academy of Medical Sciences
* South China Agricultural University
* Beijing Institute of Microbiology and Epidemiology

```


<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [Detailed scientific credit](https://nextstrain.org/ncov/2020-03-05?d=map&c=author)

これらのデータは[GISAID](https://gisaid.org)を通じて共有されています.
私たちはコントリビュータに謝意を表します.

<br>

それぞれの研究室によって共有されたシーケンスを右側に挙げています.

<!-- This is right-side text -->
```auspiceMainDisplayMarkdown

The SARS-CoV-2 genomes were generously shared by scientists at these submitting labs:

* NSW Health Pathology - Institute of Clinical Pathology and Medical Research; Westmead Hospital; University of Sydney
	* Australia/NSW01/2020
	* Australia/NSW05/2020
	* Sydney/2/2020

* Public Health Virology Laboratory
	* Australia/QLD01/2020
	* Australia/QLD02/2020
	* Australia/QLD03/2020
	* Australia/QLD04/2020

* Collaboration between the University of Melbourne at The Peter Doherty Institute for Infection and Immunity, and the Victorian Infectious Disease Reference Laboratory
	* Australia/VIC01/2020

* National Institute for Viral Disease Control & Prevention, CCDC
	* Beijing/IVDC-BJ-005/2020
	* Chongqing/IVDC-CQ-001/2020
	* Jiangsu/IVDC-JS-001/2020
	* Jiangxi/IVDC-JX-002/2020
	* Shandong/IVDC-SD-001/2020
	* Shanghai/IVDC-SH-001/2020
	* Sichuan/IVDC-SC-001/2020
	* Yunnan/IVDC-YN-003/2020

* KU Leuven, Clinical and Epidemiological Virology
	* Belgium/GHB-03021/2020

* Instituto Adolfo Lutz Interdisciplinary Procedures Center Strategic Laboratory
	* Brazil/SPBR-01/2020

* Virology Unit, Institut Pasteur du Cambodge (Sequencing done by: Jessica E Manning/Jennifer A Bohl at Malaria and Vector Research Research Laboratory, National Institute of Allergy and Infectious Diseases and Vida Ahyong from Chan-Zuckerberg Biohub)
	* Cambodia/0012/2020

* BCCDC Public Health Laboratory
	* Canada/BC_37_0-2/2020

* Technology Centre, Guangzhou Customs
	* China/IQTC01/2020
	* China/IQTC02/2020

* Key Laboratory of Human Diseases, Comparative Medicine, Institute of Laboratory Animal Science
	* China/WH-09/2020

* State Key Laboratory of Virology, Wuhan University
	* China/WHU01/2020
	* China/WHU02/2020

* Chongqing Municipal Center for Disease Control and Prevention
	* Chongqing/YC01/2020
	* Chongqing/ZX01/2020

* Respiratory Virus Unit, Microbiology Services Colindale, Public Health England
	* England/01/2020
	* England/02/2020
	* England/09c/2020

* Department of Virology, University of Helsinki and Helsinki University Hospital, Helsinki, Finland
	* Finland/1/2020

* Department of Virology Faculty of Medicine, Medicum University of Helsinki
	* Finland/FIN-25/2020

* Guangdong Provincial Center for Diseases Control and Prevention
	* Foshan/20SF207/2020
	* Foshan/20SF210/2020
	* Foshan/20SF211/2020
	* Guangdong/20SF201/2020
	* Guangzhou/20SF206/2020

* National Reference Center for Viruses of Respiratory Infections, Institut Pasteur, Paris
	* France/IDF0372-isl/2020
	* France/IDF0372/2020
	* France/IDF0373/2020
	* France/IDF0515-isl/2020
	* France/IDF0515/2020
	* France/IDF0626/2020

* Laboratoire Virpath, CIRI U111, UCBL1, INSERM, CNRS, ENS Lyon
	* France/IDF0386-islP1/2020
	* France/IDF0386-islP3/2020
	* France/IDF0571/2020

* CNR Virus des Infections Respiratoires - France SUD
	* France/RA739/2020

* Fujian Center for Disease Control and Prevention
	* Fujian/13/2020
	* Fujian/8/2020

* Charite Universitatsmedizin Berlin, Institute of Virology
	* Germany/Baden-Wuerttemberg-1/2020
	* Germany/BavPat1/2020

* Department of Microbiology, Guangdong Provincial Center for Diseases Control and Prevention
	* Guangdong/20SF012/2020
	* Guangdong/20SF013/2020
	* Guangdong/20SF014/2020
	* Guangdong/20SF025/2020
	* Guangdong/20SF028/2020
	* Guangdong/20SF040/2020

* Guangdong Provincial Center for Disease Control and Prevention
	* Guangdong/20SF174/2020

* Hangzhou Center for Disease and Control Microbiology Lab
	* Hangzhou/HZ-1/2020

* Hangzhou Center for Disease Control and Prevention
	* Hangzhou/HZCDC0001/2020

* Second Hospital of Anhui Medical University
	* Hefei/2/2020

* National Institute for Viral Disease Control & Prevention, China CDC
	* Henan/IVDC-HeN-002/2020

* School of Public Health, The University of Hon g Kong
	* HongKong/VB20026565/2020
	* HongKong/VM20001061/2020

* The University of Hong Kong
	* HongKong/VM20001988/2020
	* Nepal/61/2020

* Virology Laboratory, Scientific Department, Army Medical Center
	* Italy/CDG1/2020
	* Italy/SPL1/2020

* Laboratory of Virology, INMI Lazzaro Spallanzani IRCCS
	* Italy/INMI1-cs/2020
	* Italy/INMI1-isl/2020

* Pathogen Genomics Center, National Institute of Infectious Diseases
	* Japan/AI/I-004/2020
	* Japan/KY-V-029/2020
	* Japan/NA-20-05-1/2020
	* Japan/OS-20-07-1/2020
	* Japan/TY-WK-012/2020
	* Japan/TY-WK-501/2020
	* Japan/TY-WK-521/2020

* Takayuki Hishiki Kanagawa Prefectural Institute of Public Health, Department of Microbiology
	* Japan/Hu_DP_Kng_19-020/2020
	* Japan/Hu_DP_Kng_19-027/2020

* Jiangsu Provincial Center for Disease Control & Prevention
	* Jiangsu/JS01/2020
	* Jiangsu/JS02/2020
	* Jiangsu/JS03/2020

* Hubei Provincial Center for Disease Control and Prevention
	* Jingzhou/HBCDC-HB-01/2020
	* Tianmen/HBCDC-HB-07/2020
	* Wuhan/HBCDC-HB-01/2019
	* Wuhan/HBCDC-HB-02/2019
	* Wuhan/HBCDC-HB-02/2020
	* Wuhan/HBCDC-HB-03/2019
	* Wuhan/HBCDC-HB-03/2020
	* Wuhan/HBCDC-HB-04/2019
	* Wuhan/HBCDC-HB-04/2020
	* Wuhan/HBCDC-HB-05/2020
	* Wuhan/HBCDC-HB-06/2020

* Division of Viral Diseases, Center for Laboratory Control of Infectious Diseases, Korea Centers for Diseases Control and Prevention
	* Korea/KCDC05/2020
	* Korea/KCDC06/2020
	* Korea/KCDC07/2020
	* Korea/KCDC12/2020
	* Korea/KCDC24/2020

* Instituto de Diagnostico y Referencia Epidemiologicos (INDRE)
	* Mexico/CDMX/InDRE_01/2020

* 1. Department of Medical Sciences, Ministry of Public Health, Thailand 2. Thai Red Cross Emerging Infectious Diseases - Health Science Centre 3. Department of Disease Control, Ministry of Public Health, Thailand
	* Nonthaburi/61/2020
	* Nonthaburi/74/2020

* Li Ka Shing Faculty of Medicine, The University of Hong Kong
	* Shenzhen/HKU-SZ-002/2020
	* Shenzhen/HKU-SZ-005/2020

* Shenzhen Key Laboratory of Pathogen and Immunity, National Clinical Research Center for Infectious Disease, Shenzhen Third People's Hospital
	* Shenzhen/SZTH-001/2020
	* Shenzhen/SZTH-002/2020
	* Shenzhen/SZTH-003/2020
	* Shenzhen/SZTH-004/2020

* National Public Health Laboratory
	* Singapore/1/2020
	* Singapore/11/2020

* National Centre for Infectious Diseases, National Centre for Infectious Diseases
	* Singapore/10/2020

* Programme in Emerging Infectious Diseases, Duke-NUS Medical School
	* Singapore/2/2020
	* Singapore/3/2020
	* Singapore/4/2020
	* Singapore/5/2020
	* Singapore/6/2020

* National Public Health Laboratory, National Centre for Infectious Diseases
	* Singapore/7/2020
	* Singapore/8/2020
	* Singapore/9/2020

* Korea Centers for Disease Control & Prevention (KCDC) Center for Laboratory Control of Infectious Diseases Division of Viral Diseases
	* SouthKorea/KCDC03/2020

* Department of Clinical Diagnostics
	* SouthKorea/SNU01/2020

* Unit for Laboratory Development and Technology Transfer, Public Health Agency of Sweden
	* Sweden/01/2020

* NSW Health Pathology - Institute of Clinical Pathology and Medical Research; Centre for Infectious Diseases and Microbiology Laboratory Services; Westmead Hospital; University of Sydney
	* Sydney/3/2020

* Centers for Disease Control, R.O.C. (Taiwan)
	* Taiwan/2/2020

* Taiwan Centers for Disease Control
	* Taiwan/3/2020
	* Taiwan/4/2020

* Department of Laboratory Medicine, Lin-Kou Chang Gung Memorial Hospital, Taoyuan, Taiwan.
	* Taiwan/CGMH-CGU-01/2020

* Microbial Genomics Core Lab, National Taiwan University Centers of Genomic and Precision Medicine
	* Taiwan/NTU01/2020
	* Taiwan/NTU02/2020

* Pathogen Discovery, Respiratory Viruses Branch, Division of Viral Diseases, Centers for Disease Control and Prevention
	* USA/AZ1/2020
	* USA/CA1/2020
	* USA/CA2/2020
	* USA/CA3/2020
	* USA/CA4/2020
	* USA/CA5/2020
	* USA/CA6/2020
	* USA/CA7/2020
	* USA/CA8/2020
	* USA/CA9/2020
	* USA/IL1/2020
	* USA/IL2/2020
	* USA/MA1/2020
	* USA/TX1/2020
	* USA/WA1-A12/2020
	* USA/WA1-F6/2020
	* USA/WI1/2020

* Division of Viral Diseases, Centers for Disease Control and Prevention
	* USA/WA1/2020

* Seattle Flu Study
	* USA/WA2/2020

* National Influenza Center - National Institute of Hygiene and Epidemiology (NIHE)
	* Vietnam/VR03-38142/2020

* National Institute for Communicable Disease Control and Prevention (ICDC) Chinese Center for Disease Control and Prevention (China CDC)
	* Wuhan-Hu-1/2019

* Institute of Pathogen Biology, Chinese Academy of Medical Sciences & Peking Union Medical College
	* Wuhan/IPBCAMS-WH-01/2019
	* Wuhan/IPBCAMS-WH-02/2019
	* Wuhan/IPBCAMS-WH-03/2019
	* Wuhan/IPBCAMS-WH-04/2019
	* Wuhan/IPBCAMS-WH-05/2020

* National Institute for Viral Disease Control and Prevention, China CDC
	* Wuhan/IVDC-HB-01/2019
	* Wuhan/IVDC-HB-04/2020
	* Wuhan/IVDC-HB-05/2019

* Institute of Viral Disease Control and Prevention, China CDC
	* Wuhan/IVDC-HB-envF13-20/2020
	* Wuhan/IVDC-HB-envF13-21/2020
	* Wuhan/IVDC-HB-envF13/2020
	* Wuhan/IVDC-HB-envF54/2020

* BGI & Institute of Microbiology, Chinese Academy of Sciences & Shandong First Medical University & Shandong Academy of Medical Sciences & General Hospital of Central Theater Command of People's Liberation Army of China
	* Wuhan/WH01/2019
	* Wuhan/WH02/2019
	* Wuhan/WH03/2020
	* Wuhan/WH04/2020

* Beijing Genomics Institute (BGI)
	* Wuhan/WH05/2020

* Wuhan Institute of Virology, Chinese Academy of Sciences
	* Wuhan/WIV02/2019
	* Wuhan/WIV04/2019
	* Wuhan/WIV05/2019
	* Wuhan/WIV06/2019
	* Wuhan/WIV07/2019

* Department of Microbiology, Zhejiang Provincial Center for Disease Control and Prevention
	* Zhejiang/WZ-01/2020
	* Zhejiang/WZ-02/2020



```
