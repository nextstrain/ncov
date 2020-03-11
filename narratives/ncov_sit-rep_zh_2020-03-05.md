---
title: 新型冠状病毒（COVID-19）流行病学基因组分析（状况报告2020-03-05）
authors: "Trevor Bedford, Richard Neher, James Hadfield, Emma Hodcroft, Misja Ilcisin, Nicola Müller, Derek Zhang, Fengjun Zhang"
authorLinks: "https://nextstrain.org"
affiliations: "Fred Hutch, Seattle, USA and Biozentrum, Basel, Switzerland"
date: "2020年3月5日"
dataset: "https://nextstrain.org/ncov/2020-03-05"

abstract: "这份报告使用了在GISAID和Genbank数据库公开共享的2019新型冠状病毒（COVID-19）基因组数据来推算疫情的传播速度和方式。一旦收集到新的病毒数据，我们将陆续更新该状况报告。我们建议您使用桌面浏览器游览该网站。"

---
<!-- Translators: Only text after : in the above ^ needs to be translated -->
<!-- Comment tags like these do not need to be translated, they are only to help you! -->
<!-- Ensure that links always end in a 'letter' (. counts) If some kind of text doesn't follow them, it breaks the slide. -->
<!-- numbers can be tagged ilke this: 161</tag> - this is just for us to help find them to update! Just leave in the </tag> bit. -->

<!-- This is left-side text -->
# [报告摘要](https://nextstrain.org/ncov/2020-03-05)

<!-- This is right-side text -->
```auspiceMainDisplayMarkdown
## 报告摘要

通过使用 169</tag> 个公开共享的病毒基因组, 我们分析了COVID-19新型冠状病毒的遗传多样性，用以描述疫情在各个地区的传播情况，并且推断病毒共同祖先的日期。



我们发现：
* 新冠病毒COVID-19至少被传入意大利两次，并且都在意大利引起了社区传播。([详细分析](https://nextstrain.org/narratives/ncov/sit-rep/zh/2020-03-05?n=7))
* 在意大利发生的人际传播也输出到了其他国家———有来自六个不同国家的序列按照核酸多样性被分在一个基因簇中，数据显示这些病例是意大利输出的病例。([详细分析](https://nextstrain.org/narratives/ncov/sit-rep/zh/2020-03-05?n=9))
* 根据现有的核酸序列数据可推得，COVID-19 在大西雅图地区从一月中旬起就已经在未被检测到的情况下发生了持续传播。([详细分析](https://nextstrain.org/narratives/ncov/sit-rep/zh/2020-03-05?n=10))
* 这项分析中包括的所有病毒个体可能有一个共同祖先，其出现时间范围在在2019年11月中旬至12月中旬之间。([详细分析](https://nextstrain.org/narratives/ncov/sit-rep/zh/2020-03-05?n=11))

```


<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [什么是冠状病毒？](https://nextstrain.org/ncov/2020-03-05)

### 更多资料:

* 维基百科上 COVID-19 疫情的总结 [中文维基](https://zh.wikipedia.org/wiki/2019%E5%86%A0%E7%8B%80%E7%97%85%E6%AF%92%E7%97%85%E7%96%AB%E6%83%85) _2020-01-30_
* 美国疾病预防控制中心提供的材料 [CDC](https://www.cdc.gov/coronavirus/2019-ncov/index-Chinese.html) _2020-01-29_

<!-- This is right-side text -->
```auspiceMainDisplayMarkdown

## COVID-19 的一些其他信息


我们准备了一些值得一读的材料，让您了解COVID-19（2019冠状病毒病）及其致病病毒SARS-CoV-2（严重急性呼吸综合征冠状病毒2型）。以下材料来源于 [Nextstrain](https://nextstrain.org) 网站的英文文档，我们建议您用谷歌翻译或者类似的机器翻译服务辅助浏览。
这些信息能使您更容易理解我们在这篇报告中的分析和结论。

<div>
  <a href="https://nextstrain.org/help/coronavirus/human-CoV"><img alt="microscopy image of coronaviruses" width="100" src="https://nextstrain.org/static/ncov_narrative-76cfd610d11ef708d213a3170de9519f.png"/> 冠状病毒的背景知识 </a>

  <a href="https://nextstrain.org/help/coronavirus/SARS-CoV-2"><img alt="illustration of a coronavirus" width="100" src="http://data.nextstrain.org/img_nCoV-CDC.jpg"/> 近期 COVID-19 疫情总结 </a>

  <a href="https://nextstrain.org/help/general/how-to-read-a-tree"><img alt="cartoon of a phylogenetic tree" width="100" src="http://data.nextstrain.org/img_toy_alignment_mini.png"/> 如何解读进化树 </a>

</div>

## Nextstrain 状况描述报告

以下页面包含使用 [Nextstrain](https://nextstrain.org) 执行的分析。
滚动左侧侧栏将显示文本段落的同时，在右侧会显示相应的基因组数据可视化图。

对于一种新型的、基因组较大的RNA病毒，在疫情正在爆发的当下，这么迅速就能获得多条全基因组序列是一项了不起的成就。
而这些成果是由世界各地的科学家快速公开共享基因组数据和分析而实现的（您可阅读本报告最后一页，我们列出了详细的序列提供单位）。


```


<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [常见问题和误区](https://nextstrain.org/ncov/2020-03-05)

### 更多详细资料：

_以下材料均来自英文网站，我们建议您使用机器翻译服务辅助阅读。_

* "不要相信关于冠状病毒和艾滋病的阴谋论" [文章链接](https://massivesci.com/notes/wuhan-coronavirus-ncov-sars-mers-hiv-human-immunodeficiency-virus/) _2020-01-31_

* "毫无根据的阴谋论声称新的冠状病毒是生物工程制造的" [文章链接](https://www.factcheck.org/2020/02/baseless-conspiracy-theories-claim-new-coronavirus-was-bioengineered/) _2020-02-07_

* "不，武汉冠状病毒并没有被人为插入艾滋病病毒片段。" [文章链接](https://www.forbes.com/sites/victoriaforster/2020/02/02/no-coronavirus-was-not-bioengineered-to-put-pieces-of-hiv-in-it/#5d339e8e56cb) _2020-02-02_

* "关于冠状病毒的辟谣" [文章链接](https://factcheck.afp.com/busting-coronavirus-myths) _2020-02-19_

<!-- This is right-side text -->
```auspiceMainDisplayMarkdown

## 常见问题和误区

_以下材料均来自英文网站，我们建议您使用机器翻译服务辅助阅读。_

### 常见问题

我们知道很多人对 COVID-19 有疑问，
我们特地设立了一个指南来回答一些最常被问到的问题。[FAQ](https://nextstrain.org/help/coronavirus/FAQ):

<div>

  <a href="https://nextstrain.org/help/coronavirus/FAQ"><img alt="picture of a question mark" width="100" src="http://data.nextstrain.org/img_question-mark.jpg"/> COVID-19 常见问题</a>

</div>


### 常见误区

目前社会上流行着许多关于新型冠状病毒起源的谣言。
但在像这样的疫情爆发期间，不正确信息的传播可能会引发更多的恐慌，使人们丧失对科学家和政府的信任，从而可能不遵循建议或不采取适当的预防措施。

科学家们在下面的文章中了解释了为什么一些观点是不正确的：

<div>

  <a href="http://virological.org/t/ncovs-relationship-to-bat-coronaviruses-recombination-signals-no-snakes-no-evidence-the-2019-ncov-lineage-is-recombinant/331"><img alt="picture of a snake" width="100" src="http://data.nextstrain.org/img_snake-freeToUse.jpg"/> SARS-CoV-2 并不起源于"蛇"（专业论坛帖）</a>
  <a href="https://twitter.com/trvrb/status/1223666856923291648"><img alt="illustration of HIV" width="100" src="http://data.nextstrain.org/img_HIV-wiki.jpg"/> 'HIV基因工程'之说错在哪里（推特系列帖）</a>


</div>


```


<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [系统进化树分析](https://nextstrain.org/ncov/2020-03-05?d=tree)

我们在这里给出了由 169</tag> 株 SARS-CoV-2（即 COVID-19 致病病毒）病毒个体的公开序列构成的系统进化树。
分析的方法可在 [此GitHub仓库中找到](https://github.com/nextstrain/ncov)。

<br>

数据点的颜色表示分离出该病毒的各国地区或美国各州，x 轴则表示取样日期。

y 轴显示的是序列之间的关系，请注意它没有测量单位。

<br>

取样日期信息虽然很有用，但它并不能总是准确地显示两个序列在遗传上的相关性。

两个完全相同的序列可能具有不同的取样日期，因而在此图中，它们会看起来相距甚远。

<!-- There is NO right-side text -->


<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [系统进化树上的“差异程度”](https://nextstrain.org/ncov/2020-03-05?d=tree&m=div)

我们现在可以从“时间”（Time）视图模式转换成“差异程度”（Divergence）视图模式，使 x 轴按“差异程度”排列。

<br>

请注意之前看起来关联性不大的序列，现在可能会排在一条竖直线上。
您可以通过上下滚动在上一张幻灯片和这张幻灯片之间切换，以查看系统进化树是如何变化的。

<br>

遗传差异程度是用基因组中的核酸变化（自然突变）的数量来衡量的。
有些序列中的突变数量为零，也就是说它们与系统进化树的“根部”序列一致。
而其他病毒个体中的突变数量则有1到11个不等。

<br>

在不断变化的流行病爆发情况下，对新发现的大型RNA病毒的基因组进行测序是一项具有挑战性的任务。
因为在这些序列中观察到的核酸差异可能是由于测序误差，而不是实际的核酸突变。
特别是基因组末端的插入、缺失等差异有更高的概率是测序误差导致的，所以在本次分析中，我们没有包括这些区域的差异。

<br>

在后面的页面中，我们会根据所讨论的问题，有时显示“时间”视图中的系统进化树，有时会展示“差异程度”视图中的系统进化树。

<!-- There is NO right-side text -->


<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [系统进化树的解释](https://nextstrain.org/ncov/2020-03-05?d=tree)

我们目前有来自五个不同大陆的样本的序列。
虽然早期的病例都或多或少与武汉华南海鲜市场有着密切联系，但我们现在观察到了各种不同情况的病例，有的来自于社区传播，另一些则是从中国以外的国家输入的。

<br>

通常，如果病毒来自于从动物宿主到人的反复多次感染，病毒的种群基因组将显示出显著的核酸碱基多样性（拉萨病毒、埃博拉病毒、MERS冠状病毒和禽流感都是如此）。

然而，此次我们观察到来源于人类病例的样本能显著地归成系统进化树上的基因簇，那么基本可以用一次爆发模型来解释疫情的起源，即只发生了一次从动物传到人的事件，随后发展成人与人之间的种群内传播。

<!-- There is NO right-side text -->


<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [病毒至少被引入意大利两次，两次都可能存在社区传播](https://nextstrain.org/ncov/2020-03-05?d=tree&f_country=Italy)

我们目前有三条来自意大利的序列，其中两条来自罗马地区，一条来自意大利北部的伦巴第。

<br>

这三条序列的共同祖先出现在疫情的早期（这个祖先株在图里的位置在系统进化树根部附近，非常接近纵轴），这个证据强有力地表明病毒至少有两次被输入到意大利，并引发了社区传播。

<br>

努诺·法里亚博士团队非常详细地分析了巴西和其他全球序列，并得出“意大利北部的疫情很可能是病毒多次传入该地区的结果，而不是来自一个单一的来源”的结论。 [点击此处可以查看原文文章](http://virological.org/t/first-cases-of-coronavirus-disease-covid-19-in-brazil-south-america-2-genomes-3rd-march-2020/409)。



<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [意大利可能存在的隐藏传播路线](https://nextstrain.org/ncov/2020-03-05?d=tree&label=clade:A1a&m=div)

其中来自罗马（2020年1月29日）的两个序列联系密切，都有中国的旅行史。

<br>

我们切换到“差异程度”视图，您可以发现这两个意大利序列是相同的，而其他邻近的序列（分别来自英国、巴西、瑞士、美国和中国）实际上与意大利序列相差2-4个突变。

<br>


该基因簇中巴西的序列（SPBR-02），已知有前往伦巴第地区米兰市的旅行历史，而有报道称该簇中的瑞士病例最近也曾前往意大利。

我们对同簇中的美国病例的旅行史一无所知，但英国的病例（England/09c）已知是从中国直接输入的。

<br>

源于中国的这例英国序列，处于早期意大利罗马的病例和已知（瑞士、巴西）或可能（美国）有意大利旅行史的病例之间。

这意味着我们不应该假设较早的意大利序列和这些较新的与意大利有关的序列有直接的关系，这些最近的意大利序列也有可能来源于另一次独立输入感染的事件。

<!-- There is NO right-side text -->



<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [从意大利伦巴第传播至全球其他地域的序列](https://nextstrain.org/ncov/2020-03-05?d=tree&label=clade:A2)

来自伦巴第（Italy/CDG1/2020）的序列与其他国家的已知有意大利旅行史的序列被分在了一起（已知有可能在意大利被感染的病例分别来自墨西哥、德国、巴西和芬兰）。

<br>

图中最下方的德国的“BavPat1”序列是在疫情更早的时候从中国传入的一例病例。
它与该簇中的另一个序列的相似性很高（只有一个突变之分），暗示欧洲可能存在起源于这个早期德国病例的、未被检测到的（“不明途径的”）疫情传播路线。

<br>

这个现象也可能是病毒两次单独输入欧洲的结果，因为可能会有来自其他地方的、未被采集到的序列分布在“BavPat1”和这个大基因簇之间。
但以目前的样本量，我们不能确认是哪种情况。

<!-- There is NO right-side text -->


<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [SARS-CoV-2 可能已在西雅图地区蔓延](https://nextstrain.org/ncov/2020-03-05?label=clade:B1%20&d=tree)

现在大西雅图地区和整个美国范围内都报告了多起COVID-19病例。
这些新分离和测序得到的序列在遗传上与一月中旬同一地区分离到的一例病例密切相关。

<br>

对此，我们有两种可能的解释。
一种解释是，该病毒可能至少两次从中国的共同源头传入大西雅图地区。
而另一种解释是，病毒已经在该地区传播了一段时间，只是没有被发现。

<br>

特雷弗·贝德福德（Nextstrain的联合创始人）就这些可能性写了一篇很棒的博客文章，你可以在这里读到 [他的相关技术博客文章](https://bedford.io/blog/ncov-cryptic-transmission/)。

<br>

华盛顿州的其他新序列还告诉我们另一件事：来自大西雅图地区的序列在核酸分析中被排列在了一起。
这强烈表明该地区已存在社区传播，且 SARS-CoV-2 病毒已经传播了一段时间。

<!-- There is NO right-side text -->



<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [测定最近共同祖先的时间](https://nextstrain.org/ncov/2020-03-05?label=clade:B1%20&d=tree)

一群个体序列的最近共同祖先时间（即tMRCA）能够提示这些个体共同的祖先株最近在何时出现。
如下图所示，这一时间可能是病毒首次进入人类群体的时间，但也可能大大晚于此时间。

<div>
  <img alt="Example phylogeny where the time of the initial zoonosis is different from the most recent common ancestor of several sequenced cases" width="500" src="https://raw.githubusercontent.com/nicfel/nCov-Nicola/master/figures/zoonosis.png"/>
</div>


<!-- This is right-side text -->
```auspiceMainDisplayMarkdown

## 病毒的最近共同祖先时间

数个研究小组已经估算了最近共同祖先时间，您可以查看 [A Rambaut的报告（英文）](http://virological.org/t/phylodynamic-analysis-of-sars-cov-2-update-2020-03-06/420) 或 [T Stadler的报告（英文）](http://virological.org/t/evolutionary-epidemiological-analysis-of-93-genomes)。

所有序列的共同祖先最有可能出现在11月中旬至12月中旬之间。
这也支持目前所有已测序病毒都来源于武汉华南海鲜市场的聚集性病例的说法 （[点击此处查看专业论坛的讨论帖](http://virological.org/t/phylodynamic-analysis-of-sars-cov-2-update-2020-03-06/420)）。


<div>
  <img alt="estimate of the tMRCA using Bayesian phylogenetics" width="500" src="https://raw.githubusercontent.com/nicfel/nCov-Nicola/master/figures/beast_coal-tmrca_2020303.png"/>
</div>

```





<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [科学贡献人员](https://nextstrain.org/ncov/2020-03-05?d=map&c=author)

我们要感谢参与此次疫情的所有科学家所做的令人惊叹和及时的工作，特别是那些在中国工作的科学家。
只有通过基因组数据和元数据的快速共享，才有可能进行这样的分析。

<br>

我们也感谢[GISAID](https://gisaid.org)提供了上传和共享这些数据的平台。

<!-- Do not need to translate insitutions names -->
<!-- This is right-side text -->
```auspiceMainDisplayMarkdown

非常感谢以下单位机构和实验室慷慨分享:

* Centre for Infectious Diseases and Microbiology Laboratory Services
* Pathology Queensland
* Monash Medical Centre
* 中国疾控中心国家病毒性疾病预防控制研究所
* KU Leuven, Clinical and Epidemiological Virology
* Hospital Israelita Albert Einstein
* Virology Unit, Institut Pasteur du Cambodge.
* BCCDC Public Health Laboratory
* 重庆市永川区疾病预防控制中心
* 重庆市忠县疾病预防控制中心
* Respiratory Virus Unit, Microbiology Services Colindale, Public Health England
* Lapland Central Hospital
* HUS Diagnostiikkakeskus, Hallinto
* 广东省疾病预防控制中心；广东省公共卫生厅
* Department of Infectious and Tropical Diseases, Bichat Claude Bernard Hospital, Paris
* Sorbonne Universite, Inserm et Assistance Publique-Hopitaux de Paris (Pitie Salpetriere)
* CNR Virus des Infections Respiratoires - France SUD
* 福建省疾病预防控制中心
* State Health Office Baden-Wuerttemberg
* Charite Universitatsmedizin Berlin, Institute of Virology; Institut fur Mikrobiologie der Bundeswehr, Munich
* 广东省疾病预防控制中心；广东省公共卫生厅
* 广东省疾病预防控制中心；广东省公共卫生研究所
* 杭州市疾病预防控制中心微生物学实验室
* 杭州市疾病预防控制中心
* 安徽医科大学第二医院
* 香港卫生署
* Department of Infectious Diseases, Istituto Superiore di Sanita, Roma , Italy
* INMI Lazzaro Spallanzani IRCCS
* Department of Infectious Diseases, Istituto Superiore di Sanita, Rome, Italy
* Department of Virology III, National Institute of Infectious Diseases
* Dept. of Virology III, National Institute of Infectious Diseases
* Dept. of Pathology, National Institute of Infectious Diseases
* NHC Key laboratory of Enteric Pathogenic Microbiology, Institute of Pathogenic Microbiology
* 荆州市疾病预防控制中心
* Division of Viral Diseases, Center for Laboratory Control of Infectious Diseases, Korea Centers for Diseases Control and Prevention
* Instituto Nacional de Enfermedades Respiratorias
* National Influenza Centre, National Public Health Laboratory, Kathmandu, Nepal
* Bamrasnaradura Hospital
* 香港大学深圳医院
* 深圳市第三人民医院
* 深圳市第三人民医院国家传染病临床研究中心深圳市病原与免疫重点实验室
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
* 天门市疾病预防控制中心
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
* 武汉金银潭医院
* 武汉市中心医院
* 华中科技大学同济医学院协和医院
* CR & Wisco General Hospital
* 武汉肺科医院
* 中国医学科学院、北京协和医学院病原生物学研究所
* 中国疾控中心病毒病预防控制研究所
* 中国人民解放军中央战区总医院
* 武汉市第四医院
* 浙江省疾病预防控制中心
* 中国科学院武汉病毒研究所
* 山东第一医科大学&山东省医学科学院
* 华南农业大学
* 北京微生物学与流行病学研究所

```


<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [科学贡献详情](https://nextstrain.org/ncov/2020-03-05?d=map&c=author)

所有基因组数据已被共享至 [GISAID](https://gisaid.org)数据库。感谢各位科学家的慷慨分享。

<br>

在右边，我们列出了每个实验室共享的序列。

<!-- This is right-side text -->
```auspiceMainDisplayMarkdown

感谢下列实验室里科学家们慷慨分享SARS-CoV-2的基因组:

* 新南威尔士州健康病理学 - 临床病理和医学研究所, 韦斯特米德医院, 悉尼大学
	* Australia/NSW01/2020
	* Australia/NSW05/2020
	* Sydney/2/2020

* 公共卫生病毒学实验室， 澳大利亚
	* Australia/QLD01/2020
	* Australia/QLD02/2020
	* Australia/QLD03/2020
	* Australia/QLD04/2020

* 墨尔本大学,The Peter Doherty感染和免疫研究所, 维多利亚传染病参考实验室之间的合作
	* Australia/VIC01/2020

* 中国疾病预防控制中心 病毒病预防控制所
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

* 科学中心, 广州海关
	* China/IQTC01/2020
	* China/IQTC02/2020

* Key Laboratory of Human Diseases, Comparative Medicine, Institute of Laboratory Animal Science
	* China/WH-09/2020

* 病毒学国家重点实验室, 武汉大学
	* China/WHU01/2020
	* China/WHU02/2020

* 重庆市疾控中心
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

* 广东省疾控中心
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

* 福建省疾控中心
	* Fujian/13/2020
	* Fujian/8/2020

* Charite Universitatsmedizin Berlin, Institute of Virology
	* Germany/Baden-Wuerttemberg-1/2020
	* Germany/BavPat1/2020

* 广东省疾控中心微生物部门
	* Guangdong/20SF012/2020
	* Guangdong/20SF013/2020
	* Guangdong/20SF014/2020
	* Guangdong/20SF025/2020
	* Guangdong/20SF028/2020
	* Guangdong/20SF040/2020

* 广东省疾控中心
	* Guangdong/20SF174/2020

* 杭州疾控中心微生物实验室
	* Hangzhou/HZ-1/2020

* 杭州疾控中心
	* Hangzhou/HZCDC0001/2020

* 安徽医科大学第二医院
	* Hefei/2/2020

* 国家病毒性疾病控制与预防研究所, 中国疾控中心
	* Henan/IVDC-HeN-002/2020

* School of Public Health, 香港大学
	* HongKong/VB20026565/2020
	* HongKong/VM20001061/2020

* 香港大学
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

* 江苏省疾控中心
	* Jiangsu/JS01/2020
	* Jiangsu/JS02/2020
	* Jiangsu/JS03/2020

* 湖北省疾控中心
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

* 李嘉诚医学院, 香港大学
	* Shenzhen/HKU-SZ-002/2020
	* Shenzhen/HKU-SZ-005/2020

* 深圳市病原与免疫重点实验室, 国家传染病临床研究中心, 深圳第三人民医院
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

* 中国疾病预防控制中心(中国疾控中心)
	* Wuhan-Hu-1/2019

* 中国医学科学院、北京协和医学院病原生物学研究所
	* Wuhan/IPBCAMS-WH-01/2019
	* Wuhan/IPBCAMS-WH-02/2019
	* Wuhan/IPBCAMS-WH-03/2019
	* Wuhan/IPBCAMS-WH-04/2019
	* Wuhan/IPBCAMS-WH-05/2020

* 中国疾控中心国家病毒病预防控制研究所
	* Wuhan/IVDC-HB-01/2019
	* Wuhan/IVDC-HB-04/2020
	* Wuhan/IVDC-HB-05/2019

* 中国疾控中心病毒病预防控制研究所
	* Wuhan/IVDC-HB-envF13-20/2020
	* Wuhan/IVDC-HB-envF13-21/2020
	* Wuhan/IVDC-HB-envF13/2020
	* Wuhan/IVDC-HB-envF54/2020

* 华大基因与中国科学院微生物学研究所、山东第一医科大学、山东医学科学院、中国人民解放军中心战区总医院
	* Wuhan/WH01/2019
	* Wuhan/WH02/2019
	* Wuhan/WH03/2020
	* Wuhan/WH04/2020

* 北京基因组研究所 (BGI)
	* Wuhan/WH05/2020

* 中国科学院武汉病毒研究所
	* Wuhan/WIV02/2019
	* Wuhan/WIV04/2019
	* Wuhan/WIV05/2019
	* Wuhan/WIV06/2019
	* Wuhan/WIV07/2019

* 浙江省疾病预防控制中心微生物室
	* Zhejiang/WZ-01/2020
	* Zhejiang/WZ-02/2020



```
