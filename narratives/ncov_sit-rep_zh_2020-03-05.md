---
title: 新型冠状病毒（COVID-19）流行病学基因组分析（状况报告2020-03-05）
authors: "Trevor Bedford, Richard Neher, James Hadfield, Emma Hodcroft, Misja Ilcisin, Nicola Müller, Derek Zhang, Fengjun Zhang"
authorLinks: "https://nextstrain.org"
affiliations: "Fred Hutch, Seattle, USA and Biozentrum, Basel, Switzerland"
date: "2020年3月5日"
dataset: "https://nextstrain.org/ncov/2020-03-05"
abstract: "这份报告使用了来自GISAID和Genbank公开共享的2019新型冠状病毒（COVID-19）基因组数据来推算疫情的流传速度和方式。一旦收集到新的病毒数据，我们将陆续更新该状况报告。我们建议您使用桌面浏览器游览该网站。"
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

通过使用 169</tag> 个公开共享的 COVID-19 新型冠状病毒基因组, 我们分析了各病毒个体的遗传多样性，以推断病毒共同祖先的日期和疫情传播速度。



我们发现：
* 在意大利，新冠病毒 SARS-CoV-2 发生了至少两次的输入感染事件，并且两次事件都伴随发生了在意大利的本地区内传播。([详细分析](https://nextstrain.org/narratives/ncov/sit-rep/zh/2020-03-05?n=7))
* 来自六个不同国家的序列，按照核酸多样性被分在一个基因簇中，这些国家的病例应该是输出自意大利的病例。([详细分析](https://nextstrain.org/narratives/ncov/sit-rep/zh/2020-03-05?n=9))
* 根据现有的核酸序列数据可推得，自一月中旬以来 COVID-19 在大西雅图地区发生了一定程度的本地传播，但没有被检测到。([详细分析](https://nextstrain.org/narratives/ncov/sit-rep/zh/2020-03-05?n=10))
* 这项分析中包括的所有病毒个体核酸序列，可能有一个共同祖先，其出现时间范围在在2019年11月中旬至12月中旬之间。([详细分析](https://nextstrain.org/narratives/ncov/sit-rep/zh/2020-03-05?n=11))

```


<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [什么是冠状病毒？](https://nextstrain.org/ncov/2020-03-05)

### 更多资料:

* 维基百科上 COVID-19 疫情的总结 [中文维基](https://zh.wikipedia.org/wiki/2019%E5%86%A0%E7%8B%80%E7%97%85%E6%AF%92%E7%97%85%E7%96%AB%E6%83%85) _2020-01-30_
* 美国疾病预防控制中心提供的材料 [CDC](https://www.cdc.gov/coronavirus/2019-ncov/index-Chinese.html) _2020-01-29_
* 中国疾病预防控制中心关于[新型冠状病毒肺炎的专题网站](http://www.chinacdc.cn/jkzt/crb/zl/szkb_11803/)

<!-- This is right-side text -->
```auspiceMainDisplayMarkdown

## COVID-19 的一些其他信息

我们准备了一些值得一读的材料，让您熟悉 COVID-19 和导致其致病病毒 SARS-CoV-2。以下材料来源于 [Nextstrain](https://nextstrain.org) 网站的英文文档，我们建议您用谷歌翻译或者类似的机器翻译服务辅助浏览。
了解这些信息后，能使您更容易理解我们在这篇报告中的分析和结论。

（注：如果您想阅读关于如何解读系统进化树的中文版资料，我们在网上也找到了浅显易懂的[中文公开文档](http://www.360doc.com/content/17/1226/18/33459258_716501139.shtml)。）

<div>
  <a href="https://nextstrain.org/help/coronavirus/human-CoV"><img alt="microscopy image of coronaviruses" width="100" src="https://nextstrain.org/static/ncov_narrative-76cfd610d11ef708d213a3170de9519f.png"/> 冠状病毒的背景 </a>

  <a href="https://nextstrain.org/help/coronavirus/SARS-CoV-2"><img alt="illustration of a coronavirus" width="100" src="http://data.nextstrain.org/img_nCoV-CDC.jpg"/> 近期 COVID-19 疫情总结 </a>

  <a href="https://nextstrain.org/help/general/how-to-read-a-tree"><img alt="cartoon of a phylogenetic tree" width="100" src="http://data.nextstrain.org/img_toy_alignment_mini.png"/> 如何解读进化树 </a>

</div>

## Nextstrain 状况描述报告

以下页面包含使用 [Nextstrain](https://nextstrain.org) 执行的分析。
滚动左侧侧栏将显示文本段落的同时，在右侧会显示相应的基因组数据可视化图。

在疫情正在爆发的当下，这么迅速就获得新发现的大型RNA病毒的全基因组是一项了不起的成就。
而这些成果是由世界各地的科学家，快速公开共享基因组数据和分析而实现的（您可阅读本报告最后一页，我们列出了详细的序列提供单位）。


```


<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [常见问题和误区](https://nextstrain.org/ncov/2020-03-05)

### 更多详细资料：

_以下材料均来自英文网站，我们建议您使用机器翻译服务辅助阅读。_

* "不要相信你听到的关于冠状病毒和艾滋病的阴谋论" [文章链接](https://massivesci.com/notes/wuhan-coronavirus-ncov-sars-mers-hiv-human-immunodeficiency-virus/) _2020-01-31_

* "毫无根据的阴谋论声称新的冠状病毒是生物工程制造的" [文章链接](https://www.factcheck.org/2020/02/baseless-conspiracy-theories-claim-new-coronavirus-was-bioengineered/) _2020-02-07_

* "不，武汉冠状病毒不是经过基因改造植入艾滋病病毒片段的。" [文章链接](https://www.forbes.com/sites/victoriaforster/2020/02/02/no-coronavirus-was-not-bioengineered-to-put-pieces-of-hiv-in-it/#5d339e8e56cb) _2020-02-02_

* "冠状病毒解密" [文章链接](https://factcheck.afp.com/busting-coronavirus-myths) _2020-02-19_

<!-- This is right-side text -->
```auspiceMainDisplayMarkdown

## 常见问题和误区

_以下材料均来自英文网站，我们建议您使用机器翻译服务辅助阅读。_

### 常见问题

我们知道很多人对 COVID-19 有疑问，
我们特地建立了一个指南来回答一些最常被问到的问题。[FAQ](https://nextstrain.org/help/coronavirus/FAQ):

<div>

  <a href="https://nextstrain.org/help/coronavirus/FAQ"><img alt="picture of a question mark" width="100" src="http://data.nextstrain.org/img_question-mark.jpg"/> COVID-19 常见问题</a>

</div>


### 常见误区

目前社会上流行着许多关于新型冠状病毒起源的谣言。
但在像这样的疫情爆发期间，不正确信息的传播可能会引发更多的恐慌，并导致人们不信任科学家和政府，这意味着人们可能不遵循建议并采取适当的预防措施。

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

颜色表示各个国家内地区或美国各州的状态，x 轴则表示样本日期。

y 轴显示的是序列之间有怎样的关系，但请注意它没有测量单位。

<br>

样本日期信息虽然很有用，但它们并不总是准确地显示两个序列在基因上是如何相关的。

两个完全相同的序列可能具有不同的样本日期，因而在此图中，它们会看起来相距甚远。

<!-- There is NO right-side text -->


<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [系统进化树上的“差异程度”](https://nextstrain.org/ncov/2020-03-05?d=tree&m=div)

您可以从“时间”（Time）视图模式转换成“差异程度”（Divergence）视图模式，使 x 轴按“差异程度”排列。

<br>

由于中文没有明确的 Divergence 概念，这里译者提供一些辅助的维基资料帮助您理解系统进化树如何阅读。

网上公开的中文文档：[如何解读系统进化树](http://www.360doc.com/content/17/1226/18/33459258_716501139.shtml)

英文维基：[核酸差异程度](https://en.wikipedia.org/wiki/Genetic_divergence)

中文维基：[数学中的 Divergence 散度](https://zh.wikipedia.org/wiki/%E6%95%A3%E5%BA%A6)

<br>

请注意之前看起来关联性不大的序列，现在可能会排在一条竖直线上。
您可以通过上下滚动在上一张幻灯片和这张幻灯片之间切换，以查看系统进化树是如何变化的。

<br>

遗传差异程度是用基因组中的核酸变化（[自然突变](https://zh.wikipedia.org/wiki/%E7%AA%81%E5%8F%98)）的数量来衡量的。
有些序列可能完全没有突变或者没有核酸变化，这意味着它们都与系统进化树根部假定的祖先株序列上一致。
而其他病毒个体则有1到11个突变。

<br>

在不断变化的流行病爆发情况下，对新发现的大型RNA病毒的基因组进行测序是具有挑战性的一项任务。
因为在这些序列中观察到的核酸差异可能是由于测序误差，而不是实际的核酸突变。
特别是基因组末端的插入、缺失和其他差异有更高的概率是测序导致的，所以为了分析准确，我们提前筛选掉了部分序列差异。

<br>

在后面的页面中，我们会有时显示“时间”视图中的系统进化树，有时会展示“差异程度”视图中的系统进化树，这将取决于我们要展示的讨论和分析的结果。

<!-- There is NO right-side text -->


<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [系统进化树的解释](https://nextstrain.org/ncov/2020-03-05?d=tree)

我们目前有来自五个不同大陆的样本的序列。
虽然早期的病例都或多或少与武汉华南野生市场有着密切联系，但以我们现在观察到各种不同的病例，表明一些地区已经有了本地社区传播，另一些则是从中国以外的国家通过旅行者输入的。

<br>

通常，从动物身上集中地反复感染人类的情况下，病毒基因组将显示出显著的核酸碱基多样性（拉萨病毒、埃博拉病毒、MERS冠状病毒和禽流感都是如此）。

能观察到人类病例样本如此显著地归成系统进化树上的基因簇，基本可以用一次爆发模型来解释疫情的起源，即该疫情源于单一事件的人畜共传染，随后发展成人与人之间的种群内传播。

<!-- There is NO right-side text -->


<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [至少两次的输入感染导致了意大利的疫情爆发，两次都可能存在社区传播](https://nextstrain.org/ncov/2020-03-05?d=tree&f_country=Italy)

我们目前有三条来自意大利的序列，其中两条来自罗马地区，一条来自意大利北部的伦巴第。

<br>

这三条序列在疫情共同来源于一个疫情早期的共同祖先（即这个祖先株在图里的位置应在系统进化树根部附近，而且非常接近纵轴），这强烈地表明至少有两个病例输入到意大利，并通过社区传播遍了该地区。

<br>

努诺·法里亚博士团队非常详细地分析了巴西和其他全球序列，并得出“意大利北部的疫情很可能是多个人传入该地区的结果，而不是来自一个单一的来源”的结论。 [点击此处可以查看原文文章](http://virological.org/t/first-cases-of-coronavirus-disease-covid-19-in-brazil-south-america-2-genomes-3rd-march-2020/409)。



<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [意大利可能存在的隐藏传播路线](https://nextstrain.org/ncov/2020-03-05?d=tree&label=clade:A1a&m=div)

其中来自罗马（2020年1月29日）的两个序列联系密切，都有去过中国的旅行历史。

<br>

我们切换到“差异程度”视图，您可以发现这两个意大利序列是相同的，而其他邻近的序列（分别来自英国、巴西、瑞士、美国和中国）实际上与意大利序列相差2-4个突变。

<br>

该基因簇中巴西的序列（SPBR-02），已知有前往伦巴第米兰的旅行历史，报道称该簇中的瑞士病例最近也曾前往意大利。

我们对同簇中的美国病例的旅行史一无所知，但英国的病例（England/09c）已知是从中国直接输入的。

<br>

来自中国的这例英国序列，处于早期意大利罗马的病例和已知或可能有意大利旅行史的病例之间。

这意味着我们不应该假设较早的意大利序列和这些较新的意大利相关序列有直接的关系，而更可能是这些最近的意大利相关序列其实来源于意大利其他独立输入感染的事件。

<!-- There is NO right-side text -->



<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [从意大利伦巴第传播至全球其他地域的序列](https://nextstrain.org/ncov/2020-03-05?d=tree&label=clade:A2)

来自伦巴第（Italy/CDG1/2020）的序列与已知有意大利旅行史的序列被分在了一起（已知有可能在意大利受感染的病例分别来自墨西哥、德国、巴西和芬兰）。

<br>

图中最下方的德国的“BavPat1”序列是在疫情更早的时候从中国传入的一例病例。
它与该簇中的另一个序列的相似性很高（只有一个突变之分），暗示欧洲其他一些未被检测到的（“不明途径的”）疫情传播路线，可能起源于这个早期的德国病例。

<br>

当然，这也可能是两次单独输入欧洲的结果，因为来自其他地方的未被采集到的序列可能会分布在“BavPat1”和这个大基因簇之间。
但以目前的样本量，我们不能确认是哪种情况。

<!-- There is NO right-side text -->


<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [SARS-CoV-2 可能已在西雅图地区蔓延](https://nextstrain.org/ncov/2020-03-05?label=clade:B1%20&d=tree)

现在大西雅图地区和整个美国报告了许多例COVID-19病例。
新分离和测序得到的病例在核酸上与一月中旬同一地区分离到的一例病例密切相关。

<br>

对此，我们有两种可能的解释。
一种解释是，该病毒可能至少两次从中国的共同源头传入大西雅图地区。
而另一种解释是，病毒已经在该地区传播了一段时间，只是没有被发现。

<br>

特雷弗·贝德福德（Nextstrain的联合创始人）就这些可能性写了一篇很棒的博客文章，你可以在这里读到 [他的相关技术博客文章](https://bedford.io/blog/ncov-cryptic-transmission/)。

<br>

华盛顿最近的其他序列告诉我们另一件事：这些序列和来自大西雅图地区的序列，在核酸分析中被排列在了一起。
这强烈表明该地区已存在社区传播，且 SARS-CoV-2 病毒已经传播了一段时间。

<!-- There is NO right-side text -->



<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [测定最近共同祖先的时间](https://nextstrain.org/ncov/2020-03-05?label=clade:B1%20&d=tree)

一群个体序列的最近共同祖先时间（即tMRCA）能够提示这些个体共同的祖先株最近在何时出现。
如下图所示，这一时间可以早于病毒首次进入人类群体的时间，但也可以大大晚于此时间。

<div>
  <img alt="Example phylogeny where the time of the initial zoonosis is different from the most recent common ancestor of several sequenced cases" width="500" src="https://raw.githubusercontent.com/nicfel/nCov-Nicola/master/figures/zoonosis.png"/>
</div>


<!-- This is right-side text -->
```auspiceMainDisplayMarkdown

## 病毒的最近共同祖先时间

数个研究小组估算了最近共同祖先时间，您可以查看 [A Rambaut的报告（英文）](http://virological.org/t/phylodynamic-analysis-of-sars-cov-2-update-2020-03-06/420) 或 [T Stadler的报告（英文）](http://virological.org/t/evolutionary-epidemiological-analysis-of-93-genomes)。

所有序列的共同祖先最有可能出现在11月中旬至12月中旬之间。
这个时间点，与目前所有已测序的病毒个体，是从由武汉华南海鲜市场序列构成的最初基因簇中的，病毒群体衍生下来的说法一致[点击此处查看专业论坛的讨论帖](http://virological.org/t/phylodynamic-analysis-of-sars-cov-2-update-2020-03-06/420)。


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
