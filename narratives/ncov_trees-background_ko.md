---
title: 계통수를 읽는 법
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
  - Kyo Bin Kang
  - Minkyu Kim
  - Taehoon Ha
date: "2020년 3월 13일"
dataset: "https://nextstrain.org/ncov/2020-03-11?d=tree&legend=open&c=country"
abstract: "이 문서는 유전체 역학적 정보를 제공하는 계통수를 어떻게 읽고 해석하는지에 대해 다룹니다. 이 웹사이트는 데스크톱 브라우저에 최적화되어 있습니다."
---
<!-- Translators: Only text after : in the above ^ needs to be translated -->
<!-- Comment tags like these do not need to be translated, they are only to help you! -->
<!-- Ensure that links always end in a 'letter' (. counts) If some kind of text doesn't follow them, it breaks the slide. -->
<!-- numbers can be tagged ilke this: 161</tag> - this is just for us to help find them to update! Just leave in the </tag> bit. -->

<!-- This is left-side text -->
# [목차](https://nextstrain.org/ncov/2020-03-11?d=tree&legend=open&c=country)

* [전염병의 전파 경로는 계통수와 어떤 관련이 있는가](https://nextstrain.org/narratives/trees-background/ko?n=2)?  
* [계통수는 어떻게 읽는가](https://nextstrain.org/narratives/trees-background/ko?n=3)?  
* ["다양성" 패널은 계통수와 어떤 관련이 있는가](https://nextstrain.org/narratives/trees-background/ko?n=4)?   
* [유전적 발산 사이의 차이 측정하기](https://nextstrain.org/narratives/trees-background/ko?n=5).  
* [시간에 따른 차이 측정하기](https://nextstrain.org/narratives/trees-background/ko?n=6).  
* [바이러스 발생 시점 추정하기](https://nextstrain.org/narratives/trees-background/ko?n=7).  
* [계통수의 (색깔로 구분된) 특징들을 어떻게 해석해야 합니까](https://nextstrain.org/narratives/trees-background/ko?n=8)?  
* [지도와 계통수는 어떤 관련이 있는가](https://nextstrain.org/narratives/trees-background/ko?n=9)?  
* [추가 읽기자료: 계통수의 통계적 불확실성](https://nextstrain.org/narratives/trees-background/ko?n=10).  
* [데이터셋에 대하여](https://nextstrain.org/narratives/trees-background/ko?n=11).  

<!-- No right-side text -->


<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [전염병의 전파 경로는 계통수와 어떤 관련이 있는가?](https://nextstrain.org/ncov/2020-03-11?d=tree&p=full)
병원체는 숙주의 체내에서 빠르게 증식하여 다른 숙주로 전염됩니다. 전염병은 한 번의 감염이 한 번 이상 뒤따르는 다른 감염으로 이어질 경우에만 유행하게 됩니다.
<br><br>
병원체가 증식하고 전파되는 과정에서 병원체의 유전체는 여러 번 복제되고, 이 복제 과정에서 자연적으로 발생하는 오류로 인한 무작위 돌연변이가 유전체에 누적됩니다. 이것은 정상적인 일입니다. 이런 무작위 돌연변이를 통해 병원체의 확산을 추적하고 전파 경로와 역학에 대해 알아낼 수 있습니다.

<!-- This is right-side text -->
```auspiceMainDisplayMarkdown
# 예시
<div width="50%" margin="auto">
<p>
<img width="500px" alt="cartoon showing how transmission tree and phylogenetic tree relate" src="https://github.com/nextstrain/nextstrain.org/raw/master/static-site/content/help/01-general/figures/infection_tree_combined.png"/>
</p>
<p>
위의 그림은 간략화된 전파 계통수의 예시입니다. 각각의 원은 감염된 환자를 의미하며, 가로 직선은 환자들의 감염 기간을 의미합니다. 연결된 원들은 한 환자로부터 다른 환자로 전염이 이뤄진 경우를 가리킵니다.
<br> <br>
지금 이 그림에서 우리는 병원체 전파 경로의 전체를 보고 있습니다. 그러나 현실에서는 모든 전파 사건이 아니라 몇몇 환자(파란색으로 표시)의 경우에 대해서만 감염체 표본를 얻을 수 있습니다. 따라서 실제 전파 경로는 추적하기가 어려우며, 감염자 숫자에 대한 대략적인 추정만이 가능합니다. 하지만 유전체 염기서열을 분석하면 전파 경로를 짐작할 수 있습니다. 이 예시의 경우, 세 번의 돌연변이(작은 다이아몬드로 표시)가 전파 과정에서 일어났습니다. 같은 돌연변이를 가지고 있는 염기 서열은 상대적으로 서로에게 더 밀접한 관련이 있을 것이며, 따라서 이런 유전체 서열 상의 돌연변이를 기준으로 같은 전파 연쇄에 속하는 밀접하게 관련있는 바이러스들을 각각의 군집들로 묶을 수 있습니다.
</p>
</div>
```


<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [계통수는 어떻게 해석하는가?](https://nextstrain.org/ncov/2020-03-11)

계통수의 x축은 변화의 정도(시간 혹은 유전적 분화 - 이후에 자세히 설명합니다)를 의미합니다. y축은 개체들을 퍼뜨려 전체적 양상을 볼 수 있게 해주는 역할을 하며, 어떠한 단위나 값을 가지지 않습니다.
<br><br>
계통수의 가지 끝 부분은 표본(이전 슬라이드에서 파란색 원과 같은 경우)을 의미합니다. 가지 안쪽의 분기점은 표본을 채집하지 못했지만 이후에 파생된 모든 감염의 원천으로 추정되는 감염사례(이전 슬라이드에서 빨간색 원과 같은 경우)를 의미합니다. 이런 관계는 채집된 병원체 표본의 돌연변이 패턴을 분석하여 유추할 수 있습니다.

<!-- This is right-side text -->
```auspiceMainDisplayMarkdown
## 예시
<div width="50%" margin="auto">
<p>
<img width="700px" alt="Example phylogeny where all or only a subset of cases are included in the final phylogeny" src="https://github.com/nextstrain/nextstrain.org/raw/master/static-site/content/help/01-general/figures/toy_alignment_tree.png"/>
</p>
<p>
위의 그림 좌측에는 계통수가 그려져 있으며, 돌연변이가 색이 칠해진 동그라미로 표시되어 있습니다. 우측에는 각 표본에 해당하는 유전체 서열이 표시되어 있으며, 돌연변이가 마찬가지의 방법으로 표시되어 있습니다. 그림에서 공통의 돌연변이를 가지고 있는 표본들이 같은 그룹으로 묶여 있는 것을 볼 수 있습니다. 두 개의 서열이 A와 B의 경우처럼 세로 직선으로 연결되어 있다면, 이는 두 표본의 서열이 완전히 동일하다는 것을 의미합니다.
<br><br>
만약 어떤 서열이 C와 E의 경우처럼 한 직선에 하나만 연결되어 있다면, 이는 해당 표본에서 다른 표본에서는 발견되지 않는 고유의 돌연변이가 발견되었음을 의미합니다. 직선의 길이가 길수록, 더 많은 돌연변이가 존재합니다.
A와 B 또한 다른 서열에서 발견되지 않는 고유의 돌연변이(녹색 원)을 가지고 있지만, 이들 두 표본의 서열은 서로 동일합니다.
<br><br>
이 계통수를 바탕으로, 우리는 A와 B가, 그리고 D와 E가 각각 매우 가까운 관계임을 알 수 있습니다. 그리고 A와 B는 D와 E보다는 C에게 상대적으로 더 가깝습니다.
</p>

### Further reading  
* [계통수를 읽는 법: Arctic Network의 설명서](https://artic.network/how-to-read-a-tree.html) (영어).  
* [계통수를 읽는 법: Khan academy의 영상](https://www.khanacademy.org/science/high-school-biology/hs-evolution/hs-phylogeny/a/phylogenetic-trees) (영어).  

</div>

```


<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# ["다양성" 패널은 계통수와 어떤 관련이 있는가?](https://nextstrain.org/ncov/2020-03-11?d=tree,entropy&c=gt-ORF1b_314&legend=open)

169개의 공개된 SARS-CoV-2 (COVID-19를 유발하는 바이러스) 계통을 살펴 보도록 하겠습니다. 이전 슬라이드에서 보신 것처럼, 우리는 바이러스의 염기서열을 정리해 보았습니다. (여기서 언급된 모든 분석이 어떻게 이루어졌는지는 [Github에서](https://github.com/nextstrain/ncov) 확인할 수 있습니다.
<br><br>
오른쪽에 계통수(위)와 염기서열 상 변이(돌연변이)를 보여주는 막대 그래프(아래)를 함께 보여주고 있습니다. 만약 이런 돌연변이들이 없었다면, 계통수를 그릴 수 없었을 것입니다. 따라서 다양성 패널과 계통수는 서로 밀접한 연관성을 가지고 있습니다.
<br><br>
이 다양성 패널(아래 그림)에서 가로축은(다 합쳐서 약 3만 개 정도 되는) 바이러스 염기서열 위치를 의미합니다. 반면 세로축은 얼마나 많은 유전적 변이가 있었는지를 보여 줍니다.
<br><br>
우리는 이러 돌연변이들 중 하나를 선택하여 이에 따라 나무의 색깔을 다르게 표시했습니다. 여기서 선택된 것은 "ORF1b" 안에 들어 있는 Codon 314입니다. 이 돌연변이가 기능적 변이(즉, 생물학적 변화를 주는 것)라고 여길 마땅한 근거는 아직 없습니다. 보통 이런 종류의 돌연변이는 서열과 서열 사이의 관계를 정의하고, 계통수를 구성하는 데 사용됩니다.

<!-- There is NO right-side text -->


<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [유전적 발산 사이의 차이 측정하기](https://nextstrain.org/ncov/2020-03-11?c=num_date&d=tree&m=div)
이것은 169개의 공개된 SARS-CoV-2(COVID-19를 유발하는 바이러스) 계통을 표시한 계통수입니다.
<br><br>
가로 축은 유전적 분화를 나타내며, 이는 계통수의 시작점(즉, 바이러스의 최초 발생 시점)에 대한 유전적 변이(돌연변이)의 횟수를 의미합니다. 일부 서열은 돌연변이가 존재하지 않아 계통수의 뿌리(중앙의 시작점)의 서열과 완전히 동일합니다. 다른 바이러스들은 한 개에서 최대 열한 개의 돌연변이를 보입니다.
<br><br>
현재로서는 이 그림이 '나무'처럼 보이지 않을 수 있습니다. 많은 서열들이 동일하기 때문인데, 여기서 동일하다는 의미는 A와 B처럼 동일 수직선 상에 위치한 것을 의미합니다. (일부는 나무의 가장 왼쪽 지점에 위치해 있습니다.) 바이러스들이 몇 개의 서로 같은 돌연변이를 공유하면서도 각자의 고유한 돌연변이를 가지고 있을 경우, '가지'를 형성하며 나무를 오른쪽으로 뻗어 나가게 합니다. 마우스 커서를 가지 위에 올려서 각 가지들이 얼마나 많은 돌연변이들을 가지고 있는지 확인할 수 있습니다.


<!-- There is NO right-side text -->

<!-- ############ SLIDE BREAK ############# -->
# [시간에 따른 차이 측정하기](https://nextstrain.org/ncov/2020-03-11?c=num_date&d=tree&legend=open)
시료가 측정된 날짜에 대한 정보를 x축에 표시하여 바이러스가 어떻게 퍼져 나갔는지를 시간에 따라 시각화할 수 있습니다. 이 경우 나무의 끝부분이 위치한 x축 상의 위치가 각 바이러스 표본들이 채집된 날짜를 나타냅니다. 나무의 가지 안쪽 분기점("결측 케이스")의 날짜는 그 하위의 바이러스 표본들이 채집된 시점과 바이러스가 변이되는 속도를 기준으로 추정합니다.
<br><br>
이전에는 동일 선상에 위치했던(동일한 염기서열을 가졌던) 서열들이 시간이 지남에 따라 얼마나 다양하게 갈라지며 분화되었는지를 확인해 보십시오.
이런 현상은 바이러스가 변이되는 속도가 확산되는 속도보다 약간 느릴 때 발생합니다. 이전 슬라이드와 본 슬라이드를 위아래로 스크롤하며 나무가 어떻게 변하는지 확인해 보십시오.

<!-- There is NO right-side text -->


<!-- ############ SLIDE BREAK ############# -->
# [바이러스 발생 시점 추정하기](https://nextstrain.org/ncov/2020-03-11?c=num_date&d=tree&legend=open)
유전체학적 지식을 이용하면 언제 바이러스가 발생했는지를, 심지어 우리가 바이러스에 대해 인지하기 이전 시점이라 하더라도, 특정할 수 있습니다. 계통수 상의 각 표본과 내부 분기점에 날짜를 특정할 수 있으므로, 이를 활용하여 나무의 '뿌리' 날짜를 추정할 수 있습니다. 나무의 '뿌리'란 우리가 현재 가지고 있는 모든 SARS-CoV-2 서열들의 "가장 가까운 공통의 조상"을 의미합니다. 이는 마치 당신의 조부모가 당신과 당신의 모든 사촌들의 "가장 가까운 공통 조상"인 것과 비슷한 의미를 갖습니다.
<br><br>
오른쪽 그림에서 가장 왼쪽에 위치한 수직선 위에 마우스 커서를 놓으면, 바이러스 발생 추정 시기가 2019년 11월 중순에서 12월 중순 사이로 추정됨을 알 수 있습니다.

<!-- There is NO right-side text -->


<!-- ############ SLIDE BREAK ############# -->
<!-- This is left-side text -->
# [계통수에 (색깔로 구분된) 특징들을 어떻게 해석해야 합니까?](https://nextstrain.org/ncov/2020-03-11)
종종 계통수는 각각의 표본 수집 지역 같은 추가적인 정보를 담고 있습니다. 이러한 정보를 바탕으로, 우리는 수학적 모델을 활용하여 내부 분기점들(표본 채집이 되지 않았지만 존재헀을 것으로 추정되는 바이러스 감염 사례)의 위치를 추정할 수 있습니다. 우리는 이를 통해 바이러스가 한 지역에서 다른 지역으로 어떻게 이동하고 있는지 이해하는 데에 도움을 줍니다.
<br><br>
하지만 이런 사항을 해석할 때는 항상 주의를 기울여야 합니다. 표본 채집이나 유전체 데이터가 충분하지 않은 경우 결과가 왜곡될 수 있기 때문입니다.

<!-- This is right-side text -->
```auspiceMainDisplayMarkdown
# 예시
<div width="50%" margin="auto">
<p>
<img width="700px" alt="Illustration showing how sampling effects interpretation of viral spread" src="https://github.com/nextstrain/nextstrain.org/raw/master/static-site/content/help/01-general/figures/introductions.png"/>
</p>
<p>
왼쪽 그림은 모든 표본이 완전히 채집된 가상의 계통수이며, 주황색과 파란색으로 각각 서로 다른 두 지역을 구분하였습니다. 나무를 따라 내려가면 색상(지역)이 주황색에서 파란색으로 바뀌는 세 지점을 관찰할 수 있습니다. 이를 통해 우리는 바이러스가 주황색 지역에서 파란색 지역으로 전해질 때 세 건의 전파 사례가 있었다고 결론내릴 수 있을 것입니다.
<br><br>
하지만, 표본 채집 여부에 따라 이런 분석 결과가 바뀔 수 있습니다. 중앙 그림의 나무를 보면 왼쪽 나무에서 주황색 표본 하나를 제거한 것을 확인 할 수 있습니다. 이제는 주황색 지역에서 파랑색 지역으로 전파되는 오직 하나의 경우만이 관찰됩니다. 그리고 이 전파 사례는 현재보다 더 이전 시점에 발생한 것으로 추정되고 있습니다.
<br><br>
오른쪽 그림은 주황색 지역에서 오직 하나의 표본만 가지고 있을 경우의 계통수입니다. 이 경우 우리는 파란색에서 주황색 지역으로 전파되는 한 번의 전파만이 있었다고 해석하게 됩니다.
<br><br>
따라서 이런 추론법은 매우 유용하지만, 반드시 주의를 기울여 해석해야 합니다.
</p>
```
<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
<!-- ############ SLIDE BREAK ############# -->
<!-- This is left-side text -->
# [지도와 계통수는 어떤 관련이 있는가?](https://nextstrain.org/ncov/2020-03-11?d=tree,map&legend=closed)

여기서 보여지는 계통수는 각 표본이 채집된 위치(및 내부 분기점의 추정 위치)에 따라 다른 색으로 표현 되었습니다. ['데이터 탐색'](https://nextstrain.org/ncov)을 클릭하면 바이러스가 확산되는 과정을 추정한 애니메이션을 재생할 수 있습니다.


<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [추가 읽기자료: 계통수의 통계적 불확실성](https://nextstrain.org/ncov/2020-03-11)

앞에서 우리는 계통수의 내부 분기점이 실제 존재하는 표본이 아닌 가상으로 추정한 케이스를 의미한다고 이야기했습니다. 사실 모든 계통수는 병원체가 어떻게 진화하고 이동했는지에 대한 문제에서 어디까지나 가설일 뿐입니다. Nextstrain에서 보여지는 모든 나무들은 점 추정(역주: point estimate. 모수를 하나의 값으로 추정하는 방식)을 사용합니다. 다시 말해, 여기서 보여주는 계통수는 존재할 수 있는 여러 버전의 계통수들 중 우리가 현재 가진 데이터를 관찰할 가능성이 최대인 버전이라고 할 수 있습니다.
<br><br>
그러나, 이러한 추정에는 언제나 불확실성이 존재합니다. 일반적으로 나무에서 데이터가 풍부한 부분은 상대적으로 예측이 더 명확하지만, 데이터가 적은 부분은 상대적으로 예측이 니다.


<!-- There is NO right-side text -->
```auspiceMainDisplayMarkdown
# 계통수의 통계적 불확실성에 대한 개념도
<div width="50%" margin="auto">
<p>
<img width="700px" alt="Illustration of the uncertainty inherent in tree reconstruction" src="https://github.com/nextstrain/nextstrain.org/raw/c69bfd0750c284ff12f33682f8d82848e13d9e15/static-site/content/help/01-general/figures/hcov_densitree.png"/>
</p>
</div>
```

<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [과학적 기여자들](https://nextstrain.org/ncov/2020-03-05?d=map&c=author)

이번 유행에 관여하고 있는 모든 과학자들, 특히 중국의 과학자들의, 작업들이 대단하고 신속했다는 점에 감사를 표하고자 합니다. 이 분석들은 유전체 정보 및 메타정보의 빠른 공유가 아니었다면 불가능했을 것입니다.

<br><br>

또한, 이 데이터들이 업로드되고 공유할 수 있는 플랫폼을 제공한 [GISAID](https://gisaid.org)에게 감사의 뜻을 전합니다.

<!-- Do not need to translate insitutions names -->
<!-- This is right-side text -->
```auspiceMainDisplayMarkdown

우리는 다음 실험실들이 제공한 데이터들에 대하여 감사의 말씀을 전합니다:

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
# [과학적 기여: 상세 내역](https://nextstrain.org/ncov/2020-03-05?d=map&c=author)

이 데이터는 [GISAID](https://gisaid.org)를 통해 공유되었습니다. 그들에게 감사의 뜻을 전합니다.


<br><br>

오른쪽에는 각 염기서열을 제공한 각각의 연구실들이 나열되어 있습니다.

<!-- This is right-side text -->
```auspiceMainDisplayMarkdown

SARS-CoV-2 유전체는 아래 연구실의 과학자들에 의해 제공되었습니다:

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
