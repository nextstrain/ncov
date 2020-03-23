---
title: كيف تقرأ وتحلل شجرة التصنيف الجيني
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
  - Ahmad Saleh
  - Fadi Alnaji
translatorLinks:
  - https://twitter.com/Ahmadsa96201545
  - https://twitter.com/biofadialnaji
date: "2020 March 13"
dataset: "https://nextstrain.org/ncov/2020-03-11?d=tree&legend=open&c=country"
abstract: " يفسر هذا التقرير كيفية قراءة وتحليل اشجار التصنيف الجيني والتي تغذي علم الوبائيات الجينومي بالمعلومات. هذا الموقع الالكتروني مهيأ للعرض على متصفحات اجهزة الحاسوب المكتبية"
---
<!-- Translators: Only text after : in the above ^ needs to be translated -->
<!-- Comment tags like these do not need to be translated, they are only to help you! -->
<!-- Ensure that links always end in a 'letter' (. counts) If some kind of text doesn't follow them, it breaks the slide. -->
<!-- numbers can be tagged ilke this: 161</tag> - this is just for us to help find them to update! Just leave in the </tag> bit. -->

<!-- This is left-side text -->
# [جدول المحتويات](https://nextstrain.org/ncov/2020-03-11?d=tree&legend=open&c=country)

* [كيف ترتبط شبكات الانتشار بأشجار التصنيف الجيني؟](https://nextstrain.org/narratives/trees-background/ar?n=2).
* [كيف تقرأ الشجرة؟](https://nextstrain.org/narratives/trees-background/ar?n=3).
* [كيف يرتبط اطار "التنوع" بالشجرة؟](https://nextstrain.org/narratives/trees-background/ar?n=4).
* [قياس الاختلافات مع الافتراق الجيني؟](https://nextstrain.org/narratives/trees-background/ar?n=5).  
* [قياس الاختلافات مع الزمن](https://nextstrain.org/narratives/trees-background/ar?n=6).  
* [كيف تحدد موعد بداية الانفجار الوبائي؟](https://nextstrain.org/narratives/trees-background/ar?n=7).
* [كيف احلل الصفات (الالوان) على الشجرة؟](https://nextstrain.org/narratives/trees-background/ar?n=8).
* [كيف ترتبط الخريطة بالشجرة؟](https://nextstrain.org/narratives/trees-background/ar?n=9).
* [قراءات متقدمة: عدم اليقين بالاشجار](https://nextstrain.org/narratives/trees-background/ar?n=10).  
* [عن مجموعة البيانات](https://nextstrain.org/narratives/trees-background/ar?n=11).  

<!-- No right-side text -->


<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [كيف ترتبط شبكات الانتشار باشجار التصنيف الجيني؟](https://nextstrain.org/ncov/2020-03-11?d=tree&p=full)
تنتشر مسببات الأمراض من خلال تكاثرها السريع في عائل اول ثم تنتقل منه لعائل آخر. يتصاعد الإنتشار الوبائي فقط عندما تؤدي حالة اصابة واحدة الى اكثر من حالة عدوى بشكل متتالي.
<br><br>
كلما تكاثر مسبب المرض وانتشر فإن جينومه سيتناسخ عدة مرات مما يؤدي الى تراكم العديد من الطفرات (أخطاء أثناء عملية نسخ الجينوم) وهذا أمر طبيعي. مثل هذه الطفرات العشوائية يمكن ان تساعد في تتبع انتشار مسبب المرض وتعلمنا عن طرق وحركية الانتشار.


<!-- This is right-side text --
```auspiceMainDisplayMarkdown
# مثال
<div width="50%" margin="auto">
<p>
<img width="500px" alt="cartoon showing how transmission tree and phylogenetic tree relate" src="https://github.com/nextstrain/nextstrain.org/raw/master/static-site/content/help/01-general/figures/infection_tree_combined.png"/>
</p>
يظهر الشكل الموضح أعلاه عرض لشجرة انتشار العدوى. حيث أن كل دائرة تمثل حالة (شخص مصاب), والخطوط الافقية تشير الى طول فترة اصاباتهم. الحالات المترابطة تشير الى انتشار العدوى من شخص الى الآخر. يظهر هنا صورة كاملة لشجرة انتشار العدوى. عمليا، مع ذلك، فقط بعض الحالات تم جمعها (بالازرق) والاخرى تقديريه؛ ولذلك شجرة الانتشار ليست أكيدة بالمطلق، وما هو متوفر فقط تقدير تقريبي عن عدد الحالات.  تسلسل الجينومات يجعلنا قادرين على تحليل اجزاء من شجرة الانتشار. في هذا المثال، هناك ثلاث طفرات ( على شكل معينات صغيرة) مشار لها على الشجرة. التسلسلات التي تمتلك نفس الطفرات تعتبر متشابهة وبالتالي فإن هذه الطفرات تتيح لنا وضع العينات في  مجموعات مختلفه، كل مجموعة تحتوي سلاسل الفيروسات الأكثر تقاربا للحد الذي يمكن اعتبارها بأنها تنتمي لنفس بؤرة الانتشار.
</p>
</div>

```
<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [كيف تقرأ الشجرة؟](https://nextstrain.org/ncov/2020-03-11)

المحور السيني للشجرة يمثل درجة الاختلاف (في الزمن او التباعد الجيني-- سنشرح  ذلك لاحقاً). اما المحور الصادي فهو فقط من اجل المساعدة على فصل الاشياء بعيداً عن بعضها لنرى كل شيء بوضوح، ولذلك ليس له أي وحدات قياس.
<br><br>
نهايات فروع الشجرة تمثل العينات (وهذا يعني الحالات الزرقاء في الشريحة السابقة). العقد الداخلية للشجرة تمثل الحالات التي لم يتم جمعها، ولكن نعتقد انها مصدر جميع الحالات التي انحدرت منها ( وهذا يعني العقد الحمراء في الشريحة السابقة). هذا العلاقات تم ربطها من خلال تحليل انماط الطفرات التي تم ملاحظتها في العينات التي تم جمعها من الحالات المصابة.


<!-- This is right-side text -->
```auspiceMainDisplayMarkdown
## مثال
<div width="50%" margin="auto">
<p>
<img width="700px" alt="Example phylogeny where all or only a subset of cases are included in the final phylogeny" src="https://github.com/nextstrain/nextstrain.org/raw/master/static-site/content/help/01-general/figures/toy_alignment_tree.png"/>
</p>
<p>
نرى في الشكل الموضح اعلاه شجرة تصنيف جيني على اليسار، حيث ان الطفرات مشار لها بدوائر ملونة. علي اليمين يوجد السلاسل الجينية المقابلة لتلك الطفرات على اليسار. هنا أيضا يشار للطفرات بدوائر ملونة ايضا. نستطيع مشاهدة ان التسلسلات التي تشترك بنفس الطفرات تتجمع معا (A و B). عندما تظهر التسلسلات مرتبطة معا بخط عمودي فهذا يعني انه لا يوجد هناك فروقات فيما بينهم-- اي سلاسلهم الجينية متطابقة.
<br><br>
عند وجود سلسلة على خط طويل بمفردها (C أو E) فهذا يعني أنه يمتلك طفرة فريدة غير موجودة في السلال الأخرى.
<br><br>
كلما كان طول الخط اكبر كلما كانت الطفرات أكثر.السلاسل A و B أيضاً تمتلك طفرات (الدوائر الخضراء) غير مشتركة مع تسلسلات أخرى، لكنها متطابقة مع بعضها.
<br><br>
اعتمادا على هذه الشجرة نستطيع الإستنتاج بأن  A و B مرتبطتان جدا ببعضهما وان D و E مرتبطتان جدا ببعضهما في حين ان  A و B  اكثر ارتباطا ب C من ارتباطهما ب D و E.

</p>

### قراءات اضافية
* [How to read a tree: tutorial from Arctic Network](https://artic.network/how-to-read-a-tree.html).  
* [How to read a tree: video from Khan academy](https://www.khanacademy.org/science/high-school-biology/hs-evolution/hs-phylogeny/a/phylogenetic-trees).  

</div>

```


<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [كيف يرتبط اطار "التنوع" بالشجرة؟](https://nextstrain.org/ncov/2020-03-11?d=tree,entropy&c=gt-ORF1b_314&legend=open)

لنلق نظرة على أول ١٦٩ سلالة لفيروس سارس كو في ٢ (الفيروس المسبب ل كوفيد ١٩) والتي تم مشاركتها للعموم. كما المثال في الصفحة السابقة، قمنا بمحاذاة هذه التسلسلات الفيروسية لمقارنتها (يمكنك مشاهدة كيف تم عمل كل هذه التحليلات المذكورة هنا على موقع [on GitHub](https://github.com/nextstrain/ncov)).
<br><br>
نقوم هنا بعرض شجرة التصنيف الجيني فوق مخطط يبين التغيرات (اي الطفرات) في الجينوم (أسفل الصفحة). بدون هذه الطفرات لا يمكننا بناء الشجرة، لذلك فإن كلاهما مرتبط بالاخر.
<br><br>
في اطار "التنوع" هذا، المحور الافقي يمثل كل نيوكليوتيد على جينوم الفيروس. المحور العمودي يشير الى كمية الاختلاف (الطفرات) على كل نيوكليوتيد.
<br><br>
قمنا بتلوين الشجرة بناءاً على على واحدة من هذه الطفرات -- في هذه الحالة استخدمنا كودون رقم ٣١٤ في "اطار القراءة ١ ب". لا يوجد سبب للظن بأن هناك افضلية لان تكون هذه الطفرة طفرة وظيفية (اي تسبب اي تغير وظيفي). انها بالضبط طفرات كهذه هي التي نستخدمها لتحديد العلاقات بيت التسلسلات ونبني منها الشجرة.

<!-- There is NO right-side text -->


<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [قياس الاختلافات مع الافتراق الجيني](https://nextstrain.org/ncov/2020-03-11?c=num_date&d=tree&m=div)
هذه شجرة تصنيف جيني لاول ١٦٩ سلالة لفيروسات سارس كو في ٢ (االفيروسات المسببة ل كوفيد ١٩) والتي تم مشاركتها للعموم.
<br><br>
هنا يشير المحور الافقي للاختلاف، والذي هو عدد التغيرات (الطفرات) في الجينوم مقارنة بجذر الشجرة ( اي بداية تفشي الوباء) بعض التسلسلات يمكن ان لا تمتلك اي طفرات-- وهذا يعني انها متطابقة مع جذر (مركز) الشجرة. فيروسات اخرى تمتلك ما بين طفرة الى احد عشرة طفرة.
<br><br>
في هذه اللحظة قد لا تبدو بانها تشبه "شجرة". العديد من التسلسلات متطابقة-ولذلك جاءت فوق بعضها بخطوط عمودية مثل A و B ( بعضها على  أقصى يسار الشجرة). بعضها الاخر يمتلك طفرات فريدة او مشتركة ولذلك تتواجد على خطوط او 'فروع' تتجه نحو اليمين. يمكنك مشاهدة كم عدد الطفرات على الفرع من خلال تمرير الفأرة فوقها.

<!-- There is NO right-side text -->

<!-- ############ SLIDE BREAK ############# -->
# [قياس الاختلافات مع الزمن](https://nextstrain.org/ncov/2020-03-11?c=num_date&d=tree&legend=open)
نستطيع ايضا رؤية كيف انتشرت الفيروسات مع الزمن من خلال استخدام تاريخ جمع العينات كمحور سيني. هنا يمثل المحور السيني تاريخ جمع العينة لكل فيروس. مواقع الأطراف تبين تاريخ اخذ تلك العينات. تواريخ العقد الداخلية- - "الحالات المفقودة"-- يتم استنتاجها بناءا على متى تم جمع العينات المنحدرة منها ومعدل حدوث الطفرات للفيروس.
<br><br>
لاحظ كم عدد التسلسلات التي جاءت مسبقا فوق بعضها (مشيرة للجينومات المتطابقة) والتي أصبحت متباعدة مع الزمن.
هذا يحدث عندما يكون معدل حدوث الطفرات في الفيروس اقل بقليل من من معدل انتشاره. يمكنك التنقل بين الشريحة  السابقة والشريحة الحالية لرؤية كيف تتغير الشجرة.
<!-- There is NO right-side text -->


<!-- ############ SLIDE BREAK ############# -->
# [تحديد موعد بدء تفشي الوباء](https://nextstrain.org/ncov/2020-03-11?c=num_date&d=tree&legend=open)

نستطيع أيضاً استخدام علم الجينومات من اجل تحديد متى بدأ تفشي الوباء حتى لو حدث ذلك قبل ان ندرك حدوثه. لاننا نستطيع تحديد تواريخ كل عينة وكل عقدة في الشجرة فإننا نستطيع استخدام ذاك لاستنتاج تاريخ جذر الشجرة. وهذا يمثل " السلف الاحدث والاكثر تشابها" لكل تسلسلات فيروسات سارس كو في ٢ التي نملكها حتى الان. تماما مثل جدك أو جدتك " السلف الاحدث" لك ولكل ابناء عمومتك.
<br><br>
اذا قمت بتحريك الفأرة فوق الخط العمودي في اقصى اليسار يمكنك مشاهدة تاريخ البداية الذي تم استنتاجه لوقوع تفشي الوباء ؛ وهو ما بين وسط نوفمبر ووسط ديسمبر ٢٠١٩.

<!-- There is NO right-side text -->


<!-- ############ SLIDE BREAK ############# -->
<!-- This is left-side text -->
# [كيف احلل الصفات (الالوان) على الشجرة؟](https://nextstrain.org/ncov/2020-03-11)
تحتوي شجرة التصنيف الجيني عادة على معلومات إضافية مثل موقع جمع كل عينة. من هنا نستطيع استنتاج مواقع العقد الداخلية (وسيط مفترض، حالات لم تجمع) باستخدام نماذج رياضية. يمكن ان بساعدنا هذا على فهم كيف يتنقل الفيروس من موقع الى الاخر.
<br><br>
عمليات الاستنتاج هذه يجب ان يتم عملها بحذر حيث ان ندرة عمليات الجمع او تحديد التسلسلات يمكن ان تؤثر بشكل ملحوظ على هذه الاستنتاجات.


<!-- This is right-side text -->
```auspiceMainDisplayMarkdown
# مثال
<div width="50%" margin="auto">
<p>
<img width="700px" alt="Illustration showing how sampling effects interpretation of viral spread" src="https://github.com/nextstrain/nextstrain.org/raw/master/static-site/content/help/01-general/figures/introductions.png"/>
</p>
<p>
نعرض على اليسار شجرة تصنيف جيني لمجموعة كاملة من العينات، مع عينات من موقعين مختلفين موضحة باللونين البرتقالي والازرق. كلما انتقلنا في الشجرة من اعلى لاسفل نلاحظ تغير اللون من البرتفالي للازرق في ثلاث مواقع منفصلة. نستنتج من ذلك بانه كان هناك ثلاثة حالات ادخال مختلفة للفيروس من المنطقة البرتقالية الى المنطقة الزرقاء. لكن هذا الاستنتاج يعتمد على جمع العينات: في الشجرة الوسطى قمنا بازالة عينة برتقالية واحدة. نلاحظ الان تبدل واحد فقط من البرتقالي للازرق، مقترحا بانه كان هناك حالة ادخال واحدة الى المنطقة الزرقاء وحدثت قبل الوقت المقترح في الشجرة الأولى (اقصى اليسار). في المثال الاخير، لدينا تسلسل واحد فقط من المنطقة البرتقالية ليبدو وكأن هناك ادخال واحد من المنطقة البرتقالية للمنطقة الزرقاء. لذلك، وعلى الرغم من أهمية هذه الاستنتاجات الا انه يجب توخي الحذر في تحليلها .
</p>

```
<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
<!-- ############ SLIDE BREAK ############# -->
<!-- This is left-side text -->
# [كيف ترتبط الخارطة بالشجرة؟](https://nextstrain.org/ncov/2020-03-11?d=tree,map&legend=closed)
نعرض هنا الشجرة ملونة بموقع كل عينة (والموقع المستنتج لكل عقدة داخلية)
اذا ضغطت هنا ['استعرض البيانات'](https://nextstrain.org/ncov), يمكنك التلاعب في البيانات المرئية  لكيفية الانتشار المستنتج للفيروس خلال فترة الوباء.



<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [قراءة متقدمة: عدم اليقين بالاشجار](https://nextstrain.org/ncov/2020-03-11)
تحدثنا سابقاً عن كيف ان العقد الداخلية تمثل عينات حالات-افتراضية لم تجمع. في الحقيقة، كل الاشجار تقدم - فرضيات - عن كيفية تطور وانتقال مسبب المرض. الاشجار التي نعرضها على موقع نيكست سترين هي تقديرات دقيقة-- وهذا يعني وجود احتمالية عالية لحدوث هذه التقديرات
<br><br>
ومع ذلك فهناك دائما عدم يقين في هذه التقديرات. وبشكل عام فإن نسبة اليقين تكون اعلى في  اجزاء الشجرة التي تمثل المناطق التي تم جمع العينات منها بكثافة.  في حين ان نسبة اليقين تقل في المناطق التي كان جمع العينات فيها قليلا.

```auspiceMainDisplayMarkdown
# شرح توضيحي
<div width="50%" margin="auto">
<p>
<img width="700px" alt="Illustration of the uncertainty inherent in tree reconstruction" src="https://github.com/nextstrain/nextstrain.org/raw/c69bfd0750c284ff12f33682f8d82848e13d9e15/static-site/content/help/01-general/figures/hcov_densitree.png"/>
</p>
</div>
```

<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [الحقوق العلمية](https://nextstrain.org/ncov/2020-03-05?d=map&c=author)

نود هنا ان نقدر العمل السريع والرائع من قبل جميع العلماء الذين  عملوا على هذا الوباء، ولكن نخص اولئك العاملين في الصين.
ما كان لمثل هذه التحليلات ان تكون لولا المشاركة السريعة لبيانات الجينومات.

<br><br>

عظيم شكرنا موجه الى [GISAID](https://gisaid.org) على تزويدنا بمنصتهم والتي تم من خلالها رفع وتحليل البيانات التي تم مشاركتها.

<!-- Do not need to translate insitutions names -->
<!-- This is right-side text -->
```auspiceMainDisplayMarkdown

نحن ممتنون للبيانات الاصلية التي تم جمعها من قبل المختبرات التالية:

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
# [الحقوق العلمية التفصيلية](https://nextstrain.org/ncov/2020-03-05?d=map&c=author)

تم مشاركة هذه البيانات من خلال[GISAID](https://gisaid.org).
نقدر ونثمن عاليا مساهماتهم

<br><br>

على اليمين نقدم التسلسلات المحددة والتي تم مشاركتها من كل مختبر

<!-- This is right-side text -->
```auspiceMainDisplayMarkdown

تسلسلات جينومات سارس كو في ٢ تم مشاركتها بجهود كريمة مم علماء في المختبرات التالية:

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
