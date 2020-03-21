---
title: Πώς να ερμηνεύσετε φυλογενετικά δέντρα
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
- Sotiris Salloumis
- akifoss
translatorLinks:
- https://github.com/codergr
- https://github.com/akifoss
date: "13 Μαρτίου 2020"
dataset: "https://nextstrain.org/ncov/2020-03-11?d=tree&legend=open&c=country"
abstract: "Αυτή η παρουσίαση εξηγεί πώς να διαβάζουμε και να ερμηνεύουμε τα φυλογενετικά δέντρα που ενημερώνουν την γονιδιωματική επιδημιολογία. Αυτός ο ιστότοπος έχει βελτιστοποιηθεί για προβολή σε προγράμματα περιήγησης για επιφάνεια εργασίας."
---
<!-- Translators: Only text after : in the above ^ needs to be translated -->
<!-- Comment tags like these do not need to be translated, they are only to help you! -->
<!-- Ensure that links always end in a 'letter' (. counts) If some kind of text doesn't follow them, it breaks the slide. -->
<!-- numbers can be tagged ilke this: 161</tag> - this is just for us to help find them to update! Just leave in the </tag> bit. -->

<!-- This is left-side text -->
# [Πίνακας περιεχομένων](https://nextstrain.org/ncov/2020-03-11?d=tree&legend=open&c=country)

* [Πώς συνδέονται τα δίκτυα μεταφοράς με φυλογενετικά δέντρα](https://nextstrain.org/narratives/trees-background/el?n=2);
* [Πώς διαβάζουμε ένα δένδρο](https://nextstrain.org/narratives/trees-background/el?n=3);
* [Πώς αναφέρεται το πάνελ "ποικιλομορφία" στο δέντρο](https://nextstrain.org/narratives/trees-background/el?n=4);
* [Μετρήσεις διαφορών με γενετική απόκλιση](https://nextstrain.org/narratives/trees-background/el?n=5).
* [Μέτρηση των διαφορών με την πάροδο του χρόνου](https://nextstrain.org/narratives/trees-background/el?n=6).
* [Πώς βρίσκουμε την ημερομηνία έναρξης μιας επιδημίας](https://nextstrain.org/narratives/trees-background/el?n=7);
* [Πώς πρέπει να ερμηνεύσω τα χαρακτηριστικά (χρώματα) στο δέντρο](https://nextstrain.org/narratives/trees-background/el?n=8);
* [Πώς σχετίζεται ο χάρτης με το δέντρο](https://nextstrain.org/narratives/trees-background/el?n=9);
* [Επιπλέον ανάγνωση: Aβεβαιότητα στα δέντρα](https://nextstrain.org/narratives/trees-background/el?n=10).
* [Σχετικά με το σύνολο δεδομένων](https://nextstrain.org/narratives/trees-background/el?n=11).

<!-- No right-side text -->


<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [Πώς συνδέονται τα δίκτυα μεταφοράς με φυλογενετικά δέντρα;](https://nextstrain.org/ncov/2020-03-11?d=tree&p=full)
Τα παθογόνα διαδίδονται μέσω ταχείας αναπαραγωγής σε ένα φορέα ακολουθούμενα από μετάδοση σε άλλο φορέα. Μια επιδημία μπορεί να αυξηθεί μόνο όταν μία μόλυνση έχει ως αποτέλεσμα περισσότερες από μία επόμενες μολύνσεις.
<br><br>
Καθώς το παθογόνο αναπαράγεται και εξαπλώνεται, το γονιδίωμα του πρέπει να αναπαραχθεί πολλές φορές και τυχαίες μεταλλάξεις (λάθη αντιγραφής) θα συσσωρεύονται στο γονιδίωμα. αυτό είναι φυσιολογικό. Τέτοιες τυχαίες μεταλλάξεις μπορούν να βοηθήσουν στην παρακολούθηση της εξάπλωσης του παθογόνου παράγοντα και να μάθουμε για τις διαδρομές μετάδοσης και τη δυναμική του.

<!-- This is right-side text -->
```auspiceMainDisplayMarkdown
# Παράδειγμα
<div width="50%" margin="auto">
<p>
<img width="500px" alt="σκίτσο το οποίο παρουσιάζει πως σχετίζονται το δένδρο μετάδοση με το φυλογενετικό δέντρο" src="https://github.com/nextstrain/nextstrain.org/raw/master/static-site/content/help/01-general/figures/infection_tree_combined.png"/>
</p>
<p>
Η παραπάνω εικόνα δείχνει ένα σκίτσο δέντρου μετάδοσης. Κάθε κύκλος αντιπροσωπεύει μια περίπτωση (μολυσμένο άτομο), με οριζόντιες γραμμές που υποδηλώνουν τη διάρκεια της μόλυνσης τους. Οι συνδεδεμένες περιπτώσεις αντιπροσωπεύουν μεταδόσεις από ένα άτομο στο άλλο.
<br><br>
Εδώ βλέπουμε την πλήρη εικόνα του δέντρου μετάδοσης. Στην πράξη, όμως, δειγματίζεται μόνο ένα υποσύνολο περιπτώσεων (μπλε). το δέντρο μετάδοσης είναι άγνωστο και συνήθως είναι διαθέσιμες μόνο οι ακαθάριστες εκτιμήσεις των αριθμών υποθέσεων. Οι αλληλουχίες γονιδιώματος μας επιτρέπουν να συναγάγουμε τμήματα του δέντρου μετάδοσης. Σε αυτό το παράδειγμα, τρεις μεταλλάξεις (μικρά διαμάντια) επισημαίνονται στο δέντρο. Οι αλληλουχίες που έχουν τις ίδιες μεταλλάξεις είναι πιο στενά συνδεδεμένες, έτσι ώστε αυτές οι μεταλλάξεις μας επιτρέπουν να ομαδοποιούμε δείγματα σε συστάδες στενά συγγενών ιών που ανήκουν στις ίδιες αλυσίδες μετάδοσης.
</p>
</div>
```


<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [Πώς διαβάζω ένα δέντρο;](https://nextstrain.org/ncov/2020-03-11)

Ο άξονας x ενός δέντρου αντιπροσωπεύει το βαθμό διαφοράς (σε χρόνο ή γενετική απόκλιση - θα το συζητήσουμε στο επόμενο). Ο άξονας y απλά βοηθάει να απλωθούν τα πράγματα έτσι ώστε να μπορούμε να δούμε τα πάντα. δεν διαθέτει μονάδες μέτρησης.
<br><br>
Οι άκρες του δέντρου αντιπροσωπεύουν δείγματα (δηλαδή, μπλε κουτιά από την τελευταία διαφάνεια). Οι εσωτερικοί κόμβοι αντιπροσωπεύουν περιπτώσεις που δεν ελήφθησαν δείγματα, αλλά πιστεύουμε ότι ήταν η πηγή όλων των περιπτώσεων που προέρχονται από αυτές (δηλαδή οι κόκκινες κόμβοι από την τελευταία διαφάνεια). Αυτές οι σχέσεις συνάγονται από την ανάλυση του μοτίβου των μεταλλάξεων που παρατηρήθηκαν στις περιπτώσεις δειγματοληψίας.

<!-- This is right-side text -->
```auspiceMainDisplayMarkdown
## Ένα παράδειγμα
<div width="50%" margin="auto">
<p>
<img width="700px" alt="Παράδειγμα φυλογένεσης όπου όλα ή μόνο ένα υποσύνολο περιπτώσεων περιλαμβάνονται στην τελική φυλογενία" src="https://github.com/nextstrain/nextstrain.org/raw/master/static-site/content/help/01-general/figures/toy_alignment_tree.png"/>
</p>
<p>
Πάνω, βλέπουμε μια εικόνα με ένα φυλογενετικό δέντρο στα αριστερά, όπου οι μεταλλάξεις εμφανίζονται ως έγχρωμοι κύκλοι. Δεξιά είναι οι αντίστοιχες ακολουθίες, επίσης με μεταλλάξεις που εμφανίζονται ως έγχρωμοι κύκλοι. Μπορούμε να δούμε ότι αλληλουχίες που μοιράζονται την ίδια ομάδα μεταλλάξεων μαζί. Όταν οι ακολουθίες εμφανίζονται συνδεδεμένες με μία επίπεδη κατακόρυφη γραμμή, όπως οι Α και Β, αυτό σημαίνει ότι δεν υπάρχουν διαφορές μεταξύ τους - οι ακολουθίες τους είναι ίδιες.
<br><br>
Όταν μια ακολουθία βρίσκεται σε μια μακρά γραμμή από μόνη της, όπως το C ή το Ε, αυτό σημαίνει ότι έχει μοναδικές μεταλλάξεις που δεν βρέθηκαν σε άλλες ακολουθίες. Όσο μεγαλύτερη είναι η γραμμή, τόσο περισσότερες μεταλλάξεις. Τα Α και Β έχουν επίσης μοναδικές μεταλλάξεις (τον πράσινο κύκλο) που δεν μοιράζονται οι άλλες αλληλουχίες, αλλά είναι ταυτόσημες μεταξύ τους.
<br><br>
Με βάση αυτό το δέντρο, θα καταλήγαμε στο συμπέρασμα ότι η A και B είναι στενά συνδεδεμένες μεταξύ τους και η D και E είναι στενά συνδεδεμένες μεταξύ τους. Τα A και B σχετίζονται στενότερα με το C από ότι στο D και E.
</p>

### Επιπλέον ανάγνωση  
* [Πώς να διαβάσετε ένα δέντρο: εκπαιδευτικό υλικό από Arctic Network](https://artic.network/how-to-read-a-tree.html).  
* [Πώς να διαβάσετε ένα δέντρο: video από Khan academy](https://www.khanacademy.org/science/high-school-biology/hs-evolution/hs-phylogeny/a/phylogenetic-trees).  

</div>

```


<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [Πώς αναφέρεται το πάνελ "ποικιλομορφία" στο δέντρο;](https://nextstrain.org/ncov/2020-03-11?d=tree,entropy&c=gt-ORF1b_314&legend=open)

Ας ρίξουμε μια ματιά στα πρώτα 169 </ tag> στελέχη του SARS-CoV-2 (τον ιό που προκαλεί το COVID-19) που έχουν κοινοποιηθεί δημόσια. Όπως και στην τελευταία σελίδα, δημιουργήσαμε μια ευθυγράμμιση αυτών των ακολουθιών ιών (μπορείτε να δείτε πώς όλες οι αναλύσεις που αναφέρθηκαν εδώ έγιναν [στό GitHub](https://github.com/nextstrain/ncov)).
<br><br>
Εδώ εμφανίζουμε το φυλογενετικό δέντρο πάνω από ένα διάγραμμα ράβδων που δείχνει τη μεταβολή (δηλ. Τις μεταλλάξεις) στο γονιδίωμα. Χωρίς αυτές τις μεταλλάξεις δεν μπορούσαμε να οικοδομήσουμε το δέντρο, έτσι ώστε οι δύο είναι στενά συνδεδεμένες.
<br><br>
Σε αυτό το πάνελ "ποικιλομορφίας", ο οριζόντιος άξονας είναι κάθε θέση στο ιικό γονιδίωμα (και όλες οι τριάντα χιλιάδες περίπου αυτών!).Ο κάθετος άξονας υποδεικνύει πόση μεταβλητότητα υπάρχει σε κάθε τοποθεσία.
<br><br>
Χρωματίσαμε το δέντρο σύμφωνα με μία από αυτές τις μεταλλάξεις - στην περίπτωση αυτή κωδικόνιο 314 στο γονίδιο "ORF1b".
Δεν υπάρχει κανένας λόγος a priori να πιστεύουμε ότι αυτή η μετάλλαξη είναι μια λειτουργική μετάλλαξη (δηλ. Προσδίδοντας οποιαδήποτε βιολογική αλλαγή). Είναι ακριβώς μεταλλάξεις όπως αυτό που χρησιμοποιούμε για να καθορίσουμε τις σχέσεις μεταξύ των αλληλουχιών και να κατασκευάσουμε το δέντρο.

<!-- There is NO right-side text -->


<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [Μετρήσεις διαφορών με γενετική απόκλιση](https://nextstrain.org/ncov/2020-03-11?c=num_date&d=tree&m=div)
Αυτή είναι μια φυλογένεση των πρώτων 169 </tag> στελεχών του SARS-CoV-2 (ο ιός που προκαλεί το COVID-19) που έχουν κοινοποιηθεί δημόσια.
<br><br>
Εδώ, ο οριζόντιος άξονας υποδεικνύει την απόκλιση, που είναι ο αριθμός αλλαγών (μεταλλάξεων) στο γονιδίωμα, σε σχέση με τη ρίζα του δέντρου (δηλ. Την έναρξη της εστίας).
Ορισμένες αλληλουχίες μπορεί να έχουν μηδενικές μεταλλάξεις - δηλαδή όλες είναι ίδιες με τη ρίζα (κέντρο) του δέντρου.
Άλλοι ιοί έχουν μεταξύ ενός και έντεκα μεταλλάξεις.
<br><br>
Προς το παρόν, αυτό δεν μπορεί να μοιάζει πολύ με ένα «δέντρο». Πολλές από τις ακολουθίες είναι ίδιες - κάθονται μαζί σε κατακόρυφες γραμμές όπως οι Α και Β (μερικές βρίσκονται στο αριστερό μέρος του δέντρου).
Άλλοι έχουν μοναδικές ή κοινές μεταλλάξεις και έτσι κάθονται στις γραμμές ή «κλάδους», πηγαίνοντας δεξιά.
Μπορείτε να δείτε πόσες μεταλλάξεις έχει ένα υποκατάστημα τοποθετώντας το ποντίκι πάνω από αυτό.

<!-- There is NO right-side text -->

<!-- ############ SLIDE BREAK ############# -->
# [Μέτρηση των διαφορών με την πάροδο του χρόνου](https://nextstrain.org/ncov/2020-03-11?c=num_date&d=tree&legend=open)
Μπορούμε επίσης να απεικονίσουμε πώς ο ιός έχει εξαπλωθεί με την πάροδο του χρόνου χρησιμοποιώντας την ημερομηνία δειγματοληψίας ως τον άξονα x.
<br><br>Εδώ, ο άξονας x αντιπροσωπεύει την ημερομηνία δειγματοληψίας για κάθε ιό. Οι θέσεις των συμβουλών αντικατοπτρίζουν την ημερομηνία λήψης αυτών των δειγμάτων. Οι ημερομηνίες των εσωτερικών κόμβων - οι "περιπτώσεις που λείπουν" - συνάγονται βάσει του χρόνου δειγματοληψίας των απογόνων τους και του ρυθμού με τον οποίο μεταλλάσσεται ο ιός.
<br><br>
Παρατηρήστε πόσες αλληλουχίες που προηγουμένως κάθισαν σε μια γραμμή (που υποδεικνύουν ταυτόσημα γονιδιώματα) τώρα διαδίδονται διαχρονικά.Αυτό συμβαίνει όταν ο ρυθμός με τον οποίο μεταλλάσσεται ο ιός είναι ελαφρώς βραδύτερος από τον ρυθμό με τον οποίο εξαπλώνεται.Μπορείτε να μετακινηθείτε προς τα πάνω και προς τα κάτω μεταξύ της προηγούμενης διαφάνειας και αυτής, για να δείτε πώς αλλάζει το δέντρο.
<!-- There is NO right-side text -->


<!-- ############ SLIDE BREAK ############# -->
# [Έναρξη της επιδημίας](https://nextstrain.org/ncov/2020-03-11?c=num_date&d=tree&legend=open)

Μπορούμε επίσης να χρησιμοποιήσουμε τη γονιδιωματική για να ορίσουμε μια ημερομηνία όταν ξεκίνησε ένα ξέσπασμα, ακόμα κι αν αυτό ήταν πριν συνειδητοποιήσουμε ότι συνέβαινε.
<br><br>
Επειδή μπορούμε να ορίσουμε ημερομηνίες σε κάθε δείγμα και κόμβο στο δέντρο, μπορούμε να το χρησιμοποιήσουμε για να συμπεράνουμε την ημερομηνία της "ρίζας" του δέντρου. Αυτό αντιπροσωπεύει τον "πιο πρόσφατο κοινό πρόγονο" όλων των ακολουθιών SARS-CoV-2 που έχουμε μέχρι στιγμής. Π.χ. οι παππούδες σας είναι οι "πιο πρόσφατοι κοινοί πρόγονοι" σας και όλων των πρώτων ξαδέρφων σας.
<br><br>
Αν τοποθετήσετε το κέρσορα του ποντικιου  πάνω από την αριστερή κάθετη γραμμή, μπορείτε να δείτε ότι η υπογεγραμμένη ημερομηνία έναρξης είναι μεταξύ των μέσων Νοεμβρίου και τα μέσα Δεκεμβρίου του 2019 για αυτή τη συγκεκριμένη εστία.

<!-- There is NO right-side text -->


<!-- ############ SLIDE BREAK ############# -->
<!-- This is left-side text -->
# [Πώς πρέπει να ερμηνεύσω τα χαρακτηριστικά (χρώματα) στο δέντρο;](https://nextstrain.org/ncov/2020-03-11)
Φυλογενετικά δέντρα συχνά περιέχουν πρόσθετες πληροφορίες, όπως η θέση κάθε συλλογής δειγμάτων.
<br><br>
Από αυτό, μπορούμε να συμπεράνουμε τις θέσεις των εσωτερικών κόμβων (υποθετικές ενδιάμεσες, μη δειγματοληπτικές περιπτώσεις) χρησιμοποιώντας μαθηματικά μοντέλα. Αυτό μπορεί να μας βοηθήσει να καταλάβουμε πώς μετακινείται ο ιός από τη μια θέση στην άλλη.
<br><br>
Η ερμηνεία αυτών θα πρέπει, ωστόσο, να γίνεται με προσοχή, καθώς η δειγματοληψία και η αλληλουχία ή η έλλειψή τους μπορεί να επηρεάσει σημαντικά την ερμηνεία.

<!-- This is right-side text -->
```auspiceMainDisplayMarkdown
# Ενα παράδειγμα
<div width="50%" margin="auto">
<p>
<img width="700px" alt="Εικόνα που δείχνει πώς ερμηνεύεται η δειγματοληπτική επίδραση της εξάπλωσης του ιού" src="https://github.com/nextstrain/nextstrain.org/raw/master/static-site/content/help/01-general/figures/introductions.png"/>
</p>
<p>
Αριστερά, παρουσιάζουμε ένα πλήρως δειγματοληπτό φυλογενετικό δέντρο, με δείγματα από δύο διαφορετικές τοποθεσίες που σημειώνονται με πορτοκαλί και μπλε χρώμα. Καθώς περπατάμε στο δέντρο παρατηρούμε τρεις περιπτώσεις όπου το χρώμα (θέση) αλλάζει από πορτοκαλί σε μπλε. Από αυτό, καταλήγουμε στο συμπέρασμα ότι υπήρχαν τρεις διαφορετικές εισαγωγές από την πορτοκαλή τοποθεσία έως την μπλε θέση.

<br><br>

Αλλά αυτή η ερμηνεία βασίζεται στη δειγματοληψία: στο μεσαίο δέντρο, έχουμε αφαιρέσει ένα πορτοκαλί δείγμα. Παρατηρούμε τώρα μόνο ένα διακόπτη από πορτοκαλί σε μπλε, υποδηλώνοντας ότι υπήρχε μόνο μία εισαγωγή στο μπλε που συνέβη πολύ νωρίτερα.
<br><br>
Στο τελευταίο παράδειγμα, έχουμε μόνο μία ακολουθία από πορτοκαλί, που θα μπορούσε να μας οδηγήσει στο να σκεφτούμε ότι υπήρχε μία εισαγωγή από πορτοκαλί σε μπλε.

<br><br>

Έτσι, ενώ αυτά τα συμπεράσματα μπορούν να είναι ανεκτίμητα, πρέπει επίσης να ερμηνεύονται με προσοχή.
</p>
```
<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
<!-- ############ SLIDE BREAK ############# -->
<!-- This is left-side text -->
# [Πώς σχετίζεται ο χάρτης με το δέντρο;](https://nextstrain.org/ncov/2020-03-11?d=tree,map&legend=closed)

Εδώ, δείχνουμε το δέντρο χρωματισμένο από τη θέση κάθε δείγματος (και τη συναγόμενη θέση για κάθε εσωτερικό κόμβο).
Εάν κάνετε κλικ στο ['Εξερευνήστε τα δεδομένα'](https://nextstrain.org/ncov), μπορείτε να παίξετε ένα κινούμενο σχέδιο της υποτιθέμενης εξάπλωσης του ιού κατά τη διάρκεια της επιδημίας.


<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [Προηγμένη ανάγνωση: αβεβαιότητα στα δέντρα](https://nextstrain.org/ncov/2020-03-11)
Νωρίτερα, μιλήσαμε για τον τρόπο με τον οποίο οι εσωτερικοί κόμβοι αντιπροσωπεύουν _υποθετικές_ περιπτώσεις χωρίς δείγμα. Στην πραγματικότητα, όλα τα δέντρα αντιπροσωπεύουν _υποθέσεις_ για τον τρόπο με τον οποίο ένα παθογόνο έχει εξελιχθεί και μετακινηθεί με την πάροδο του χρόνου. Τα δέντρα που παρουσιάζουμε στο Nextstrain είναι εκτιμήσεις σημείων -- δηλαδή την έκδοση αυτής της ιστορίας που μεγιστοποιεί την πιθανότητα παρατήρησης των δεδομένων που κάνουμε.
<br><br>
Ωστόσο, στις εκτιμήσεις αυτές υπάρχει πάντα αβεβαιότητα. Σε γενικές γραμμές, τμήματα του δέντρου που είναι πυκνά δειγματοληψμένα είναι πιο σίγουρα. οι περιοχές που έχουν αραιά δειγματοληψία είναι λιγότερο σίγουρες.

```auspiceMainDisplayMarkdown
# Μια απεικόνιση
<div width="50%" margin="auto">
<p>
<img width="700px" alt="Απεικόνιση της αβεβαιότητας που είναι εγγενής στην ανακατασκευή δέντρων" src="https://github.com/nextstrain/nextstrain.org/raw/c69bfd0750c284ff12f33682f8d82848e13d9e15/static-site/content/help/01-general/figures/hcov_densitree.png"/>
</p>
</div>
```

<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text -->
# [Επιστημονική πίστωση](https://nextstrain.org/ncov/2020-03-05?d=map&c=author)

Θα θέλαμε να αναγνωρίσουμε την εκπληκτική και έγκαιρη εργασία όλων των επιστημόνων που εμπλέκονται σε αυτό το ξέσπασμα, αλλά κυρίως εκείνων που εργάζονται στην Κίνα.
Μόνο μέσω της ταχείας ανταλλαγής γονιδιωματικών δεδομένων και μεταδεδομένων είναι δυνατές αναλύσεις όπως αυτές.

<br><br>

Ευχαριστούμε επίσης το [GISAID] (https://gisaid.org) για την παροχή της πλατφόρμας μέσω της οποίας μπορούν να μεταφορτωθούν και να μοιραστούν αυτά τα δεδομένα.

<!-- Do not need to translate insitutions names -->
<!-- This is right-side text -->
```auspiceMainDisplayMarkdown

Είμαστε ευγνώμονες για τα δεδομένα που συγκεντρώθηκαν από αυτά τα εργαστήρια που προέρχονταν:

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
# [Λεπτομερής επιστημονική πίστωση](https://nextstrain.org/ncov/2020-03-05?d=map&c=author)

Αυτά τα δεδομένα μοιράστηκαν μέσω[GISAID](https://gisaid.org).
Ευχαριστούμε πολύ για τη συμβολή τους.

<br><br>

Στα δεξιά δίνουμε συγκεκριμένες ακολουθίες που μοιράζονται σε κάθε εργαστήριο.

<!-- This is right-side text -->
```auspiceMainDisplayMarkdown

Τα γονιδιώματα SARS-CoV-2 μοιράστηκαν γενναιόδωρα από επιστήμονες σε αυτά τα εργαστήρια υποβολής:

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
