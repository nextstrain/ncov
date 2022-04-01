Clade Naming & Definitions
==========================

The nomenclature used by Nextstrain to designate clades for SARS-CoV-2 is driven by the following objectives:

-  label genetically well defined clades that have reached significant frequency and geographic spread,
-  allow for transient clade designations that are elevated to major clades if they persist and rise in frequency,
-  provide memorable but informative names,
-  gracefully handle clade naming in the upcoming years as SARS-CoV-2 becomes a seasonal virus.

.. contents:: Table of Contents
   :local:

Major clades
------------

Definition
~~~~~~~~~~

We name a new major clade when it reaches a frequency of 20% globally at any time point. When calculating these frequencies, care has to be taken to achieve approximately even sampling of sequences in time and space since sequencing effort varies strongly between countries. A clade name consists of the year it emerged and the next available letter in the alphabet. A new clade should be at least 2 mutations away from its parent major clade.

Naming
~~~~~~

We name major clades by the year they are estimated to have emerged and a letter, e.g. 19A, 19B, 20A. The yearly reset of letters will ensure that we don't progress too far into the alphabet, while the year-prefix provides immediate context on the origin of the clade that will become increasingly important going forward. These are meant as major genetic groupings and not intended to completely resolve genetic diversity.

The hierarchical structure of clades is sometimes of interest. Here, the "derivation" of a major clade can be labeled with the familiar "." notation as in 19A.20A.20C for the major clade 20C.

Subclades
---------

Within these major clades, we subclades, which we will label by their parent clade and the nucleotide mutation(s) that defines them (ex: 19A/28688C). It should be noted however, that these mutations are only meaningful in that they define the clade. Once a subclade reaches (soft) criteria on frequency, spread, and genetic distinctiveness, it will be renamed to a major clade (hypothetically 19A/28688C to 20D).

Current Clades
--------------

+-----------------+--------------------------------------------+-------------------------+-------------------------+
| Clade           | Primary Countries                          | Mutations               | Max Frequency           |
+=================+============================================+=========================+=========================+
| 19A             | Asia: China/Thailand                       | Root clade              | 65-47% Globally in Jan  |
+-----------------+--------------------------------------------+-------------------------+-------------------------+
| 19B             | Asia: China                                | C8782T T28144C          | 28-33% Globally in Jan  |
+-----------------+--------------------------------------------+-------------------------+-------------------------+
| 20A             | N America/Europe/Asia: USA, Belgium, India | C14408T A23403G         | 41-46% Globally Apr-May |
+-----------------+--------------------------------------------+-------------------------+-------------------------+
| 20B             | Europe: UK, Belgium, Sweden                | G28881A G28882A G28883C | 19-20% Globally Mar-Apr |
+-----------------+--------------------------------------------+-------------------------+-------------------------+
| 20C             | N America: USA                             | C1059T G25563T          | 19-21% Globally Apr     |
+-----------------+--------------------------------------------+-------------------------+-------------------------+

You can view the current clades on the Global SARS-CoV-2 Nextstrain tree `here <https://nextstrain.org/ncov/global?branchLabel=clade&c=clade_membership>`__.

Identifying Nextstrain Clades
-----------------------------

To make it easy for users to identify the Nextstrain clade of their own sequences, we provide a clade assignment tool at `clades.nextstrain.org <https://clades.nextstrain.org/>`__. In addition to assigning clades, this tool will call mutations in your sequences relative to the reference and performs some basic QC.

You can also use the `simple python script <https://github.com/nextstrain/ncov/blob/master/scripts/assign_clades.py>`__ to assign appropriate clades to sequences in a fasta file. This script is part of the ``ncov`` GitHub repository, but does not require running any other part of the workflow. However, ``augur`` :doc:`must be installed <augur:installation/installation>` to run the script.

Note when running this script you can supply ``--sequences`` if your sequences require aligning first. If you already have aligned your sequences to the ``ncov`` repository reference (for example, from running the repository), you can supply ``--alignment``. If you supply sequences that are not aligned to the ``ncov`` reference, you may get bad results!
