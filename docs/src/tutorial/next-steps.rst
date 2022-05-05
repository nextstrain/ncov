==========
Next steps
==========

Congratulations! You have completed all of the tutorials for the ncov workflow. Read on for some next steps.

.. contents:: Table of Contents
   :local:

.. _create-analysis-directory:

Create your own analysis directory
==================================

On a web browser:

1. `Sign up for a GitHub account <https://github.com/signup>`__ if you do not already have one.
2. Create a repository from the ``ncov-tutorial`` template repository:

   1. Go to https://github.com/nextstrain/ncov-tutorial.
   2. Click **Use this template**.
   3. Give your repository a name. We recommend ``my-ncov-analyses`` and will use that name in the following steps.
   4. Click **Create repository from template**.

In a command prompt:

1. Go to the ``ncov/`` directory.
2. Clone your new repository, replacing ``<username>`` with your own username:

   .. code:: text

      git clone https://github.com/<username>/my-ncov-analyses

3. Read the next section to learn how to modify ``genomic-surveillance.yaml``.

Modify the genomic surveillance workflow configuration
======================================================

Instead of an Idaho-focused workflow config, you can provide your own data for the ``custom_data`` input. Follow the same steps in the tutorial for GISAID download but select your own set of sequences and rename your ``metadata.tsv`` and ``sequences.fasta`` files accordingly.

   .. note::

      Workflow run time increases with the number of sequences, and the GISAID web interface has a maximum of 5,000 sequences per download.

Then, use the following steps to customize names, titles, and context:

1. Change the ``custom_data`` input filenames from ``idaho.metadata.tsv`` and ``idaho.sequences.fasta`` to your own.
2. Change the regional input dataset from North America to an appropriate region for your custom focal data. :doc:`See the complete list of available URLs <../reference/remote_inputs>`.
3. Rename the output dataset from ``idaho`` to your own. Note the name restrictions.
4. Reword the output dataset title to your own.
5. Rename the subsampling scheme from ``idaho_scheme`` to your own. Note the name restrictions.
6. For each sample, increase the ``max_sequences`` to your own.
7. Rename the ``usa_context`` sample and update the ``query`` accordingly.

.. warning::

   File paths in the :term:`config files <config file>` must start with the :term:`analysis directory`. For example, in the tutorial:

   .. code:: yaml

      auspice_config: ncov-tutorial/auspice-config-custom-data.json

   Now that you have created your own analysis directory, this must be modified, e.g.

   .. code:: yaml

      auspice_config: my-ncov-analyses/auspice-config-custom-data.json

Additional resources
====================

- Learn more about genomic epidemiology:

   - `An applied genomic epidemiological handbook <https://alliblk.github.io/genepi-book/intro.html>`__ by Allison Black and Gytis Dudas
   - `Genomic Epidemiology Seminar Series <https://czgenepi.org/resources>`__ by Chan Zuckerberg Initiative Genomic Epidemiology (CZ GEN EPI)
   - `COVID-19 Genomic Epidemiology Toolkit <https://www.cdc.gov/amd/training/covid-19-gen-epi-toolkit.html>`__ by Centers for Disease Control and Prevention (CDC)

- :doc:`Review all possible options to configure your SARS-CoV-2 analyses with Nextstrain <../reference/workflow-config-file>`.
- Watch `this 1-hour video overview <https://youtu.be/m4_F2tG58Pc>`__ by Heather Blankenship on how to deploy Nextstrain for a Public Health lab.
