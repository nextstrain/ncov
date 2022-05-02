Setup and installation
======================

The following steps will prepare you to run complete analyses of SARS-CoV-2 data by installing required software and running a simple example workflow.

.. contents:: Table of Contents
   :local:

Install Nextstrain components
--------------------------------

:doc:`Follow instructions to install Nextstrain components <docs.nextstrain.org:install>`.

.. note::

   If using the :term:`native runtime <docs.nextstrain.org:runtime>`, install these additional packages necessary to run the ncov workflow. Make sure to activate the correct conda environment.

   .. code:: bash

      mamba install -c conda-forge -c bioconda \
        epiweeks nextclade nextalign pangolin pangolearn

Download the ncov workflow
-----------------------------

Use Git to download a copy of the ncov repository containing the workflow and this tutorial.

.. code:: bash

   git clone https://github.com/nextstrain/ncov.git
