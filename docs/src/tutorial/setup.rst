Setup and installation
======================

The following steps will prepare you to run complete analyses of SARS-CoV-2 data by installing required software and running a simple example workflow.

.. contents:: Table of Contents
   :local:

1. Install Nextstrain components
--------------------------------

Follow instructions to install Nextstrain components :doc:`here <docs.nextstrain.org:install>`.

.. note::

   If using the :term:`native runtime <docs.nextstrain.org:runtime>`, install these additional packages necessary to run the ncov workflow. Make sure the correct conda environment is activated.

   .. code:: bash

      mamba install -c conda-forge -c bioconda epiweeks nextclade nextalign pangolin pangolearn

2. Download the ncov workflow
-----------------------------

Download the workflow
~~~~~~~~~~~~~~~~~~~~~

Use Git to download a copy of the ncov repository containing the workflow and this tutorial.

.. code:: bash

   git clone https://github.com/nextstrain/ncov.git
