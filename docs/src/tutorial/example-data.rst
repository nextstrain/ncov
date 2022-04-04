Run the workflow using example data
===================================

The aim of this first tutorial is to introduce our SARS-CoV-2 workflow.
To do this, we will run the workflow using a small set of reference data which we provide.
This tutorial leads on to subsequent tutorials where we will walk through more complex scenarios.

.. contents:: Table of Contents
   :local:

Prerequisites
-------------

1. :doc:`setup`. These instructions will install all of the software you need to complete this tutorial and others.

Setup
-----

1. Activate the ``nextstrain`` conda environment:

   .. code:: text

      conda activate nextstrain

2. Change directory to the ``ncov`` directory:

   .. code:: text

      cd ncov

3. Download the example tutorial repository into a new directory ``my-analyses/``:

   .. code:: text

      git clone https://github.com/nextstrain/ncov-tutorial my-analyses

Run the workflow
----------------

From within the ``ncov/`` directory, run the ``ncov`` workflow using a configuration file provided in the tutorial directory:

.. code:: text

   nextstrain build . --cores all --configfile my-analyses/example-data.yaml

Break down the command
~~~~~~~~~~~~~~~~~~~~~~

The workflow can take several minutes to run. While it is running, you can learn about the parts of this command:

- ``nextstrain build .``
   - This tells the :term:`docs.nextstrain.org:Nextstrain CLI` to :term:`build <docs.nextstrain.org:build (verb)>` the workflow from ``.``, the current directory. All subsequent command-line parameters are passed to the workflow manager, Snakemake.
- ``--cores all``
   - This required Snakemake parameter specifies the number of CPU cores to use (`more info <https://snakemake.readthedocs.io/en/stable/executing/cli.html>`_).
- ``--configfile my-analyses/example-data.yaml``
   - ``--configfile`` is another Snakemake parameter used to configure the ncov workflow.
   - ``my-analyses/example-data.yaml`` is a YAML file which provides custom workflow configuration including inputs and outputs. Contents with comments excluded:

      .. code-block:: yaml

         inputs:
           - name: reference_data
             metadata: https://data.nextstrain.org/files/ncov/open/reference/metadata.tsv.xz
             sequences: https://data.nextstrain.org/files/ncov/open/reference/sequences.fasta.xz

         refine:
           root: "Wuhan-Hu-1/2019"

      This provides the workflow with one input named ``reference_data``, which is a small dataset maintained by the Nextstrain team. The metadata and sequences files are downloaded directly from the associated URLs. `See the complete list of SARS-CoV-2 datasets we provide through data.nextstrain.org <https://docs.nextstrain.org/projects/ncov/en/latest/reference/remote_inputs.html>`_.

      The ``refine`` entry specifies the root sequence for the example GenBank data.

      For more information, visit `the complete configuration guide <../reference/configuration.html>`_.

The workflow output produces a new directory ``auspice/`` containing a file ``ncov_default-build.json``, which will be visualized in the following section. The workflow also produces intermediate files in a new ``results/`` directory.

Visualize the results
---------------------

Run this command to start the :term:`docs.nextstrain.org:Auspice` server, providing ``auspice/`` as the directory containing output dataset files:

.. code:: text

   nextstrain view auspice/

Navigate to ``http://127.0.0.1:4000/ncov/default-build``. The resulting :term:`docs.nextstrain.org:dataset` should show a phylogeny of ~200 sequences:

.. figure:: ../images/dataset-example-data.png
   :alt: Phylogenetic tree from the "example data" tutorial as visualized in Auspice

To stop the server, press :kbd:`Control-C` on your keyboard.

.. note::

   You can also view the results by dragging the file ``auspice/ncov_default-build.json`` onto `auspice.us <https://auspice.us>`__.
