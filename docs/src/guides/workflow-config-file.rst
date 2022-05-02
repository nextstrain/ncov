Workflow config file guide
==========================

This is a guide for common use cases of a :term:`workflow config file <config file>`. For a detailed reference, see :doc:`../reference/workflow-config-file`.

.. contents:: Table of Contents
   :local:

Changing parameters
-------------------

You can configure most steps of the workflow by specifying values in a workflow config file. We've provided reasonable defaults in ``defaults/parameters.yaml``; these are the same values the Nextstrain team uses for our analyses.

We encourage you to take a few minutes to skim through `the default workflow config file <https://github.com/nextstrain/ncov/blob/master/defaults/parameters.yaml>`__. Although these default values should be fine for most users, it's helpful to get a sense for what options are available.

If you'd like to tweak the parameterization, **you can override any of these values by specifying them in the workflow config file. Any values not overridden in this way will fall back to the default values.** Keeping build-specific parameters separate this way prevents mixups of settings between runs, and gives you a cleaner file to work with (rather than having to wrestle the *entire* default workflow config file).

Adding custom rules
-------------------

Insert your own custom Snakemake rules into the workflow without modifying the existing Snakemake files, by defining a list of ``custom_rules`` in your workflow config file. Each entry in the ``custom_rules`` list should be a path to a valid Snakemake file (e.g., ``my_rules.smk``). The workflow will detect these custom rules and include them after all other rules have been defined.

To modify rules, create a new :term:`customization file` named ``my_rules.smk`` and add the ``custom_rules`` section in a workflow config file:

.. code:: yaml

   custom_rules:
     - my-ncov-analyses/my_rules.smk

As an example, the Nextstrain team's workflow defines custom export rules that modify the default Auspice JSONs. These rules are defined in the workflow config file as follows:

.. code:: yaml

   custom_rules:
     - workflow/snakemake_rules/export_for_nextstrain.smk

Adding a new place
------------------

Places are defined as one of:

- ``region`` (e.g., ``North America``, ``Asia``)
- ``country``
- ``division`` (i.e., state, province, or canton)
- ``location`` (i.e., a county or city within a division)

To define a new place, you'll need to specify its GPS coordinates and a color.

1. Add a line to ``defaults/lat_longs.tsv``. This file is separated into sections for each geographic resolution. This looks like:

::

   # resolution    place   latitude    longitude
   location    Abondant    48.790785   1.420178

..

   Note: keep in mind that ``0.0`` longitude is the prime meridian; to specify something in the Western hemisphere, you'll need to enter a *negative* value for longitude. Similarly, to specify something in the Southern hemisphere, you'll need to enter a *negative* value for latitude

2. Add an entry to ``color_ordering.tsv`` such that your newly-defined place is next to geographically nearby places in the list.

Subsampling
-----------

Basic subsampling
~~~~~~~~~~~~~~~~~

Reasonable defaults are pre-defined. You can find a :ref:`list of them here <configuration-subsampling>`.

Custom subsampling schemes
~~~~~~~~~~~~~~~~~~~~~~~~~~

We implement hierarchical subsampling by producing multiple samples at different geographic scales and merge these samples into one file for further analysis. A build can specify any number of such samples which can be flexibly restricted to particular meta data fields and subsampled from groups with particular properties. When specifying subsampling in this way, we'll first take sequences from the 'focal' area, and the select samples from other geographical areas. Read further for information on how we select these samples. Here, we'll look at `the advanced example config file <https://github.com/nextstrain/ncov-tutorial/blob/main/examples/example_advanced_customization/builds.yaml>`__ to explain some of the options.

When specifying how many sequences you want in a subsampling level (for example, from a country or a region), you can do this using either ``seq_per_group`` or ``max_sequences`` - these work with the ``group_by`` argument. For example, ``switzerland`` subsampling rules in the advanced example looks like this:

.. code:: yaml

   switzerland:
     # Focal samples for country
     country:
       group_by: "division year month"
       max_sequences: 1500
       exclude: "--exclude-where 'country!={country}'"
     # Contextual samples from country's region
     region:
       group_by: "country year month"
       seq_per_group: 20
       exclude: "--exclude-where 'country={country}' 'region!={region}'"
       priorities:
         type: "proximity"
         focus: "country"
     # Contextual samples from the rest of the world,
     # excluding the current region to avoid resampling.
     global:
       group_by: "country year month"
       seq_per_group: 10
       exclude: "--exclude-where 'region={region}'"
       priorities:
         type: "proximity"
         focus: "country"

For ``country``-level sampling above, we specify that we want a maximum of 1,500 sequences from the country in question (here, Switzerland). Since we set ``group_by`` to “division year month”, all the Swiss sequences will be divided into groups by their division, month, and year of sampling, and the code will try to equally sample from each group to reach 1,500 sequences total.

Alternatively, in the ``region``-level sampling, we set ``seq_per_group`` to 20. This means that all the sequences from Europe (excluding Switzerland) will be divided into groups by their sampling country, month, and year (as defined by ``group_by``), and then 20 sequences will taken from each group (if there are fewer than 20 in any given group, all of the samples from that group will be taken).

Now we'll look at a subsampling scheme which defines a multi-``canton`` build. Cantons are regional divisions in Switzerland - below 'country,' but above 'location' (often city-level). In the advanced example, we'd like to be able to specify a set of neighboring 'cantons' and do focal sampling there, with contextual samples from elsewhere in the country, other countries in the region, and other regions in the world.

For cantons this looks like this:

.. code:: yaml

   # This build will take from 3 cantons - we have a sample rule for each,
   # rather than just one division that's focal build
   lac-leman:
     # focal samples
     geneva:
       group_by: "year month"
       seq_per_group: 300
       exclude: "--exclude-where 'division!=geneva'"
     vaud:
       group_by: "year month"
       seq_per_group: 300
       exclude: "--exclude-where 'division!=vaud'"
     valais:
       group_by: "year month"
       seq_per_group: 300
       exclude: "--exclude-where 'division!=valais'"

     # Contextual samples from the country
     country:
       group_by: "division year month"
       seq_per_group: 20
       exclude: "--exclude-where 'country!=switzerland'"

     # Contextual samples from division's region
     region:
       group_by: "country year month"
       seq_per_group: 10
       exclude: "--exclude-where 'region!=europe'"
       priorities:
         type: "proximity"
         focus: "country"
     # Contextual samples from the rest of the world, excluding the current
     # division to avoid resampling.
     global:
       group_by: "country year month"
       seq_per_group: 5
       exclude: "--exclude-where 'region=europe'"
       priorities:
         type: "proximity"
         focus: "country"

All entries above canton level (the 'contextual' samples) specify priorities. Currently, we have only implemented one type of priority called ``proximity``. It attempts to selected sequences as close as possible to the focal samples specified as ``focus: division``. The argument of the latter has to match the name of one of the other subsamples.

In addition to the ``exclude`` filter, you can also specify strains to keep by providing a ``query``. The ``query`` field uses augur filter's ``--query`` argument (introduced in version 8.0.0) and supports `pandas-style logical operators <https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#indexing-query>`__. For example, the following exclusionary filter,

.. code:: yaml

   exclude: "--exclude-where 'region!={region}' 'country!={country}' 'division!={division}'"

can also be written as an inclusionary filter like so:

.. code:: yaml

   query: --query "(region == {region}) & (country == {country}) & (division == '{division}')"

If you need parameters in a way that isn't represented by the config file, `create a new issue in the ncov repository <https://github.com/nextstrain/ncov/issues/new>`__ to let us know.

Ancestral trait reconstruction
------------------------------

Trait reconstruction is the process by which augur infers the most likely metadata value of an internal node. For example, if an internal node (which always represents a hypothesized, ancestral virus / case) has 3 descendants, all of which were isolated in Washington State, we might infer that the ancestor was most likely also circulating in Washington State (see :doc:`../visualization/interpretation` for more).

For each build, you can specify which categorical metadata fields to use for trait reconstruction.

.. raw:: html

   <!-- TODO: can someone please check this section for me? the existing docs were unclear to me -->

To specify this on a per-build basis, add a block like the following to your workflow config file:

.. code:: yaml

   traits:
     my_north_america_build: ### build name
       sampling_bias_correction: 2.5
       columns: ["country", "division"] ### traits to reconstruct; must match column names in metadata.tsv

Labeling clades
---------------

We assign clade labels according to :doc:`this schema <../reference/naming_clades>`.

Because the exact topology of the tree will vary across runs, clades are defined based on their unique mutations. These are specified in ``defaults/clades.tsv`` like so:

::

   # clade gene    site    alt

   A1a ORF3a   251 V
   A1a ORF1a   3606    F
