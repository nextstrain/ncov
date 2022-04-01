Nextstrain overview
===================

Nextstrain has two main parts:

- :term:`docs.nextstrain.org:Augur` **performs the bioinformatic analyses** required to produce a tree, map, and other inferences from your input data.
- The outputs of Augur form the input for :term:`docs.nextstrain.org:Auspice`, **which provides the visualizations** you see on Nextstrain.org

You can find more information about how these tools fit together :doc:`here <docs.nextstrain.org:learn/parts>`. We'll come back to Auspice when we get to the :doc:`visualization <../visualization/sharing>` section.

First, let's take a look at how Augur works.

How bioinformatic analyses are managed
--------------------------------------

At its core, Augur is a collection of Python scripts, each of which handles one step in the bioinformatic analyses necessary for visualization with Auspice.

As you might imagine, keeping track of the input and output files from each step individually can get very confusing, very quickly. So, **to manage all of these steps, we use a workflow manager called Snakemake**.

.. note::

   There are many other workflow managers out there, such as Nextflow. While we fully encourage you to use whichever workflow tools you prefer, we only provide support and maintenance for Snakemake.

Snakemake is an incredibly powerful workflow manager with many complex features. For our purposes, though, we only need to understand a few things:

-  **Each step in a workflow is called a "rule."** The inputs, outputs, and shell commands for each step/rule are defined in a ``.smk`` file.
-  Each rule has a number of **parameters, which are specified in a ``.yaml`` file**.
-  Each rule produces **output (called a "dependency") which may be used as input to other rules**.

Overview of a Nextstrain build
------------------------------

Below is an illustration of each step in a standard :term:`Nextstrain build <docs.nextstrain.org:build>`. Dependencies (output files from one step that act as input to the next) are indicated by grey arrows. Input files which must be provided are indicated with red outlines. As you can see in yellow, the final output is a JSON file for visualization in Auspice.

Required input files (e.g. the sequence data generated in the `data preparation section <../guides/data-prep>`__, or other files which are part of this repo) are indicated with red outlines. We'll walk through each of these in detail in the next section.

.. figure:: ../images/basic_nextstrain_build.png
   :alt: nextstrain_build

Running multiple builds
-----------------------

It is common practice to run several related builds. For example, to run one analysis on just your data and another analysis that incorporates background / contextual sequences, you could configure two different builds.

The ncov workflow facilitates this through the ``builds`` section in a :term:`workflow config file <config file>`. This is covered in more detail in the :doc:`genomic surveillance tutorial <../tutorial/genomic-surveillance>`.

We encourage you to take a look at `main_workflow.smk <https://github.com/nextstrain/ncov/blob/master/workflow/snakemake_rules/main_workflow.smk>`__ to see what each rule is doing in more detail.

.. note::

   Not all of the rules included are essential, or may even be desirable for your analysis. Your workflow may be able to be made a lot simpler, depending on your goals.
