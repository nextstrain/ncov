Setup and installation
======================

The following steps will prepare you to run complete analyses of SARS-CoV-2 data by installing required software and running a simple example workflow.

1. Setup your Nextstrain environment with Conda
-----------------------------------------------

The following instructions use `Conda <https://docs.conda.io/en/latest/>`__ to install the tools you'll need for this tutorial. Conda is a package and environment management system that allows you to install Python and other software into controlled environments without disrupting other software you have installed (e.g., on your computer, your shared cluster, etc.).

.. note::

   If you use Microsoft Windows, install the Windows Subsystem for Linux (WSL). Follow instructions to open a new WSL window for your Linux distribution and then run the following commands.

Create the Nextstrain environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`Install Miniconda with Python 3 for your operating system <https://docs.conda.io/en/latest/miniconda.html>`__, update Conda to the latest version, and install Mamba for faster installation.

.. code:: bash

   conda update -n base conda
   conda install -n base -c conda-forge mamba

Create a Conda environment named ``nextstrain``. This command will install `Git <https://git-scm.com/>`__, `Snakemake <https://snakemake.readthedocs.io/en/stable/>`__, and `Nextstrain's full toolkit <https://docs.nextstrain.org/en/latest/install-nextstrain.html>`__, the tools you'll need to work through this tutorial.

.. code:: bash

   mamba create -n nextstrain -c conda-forge -c bioconda \
     augur auspice nextstrain-cli nextalign snakemake awscli git pip

Confirm that the installation worked.

.. code:: bash

   conda activate nextstrain
   nextstrain check-setup --set-default

The final output from the last command should look like this:

.. code:: bash

   Setting default environment to native.

.. warning::

   If the final output from the last command is ``Setting default environment to docker.``, close Docker and run the last command again. The Nextstrain CLI will prefer the Docker environment when it is available, but this tutorial assumes you are running commands in the native environment.

Update the Nextstrain environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Update your environment regularly to the latest versions of these tools.

.. code:: bash

   # Update Conda and Mamba.
   mamba update -n base conda mamba

   # Update tools in the Nextstrain environment.
   conda activate nextstrain
   mamba update --all -c conda-forge -c bioconda

2. Download the ncov workflow
-----------------------------

Download the workflow
~~~~~~~~~~~~~~~~~~~~~

Use Git to download a copy of the ncov repository containing the workflow and this tutorial.

.. code:: bash

   git clone https://github.com/nextstrain/ncov.git
   cd ncov

Alternately, `download a compressed copy of the ncov repository <https://github.com/nextstrain/ncov/archive/refs/heads/master.zip>`__ called ``ncov-master.zip``. Open this file to decompress it and create a directory called ``ncov-master/`` with the contents of the workflow in it. Navigate to this directory from the command line.

Update the workflow
~~~~~~~~~~~~~~~~~~~

We update the official workflow regularly with:

-  `curated metadata including latitudes/longitudes, clade annotations, and low quality sequences <https://github.com/nextstrain/ncov/commits/master>`__
-  bug fixes
-  `new features <../reference/change_log>`__

Update your local copy of the workflow, to benefit from these changes.

.. code:: bash

   # Download and apply changes from the Nextstrain team.
   # This only works if there is no conflict with your local repository.
   git pull --ff-only origin master

   # OR:

   # Alternately, download and apply changes from the Nextstrain team
   # and then replay your local changes on top of those incoming changes.
   git pull --rebase origin master

Alternately, download a specific version of the workflow that you know works for you. We create new `releases of the workflow <https://github.com/nextstrain/ncov/releases/>`__ any time we introduce breaking changes, so you can choose when to update based on `what has changed <../reference/change_log>`__.

.. code:: bash

   # Download version 7 (v7) of the workflow.
   curl -OL https://github.com/nextstrain/ncov/archive/refs/tags/v7.zip

   # Uncompress the workflow.
   unzip v7.zip

   # Change into the workflow's directory.
   cd ncov-7/

3. Run a basic analysis with example data
-----------------------------------------

Run a basic workflow with example data, to confirm that your Nextstrain environment is properly configured.

.. code:: bash

   nextstrain build . --cores 4 --use-conda \
     --configfile ./my_profiles/getting_started/builds.yaml

The ``nextstrain build`` command runs a pathogen analysis or “`build <https://docs.nextstrain.org/projects/augur/en/stable/faq/what-is-a-build.html>`__” defined by a Snakemake workflow in a specific directory. Since our ``Snakefile`` is in the current directory, we specify the build directory as ``.``. All other arguments pass through to Snakemake.

The ``--use-conda`` flag `instructs Snakemake to maintain a workflow-specific Conda environment <https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#integrated-package-management>`__ defined in ``workflow/envs/nextstrain.yaml``. Snakemake will activate this environment before it runs each rule, ensuring each rule runs with the necessary software. The first time you use this flag, Snakemake will download and install all Nextstrain tools required to run the workflow. This can take a little time. The next time you use this flag, Snakemake will detect the existing environment and use it immediately.

The ``getting_started`` build produces a minimal global phylogeny for visualization in Auspice. This workflow should complete in about 10 minutes on a MacBook Pro (2.7 GHz Intel Core i5) with four cores, including 5 minutes to create the Conda environment and 5 minutes to run the workflow.

4. Visualize the phylogeny for example data
-------------------------------------------

`Open http://auspice.us <http://auspice.us>`__ in your browser. Drag and drop the JSON file ``auspice/ncov_global.json`` anywhere on the http://auspice.us landing page, to visualize the resulting phylogeny. The resulting phylogeny should look something like this.

.. figure:: ../images/getting-started-tree.png
   :alt: Phylogenetic tree from the “getting started” build as visualized in Auspice

   Phylogenetic tree from the “getting started” build as visualized in Auspice
