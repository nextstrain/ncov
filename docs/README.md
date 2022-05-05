# ncov 'Read The Docs' Documentation.


## Building the docs

Build dependencies are managed with [Conda](https://conda.io).
Install them
into an isolated environment named `ncov-docs` with:

    mamba env create -f=conda.yml

Enter the environment with:

    mamba activate ncov-docs

You can now build the documentation with:

    make html

which invokes Sphinx to build static HTML pages in `build/html/`.
You can view them by running:

    open build/html/index.html

You can clean the build directory for a fresh start with:

    make clean

Leave the environment with:

    conda deactivate
