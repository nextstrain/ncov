# ncov 'Read The Docs' Documentation.


## Building the docs

Build dependencies are managed with [Conda](https://conda.io).
Install them
into an isolated environment named `ncov-docs` with:

    conda env create -f=conda.yml

Enter the environment with:

    conda activate ncov-docs

You can now build the documentation with:

    make html

which invokes Sphinx to build static HTML pages in `build/html/`.
You can view them by running:

    open build/html/index.html


To monitor the source files for changes and automatically rebuild as necessary,
run:

    make livehtml

and then open <http://localhost:8000>.  Pages open in the browser will
automatically refresh when they're rebuilt.

You can clean the build directory for a fresh start with:

    make clean

Leave the environment with:

    conda deactivate
