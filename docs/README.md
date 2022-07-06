# ncov 'Read The Docs' Documentation.


## Building the docs

Build dependencies are managed with Pip.
Install them into an isolated environment under `docs/env/` with:

    python3 -m venv env/
    ./env/bin/python3 -m pip install --upgrade pip setuptools wheel
    ./env/bin/python3 -m pip install -r requirements.txt

Enter the environment with:

    source ./env/bin/activate

You can now build the documentation with:

    make html

which invokes Sphinx to build static HTML pages in `build/html/`.
You can view them by running:

    open build/html/index.html

You can clean the build directory for a fresh start with:

    make clean

Leave the environment with:

    deactivate
