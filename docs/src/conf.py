# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
sys.path.insert(0, os.path.abspath('..'))


# -- Project information -----------------------------------------------------

from datetime import date
import subprocess

def git_authors():
    result = subprocess.run(
        ["git", "shortlog", "--summary", "HEAD"],
        stdout = subprocess.PIPE,
        check  = True)

    names = [
        line.strip().split("\t")[1]
            for line in result.stdout.decode("utf-8").splitlines()
    ]

    return names

def prose_list(items):
    if not items:
        return ""
    if len(items) == 1:
        return items[0]
    elif len(items) == 2:
        return " and ".join(items)
    else:
        return ", ".join([*items[0:-1], "and " + items[-1]])

project = 'SARS-CoV-2 Workflow'
#TODO copyright date change? copied from augur
copyright = '2014–%d Trevor Bedford and Richard Neher' % (date.today().year)
author = prose_list(git_authors())


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['recommonmark', 'sphinx.ext.autodoc', 'sphinxarg.ext', 'sphinx.ext.napoleon', 'sphinx_markdown_tables', 'sphinx.ext.intersphinx']

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# A string of reStructuredText that will be included at the end of every source
# file that is read. This is a possible place to add substitutions that should
# be available in every file.
rst_epilog = f"""
.. |authors| replace:: {author}
"""


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'nextstrain-sphinx-theme'

html_theme_options = {
    'logo_only': False, # if True, don't display project name at top of the sidebar
    'collapse_navigation': False, # if True, no [+] icons in sidebar
    'titles_only': True, # if True, page subheadings not included in nav
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# These paths are either relative to html_static_path
# or fully qualified paths (eg. https://...)
html_css_files = [
    'css/custom.css',
]

# -- Cross-project references ------------------------------------------------

intersphinx_mapping = {
    'docs.nextstrain.org': ('https://docs.nextstrain.org/en/latest/', None),
}
