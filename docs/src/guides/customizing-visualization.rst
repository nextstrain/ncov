Customizing visualization
=========================

Visualization options can be configured in either a :term:`workflow config file<config file>` or a :term:`Auspice config file`, depending on the option.

.. contents:: Table of Contents
   :local:

Options in the workflow config file
-----------------------------------

These options can be coded into the workflow config file directly without requiring a custom Auspice config file.

Custom color schemes
~~~~~~~~~~~~~~~~~~~~

To specify a custom color scale:

1. Add a ``colors.tsv`` file, where each line is a tab-delimited list of a metadata column name; a metadata value; and a corresponding hex code. Example:

   ::

      country Russia  #5E1D9D
      country Serbia  #4D22AD
      country Europe  #4530BB
      ...

2. Update your workflow config file with a reference:

   .. code:: yaml

      files:
        colors: "my-ncov-analyses/colors.tsv"

Changing the dataset description
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The dataset description, which appears below the visualizations, is read from a file which is specified in the workflow config file. Per-build description can be set by specifying them in the workflow config file:

.. code:: yaml

   builds:
     north-america: # name of the build
       description: my-ncov-analyses/north-america-description.md

If that is not provided, then a per-run description is used, also specified in the workflow config file:

.. code:: yaml

   files:
     description: my-ncov-analyses/my_description.md

Options in the Auspice config file
----------------------------------

These options require creating an Auspice config file, used to configure :term:`docs.nextstrain.org:Auspice`. It is specified in the workflow config file using the ``auspice_config`` entry. Example:

.. code:: yaml

   auspice_config: ncov-tutorial/auspice-config-custom-data.json

This overrides the default Auspice config file, ``defaults/auspice_config.json``.

Adding custom metadata fields to color by
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. Add a :doc:`valid metadata column <./data-prep/local-data>` to your ``metadata.tsv``
2. Add an entry to the ``colorings`` block of the Auspice config file:

   .. code:: json

      "colorings": [
      {
         "key": "location",
         "title": "Location",
         "type": "categorical"
      },
      {
         "key": "metadata_column_name",
         "title": "Display name for interface",
         "type": "<categorical/continuous>"
      }
      ]

Choosing defaults
~~~~~~~~~~~~~~~~~

You can specify the default view in the ``display_defaults`` block of an Auspice config file:

.. code:: json

   "display_defaults": {
     "color_by": "division",
     "distance_measure": "num_date",
     "geo_resolution": "division",
     "map_triplicate": true,
     "branch_label": "none"
   },

Choosing panels to display
~~~~~~~~~~~~~~~~~~~~~~~~~~

Similarly, you can choose which panels to enable in the ``panels`` block:

.. code:: json

   "panels": [
     "tree",
     "map",
     "entropy"
   ]
