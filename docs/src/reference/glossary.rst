========
Glossary
========

.. glossary::

   analysis directory

      The folder within ``ncov/`` where :term:`customization files <customization file>` live. Previously this was ``my_profiles/`` but we now allow any name of choice, and provide `ncov-tutorial <https://github.com/nextstrain/ncov-tutorial>`__ as a starter template.

   Auspice config file
      also ``auspice_config.json``

      A JSON file used to configure visualization in :term:`docs.nextstrain.org:Auspice`.

   config file
      also *workflow config file*, *workflow configuration file*, ``builds.yaml``

      A YAML file used to `configure <https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html#standard-configuration>`__ the :term:`Snakemake` workflow (via the ``--configfile`` option). Appends to and overrides default configuration in ``defaults/parameters.yaml``. For the :term:`ncov workflow`, this file must follow a :doc:`specific format <workflow-config-file>`.

   customization file

      A file used to customize the :term:`ncov workflow`.

      Examples: :term:`Auspice config file`, :term:`workflow config file<config file>`, :term:`default files`

   default files

      Default :term:`customization files <customization file>` provided in ``ncov/defaults/``.

   ncov workflow
      also *SARS-CoV-2 workflow*

      The workflow used to automate execution of :term:`builds<docs.nextstrain.org:build>`. Implemented in :term:`Snakemake`.

   Snakemake

      The workflow manager used in the :term:`ncov workflow`.
