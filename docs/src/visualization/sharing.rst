Visualizing and sharing results
===============================

`Nextstrain.org <https://www.nextstrain.org/ncov>`__ uses Auspice to visualize JSON files that are created by Augur. While this is the most visible example, there are many other ways to use Auspice to visualize your results.

We'll walk through each of these in detail, roughly in order from the simplest and most private to the slightly more complex and publicly shareable.

If none of these options meet your needs, please `get in touch <https://discussion.nextstrain.org/>`__!

--------------

First: a note on sensitive or private metadata
----------------------------------------------

Below, we describe the privacy considerations for each of the available options.

**One approach to handling sensitive metadata is to simply keep it entirely separate. With any of these options, you can also choose to drag-and-drop a separate TSV file onto the Auspice visualization.**

Doing so will enable you to color by any of the values in this extra metadata file, but none of that data ever leaves your local computer.

How to view visualize private metadata
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

1. Create a TSV file with a ``strain`` or ``name`` column that matches all the samples in your dataset
2. Add your sensitive metadata to the remaining columns
3. On your computer, drag and drop the file onto the browser window where Auspice is visualizing your JSON

*For more help formatting this metadata file, including how to do so using Excel,* :doc:`see this page <../guides/data-prep/index>`.

--------------

Option 1: Drag-and-drop web-based visualization
-----------------------------------------------

-  **Quickstart**: Drag-and-drop the file from ``./auspice/sarscov2_global.json`` onto the page at https://auspice.us.
-  **Advantages:** Quick, no-setup viewing of results, including sensitive data.
-  **Limitations:** Requires separate management of JSON file sharing and version control. Sharing a specific view via a URL isn't possible with this method.

How to view
^^^^^^^^^^^

1. Navigate to https://auspice.us
2. Drag the output JSON file from ``./auspice/<buildname>.json`` onto the page
3. [Optional] drag and drop a TSV with additional or private metadata onto the page (see above)

How to share
^^^^^^^^^^^^

Share the JSON file and instructions directly.

Privacy and security
^^^^^^^^^^^^^^^^^^^^

When your browser connects to auspice.us, it downloads from the server a version of the Auspice code which runs solely on your computer, within your browser. Then, when you drag a file onto the page, that code processes the data in your browser and displays it to you without ever sending it back to the auspice.us server. All the heavy bioinformatics computations were already performed and stored in the file you provide, which is what lets everything work quickly just on your computer.

Option 2: Nextstrain community pages
------------------------------------

-  **Example:** `CZBiohub's California COVID Tracker <https://nextstrain.org/community/czbiohub/covidtracker/>`__
-  **Advantages:** Fully featured, plug-and-play visualization of any JSON file hosted on Github.
-  **Limitations:** Only available for publicly viewable JSON files in public repositories.

How to get started
^^^^^^^^^^^^^^^^^^

Quickstart:

1. Put your JSON in a github repository like so: ``myGithubOrganization/myRepository/auspice/<myBuildName>.json``
2. Navigate to ``https://nextstrain.org/community/myGithubOrganization/myRepository/myBuildName`` 3. [Optional] Drag and drop a TSV with additional or private metadata onto the page (see above)

Check out our `full guide to community pages here <https://nextstrain.org/docs/contributing/community-builds>`__.

.. _privacy-and-security-1:

Privacy and security
^^^^^^^^^^^^^^^^^^^^

Community builds are visible to anyone with the URL.

Option 3: Local viewing on your computer with Auspice
-----------------------------------------------------

-  **Quickstart**: ``ncov$ auspice view``
-  **Advantages:** Offline, entirely local viewing of results, including sensitive data.
-  **Limitations:** Requires collaborators to install Auspice locally. Requires separate management of JSON file sharing and version control. Sharing a specific view via a URL isn't possible with this method.

.. _how-to-view-1:

How to view
^^^^^^^^^^^

1. Follow the instructions `here <https://nextstrain.github.io/auspice/introduction/install>`__ to install Auspice on your computer.
2. Make sure the JSON you'd like to visualize is in ``./auspice/<mybuildname>.json``; alternatively, pass the ``--datasetDir`` flag to specify another directory.
3. Run ``auspice view`` and select the build of interest.
4. [Optional] drag and drop a TSV with additional or private metadata onto the page (see above)

.. _how-to-share-1:

How to share
^^^^^^^^^^^^

Share the JSON file and instructions directly.

.. _privacy-and-security-2:

Privacy and security
^^^^^^^^^^^^^^^^^^^^

When running locally, both the server and the client run on your computer; no internet connection is required. No data ever leaves your local machine.

Option 4: Sharing with Nextstrain Groups
----------------------------------------

-  **Example:** https://nextstrain.org/groups/blab/
-  **Advantages:** Web-based viewing of results with full authentication / login controls; accommodates both public and private datasets. Sharing a specific view via URL is possible with this method.
-  **Limitations:** Setup is slightly more involved, but we're ready to help!

.. _how-to-get-started-1:

How to get started
^^^^^^^^^^^^^^^^^^

Nextstrain Groups are a new feature; if you'd like to use this option, please `get in touch <mailto:hello@nextstrain.org>`__ and we'll help you get started right away!

.. _privacy-and-security-3:

Privacy and security
^^^^^^^^^^^^^^^^^^^^

With Nextstrain Groups, you can choose whether each dataset is publicly viewable or private to only other users in your group. Data is hosted in an AWS S3 bucket under your control, and is not shared with the Nextstrain team or anyone else.

Option 5: Deploying your own Auspice server
-------------------------------------------

-  **Advantages:** Fully-featured Auspice instance, natively hosted on your own domain.
-  **Limitations:** More technically involved, especially if user authentication is required.

.. _how-to-get-started-2:

How to get started
^^^^^^^^^^^^^^^^^^

`See our guide here <https://nextstrain.github.io/auspice/server/introduction>`__

.. _privacy-and-security-4:

Privacy and security
^^^^^^^^^^^^^^^^^^^^

Independently hosted Auspice servers can be configured with any security protocols necessary.
