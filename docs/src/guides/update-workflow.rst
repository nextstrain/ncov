Update the workflow
===================

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
