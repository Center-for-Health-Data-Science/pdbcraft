Installing in a Python virtual environment with ``conda``
=========================================================

This section guides you in installing the pdbcraft package in a ``conda`` environment.

Step 1 - Install ``conda``
--------------------------

Go `here <https://docs.conda.io/en/latest/miniconda.html>`_ for detailed instructions on how to install ``conda``.

Installing ``miniconda`` rather than the full ``anaconda`` package is advised.

Once ``conda`` is installed on your system, you can create a virtual environment.

Step 2 - Create the ``conda`` environment
-----------------------------------------

You can now create your new ``conda`` environment.

.. code-block:: shell
    
    conda env create --name pdbcraft-env python=3.11

We can specify the name of the environment with the ``--name`` option, and the Python interpreter used in the environment with the ``python=`` option.

Remember that pdbcraft needs Python 3.11 or higher.

Step 3 - Activate the environment
---------------------------------

You can activate the ``conda`` environment using the ``conda activate`` command.

.. code-block:: shell
    
    conda activate pdbcraft-env

Step 4 - Get pdbcraft
---------------------

Clone the source code from its GitHub repository within a directory of your choice and enter the local copy of the repository.

.. code-block:: shell

    git clone https://github.com/Center-for-Health-Data-Science/pdbcraft.git

If the ``git`` command is unavailable, you can download the repository content as a ZIP file from the pdbcraft GitHub repository web page and unzip it.

Step 5 - Install pdbcraft
-------------------------

You can now install pdbcraft using ``pip``:

.. code-block:: shell

    pip install ./pdbcraft

pdbcraft should now be installed.

Every time you need to run pdbcraft after opening a new shell, just run step 3 beforehand.