Installing in a Python virtual environment with ``virtualenv``
==============================================================

This section guides you in installing the pdbcraft package in a virtual environment, meaning an instance of Python that is isolated from your system.

Step 1 - Install ``virtualenv``
-------------------------------

First, check if the ``virtualenv`` Python package is installed in your system. This can be done by verifying whether the ``virtualenv`` command is available.

You can install the ``virtualenv`` package for your local user using ``pip``:

.. code-block:: shell

    pip install --user virtualenv

If the installation is successful, the ``virtualenv`` command will be available.

Step 2 - Create the virtual environment
---------------------------------------

Create your virtual environment in a directory of your choice (in this case, it will be ``./pdbcraft-env``):

.. code-block:: shell

    virtualenv -p /usr/bin/python3.11 pdbcraft-env

You should replace the value of option ``-p`` with to the location of the Python interpreter you want to use inside the virtual environment.

Remember that pdbcraft needs Python 3.11 or higher.

Step 3 - Activate the environment
---------------------------------

Activate the environment:

.. code-block:: shell

    source pdbcraft-env/bin/activate

Step 4 - Get pdbcraft
---------------------

First, move to the directory where you wanto to install pdbcraft. Do not install it inside the ``pdbcraft-env`` directory. 

Clone the source code from its GitHub repository within a directory of your choice and enter the local copy of the repository.

.. code-block:: shell

    git clone https://github.com/Center-for-Health-Data-Science/pdbcraft.git

If the ``git`` command is unavailable, you can download the repository content as a ZIP file from the web page of pdbcraft's GitHub repository and unzip it.

Step 5 - Install pdbcraft
-------------------------

You can now install pdbcraft using ``pip``:

.. code-block:: shell

    pip install ./pdbcraft

pdbcraft should now be installed.

Every time you need to run pdbcraft after opening a new shell, just run step 3 beforehand.
