.. _install:

Installation
============

Dependencies
------------

These packages are required:

 * `NumPy`_
 * `Matplotlib`_

Manual installation
-------------------

Manual installation is performed with the command:

.. code-block:: bash

    python setup.py install --prefix=<prefix>
    # or
    python setup.py install --home=<my-python-home>

One may also wish to set the following environment variables

.. code-block:: bash

    export PYTHONPATH=<my-python-home>/lib/<python-version>/site-packages
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:<my-python-home>/lib
