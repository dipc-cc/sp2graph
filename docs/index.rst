.. sp2graph documentation master file, created by
   sphinx-quickstart on Fri Feb 23 13:48:57 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


|license|_


Welcome to sp2graph documentation!
======================================

The goal of **sp2graph** project is to provide tools for analyzing the bond-order of *sp2*-hybridized carbon nanostructures.
This is particular interesting for the interpretation of recent experiments performed with scanning probe microscopy on graphene-based nanostructures synthesized on metallic surfaces.

**sp2graph** was initiated by Pedro Brandimarte and Thomas Frederiksen at DIPC in June 2018.

The latest release can obtained `here <releases_>`_ and the development version through:

.. code-block:: bash

    git clone https://github.com/dipc-cc/sp2graph.git


Objectives
----------

Given the *xyz* geometry of a planar carbon-based structure (*sp2*-hybridized), *sp2graph* aims to provide:

1. all possible Kekulé representations
2. the most stable Clar structure according to Clar's &pi;-sextet theory
3. allow for user defined initial constrains, by imposing:

 * single or double bonds at specified connections
 * radicals at specific sites

4. estimation of the most stable Kekulé structure through simple nearest-neighbor tight-binding model


Contributions, issues and bugs
------------------------------

Contributions are highly appreciated!

If you find any bugs please form a `bug report/issue <issue_>`_.

If you have a fix please consider adding a `pull request <pr_>`_.


Funding
-------
Financial support from Spanish AEI (FIS2017-83780-P "GRANAS") is acknowledged.


Indices
-------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: User guide

   install
   examples

.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: API documentation

   api-gen/sp2graph


.. |license| image:: https://img.shields.io/badge/License-LGPL%20v3-blue.svg
.. _license: https://www.gnu.org/licenses/lgpl-3.0
