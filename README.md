[![License: LGPL v3](https://img.shields.io/badge/License-LGPL%20v3-blue.svg)](https://www.gnu.org/licenses/lgpl-3.0)

![sp2graph logo](/doc/images/sp2graph_logo.png "sp2graph code")

The goal of _sp2graph_ project is to provide tools for analyzing the bond-order of *sp<sup>2</sup>*-hybridized carbon nanostructures. This is particular interesting for the interpretation of recent experiments performed with scanning tunneling microscopy on graphene-based nanostructures synthesized on metallic surfaces.

## Objectives ##
Given the `xyz` geometry of a planar carbon-based structure (*sp<sup>2</sup>*-hybridized), __sp2graph__ aims to provide:

   - all possible Kekulé representations
   - the most stable structure according to Clar's &pi;-sextet theory
   - allow for user defined initial constrains, by imposing:
      - single or double bonds at specified connections
      - radicals at specific sites
   - estimation of the most stable structure through first-neighbors tight-binding model

## Introduction ##


![anthracene resonant structures](/doc/images/anthracene.png "anthracene")  
**anthrace:** *Kekulé resonance structures with the corresponding Clar sextets*.

![phenanthrene resonant structures](/doc/images/phenanthrene.png "phenanthrene")  
**phenanthrene:** *Kekulé resonance structures with the corresponding Clar sextets*

## Methodology ##


## Installation ##
For __sp2graph__ installation the following packages are required:

   - numpy

Manual installation of  is performed with the command:

    python setup.py install --prefix=<prefix>
    # or
    python setup.py install --home=<my-python-home>

One may also wish to set the following environment variables:

    export PYTHONPATH=<my-python-home>/lib/<python-version>/site-packages
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:<my-python-home>/lib

## Contributions, issues and bugs ##
Contributions are highly appreciated.

If you find any bug please form a [bug report/issue][issues]

If you have a fix please consider adding a [pull request][pulls].

## License ##
The __sp2graph__ license is [LGPL][lgpl], please see the LICENSE file.

<!---
Links to external and internal sites.
-->
[lgpl]: http://www.gnu.org/licenses/lgpl.html
[issues]: https://github.com/dipc-cc/sp2graph/issues
[pulls]: https://github.com/dipc-cc/sp2graph/pulls
