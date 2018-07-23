[![License: LGPL v3](https://img.shields.io/badge/License-LGPL%20v3-blue.svg)](https://www.gnu.org/licenses/lgpl-3.0)

![GitHub Logo][logo]

The goal of _sp2graph_ project is to provide tools for analyzing the bond-order of sp2-hybridized carbon nanostructures. This is particular interesting for the interpretation of recent experiments performed with scanning tunneling microscopy on graphene-based nanostructures synthesized on metallic surfaces.

## Objectives ##
Given the `xyz` geometry of a planar carbon-based structure (sp2 hybridized), _sp2graph_ aims to provide:

   - all possible Kekulé representations
   - the most stable structure according to Clar's sextet theory
   - allow for user defined constrains, such as to impose single or double bonds at specified connections and/or radicals
   - estimation of the most stable structure through first-neighbors tight-binding model

## Introduction ##

![anthracene](/doc/images/anthracene12_rot.png)

![phenanthrene](/doc/images/phenanthrene12_rot.png)

## Methodology ##


## Installation ##
For _sp2graph_ installation the following packages are required:

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

If you find any bugs please form a [bug report/issue][issues]

If you have a fix please consider adding a [pull request][pulls].

## License ##
The _sp2graph_ license is [LGPL][lgpl], please see the LICENSE file.

<!---
Links to external and internal sites.
-->
[logo]: /doc/images/sp2graph_logo.png
[lgpl]: http://www.gnu.org/licenses/lgpl.html
[issues]: https://github.com/dipc-cc/sp2graph/issues
[pulls]: https://github.com/dipc-cc/sp2graph/pulls
