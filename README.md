[![License: LGPL v3](https://img.shields.io/badge/License-LGPL%20v3-blue.svg)](https://www.gnu.org/licenses/lgpl-3.0)

![GitHub Logo][logo]

# sp2graph #

The goal of sp2praph project is to provide tools for analysing the bond-order of sp2 carbon nanostructures.

## Objectives ##
Given a planar carbon-based structure (sp2 hydridized), sp2graph aims to provide:

   - all possible Kekulé representations
   - the most stable structure according to Clar's sextet theory
   - allow for user defined constrains, such as to impose single or double bonds at specified connections and/or radicals

## Methodology ##



## Dependencies ##
Before installation of Inelastica the following packages are required
   - numpy

## Installation ##
Manual installation of  is performed with the command

    python setup.py install --prefix=<prefix>
    # or
    python setup.py install --home=<my-python-home>

One may also wish to set the following environment variables

    export PYTHONPATH=<my-python-home>/lib/<python-version>/site-packages
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:<my-python-home>/lib

## Contributions, issues and bugs ##
Contributions are highly appreciated.

If you find any bugs please form a [bug report/issue][issues]

If you have a fix please consider adding a [pull request][pulls].

## License ##
The sp2graph license is [LGPL][lgpl], please see the LICENSE file.

<!---
Links to external and internal sites.
-->
[logo]: https://github.com/dipc-cc/sp2graph/tree/master/doc/images/sp2graph_logo.png
[lgpl]: http://www.gnu.org/licenses/lgpl.html
[issues]: https://github.com/dipc-cc/sp2graph/issues
[pulls]: https://github.com/dipc-cc/sp2graph/pulls
