#!/usr/bin/env python

import os
from distutils.core import setup

# Create list of all sub-directories with
#   __init__.py files...
packages = []
for subdir, dirs, files in os.walk('sp2graph'):
    if '__init__.py' in files:
        packages.append(subdir.replace(os.sep, '.'))

# Main setup of python modules
setup(name='sp2graph',
      description='Bond-order graph analysis of sp2 carbon nanostructures',
      url='https://github.com/dipc-cc/sp2graph',
      license='GPL-3.0',
      packages=packages)
