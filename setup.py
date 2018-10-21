#!/usr/bin/env python

import os
from distutils.core import setup

# Create list of all sub-directories with
#   __init__.py files...
packages = []
for subdir, dirs, files in os.walk('sp2graph'):
    if '__init__.py' in files:
        packages.append(subdir.replace(os.sep, '.'))


def git_version():
    # Default release info
    MAJOR = 0
    MINOR = 0
    MICRO = 1
    VERSION = [MAJOR, MINOR, MICRO]
    # Git revision prior to release:
    GIT_REVISION = "0434304ad39d4101d5c8761f1a485dfacd8b59d8"
    GIT_LABEL = '.'.join(map(str, [MAJOR, MINOR, MICRO]))

    def _minimal_ext_cmd(cmd):
        # construct minimal environment
        env = {}
        for k in ['SYSTEMROOT', 'PATH']:
            v = os.environ.get(k)
            if v is not None:
                env[k] = v
        # LANGUAGE is used on win32
        env['LANGUAGE'] = 'C'
        env['LANG'] = 'C'
        env['LC_ALL'] = 'C'
        out = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE, env=env).communicate()[0]
        return out.strip().decode('ascii')

    current_path = os.path.dirname(os.path.realpath(__file__))

    try:
        # Get top-level directory
        git_dir = _minimal_ext_cmd(['git', 'rev-parse', '--show-toplevel'])
        # Assert that the git-directory is consistent with this setup.py script
        if git_dir != current_path:
            raise ValueError('Not executing the top-setup.py script')

        # Get latest revision tag
        rev = _minimal_ext_cmd(['git', 'rev-parse', 'HEAD'])
        if len(rev) > 7:
            GIT_REVISION = rev
        # Get latest tag
        tag = _minimal_ext_cmd(['git', 'describe', '--abbrev=0', '--tags'])
        if len(tag) > 4:
            VERSION = tag[1:].split('.')
        # Get complete "git describe" string
        label = _minimal_ext_cmd(['git', 'describe', '--tags'])
        if len(label) > 7:
            GIT_LABEL = label
        # Get number of commits since tag
        count = _minimal_ext_cmd(['git', 'rev-list', tag + '..', '--count'])
        if len(count) == 0:
            count = '1'

    except Exception as e:
        count = '0'

    return GIT_REVISION, VERSION, int(count), GIT_LABEL


def write_version(filename='sp2graph/info.py'):
    version_str = """# This file is automatically generated from sp2graph setup.py
# Git information (specific commit, etc.)
git_revision = '{git}'
git_revision_short = git_revision[:7]
git_count = {count}

# Version information
major   = {version[0]}
minor   = {version[1]}
micro   = {version[2]}

# Release tag
release = 'v'+'.'.join(map(str,[major, minor, micro]))

# Version (release + count)
version = release
if git_count > 0:
    # Add git-revision to the version string
    version += '+' + str(git_count)

# Extensive version description
label   = '{description}'
"""
    # If we are in git we try and fetch the git version as well
    GIT_REV, GIT_VER, GIT_COUNT, GIT_LAB = git_version()
    with open(filename, 'w') as fh:
        fh.write(version_str.format(version=GIT_VER,
                                    count=GIT_COUNT,
                                    git=GIT_REV,
                                    description=GIT_LAB))

write_version()

# Main setup of python modules
setup(name='sp2graph',
      description='Bond-order graph analysis of sp2 carbon nanostructures',
      url='https://github.com/dipc-cc/sp2graph',
      license='GPL-3.0',
      packages=packages)
