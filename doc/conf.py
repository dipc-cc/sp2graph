# -*- coding: utf-8 -*-
#
# sp2graph configuration file for the Sphinx documentation builder,
# initially created by sphinx-quickstart on Feb 23 2018.
#
# This file does only contain a selection of the most common options. For a
# full list see the documentation:
# http://www.sphinx-doc.org/en/master/config.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
_this_dir = os.path.dirname(__file__)
sys.path.insert(0, os.path.abspath(_this_dir))
sys.path.insert(0, os.path.abspath(os.path.dirname(_this_dir)))
#sys.path.insert(0, os.path.abspath('..'))
#sys.path.insert(0, os.path.abspath('../sp2graph'))


# -- Project information -----------------------------------------------------

project = u'sp2graph'
author = u'Pedro Brandimarte and Thomas Frederiksen'
copyright = author + u', 2018'

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.
# The short X.Y version
#version = u'1.3'
# The full version, including alpha/beta/rc tags
#release = u'1.3.0'

import sp2graph.info as info
release = info.release # release tag
version = info.release


# The default language to highlight source code in.
highlight_language = 'python'

# -- General configuration ---------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
#
# needs_sphinx = '1.0'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.doctest',
    'sphinx.ext.intersphinx',
    'sphinx.ext.todo',
    'sphinx.ext.coverage',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
    'sphinx.ext.autosummary',
    'numpydoc',
    'sphinx.ext.extlinks',
    'sphinx.ext.inheritance_diagram',
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#
# source_suffix = ['.rst', '.md']
source_suffix = '.rst'

# The master toctree document.
master_doc = 'index'

# Insert the links into the epilog (globally)
# This means that every document has access to the links
rst_epilog = ''.join(open('epilog.dummy').readlines())

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#
# This is also used if you do content translation via gettext catalogs.
# Usually you set "language" from the command line for these cases.
language = u'en'

# Add __init__ classes to the documentation
autoclass_content = 'class'
autodoc_default_flags = ['members', 'undoc-members',
                         'inherited-members', 'no-show-inheritance']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path .
exclude_patterns = ['build', '**/setupegg.py', '**/setup.rst', '**/tests', '**.ipynb_checkpoints', '**/info.py']

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# A list of ignored prefixes for module index sorting.
modindex_common_prefix = ['sp2graph.']

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
#
# html_theme_options = {}

# The name for this set of Sphinx documents.  If None, it defaults to
# "<project> v<release> documentation".
html_title = 'sp2graph %s documentation' %release

# A shorter title for the navigation bar.  Default is the same as html_title.
html_short_title = 'sp2graph'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
if os.path.exists('_static'):
    html_static_path = ['_static']
else:
    html_static_path = []

# Custom sidebar templates, must be a dictionary that maps document names
# to template names.
#
# The default sidebars (for documents that don't match any pattern) are
# defined by theme itself.  Builtin themes are using these templates by
# default: ``['localtoc.html', 'relations.html', 'sourcelink.html',
# 'searchbox.html']``.
#
# html_sidebars = {}

# Project logo.
html_logo = 'images/sp2graph_logo.png'

# -- Options for HTMLHelp output ---------------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = 'sp2graph'


# -- Options for LaTeX output ------------------------------------------------

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    #
    'papersize': 'a4paper',

    # The font size ('10pt', '11pt' or '12pt').
    #
    'pointsize': '10pt',

    # Additional stuff for the LaTeX preamble.
    #
    # 'preamble': '',

    # Latex figure (float) alignment
    #
    'figure_align': 'htbp',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    (master_doc, 'sp2graph.tex', u'sp2graph Documentation',
     u'Pedro Brandimarte and Thomas Frederiksen', 'manual'),
]


# -- Options for manual page output ------------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    (master_doc, 'sp2graph', u'sp2graph Documentation',
     [author], 1)
]


# -- Options for Texinfo output ----------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (master_doc, 'sp2graph', u'sp2graph Documentation', author, None,
     'Python package for tools for analyzing the bond-order of carbon nanostructures.',
     'Miscellaneous'),
]


# -- Extension configuration -------------------------------------------------

# -- Options for intersphinx extension ---------------------------------------

# Example configuration for intersphinx: refer to the Python standard library.
intersphinx_mapping = {
    'python': ('https://docs.python.org/', None),
    'numpy': ('https://docs.scipy.org/doc/numpy', None),
    'scipy': ('https://docs.scipy.org/doc/scipy/reference', None),
}

# -- Options for todo extension ----------------------------------------------

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = True

# These two options should solve the "toctree contains reference to nonexisting document"
# problem.
# See here: numpydoc #69
class_members_toctree = False
# If this is false we do not have double method sections
numpydoc_show_class_members = False
