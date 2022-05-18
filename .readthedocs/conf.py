# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# -- pyrosetta install -------------------------------------------------------
import os

# import pyrosetta_help
# pyrosetta_help.install_pyrosetta(username=os.environ.get('PYROSETTA_USERNAME'),
#                                  password=os.environ.get('PYROSETTA_PASSWORD'),
#                                  )


from unittest.mock import MagicMock
import sys

sys.modules['pyrosetta'] = MagicMock()



# -- Project information -----------------------------------------------------

project = 'rdkit-to-params'
copyright = '2022, University of Oxford'
github_username = 'matteoferla'
github_repository = 'rdkit_to_params'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'readthedocs_ext.readthedocs',
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode',
    'sphinx.ext.todo',
    'sphinx_toolbox.more_autodoc',
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#
# This is also used if you do content translation via gettext catalogs.
# Usually you set "language" from the command line for these cases.
language = 'en'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------

html_theme = 'sphinx_rtd_theme'
html_static_path = [] #'_static'


# -- Extension configuration -------------------------------------------------

always_document_param_types = True
typehints_defaults = 'braces'
todo_include_todos = True

# --- add md files ---------------------------------------------------------

import m2r2  # noqa
import os, re

repo_base_path = os.path.abspath("../../")

def convert_write(markdown_filename, srt_filename):
    # unlike Fragmenstein there are no images to convert
    # so we can just copy the file
    with open(markdown_filename) as fh:
        markdown_block = fh.read()
    #markdown_block = re.sub(r'\[(?P<label>.*?)\]\((?P<link>.*?)\)', fix_md_link, markdown_block)
    rst_block = m2r2.convert(markdown_block)
    with open(srt_filename, 'w') as fh:
        fh.write(rst_block)

convert_write(os.path.join(repo_base_path, 'README.md'), 'introduction.rst')
convert_write(os.path.join(repo_base_path, 'database_entries.md'), 'database_entries.rst')
convert_write(os.path.join(repo_base_path, 'atom_types.md'), 'atom_types.rst')