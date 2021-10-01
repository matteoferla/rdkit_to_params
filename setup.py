########################################################################################################################
__doc__ = 'README.md'
__author__ = "Matteo Ferla"
__url__ = "https://github.com/matteoferla/rdkit_to_params"
__email__ = "matteo.ferla@gmail.com"
__date__ = "29 January 2021 A.D."
__license__ = "MIT"
__version__ = "1.1.13"
__citation__ = "https://advances.sciencemag.org/content/7/16/eabf8711 (Fragmenstein)"

# ---------- imports  --------------------------------------------------------------------------------------------------
# remember it's `python setup.py sdist` and `python -m twine upload dist/rdkit_to_params-1.0.5.tar.gz`

from setuptools import setup, find_packages
from warnings import warn
from importlib import util
import os

this_directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()
description = 'Create or modify Rosetta params files (topology files) from scratch, RDKit mols or another params file.'

# ---------- Non pip modules  ------------------------------------------------------------------------------------------

if not util.find_spec('rdkit'):
    warn('Albeit optional, a lot of this code relies on rdkit which cannot be pip installed.' +
               'To install try either ' +
               'conda install -c conda-forge rdkit or ' +
               'sudo apt-get/brew install python3-rdkit or visit rdkit documentation.' +
         'or use the new pip install rdkit-pypi'
         )

if not util.find_spec('pyrosetta'):
    warn('The minimisation part of this code uses pyrosetta, which has to be downloaded from ' +
         'the Rosetta software site due to licencing. Without it only the classes Monster and Rectifier will work')

# ---------- setuptools.setup ------------------------------------------------------------------------------------------

setup(
    name='rdkit_to_params',
    version=__version__,
    python_requires='>3.6',
    packages=find_packages(),
    url=__url__,
    license=__license__,
    author=__author__,
    author_email=__email__,
    description=description,
    long_description=long_description,
    long_description_content_type='text/markdown',
    test_suite='tests',
    install_requires=[],  # none pip
    extras_require={},  # none pip
    entry_points={
            'console_scripts': ['smiles-to-params=rdkit_to_params.cli:smiles'],
        }
)
