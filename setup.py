########################################################################################################################
__doc__ = \
    """
    Pip install file. All the files have this block here, so in PyCharm I can change them en-mass.
    """

__author__ = "Matteo Ferla. [Github](https://github.com/matteoferla)"
__email__ = "matteo.ferla@gmail.com"
__date__ = "10 July 2020 A.D."
__license__ = "MIT"
__version__ = "1.1.2"
__citation__ = "None."

########################################################################################################################


import setuptools
import os

this_directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

print('rdkit and pyrosetta are optional... but are the main reason this was written.')
print('The former is conda, the latter is downloaded from RosettaCommons.')

setuptools.setup(
                name='rdkit_to_params',
                version=__version__,
                python_requires='>3.6',
                packages=setuptools.find_packages(),
                url='https://github.com/matteoferla/rdkit_to_params',
                license='MIT',
                author=__author__,
                author_email=__email__,
                description='Create or modify Rosetta params files (topology files) from scratch, RDKit mols or another params file.',
                long_description=long_description,
                long_description_content_type='text/markdown',
                test_suite='tests',
                install_requires=[
                ],
                extras_require = {
                        'rdkit_feature':  ["rdkit"],
                        'testing_feature': ["pyrosetta"]
                    }
)
