from setuptools import setup
import os

this_directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

print('rdkit and pyrosetta are optional... but are the main reason this was written.')
print('The former is conda, the latter is downloaded from RosettaCommons.')

setup(
    name='rdkit_to_params',
    version='1.0.2',
    python_requires='>3.6',
    packages=['rdkit_to_params'],
    url='https://github.com/matteoferla/rdkit_to_params',
    license='MIT',
    author='Matteo Ferla',
    author_email='matteo.ferla@gmail.com',
    description='Create or modify Rosetta params files (topology files) from scratch, RDKit mols or another params file.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    install_requires=[
    ]
)
