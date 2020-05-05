from setuptools import setup

print('rdkit and pyrosetta are optional... but are the main reason this was written.')
print('The former is conda, the latter is downloaded from RosettaCommons.')

setup(
    name='rdkit_to_params',
    version='1',
    python_requires='>3.6',
    packages=['rdkit_to_params'],
    url='https://github.com/matteoferla/rdkit_to_params',
    license='MIT',
    author='Matteo Ferla',
    author_email='matteo.ferla@gmail.com',
    description='Create or modify Rosetta params files (topology files) from scratch, RDKit mols or another params file.',
    install_requires=[
        #'numpy' is required for rdkit part and its already a dependency of it.
    ]
)
