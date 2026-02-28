"""Test suite for rdkit_to_params package."""
from pathlib import Path

import pyrosetta
from rdkit import Chem
from rdkit.Chem import AllChem

from rdkit_to_params import DummyMasker, Params

# Fixtures are defined in conftest.py

def test_load(official_phe_params):
    p = Params.load(str(official_phe_params))
    assert p.NAME == 'PHE'

def test_round(official_phe_params):
    p = Params.load(str(official_phe_params))
    p.IO_STRING[0].name3 = 'PHX'
    p.IO_STRING[0].name1 = 'Z'
    p.AA = 'UNK'  # If it's not one of the twenty (plus extras), UNK!
    del p.ROTAMER_AA[0]
    p.rename_atom(' CB ', ' CX ')  # this renames
    p.dump('fake.params')
    pose = p.test()
    buffer = pyrosetta.rosetta.std.stringbuf()
    pose.dump_pdb(pyrosetta.rosetta.std.ostream(buffer))
    pdbblock = buffer.str()
    assert Chem.MolFromPDBBlock(pdbblock) is not None

def test_aa_smiles():
    p = Params.from_smiles('*C(=O)C(Cc1ccccc1)N*',  # recognised as amino acid.
                           name='PHX',  # optional.
                           # atomnames={4: 'CZ'}  # optional, rando atom name for CB
                           )
    assert p.is_aminoacid(), 'Should be an amino acid.'

def test_renames():
    p = Params.from_smiles('CC(ONC)O', atomnames=['CA', 'CB', 'OX', 'ON', 'CX', 'CG'])
    p.rename_atom_by_name('CA', 'CZ')
    assert p.mol.GetAtomWithIdx(0).GetPDBResidueInfo().GetName() == ' CZ '
    assert p.ATOM[0].name == ' CZ '
    p.test()

def test_covalent_from_mol():
    mol = Chem.MolFromSmiles('CC(=O)C(Cc1ccccc1)N*')
    p = Params.from_mol(mol, name='PHX')

def test_dummy_masker():
    mol = Chem.MolFromSmiles('*CCC')
    mol = AllChem.AddHs(mol)
    with DummyMasker(mol, placekeeper_zahl=8):
        assert mol.GetAtomWithIdx(0).GetSymbol() == 'O'
        AllChem.ComputeGasteigerCharges(mol, throwOnParamFailure=True)
        AllChem.EmbedMolecule(mol)
    assert mol.GetAtomWithIdx(0).GetSymbol() == '*'
    assert mol.GetAtomWithIdx(1).GetDoubleProp('_GasteigerCharge') > -1

