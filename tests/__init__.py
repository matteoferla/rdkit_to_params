import unittest
from warnings import warn
import pyrosetta
pyrosetta.init(extra_options='-mute all')  # required for test
from rdkit_to_params import Params
from rdkit import Chem

class ParamTestCases(unittest.TestCase):
    def test_load(self):
        p = Params.load('../example/official_PHE.params')
        self.assertEqual(p.NAME, 'PHE')

    def test_round(self):
        p = Params.load('../example/official_PHE.params')
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
        self.assertIsNotNone(Chem.MolFromPDBBlock(pdbblock))

    def test_smiles(self):
        p = Params.from_smiles('*C(=O)C(Cc1ccccc1)[NH]*',  # recognised as amino acid.
                               name='PHX',  # optional.
                               atomnames={4: 'CZ'}  # optional, rando atom name for CB
                               )
        self.assertTrue(p.is_aminoacid())

    def test_renames(self):
        p = Params.from_smiles('CC(ONC)O', atomnames=['CA', 'CB', 'OX', 'ON', 'CX', 'CG'])
        p.rename_atom_by_name('CA', 'CZ')
        self.assertEqual(p.mol.GetAtomWithIdx(0).GetPDBResidueInfo().GetName(), ' CZ ')
        self.assertEqual(p.ATOM[0].name, ' CZ ')
        p.test()

if __name__ == '__main__':
    unittest.main()
