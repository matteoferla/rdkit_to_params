from rdkit_to_params import Params, DummyMasker
from rdkit import Chem
from rdkit.Chem import AllChem
import pyrosetta
import pyrosetta.rosetta.protocols as prp
from pyrosetta.rosetta.std import ostringstream

def test_covalent():

    mol = Chem.MolFromSmiles('*Nc1ccccc1C(=O)[O-]')
    mol = AllChem.AddHs(mol)
    with DummyMasker(mol, placekeeper_zahl=8):
        assert mol.GetAtomWithIdx(0).GetSymbol() == 'O'
        AllChem.ComputeGasteigerCharges(mol, throwOnParamFailure=True)
        AllChem.EmbedMultipleConfs(mol, 2)

    assert mol.GetAtomWithIdx(0).GetSymbol() == '*'
    assert mol.GetAtomWithIdx(1).GetDoubleProp('_GasteigerCharge') > -1

    params = Params.from_mol(mol, name='LIG')
    params.dump_pdb_conf(filename='conf.pdb')
    params.comments.append('foo')
    pose = params.to_pose()
    cys = pyrosetta.pose_from_sequence('C')
    pose.append_pose_by_jump(cys, 1)

    pyrosetta.rosetta.core.util.add_covalent_linkage(pose=pose,
                                                     resA_pos=1, resB_pos=2,
                                                     resA_At=pose.residue(1).atom_index('N1'),
                                                     resB_At=pose.residue(2).atom_index('SG'),
                                                     remove_hydrogens=True)
    scorefxn = pyrosetta.create_score_function('ref2015_cart')
    relax = prp.relax.FastRelax(scorefxn,15)
    movemap = pyrosetta.MoveMap()
    movemap.set_chi(True)
    movemap.set_bb(True)
    movemap.set_jump(True)
    relax.cartesian(True)
    relax.minimize_bond_angles(True)
    relax.minimize_bond_lengths(True)
    relax.set_movemap(movemap)
    relax.apply(pose)
    assert scorefxn(pose) < 50.
    buffer = ostringstream()
    pose.dump_pdb(buffer)
    pdb_string = buffer.str()
    roundtrip = Chem.MolFromPDBBlock(pdb_string.split('# All')[0],
                                     proximityBonding=True)