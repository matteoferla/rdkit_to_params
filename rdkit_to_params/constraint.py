########################################################################################################################
__doc__ = \
    """
This contains ``make_constraint`` which is creates a constraint file.
It is completely independent and different in style because it was different.
    """

__author__ = "Matteo Ferla. [Github](https://github.com/matteoferla)"
__email__ = "matteo.ferla@gmail.com"
__date__ = "2020 A.D."
__license__ = "MIT"
__version__ = "1"
__citation__ = "None."

########################################################################################################################

from typing import Tuple, List, Union
from rdkit import Chem
from rdkit.Chem import AllChem

class Constraints:
    def __init__(self, smiles: Tuple[str, str],
                       names: List[str],
                       covalent_res: Union[str, int],
                       target_res: Union[str, int]):
        """
        Give smiles of the two sides (with *) and the names, residue numbers convert.
        It requires 2 atoms on each side in addition to the attachment point.
        Note that the atom names are stored in a non-standard way out of laziness. in the Atom prop '_AtomName'.
        The instance has the following attributes

        * ``smiles``: stored input smiles string
        *``names``: stored input list of names
        * ``covalent_res``: stored input liagnd residue
        * ``target_res``: stored input protein residue
        * ``cov_template``: Chem.Mol from first smiles.
        * ``target_template``: Chem.Mol from first smiles
        * ``combo``: combined templates
        * ``atom_pair_constraint``: AtomPair
        * ``angle_constraint``: Angle . NB. this is two lines.
        * ``dihedral_constaint``: dihedral

        Class methods:
        * ``assign_names(mol, list)`` classmethods that assigns names to a mol in place.
        * ``join_by_dummy(molA, molB)`` classmethods that returns a joined molecule

        :param smiles: a tuple/list of two string. The first is the ligand, the second is the peptide.
        :param names: a list of atom names. The '*' will need a name -but will be ignored-, but not the H.
        :param covalent_res: ligand residue in pose or PDB format (12 vs. 12A)
        :param target_res:  peptide residue in pose or PDB format (12 vs. 12A)
        """
        # store inputs
        self.smiles = smiles
        self.names = names
        self.covalent_res = covalent_res
        self.target_res = target_res
        # Convert the smiles
        self.cov_template = Chem.MolFromSmiles(smiles[0])
        self.target_template = Chem.MolFromSmiles(smiles[1])
        # check all is good
        assert self.cov_template.GetNumAtoms() + self.target_template.GetNumAtoms() == len(names), \
            f'{len(names)} were provided but the SMILES have {self.cov_template.GetNumAtoms()}+{self.target_template.GetNumAtoms()} atoms.'
        # assign names
        self.assign_names(self.cov_template, self.names[:self.cov_template.GetNumAtoms()])
        self.assign_names(self.target_template, self.names[self.cov_template.GetNumAtoms():])
        cov_con_name = self.get_conn(self.cov_template).GetProp('_AtomName')
        target_con_name = self.get_conn(self.target_template).GetProp('_AtomName')
        # combine and embed
        self.combo = self.join_by_dummy(self.cov_template, self.target_template)
        conf = self._get_conformer(self.combo)
        # get key atoms
        cov_con = self.get_atom(self.combo, cov_con_name)
        cov_fore = cov_con.GetNeighbors()[0]
        cov_fore_name = cov_fore.GetProp('_AtomName')
        target_con = self.get_atom(self.combo, target_con_name)
        target_fore = target_con.GetNeighbors()[0]
        target_fore_name = target_fore.GetProp('_AtomName')
        # do maths
        ## Note: constraint is in Radian not Degree...
        dist = Chem.rdMolTransforms.GetBondLength(conf, cov_con.GetIdx(), target_con.GetIdx())
        angle_target = Chem.rdMolTransforms.GetAngleRad(conf, target_fore.GetIdx(), target_con.GetIdx(), cov_con.GetIdx())
        angle_covalent = Chem.rdMolTransforms.GetAngleRad(conf, target_con.GetIdx(), cov_con.GetIdx(), cov_fore.GetIdx())
        dihedral = Chem.rdMolTransforms.GetDihedralRad(conf, target_fore.GetIdx(), target_con.GetIdx(), cov_con.GetIdx(),
                                                       cov_fore.GetIdx())
        self.atom_pair_constraint = f'AtomPair {target_con_name} {target_res} {cov_con_name} {covalent_res} HARMONIC {dist:.2f} 0.2\n'
        self.angle_constraint = f'Angle {target_fore_name} {target_res} {target_con_name} {target_res} {cov_con_name} {covalent_res} HARMONIC {angle_target:.2f} 0.35\n' +\
                                f'Angle {target_con_name} {target_res} {cov_con_name} {covalent_res} {cov_fore_name} {covalent_res} HARMONIC {angle_covalent:.2f} 0.35\n'
        self.dihedral_constaint = f'Dihedral {target_fore_name} {target_res} {target_con_name} {target_res} {cov_con_name} {covalent_res} {cov_fore_name} {covalent_res} CIRCULARHARMONIC {dihedral:.2f} 0.35\n'

    @classmethod
    def assign_names(cls, mol: Chem.Mol, names: List[str]) -> None:
        """
        Stores names of atoms as given in the list.
        totally non-standard way. PDBInfo is correct way. But too much effort.
        """
        for name, atom in zip(names, mol.GetAtoms()):
            atom.SetProp('_AtomName', name)  # Nonstandard. do not copy

    @classmethod
    def get_conn(cls, mol: Chem.Mol) -> Chem.Atom:
        """
        Get connecting atom of mol.
        """
        for atom in mol.GetAtoms():
            if atom.GetSymbol() == '*':
                return atom.GetNeighbors()[0]
        else:
            raise ValueError('Dummy atom not found')

    @classmethod
    def join_by_dummy(cls, a: Chem.Mol, b: Chem.Mol) -> Chem.Mol:
        # So I was worried that joining by the connect neighbour and then deleting the dummy
        # may cause issues of valence. So I did it this more convoluted way.
        # but did not check first... and this approach leads to sanitisation...
        for atom in a.GetAtoms():
            atom.SetDoubleProp('_originalIdx', atom.GetIdx())
        conn = cls.get_conn(a)
        mod_a = Chem.DeleteSubstructs(a, Chem.MolFromSmiles('*'))
        i = cls._get_new_index(mod_a, conn.GetIdx())
        combo = Chem.ReplaceSubstructs(b, Chem.MolFromSmiles('*'), mod_a, replacementConnectionPoint=i)[0]
        combo.UpdatePropertyCache()
        Chem.SanitizeMol(combo)
        return combo

    def _get_conformer(self, mol):
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        AllChem.MMFFOptimizeMolecule(mol)
        return mol.GetConformer()


    def __str__(self):
        # make
        constraints = [self.atom_pair_constraint,
                       self.angle_constraint,
                       self.dihedral_constaint]
        return ''.join(constraints)

    @classmethod
    def _get_new_index(self, mol: Chem.Mol, index: int) -> int:
        for atom in mol.GetAtoms():
            if atom.GetDoubleProp('_originalIdx') == index:
                return atom.GetIdx()

    def get_atom(self, mol: Chem.Mol, name: str) -> Chem.Atom:
        for atom in mol.GetAtoms():
            if atom.HasProp('_AtomName') and name == atom.GetProp('_AtomName'):  # Nonstandard. do not copy
                return atom
        else:
            raise ValueError(f'Atom {name} not found')


if __name__ == '__main__':
    c = Constraints(smiles=('*C(=N)', '*SC'),
                          names= ['*', 'CX', 'NY', '*', 'SG', 'CB'],
                          covalent_res= '1B',
                          target_res='145A')
    print(c)
