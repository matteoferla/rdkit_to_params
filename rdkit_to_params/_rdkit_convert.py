from __future__ import annotations

########################################################################################################################
__doc__ = \
    """
There are probably loads of mistakes here ---or more correctly
in ``_RDKitPrepMixin()._add_rtypes`` or ``_RDKitPrepMixin()._add_genrtypes``.
    """

__author__ = "Matteo Ferla. [Github](https://github.com/matteoferla)"
__email__ = "matteo.ferla@gmail.com"
__date__ = "2020 A.D."
__license__ = "MIT"
__version__ = "1"
__citation__ = "None."

########################################################################################################################

from ._rdkit_prep import _RDKitPrepMixin

from typing import List, Dict, Union, Optional
from collections import defaultdict, deque, namedtuple
from rdkit import Chem
from rdkit.Chem import AllChem
from warnings import warn

import numpy as np

class _RDKitCovertMixin(_RDKitPrepMixin):
    _Measure = namedtuple('Measure', ['distance','angle', 'torsion'])

    @classmethod
    def from_mol(cls, mol: Chem.Mol, generic: bool = False, name:Optional[str]=None) -> _RDKitCovertMixin:
        """
        :param mol: Rdkit molecule, with explicit protons and all.
        :type mol: Chem.Mol
        :param generic: convert mol with generic or classic AtomTypes?
        :type generic: bool
        :param name: 3-letter name
        :type name: str
        :rtype: instance
        """
        self = cls.load_mol(mol, generic, name) # stores and calls .fix_mol()
        self.convert_mol()
        return self

    @classmethod
    def from_smiles(cls, smiles: str, name='LIG', generic:bool=False) -> _RDKitPrepMixin:
        mol = Chem.MolFromSmiles(smiles)
        mol.SetProp('_Name', name)
        mol = AllChem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        AllChem.MMFFOptimizeMolecule(mol)
        return cls.from_mol(mol, generic, name)

    @classmethod
    def from_smiles_w_pdbfile(cls, pdb_file:str, smiles:str, generic: bool = False, name='LIG', proximityBonding:bool=True):
        """
        Assumes there is only one residue of the resn/name3

        :param pdb_file:
        :param generic: convert mol with generic or classic AtomTypes?
        :type generic: bool
        :param name: 3-letter name
        :type name: str
        :param proximityBonding: rdkit option for Chem.MolFromPDBFile. Did the author of the pdb **not** add CONECT?
        :type proximityBonding: bool
        :rtype: instance
        """
        pdb = Chem.MolFromPDBFile(pdb_file, removeHs=False, proximityBonding=proximityBonding)
        dodgy = Chem.SplitMolByPDBResidues(pdb, whiteList=[name])[name]
        good = Chem.MolFromSmiles(smiles)
        good.SetProp('_Name', name)
        good = AllChem.AddHs(good)
        AllChem.EmbedMolecule(good)
        AllChem.MMFFOptimizeMolecule(good)
        self = cls.load_mol(good, generic=generic, name=name)
        self.rename_from_template(dodgy)
        self.convert_mol()
        return self

    def convert_mol(self):
        # NAME
        if self.mol.HasProp("_Name"):
            title = self.mol.GetProp("_Name").strip()
        else:
            title = 'Unnamed molecule'
        if len(title) == 3:
            self.NAME = title
        else:
            self.comments.append(title)
            self.NAME = self._get_resn_from_PDBInfo()
        # SMILES
        self.comments.append(Chem.MolToSmiles(self.mol))
        self.TYPE.append('LIGAND')
        # ATOM & CONNECT
        for atom in self.mol.GetAtoms():
            self._parse_atom(atom)
        # CHI
        self._parse_rotatables()
        # BOND
        for bond in self.mol.GetBonds():
            self._parse_bond(bond)
        # ICOOR
        self._parse_icoors()
        # NBR
        self._find_centroid()

    def _parse_icoors(self) -> None:
        undescribed = deque(self.mol.GetAtoms())
        # TODO add override to specific first atom?
        described = []
        conf = self.mol.GetConformer()
        aforementioned = []
        rotation_count = 0
        while undescribed:
            # prevent overcycling.
            if rotation_count > self.mol.GetNumAtoms():
                d = [f'{a.GetIdx()}: {self._get_PDBInfo_atomname(a)}' for a in described]
                raise StopIteration(f'Too many cycles...{d}')
            atom = undescribed[0]
            if atom.GetSymbol() == '*' and len(described) < 3:
                warn(f'DEBUG. Dummy not allowed in first three lines...')
                undescribed.rotate(-1)
                rotation_count += 1
                continue # This cannot be in the first three lines, whereas the R group is first in a canonical SMILES.
            elif len(described) == 0: #First line
                # I assume in some cases Virtual atoms are called for...
                for parent in self._get_unseen_neighbors(atom, [], True):
                    siblings = self._get_unseen_neighbors(parent, [atom], True)
                    if len(siblings) != 0:
                        break
                else:
                    warn(f'DEBUG. First row: No siblings for {self._get_PDBInfo_atomname(atom)}!')
                    undescribed.rotate(-1)
                    continue
                atoms = [atom, atom, parent, siblings[0]]
                aforementioned = [atom, parent, siblings[0]]
            elif len(described) == 1: #Second line
                parent = described[0]
                sibling = [a for a in aforementioned if a.GetIdx() not in (parent.GetIdx(), atom.GetIdx())][0]
                atoms = [atom, parent, atom, sibling]
            elif len(described) == 2: #Third line
                parent = described[1]
                sibling = described[0]
                atoms = [atom, parent, sibling, atom]
            else:
                d_idx = [d.GetIdx() for d in described]
                d_neighs = [n for n in self._get_unseen_neighbors(atom, [], False) if n.GetIdx() in d_idx]
                if len(d_neighs) == 0:
                    warn(f'DEBUG. Other row: No parent for {self._get_PDBInfo_atomname(atom)}!')
                    undescribed.rotate(-1)
                    rotation_count += 1
                    continue
                parent = d_neighs[0]
                d_neighs = [n for n in self._get_unseen_neighbors(parent, [atom], False) if n.GetIdx() in d_idx]
                if len(d_neighs) == 0:
                    undescribed.rotate(-1)
                    warn(f'DEBUG. Other row: No siblings for  {self._get_PDBInfo_atomname(atom)}!')
                    rotation_count += 1
                    continue
                elif len(d_neighs) == 1:
                    sibling = d_neighs[0]
                    d_neighs = [n for n in self._get_unseen_neighbors(sibling, [parent, atom], False) if n.GetIdx() in d_idx]
                    if len(d_neighs) == 0:
                        undescribed.rotate(-1)
                        warn(f'DEBUG. Other row: No cousins for  {self._get_PDBInfo_atomname(atom)}')
                        rotation_count += 1
                        continue
                    else:
                        atoms = [atom, parent, sibling, d_neighs[0]]

                else:
                    atoms = [atom, parent, d_neighs[0], d_neighs[1]]
            undescribed.remove(atom)
            described.append(atom)
            m = self._get_measurements(conf, *atoms)
            self.ICOOR_INTERNAL.append(dict(child=self._get_PDBInfo_atomname(atoms[0]),
                                            phi=m.torsion,
                                            theta= m.angle,
                                            distance=m.distance,
                                            parent=self._get_PDBInfo_atomname(atoms[1]),
                                            second_parent=self._get_PDBInfo_atomname(atoms[2]),
                                            third_parent=self._get_PDBInfo_atomname(atoms[3])))

    def _parse_atom(self, atom: Chem.Atom) -> None:
        if atom.GetSymbol() == '*':
            neighbor = atom.GetNeighbors()[0]
            n_name = self._get_PDBInfo_atomname(neighbor)
            self.CONNECT.append([n_name, len(self.CONNECT) + 1, 'CONNECT'])
        else:
            d = self._get_atom_descriptors(atom)
            self.ATOM.append(d)

    def _parse_rotatables(self):
        """
        I am not sure if this is right...
        I do not know if a C(=O)-C=O is rotatable (cis-trans).

        :return:
        """
        def is_single(a, b):
            bond = self.mol.GetBondBetweenAtoms(a.GetIdx(), b.GetIdx())
            return bond.GetBondType() == Chem.BondType.SINGLE

        for atom in self.mol.GetAtoms():
            if atom.GetSymbol() == 'H':
                continue #PROTON_CHI 3 SAMPLES 2 0 180 EXTRA 1 20 thing...
            elif atom.GetSymbol() == '*':
                continue
            for neighbor in self._get_unseen_neighbors(atom, []):
                if atom.GetIdx() > neighbor.GetIdx():
                    continue
                elif is_single(atom, neighbor):
                    for grandneighbor in self._get_unseen_neighbors(neighbor, [atom]):
                        if is_single(grandneighbor, neighbor):
                            ggneighbors = self._get_unseen_neighbors(grandneighbor, [atom, neighbor])
                            if ggneighbors: # 3-4 bond does not need to be single, nor does it matter which is the 4th atom
                                self.CHI.append(dict(index=len(self.CHI) + 1,
                                                first=self._get_PDBInfo_atomname(atom),
                                                second=self._get_PDBInfo_atomname(neighbor),
                                                third=self._get_PDBInfo_atomname(grandneighbor),
                                                fourth=self._get_PDBInfo_atomname(ggneighbors[0])))

    def _parse_bond(self, bond: Chem.Bond) -> None:
        if any([atom.GetSymbol() == '*' for atom in (bond.GetBeginAtom(), bond.GetEndAtom())]):
            return None  # CONNECT.
        if bond.GetBondTypeAsDouble() == 1.5:
            order = 4
        else:
            order = int(bond.GetBondTypeAsDouble())
        self.BOND.append([self._get_PDBInfo_atomname(bond.GetBeginAtom()),
                   self._get_PDBInfo_atomname(bond.GetEndAtom()),
                   order])

    def _get_measurements(self, conf: Chem.Conformer, a:Chem.Atom, b:Chem.Atom, c:Chem.Atom, d:Chem.Atom):
        ai = a.GetIdx()
        bi = b.GetIdx()
        ci = c.GetIdx()
        di = d.GetIdx()
        dist = 0
        angle = 0
        tor = 0
        try:
            dist = Chem.rdMolTransforms.GetBondLength(conf, ai, bi)
            angle = 180 - Chem.rdMolTransforms.GetAngleDeg(conf, ai, bi, ci)
            tor = Chem.rdMolTransforms.GetDihedralDeg(conf, ai, bi, ci, di)
        except ValueError:
            pass
        if str(tor) == 'nan': #quicker than isnan.
            tor = 0
        return self._Measure(distance=dist,
                             angle= angle,
                             torsion= tor)

    def _get_atom_descriptors(self, atom: Chem.Atom) -> dict:
        return {'name': self._get_PDBInfo_atomname(atom),
               'rtype': atom.GetProp('_rType'),
               'mtype': ' X  ',
               'partial': atom.GetDoubleProp('_GasteigerCharge')}

    def _get_nondummy_neighbors(self, atom) -> List[str]:
            # returns list of names!
            return [self._get_PDBInfo_atomname(neighbor) for neighbor in atom.GetNeighbors() if neighbor.GetSymbol() != '*']

    def _get_unseen_neighbors(self, atom:Chem.Atom, seen: List[Chem.Atom], nondummy:bool=True):
        neighbors = atom.GetNeighbors()
        if nondummy:
            neighbors = [neighbor for neighbor in neighbors if neighbor.GetSymbol() != '*']
        seenIdx={a.GetIdx() for a in seen}
        return [neighbor for neighbor in neighbors if neighbor.GetIdx() not in seenIdx]

    def _find_centroid(self):
        a = Chem.Get3DDistanceMatrix(self.mol)
        n = np.linalg.norm(a, axis=0)  # Frobenius norm
        i = int(np.argmin(n))
        s = np.max(a, axis=0)[i]
        self.NBR_ATOM.append(self.mol.GetAtomWithIdx(i).GetPDBResidueInfo().GetName())
        self.NBR_RADIUS.append(str(s))

    def polish_mol(self):
        """
        The mol may be inconsistent.
        :return:
        """
        name = self.NAME[:3]
        index = 1
        for atom in self.mol.GetAtoms():
            atom.GetPDBResidueInfo().SetResidueName(name)
            atom.GetPDBResidueInfo().SetResidueIdx(index)




