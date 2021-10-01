from __future__ import annotations
########################################################################################################################
__doc__ = \
    """
    Various entry points.
    """

from ..version import *

########################################################################################################################

from rdkit import Chem
from rdkit.Chem import AllChem
from typing import *
from ._rdkit_convert import _RDKitCovertMixin
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers

import warnings

class _RDKitInitMixin(_RDKitCovertMixin):

    @classmethod
    def from_mol(cls, mol: Chem.Mol,
                 name: Optional[str] = None,
                 generic: bool = False,
                 atomnames: Optional[Dict[int, str]] = None) -> _RDKitInitMixin:
        """
        :param mol: Rdkit molecule, with explicit protons and all.
        :type mol: Chem.Mol
        :param name: 3-letter name
        :type name: str
        :param generic: convert mol with generic or classic AtomTypes?
        :type generic: bool
        :param atomnames: optional dictionary to set names.
        :rtype: instance
        """
        cls.log.debug('`from_mol` called...')
        self = cls.load_mol(mol, generic, name)  # stores and calls .fix_mol() method in prep file.
        if atomnames is None:
            pass
        else:
            self.rename(atomnames)
            self.polish_mol()
        self.convert_mol()
        return self

    @classmethod
    def from_smiles(cls, smiles: str, name='LIG', generic: bool = False,
                    atomnames: Optional[Dict[int, str]] = None) -> _RDKitInitMixin:
        """
        Make a Params instance from a smiles.

        :param smiles: SMILES to use.
        :param name: name3/resn
        :param generic: use generic atoms types
        :param atomnames: optional dictionary to set names.
        :return:
        """
        cls.log.debug('`from_smiles` called...')
        self = cls.load_smiles(smiles=smiles, name=name, generic=generic)
        if atomnames is None:
            pass
        else:
            self.rename(atomnames)
            self.polish_mol()
        self.convert_mol()
        return self

    @classmethod
    def from_smiles_w_pdbblock(cls, pdb_block: str, smiles: str, generic: bool = False, name='LIG',
                              proximityBonding: bool = True):
        """
        Assumes there is only one residue of the resn/name3

        :param pdb_block:
        :param generic: convert mol with generic or classic AtomTypes?
        :type generic: bool
        :param name: 3-letter name
        :type name: str
        :param proximityBonding: rdkit option for Chem.MolFromPDBFile. Did the author of the pdb **not** add CONECT?
        :type proximityBonding: bool
        :rtype: instance
        """
        cls.log.debug('`from_smiles_w_pdbblock` called...')
        pdb = Chem.MolFromPDBBlock(pdb_block, removeHs=True, proximityBonding=proximityBonding)
        # removeHs=True is no longer an option because it causes issues with the template step.
        return cls._from_smiles_w_pdb(pdb, smiles, generic, name)


    @classmethod
    def from_smiles_w_pdbfile(cls, pdb_file: str, smiles: str, generic: bool = False, name='LIG',
                              proximityBonding: bool = True):
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
        cls.log.debug('`from_smiles_w_pdbfile` called...')
        pdb = Chem.MolFromPDBFile(pdb_file, removeHs=True, proximityBonding=proximityBonding)
        # removeHs=True is no longer an option because it causes issues with the template step.
        return cls._from_smiles_w_pdb(pdb, smiles, generic, name)

    @classmethod
    def _from_smiles_w_pdb(cls, pdb: Chem.Mol, smiles, generic, name):
        dodgy = Chem.SplitMolByPDBResidues(pdb, whiteList=[name])[name]
        AllChem.SanitizeMol(dodgy)
        good = Chem.MolFromSmiles(smiles) # TODO switch to DummyMasker
        good.SetProp('_Name', name)
        dummies = []
        for atom in good.GetAtoms():
            if atom.GetSymbol() == '*':
                atom.SetAtomicNum(9)
                dummies.append(atom.GetIdx())
        Chem.SanitizeMol(good)
        good = AllChem.AddHs(good)
        AllChem.EmbedMolecule(good)
        AllChem.ComputeGasteigerCharges(good)
        AllChem.MMFFOptimizeMolecule(good)
        for d in dummies:
            good.GetAtomWithIdx(d).SetAtomicNum(0)
        self = cls.load_mol(good, generic=generic, name=name)
        self.move_aside()
        self.rename_from_template(dodgy)
        self.move_back()
        self.convert_mol()
        #####
        #warnings.warn('CHI DISABLED. - has issues with this mode')  # todo correct this issue!
        self.CHI.data = [] # !!!!
        return self

    @staticmethod
    def split_stereoisomers(mol: Chem.Mol) -> Dict[str, List[Chem.Mol]]:
        """
        Utility to split stereoisomers of amino acids into L/D

        :param mol:
        :return: dict of keys 'levo', 'dextro', 'neither', 'both'
        """
        # prep
        levo = Chem.MolFromSmiles('C[C@@H](C(=O))N')
        dextro = Chem.MolFromSmiles('C[C@H](C(=O))N')
        splits = {'levo': [], 'dextro': [], 'neither': [], 'both': []}
        for isomer in EnumerateStereoisomers(mol):
            is_levo = isomer.HasSubstructMatch(levo, useChirality=True)
            is_dextro = isomer.HasSubstructMatch(dextro, useChirality=True)
            if is_levo and is_dextro:  # two groups.
                box = 'both'
            elif is_levo:  # L
                box = 'levo'
            elif is_dextro:  # D
                box = 'dextro'
            else:  # gly or not amino acid
                box = 'neither'
            splits[box].append(isomer)
        return splits
