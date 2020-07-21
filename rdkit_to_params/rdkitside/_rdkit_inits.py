from __future__ import annotations
########################################################################################################################
__doc__ = \
    """
    Various entry points.
    """

__author__ = "Matteo Ferla. [Github](https://github.com/matteoferla)"
__email__ = "matteo.ferla@gmail.com"
__date__ = "10 July 2020 A.D."
__license__ = "MIT"
__version__ = "1.1.0"
__citation__ = "None."

########################################################################################################################

from rdkit import Chem
from rdkit.Chem import AllChem
from typing import Optional, Dict, List
from ._rdkit_convert import _RDKitCovertMixin

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
        pdb = Chem.MolFromPDBBlock(pdb_block, removeHs=False, proximityBonding=proximityBonding)
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
        pdb = Chem.MolFromPDBFile(pdb_file, removeHs=False, proximityBonding=proximityBonding)
        return cls._from_smiles_w_pdb(pdb, smiles, generic, name)

    @classmethod
    def _from_smiles_w_pdb(cls, pdb: Chem.Mol, smiles, generic, name):
        dodgy = Chem.SplitMolByPDBResidues(pdb, whiteList=[name])[name]
        AllChem.SanitizeMol(dodgy)
        good = Chem.MolFromSmiles(smiles)
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
        self.convert_mol()
        #####
        warnings.warn('CHI DISABLED. - has issues with this mode')  # todo correct this issue!
        self.CHI = [] # !!!!
        return self