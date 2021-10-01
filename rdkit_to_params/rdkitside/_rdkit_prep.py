from __future__ import annotations

from rdkit_to_params.entries import Entries

import os

########################################################################################################################
__doc__ = \
    """
The main class here is ``_RDKitPrepMixin``, which adds the various pre checks.
It does not rely on any ``Params`` entry stuff. So can be used by itself for testing.

    """

from ..version import *

########################################################################################################################

from typing import List, Dict, Union
from collections import defaultdict, deque, namedtuple
from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS
from warnings import warn
import re, string
from typing import Optional
from rdkit_to_params.entries import CONNECTEntry

from ._rdkit_rename import _RDKitRenameMixin

class _RDKitPrepMixin(_RDKitRenameMixin):

    def __init__(self):
        self.log.critical('This is init should not have been called! This exists to stop the IDE from getting angry/ for debugging.')
        self.NAME = 'LIG'
        self.TYPE = Entries.from_name('TYPE')
        self.mol = None
        self.generic = False
        self._rtype = []

    @classmethod
    def load_mol(cls, mol: Chem.Mol, generic:bool=False, name:Optional[str]=None) -> _RDKitPrepMixin:
        """
        A fully prepared molecule with optional dummy atoms to be coverted into a Params object

        :param mol: fully prepared molecule with optional dummy atoms
        :param generic: generic or classic atom types
        :param name: 3 letter name
        :return:
        """
        self = cls()
        self.log.debug('`load_mol` called...')
        self.mol = mol
        self.generic = generic
        self.TYPE.append('LIGAND')
        if name is None:
            self.log.warning('Residue `name` not specified defaulting to "LIG"')
            name = 'LIG'
        self.NAME = name
        self.mol.SetProp('_Name', name)
        self.fix_mol()
        # conversion elsewhere
        return self

    @classmethod
    def load_smiles(cls, smiles:str, generic:bool=False, name:Optional[str]=None) -> _RDKitPrepMixin:
        """
        A SMILES with optional dummy atoms to be coverted into a Params object

        :param smiles: can contain multiple *
        :param generic: generic or classic atom types
        :param name: 3 letter name
        :return:
        """
        cls.log.debug('`load_smiles` called...')
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError('The SMILES string could not be converted')
        mol.SetProp('_Name', name)
        mol = Chem.AddHs(mol)
        # cannot embed more than one dummy
        changed = []
        for atom in mol.GetAtoms():
            if atom.GetSymbol() == '*':
                atom.SetAtomicNum(6)
                atom.SetHybridization(Chem.HybridizationType.SP3)
                changed.append(atom.GetIdx())
        AllChem.EmbedMolecule(mol)
        AllChem.MMFFOptimizeMolecule(mol)
        AllChem.ComputeGasteigerCharges(mol)
        for atomIdx in changed:
            mol.GetAtomWithIdx(atomIdx).SetAtomicNum(0)
        # operate upon!
        return cls.load_mol(mol, generic, name)


    def fix_mol(self):
        self.log.debug('`fix_mol` called...')
        # partial charges.
        if not self.mol.GetAtomWithIdx(0).HasProp('_GasteigerCharge'):
            self.log.debug('... Adding Gasteiger')
            self._add_partial_charges()
        self._fix_atom_names()
        if self.generic is True:
            self._add_genrtypes()
        else:
            self._add_rtypes()

    def dump_pdb(self, filename: str, overwrite=False, stripped=True) -> None:
        """
        Write first conformer to a PDB file.

        :param filename: name of file.
        :param overwrite: error if exists?
        :param stripped: remove * atoms?
        :return: None
        """
        self.log.debug(f'Writing PDB to file {filename}')
        mol = self._prep_dump_pdb(filename, overwrite, stripped)
        Chem.MolToPDBFile(mol, filename)

    def dump_pdb_conf(self, filename: str, overwrite=False, stripped=True) -> int:
        """
        Write conformers to a PDB file.

        :param filename: name of file.
        :param overwrite: error if exists?
        :param stripped: remove * atoms?
        :return: number of conf written
        """
        self.log.debug(f'Writing PDB conformers to file {filename}')
        mol = self._prep_dump_pdb(filename, overwrite, stripped)
        w = Chem.PDBWriter(filename)
        for cid in range(mol.GetNumConformers()):
            w.write(mol, confId=cid)
        n = w.NumMols()
        w.close()
        return n

    def _prep_dump_pdb(self, filename: str, overwrite=False, stripped=True) -> Chem.Mol:
        if os.path.exists(filename) and overwrite:
            raise FileExistsError('The file already exists')
        if stripped:
            return self.dummyless
        else:
            return self.mol

    def dumps_pdb(self, stripped=True) -> str:
        if stripped:
            Chem.MolToPDBBlock(self.dummyless)
        else:
            Chem.MolToPDBBlock(self.mol)

    @property
    def dummyless(self):
        return Chem.DeleteSubstructs(self.mol, Chem.MolFromSmiles('*'))

    def _add_rtypes(self) -> None:
        """
        Add Rosetta Atom types to each atom.
        Mostly a guess...
        """

        # These are the elements with a single type.
        # Some have oxidation states '2p', but only one...
        element_to_type = {'B': 'Bsp2', 'F': 'F', 'CL': 'Cl', 'BR': 'Br', 'I': 'I', 'ZN': 'Zn2p', 'CO': 'Co2p',
                           'CU': 'Cu2p', 'MG': 'Mg2p', 'CA': 'Ca2p', 'Si': 'Si', 'NA': 'Na1p', 'K': 'K1p', 'HE': 'He',
                           'LI': 'Li', 'BE': 'Be', 'NE': 'Ne', 'AL': 'Al', 'AR': 'Ar', 'SC': 'Sc', 'TI': 'Ti', 'V': 'V',
                           'CR': 'Cr', 'MN': 'Mn', 'NI': 'Ni', 'GA': 'Ga', 'GE': 'Ge', 'AS': 'As', 'SE': 'Se',
                           'KR': 'Kr', 'RB': 'Rb', 'SR': 'Sr', 'Y': 'Y', 'ZR': 'Zr', 'NB': 'Nb', 'MO': 'Mo', 'TC': 'Tc',
                           'RU': 'Ru', 'RH': 'Rh', 'PD': 'Pd', 'AG': 'Ag', 'CD': 'Cd', 'IN': 'In', 'SN': 'Sn',
                           'SB': 'Sb', 'TE': 'Te', 'XE': 'Xe', 'CS': 'Cs', 'BA': 'Ba', 'LA': 'La', 'CE': 'Ce',
                           'PR': 'Pr', 'ND': 'Nd', 'PM': 'Pm', 'SM': 'Sm', 'EU': 'Eu', 'GD': 'Gd', 'TB': 'Tb',
                           'DY': 'Dy', 'HO': 'Ho', 'ER': 'Er', 'TM': 'Tm', 'YB': 'Yb', 'LU': 'Lu', 'HF': 'Hf',
                           'TA': 'Ta', 'W': 'W', 'RE': 'Re', 'OS': 'Os', 'IR': 'Ir', 'PT': 'Pt', 'AU': 'Au', 'HG': 'Hg',
                           'TL': 'Tl', 'PB': 'Pb', 'BI': 'Bi', 'PO': 'Po', 'AT': 'At', 'RN': 'Rn', 'FR': 'Fr',
                           'RA': 'Ra', 'AC': 'Ac', 'TH': 'Th', 'PA': 'Pa', 'U': 'U', 'NP': 'Np', 'PU': 'Pu', 'AM': 'Am',
                           'CM': 'Cm', 'BK': 'Bk', 'CF': 'Cf', 'ES': 'Es', 'FM': 'Fm', 'MD': 'Md', 'NO': 'No',
                           'LR': 'Lr'}

        groups = [
                  {'name': 'silane', 'SMARTS': '[Si]~O', 'types': ['Si', 'OSi']},
                  {'name': 'phosphoric', 'SMARTS': 'P[OH]', 'types': ['Pha', 'OHha', 'Hha']},
                  {'name': 'phosphate', 'SMARTS': 'P~O', 'types': ['Pha', 'OPha']},
                  {'name': 'free carbonic acid?', 'SMARTS': 'C(=O)(O)O', 'types': ['CO3', 'OC3', 'OC3']}, # no idea.
                  {'name': 'carboxylate', 'SMARTS': 'C(=O)[O-]', 'types': ['COO', 'OOC', 'OOC']},
                  {'name': 'carboxylic', 'SMARTS': 'C(=O)[OH]', 'types': ['COO', 'OOC', 'OH']},
                  {'name': 'ester', 'SMARTS': 'C(=O)OC', 'types': ['COO', 'Oet2', 'Oet3', None]},
                  {'name': 'ester?', 'SMARTS': 'C(=O)O', 'types': ['COO', 'OOC', 'OOC']},
                  {'name': '3amide', 'SMARTS': 'C(=O)N(C)C', 'types': ['CNH2', 'ONH2', 'Npro', None, None]},
                  {'name': '3amide', 'SMARTS': 'C(=O)NC', 'types': ['CNH2', 'ONH2', 'Nbb', None]},
                  {'name': 'amide', 'SMARTS': 'C(=O)N', 'types': ['CNH2', 'ONH2', 'NH2O']},
                  {'name': 'keto-aldehyde', 'SMARTS': 'C(=O)', 'types': [None, 'OOC']}, ## ???
                  {'name': '2amine', 'SMARTS': 'CNC', 'types': [None, 'Nbb', None]}, # actually Nbb is a amide.
                  {'name': 'azo', 'SMARTS': 'NN', 'types': ['NtrR', 'NtrR']},
                  {'name': 'azo', 'SMARTS': 'nn', 'types': ['NtrR', 'NtrR']},
                  {'name': 'aramid?', 'SMARTS': 'nc(o)c', 'types': ['NtrR', 'aroC', 'Oaro', None]}, # unsure how it's written
                  {'name': 'nitrile', 'SMARTS': 'C#N', 'types': ['COO', 'NtrR']}, ## the old one uses Nhis??!
                  {'name': 'imino', 'SMARTS': 'C=N', 'types': [None, 'Nhis']},
                  {'name': 'nitro', 'SMARTS': '[N+](=O)[O-]', 'types': ['Nhis','OOC', 'OOC']},
                  {'name': 'nitro_aro', 'SMARTS': 'n(=o)o', 'types': ['Nhis','OOC', 'OOC']},
                  {'name': 'furan', 'SMARTS': 'coc', 'types': ['aroC', 'Oaro', 'aroC']},
                  {'name': 'ether', 'SMARTS': 'COC', 'types': [None, 'Oet3', None]},
                  {'name': 'hydroxyl', 'SMARTS': '[OH]', 'types': ['OH']},
                  {'name': 'guanidium', 'SMARTS': 'NC(=N)N', 'types': ['Narg', 'aroC', 'Narg', 'Narg']},
                  ]
        for group in groups:
            template = Chem.MolFromSmarts(group['SMARTS'])
            template = Chem.AddHs(template, explicitOnly=True)
            types = group['types']
            for match in self.mol.GetSubstructMatches(template):
                for i, n in enumerate(types):
                    j = match[i]
                    atom = self.mol.GetAtomWithIdx(j)
                    if n is not None and not (atom.HasProp('_rType') and atom.GetProp('_rType').strip()):
                        atom.SetProp('_rType', n)
        # Simple atoms
        for atom in self.mol.GetAtoms():
            symbol = atom.GetSymbol()
            Hs = [a for a in atom.GetNeighbors() if a.GetSymbol() == 'H']
            if atom.HasProp('_rType') and atom.GetProp('_rType').strip():
                pass
            elif symbol == '*':
                atom.SetProp('_rType', 'VIRT')
            elif symbol == 'C':
                if atom.GetIsAromatic():
                    atom.SetProp('_rType', 'aroC')
                    for n in Hs:
                        n.SetProp('_rType', 'Haro')
                else:
                    atom.SetProp('_rType', 'CH' + str(len(Hs)))
            elif symbol == 'N':
                if atom.GetIsAromatic() and len(Hs) == 0:
                    atom.SetProp('_rType', 'Nhis')
                elif atom.GetIsAromatic():  # could also be NtrR...
                    atom.SetProp('_rType', 'Ntrp')
                elif atom.GetHybridization() == Chem.HybridizationType.SP3 and len(Hs) < 3:
                    atom.SetProp('_rType', 'Npro')  # Nbb is a SP2
                elif atom.GetHybridization() == Chem.HybridizationType.SP3:
                    atom.SetProp('_rType', 'Nlys')
                elif atom.GetHybridization() == Chem.HybridizationType.SP2:
                    atom.SetProp('_rType', 'NH2O')  # or Narg if the orbitals are funky...
                else:
                    raise ValueError(f'No idea what this nitrogen {atom.GetHybridization()} is')
            elif symbol == 'O':
                if atom.GetIsAromatic():
                    atom.SetProp('_rType', 'Oaro')
                elif atom.GetHybridization() == Chem.HybridizationType.SP2:
                    atom.SetProp('_rType', 'Oet2')
                elif atom.GetHybridization() == Chem.HybridizationType.SP3:
                    atom.SetProp('_rType', 'Oet3')
                else:
                    raise ValueError(f'No idea what this oxygen {atom.GetHybridization()} is')
            elif symbol == 'H':
                n = atom.GetNeighbors()[0]
                if n.GetSymbol() == 'C' and n.GetIsAromatic():
                    atom.SetProp('_rType', 'Haro')
                elif n.GetSymbol() == 'C':
                    atom.SetProp('_rType', 'Hapo')
                elif n.GetSymbol() == 'S':
                    atom.SetProp('_rType', 'HS')
                else:
                    atom.SetProp('_rType', 'Hpol')
            elif symbol == 'S':
                if len(Hs) == 1:
                    atom.SetProp('_rType', 'SH1')
                else:
                    atom.SetProp('_rType', 'S')
            elif symbol == 'P':
                atom.SetProp('_rType', 'Phos')
            elif symbol.upper() in element_to_type:
                atom.SetProp('_rType', element_to_type[symbol.upper()])
            else:
                self.log.warning(f'No idea what this {atom.GetSymbol()} {atom.GetHybridization()} is. assigning it REPLS')
                atom.SetProp('_rType', 'REPLS')

    def _add_genrtypes(self) -> None:
        """
        Add Rosetta Atom types to each atom.
        """
        groups = [{'name': 'carboxylate', 'SMARTS': 'C(=O)[O-]', 'types': ['COO', 'OOC', 'OOC']},
                  {'name': 'carboxylic', 'SMARTS': 'C(=O)[OH]', 'types': ['COO', 'OOC', 'OH']},
                  {'name': 'ester', 'SMARTS': 'C(=O)OC', 'types': ['COO', 'Oal', 'Oet', None]},
                  {'name': 'ester', 'SMARTS': 'C(=O)O', 'types': ['COO', 'OOC', 'OOC']},
                  {'name': '3amide', 'SMARTS': 'C(=O)[NH0]', 'types': ['CNH2', 'ONH2', 'Nad3']},
                  {'name': 'amide', 'SMARTS': 'C(=O)N', 'types': ['CNH2', 'ONH2', 'NH2O']},
                  {'name': 'keto-aldehyde', 'SMARTS': 'C(=O)', 'types': [None, 'Oal']},
                  {'name': '2amine', 'SMARTS': 'CNC', 'types': [None, 'Nam2', None]},
                  {'name': 'azo', 'SMARTS': 'NN', 'types': ['Nad3', 'Nad3']},
                  {'name': 'nitrile', 'SMARTS': 'C#N', 'types': ['CTp', 'NG1']},
                  {'name': 'imino', 'SMARTS': 'C=[NH]', 'types': [None, 'Nin']},
                  {'name': '2imino', 'SMARTS': 'C=[NH0]', 'types': [None, 'Nim']},
                  {'name': '2amine', 'SMARTS': 'CNC', 'types': [None, 'Nam', None]},
                  {'name': 'nitro', 'SMARTS': '[N+](=O)[O-]', 'types': ['NGb', 'Ont', 'Ont']},
                  {'name': 'nitro_aro', 'SMARTS': 'n(=o)o', 'types': ['NGb', 'Ont', 'Ont']},
                  {'name': 'furan', 'SMARTS': 'coc', 'types': ['aroC', 'Ofu', 'aroC']},
                  {'name': 'ether', 'SMARTS': 'COC', 'types': [None, 'Oet', None]},
                  {'name': 'hydroxyl', 'SMARTS': '[OH]', 'types': ['OH']},
                  {'name': 'guanidium', 'SMARTS': 'NC(=N)N', 'types': ['Narg', 'aroC', 'Narg', 'Narg']}
                  ]
        for group in groups:
            template = Chem.MolFromSmarts(group['SMARTS'])
            types = group['types']
            for match in self.mol.GetSubstructMatches(template):
                for i, n in enumerate(types):
                    j = match[i]
                    atom = self.mol.GetAtomWithIdx(j)
                    if n is not None and not (atom.HasProp('_rType') and atom.GetProp('_rType').strip()):
                        atom.SetProp('_rType', n)
        # Generic atoms
        for atom in self.mol.GetAtoms():
            symbol = atom.GetSymbol()
            Hs = [a for a in atom.GetNeighbors() if a.GetSymbol() == 'H']
            if atom.HasProp('_rType') and atom.GetProp('_rType').strip():
                pass
            elif symbol == '*':
                atom.SetProp('_rType', 'VIRT')
            elif symbol == 'C':
                if atom.GetIsAromatic():
                    atom.SetProp('_rType', 'aroC')
                    for n in Hs:
                        n.SetProp('_rType', 'Haro')
                elif atom.GetHybridization() == Chem.HybridizationType.SP3:
                    if len(Hs) > 0:
                        atom.SetProp('_rType', 'CH' + str(len(Hs)))
                    else:
                        atom.SetProp('_rType', 'CS')
                elif atom.GetHybridization() == Chem.HybridizationType.SP2:
                    if len(Hs) > 0:
                        atom.SetProp('_rType', 'CD' + str(len(Hs)))
                    else:
                        atom.SetProp('_rType', 'CD')
                elif atom.GetHybridization() == Chem.HybridizationType.SP:
                    atom.SetProp('_rType', 'CT')
                else:
                    raise ValueError(f'No idea what this carbon {atom.GetHybridization()} is')
            elif symbol == 'N':
                if atom.GetIsAromatic():
                    atom.SetProp('_rType', 'NGb')
                elif atom.GetHybridization() == Chem.HybridizationType.SP3 and len(Hs) == 0:
                    atom.SetProp('_rType', 'NG3')
                elif atom.GetHybridization() == Chem.HybridizationType.SP3:
                    atom.SetProp('_rType', 'Nam')
                elif atom.GetHybridization() == Chem.HybridizationType.SP2 and len(Hs) == 0:
                    atom.SetProp('_rType', 'NG2')
                elif atom.GetHybridization() == Chem.HybridizationType.SP2 and len(Hs) == 1:
                    atom.SetProp('_rType', 'NG21')
                elif atom.GetHybridization() == Chem.HybridizationType.SP2 and len(Hs) == 2:
                    atom.SetProp('_rType', 'NG22')
                elif atom.GetHybridization() == Chem.HybridizationType.SP and len(Hs) == 0:
                    atom.SetProp('_rType', 'NG1')
                else:
                    raise ValueError(f'No idea what this nitrogen {atom.GetHybridization()} is')
            elif symbol == 'O':
                if atom.GetIsAromatic():
                    atom.SetProp('_rType', 'Oaro')
                elif atom.GetHybridization() == Chem.HybridizationType.SP2 and len(Hs) == 0:
                    atom.SetProp('_rType', 'OG2')
                elif atom.GetHybridization() == Chem.HybridizationType.SP3 and len(Hs) == 0:
                    atom.SetProp('_rType', 'OG3')
                elif atom.GetHybridization() == Chem.HybridizationType.SP3 and len(Hs) == 1:
                    atom.SetProp('_rType', 'OG31')
                else:
                    raise ValueError(f'No idea what this oxygen {atom.GetHybridization()} is')
            elif symbol == 'H':
                n = atom.GetNeighbors()[0]
                if n.GetSymbol() == 'C' and n.GetIsAromatic():
                    atom.SetProp('_rType', 'Haro')
                elif n.GetSymbol() == 'C':
                    atom.SetProp('_rType', 'Hapo')
                elif n.GetSymbol() == 'S':
                    atom.SetProp('_rType', 'HS')
                elif n.GetSymbol() == 'N':
                    atom.SetProp('_rType', 'HN')
                elif n.GetSymbol() == 'O':
                    atom.SetProp('_rType', 'HO')
                else:
                    atom.SetProp('_rType', 'HG')
            elif symbol == 'S':
                if len(Hs) == 1:
                    atom.SetProp('_rType', 'Sth')
                elif len(Hs) == 0:
                    atom.SetProp('_rType', 'Ssl')
                elif atom.GetHybridization() == Chem.HybridizationType.SP3 and atom.GetExplicitValence() == 6:
                    atom.SetProp('_rType', 'SG5')
                elif atom.GetHybridization() == Chem.HybridizationType.SP3:
                    atom.SetProp('_rType', 'SG3')
                elif atom.GetHybridization() == Chem.HybridizationType.SP2:
                    atom.SetProp('_rType', 'SG2')
                else:
                    raise ValueError(f'No idea what this S {atom.GetHybridization()} is')
            elif symbol == 'P':
                if atom.GetHybridization() == Chem.HybridizationType.SP3 and atom.GetExplicitValence() == 6:
                    atom.SetProp('_rType', 'PG5')
                elif atom.GetHybridization() == Chem.HybridizationType.SP3:
                    atom.SetProp('_rType', 'PG3')
                else:
                    raise ValueError(f'No idea what this S {atom.GetHybridization()} is')
            elif symbol in ('F', 'Cl', 'Br', 'I'):
                n = atom.GetNeighbors()[0]
                if n.GetIsAromatic():
                    atom.SetProp('_rType', symbol + 'R')
                else:
                    atom.SetProp('_rType', symbol)
            else:
                raise ValueError(f'No idea what this {atom.GetSymbol()} {atom.GetHybridization()} is')

    def _aminoacid_override(self, elemental):
        aa = Chem.MolFromSmiles('*NCC(~O)*')
        if self.mol.HasSubstructMatch(aa):
            self.log.info('Ligand detected to be polymer!')
            self.TYPE.append('POLYMER')
            aa_map = dict([('LOWER', None),
                           (' N  ', 'Nbb'),
                           (' CA ', 'CAbb'),
                           (' C  ', 'CObb'),
                           (' O  ', 'OCbb'),
                           ('UPPER', None)
                           ])
            self.rename_by_substructure(aa, list(aa_map.keys()))
            # fix amine and H for secondary (proline)
            amine = self.get_atom_by_name(' N  ')
            hs = [n for n in amine.GetNeighbors() if n.GetSymbol() == 'H']
            if len(hs):
                self._set_PDBInfo_atomname(hs[0], ' H  ', overwrite=True)
                aa_map[' H  '] = 'HNbb'
            else:
                aa_map[' N  '] = 'Npro'
            ca = self.get_atom_by_name(' CA ')
            for neigh in ca.GetNeighbors():
                if neigh.GetSymbol() == 'H':
                    self.rename_atom(neigh, ' HA ')
            if self.greekification:
                self.greekify()
            # add rtypes
            self.retype_by_name(aa_map)
            # change conn
            elemental['CONN'] = 2  # when LOWER and UPPER exist the CONN is CONN3.
            # add properties
            self.PROPERTIES.append('PROTEIN')
            self.PROPERTIES.append('ALPHA_AA')
            self.PROPERTIES.append('L_AA')
            # HYDROPHOBIC ALIPHATIC
            self.FIRST_SIDECHAIN_ATOM.append('CB')
            self.BACKBONE_AA.append('ALA')

    def greekify(self, ascii=True):
        """
        Converts the atom names into relative names, i.e. using the Greek alphabet.

        Its name is terrible, but ``hellenify`` or ``translitterate`` are too confusing.

        If ascii is True, these will be Greek-like letters else
        they will be actual Greek Unicode letters if it is not -- A really terrible idea.
        """
        if ascii:
            greek = list('ABGDEZHTIKLMNXOPRS')
        else:
            greek = list('ΑΒΓΔΕΖΗΘΙΚΛΜΝΞΟΠΡΣΤΥΦΧΨΩ')

        greekdex = defaultdict(list)
        ca = self.get_atom_by_name('CA')
        for atom in self.mol.GetAtoms():
            is_backbone = (atom.GetPDBResidueInfo() is not None and
                    atom.GetPDBResidueInfo().GetName().strip() in ('LOWER', 'UPPER', 'N', 'CA', 'C', 'H', 'HA', 'O', 'OXT'))
            if atom.GetSymbol() != 'H' and not is_backbone:
                n = len(Chem.GetShortestPath(self.mol, ca.GetIdx(), atom.GetIdx())) - 1
                greekdex[n].append(atom)
        for k in greekdex:
            if len(greek) <= k:
                pass # finished the greek!
            elif len(greekdex[k]) == 0:
                pass # impossible tho
            elif len(greekdex[k]) == 1:
                name = f'{greekdex[k][0].GetSymbol(): >2}{greek[k]} '
                self.rename_atom(greekdex[k][0], name)
            elif len(greekdex[k]) < 36:
                l = list(string.digits+string.ascii_uppercase)[1:]
                for i, atom in enumerate(greekdex[k]):
                    name = f'{atom.GetSymbol(): >2}{greek[k]}{l[i]}'
                    self.rename_atom(greekdex[k][i], name)
            else:
                pass  # no renaming. This is insane corner case. A 36 HA AA is madness.

    def _fix_atom_names(self):
        elemental = defaultdict(int)
        seen = []
        # Amino acid overwrite.
        self._aminoacid_override(elemental)
        for i in range(self.mol.GetNumAtoms()):
            atom = self.mol.GetAtomWithIdx(i)
            el = atom.GetSymbol().upper()
            if el == '*':
                el = 'CONN'
            elemental[el] += 1  # compatible mol_to_params.py
            lamename = el + str(elemental[el])
            lamename = self.pad_name(lamename, atom)
            while lamename in seen:
                elemental[el] += 1
                lamename = el + str(elemental[el])
            name = self._set_PDBInfo_atomname(atom, lamename, overwrite=False)
            if name in seen:
                self.log.warning(f'Name clash {name}, second one now called {lamename}')
                atom.GetPDBResidueInfo().SetName(lamename)
                seen.append(lamename)
            else:
                seen.append(name)

    def _add_partial_charges_OLD(self):
        """
        This is pointless convoluted.
        :return:
        """
        raise DeprecationWarning
        demovalence = {1: 1, 2: 8, 3: 7, 4: 6, 6: 15}
        # Not using Chem PeriodicTable because this was meant to be a quickfix...
        mol = Chem.Mol(self.mol)
        AllChem.ComputeGasteigerCharges(mol, throwOnParamFailure=False)
        tries = 0
        while True:
            if tries > 5:
                AllChem.ComputeGasteigerCharges(mol, throwOnParamFailure=False)
            try:
                AllChem.ComputeGasteigerCharges(mol, throwOnParamFailure=True)
            except ValueError as err:
                warn(f'{err.__class__.__name__}: {err}')
                dodgy = re.search('parameters for Element\: ([*\w]+)', str(err)).group(1)
                for atom in mol.GetAtoms():
                    if atom.GetSymbol() == dodgy:
                        valence = atom.GetExplicitValence()
                        atom.SetAtomicNum(demovalence[valence])
                tries += 1
            else:
                break
        for i in range(self.mol.GetNumAtoms()):
            gc = mol.GetAtomWithIdx(i).GetDoubleProp('_GasteigerCharge')
            if str(gc) != 'nan':
                self.mol.GetAtomWithIdx(i).SetDoubleProp('_GasteigerCharge', gc)
            else:
                self.mol.GetAtomWithIdx(i).SetDoubleProp('_GasteigerCharge', 0.)

    def _add_partial_charges(self):
        changed = []
        for atom in self.mol.GetAtoms():
            if atom.GetSymbol() == '*':
                atom.SetAtomicNum(6)
                atom.SetHybridization(Chem.HybridizationType.SP3)
                changed.append(atom.GetIdx())
        Chem.SanitizeMol(self.mol)
        AllChem.ComputeGasteigerCharges(self.mol, throwOnParamFailure=False)
        for i, atom in enumerate(self.mol.GetAtoms()):
            if i in changed:
                atom.SetAtomicNum(0)
            gc = atom.GetDoubleProp('_GasteigerCharge')
            if str(gc) != 'nan':
                atom.SetDoubleProp('_GasteigerCharge', 0.)

    def _get_resn_from_PDBInfo(self):
        """
        Gets the residue name for PDB info.

        :return:
        """
        infos = [atom.GetPDBResidueInfo() for atom in self.mol.GetAtoms()]
        names = [i.GetResidueName() for i in infos if i is not None]
        if names:
            return names[0]
        else:
            return 'LIG'

    @classmethod
    def pad_name(self, name, atom=None):
        if name in ('CONN1','CONN2','CONN3','CONN4', 'LOWER', 'UPPER'):
            return name
        elif len(name) == 4:
            return name
        elif len(name) > 4:
            return name[:4]
        elif name[0] == ' ':
            return name.ljust(4)
        elif len(name) < 4 and (atom is None or len(atom.GetSymbol()) == 1):
            return ' ' + name.ljust(3)
        else:
            return name.ljust(4)

    def move_aside(self):
        """
        Changes the names of the atoms to not clash.

        :return:
        """
        for i in range(self.mol.GetNumAtoms()):
            atom = self.mol.GetAtomWithIdx(i)
            atomname = self._get_PDBInfo_atomname(atom, throw=False)
            if atomname:
                atom.SetProp('_OriginalName', atomname)
            if atom.GetSymbol() != '*':
                self._set_PDBInfo_atomname(atom, f'XX{i: <2}', overwrite=True)

    def move_back(self):
        """
        Removes the ugly XX!

        :return:
        """
        atomnames = []
        for atom in self.mol.GetAtoms():
            atomname = self._get_PDBInfo_atomname(atom, throw=False)
            atomnames.append(atomname)
        # find bad ones!
        for atom in self.mol.GetAtoms():
            atomname = self._get_PDBInfo_atomname(atom, throw=False)
            if 'XX' in atomname:
                original = atom.GetProp('_OriginalName')
                if original not in atomnames:
                    self._set_PDBInfo_atomname(atom, original, overwrite=True)
                else:
                    for i in range(1, 100):
                        candidate = f'{atom.GetSymbol: >2}{i: <2}'
                        if candidate not in atomnames:
                            self._set_PDBInfo_atomname(atom, candidate, overwrite=True)
                            break
                    else:
                        self.log.warning(f'Could not find a decent atomname for {atom.GetIdx()}')

    def add_Hs(self):
        """
        Add Hs before convert_mol step!
        :return:
        """
        self.log.debug('Adding hydrogens')
        self.mol = AllChem.AddHs(self.mol)
        changed = []
        Chem.SanitizeMol(self.mol)
        for atom in self.mol.GetAtoms():  # TODO switch to DummyMasker
            if atom.GetSymbol() == '*':
                atom.SetAtomicNum(6)
                atom.SetHybridization(Chem.HybridizationType.SP3)
                changed.append(atom.GetIdx())
        AllChem.EmbedMolecule(self.mol, useRandomCoords=True)
        AllChem.MMFFOptimizeMolecule(self.mol)
        AllChem.ComputeGasteigerCharges(self.mol, throwOnParamFailure=False)
        for i, atom in enumerate(self.mol.GetAtoms()):
            if i in changed:
                atom.SetAtomicNum(0)
        self.fix_mol()


