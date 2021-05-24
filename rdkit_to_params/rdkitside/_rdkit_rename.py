from rdkit import Chem
from rdkit.Chem import rdFMCS, AllChem

from typing import Sequence, List, Dict, Optional, Any, Union
from warnings import warn
import re, logging


########################################################################################################################
__doc__ = \
    """
The main class here is ``_RDKitRenameMixin``, which adds the various atom renaming methods.

    """
from ..version import *

########################################################################################################################

class _RDKitRenameMixin:
    log = logging.getLogger(__name__)

    # ============= overridden =========================================================================================

    def __init__(self):
        # will be overridden. Just here for typehinting...
        self.log.critical('WRONG INIT CALLED')
        self.mol = Chem.Mol()
        self.NAME = ''
        self.is_aminoacid = lambda: True # main class
        self.polish = lambda: None # covert class. makes sure all atoms have the same rdkit residue info
        self.rename_atom = lambda oldname, newname, overwrite=False: None
        self.pad_name = lambda name: name

    # ============= set and get ========================================================================================

    def _get_PDBInfo_atomname(self, atom, throw=True) -> str:
        info = atom.GetPDBResidueInfo()
        if info is not None:
            return info.GetName()
        elif throw:
            raise ValueError('Atoms changed but `fix_mol` was not called.')
        else:
            return ''

    def _set_PDBInfo_atomname(self, atom, name, overwrite=False):
        info = atom.GetPDBResidueInfo()
        # if has none, it is None
        if info is None:
            isHeteroAtom = not self.is_aminoacid()
            atom.SetMonomerInfo(Chem.AtomPDBResidueInfo(atomName=self.pad_name(name),
                                                        serialNumber=atom.GetIdx(),
                                                        residueName=self.NAME,
                                                        isHeteroAtom=isHeteroAtom))
            return name
        elif info.GetName() == name:
            return name
        elif overwrite:
            info.SetName(name)
            return name
        else:
            return info.GetName()


    def get_atom_by_name(self, name):
        for atom in self.mol.GetAtoms():
            info = atom.GetPDBResidueInfo()
            if info and info.GetName().strip() == name.strip():
                return atom
        else:
            raise ValueError(f'Name {name} not found.')

    # ============= retype method ======================================================================================

    def retype_by_name(self, mapping: Dict[str, str]):
        """
        Renames the atom names in ``self.mol`` by changing key to value of the mapping dictionary.

        :param mapping: old name to new name
        :return:
        """
        for name, rtype in mapping.items():
            if rtype is None:
                continue
            atom = self.get_atom_by_name(name)
            if atom is not None:
                atom.SetProp('_rType', rtype)
            else:
                self.log.warning(f'Could not find "{name}"')
        for entry in self.ATOM:
            if entry.name.strip() == name.strip():
                entry.rtype = rtype
                break
        else:
            pass

    # ============= rename methods =====================================================================================

    def rename(self, atomnames: Union[None, dict, list, str, Chem.Mol]) -> None:
        """
        Rename options for atomnames:

        * None: nothing
        * Dict: ``rename_from_dict``
        * List: ``rename_from_list``
        * str in '1:XX, 3:YY' format: ``rename_from_dict`` via ``rename_from_str``
        * str in 'XX,YY' format: ``rename_from_list`` via ``rename_from_str``
        * Chem.Mol: ``rename_from_template``

        :param atomnames: renaming values
        :return:
        """
        if atomnames is None:
            pass
        elif isinstance(atomnames, dict):
            self.rename_from_dict(atomnames)
        elif isinstance(atomnames, list):
            self.rename_from_list(atomnames)
        elif isinstance(atomnames, str):
            self.rename_from_str(atomnames)
        elif isinstance(atomnames, Chem.Mol):
            self.rename_from_template(atomnames)
        else:
            raise ValueError(f'atomnames is not None, dict, list, str or Chem.Mol but {type(atomnames)}')


    def rename_by_substructure(self, substructure: Chem.Mol, atomnames: Sequence[str]) -> List[str]:
        """
        Assigns to the atoms in self.mol the names based on the backbone template and the names variable.
        See ``_fix_atom_names`` for example usage.
        Changes also params.

        :param substructure:
        :param atomnames: the list oof new names. Falsey names will not be set.
        :return: the list of the old names.. why? May change in future.
        """
        assert self.mol.HasSubstructMatch(substructure), 'Bad backbone match'
        originals = []
        substructure = Chem.AddHs(substructure, explicitOnly=True)
        for name, idx in zip(atomnames, self.mol.GetSubstructMatch(substructure)):
            atom = self.mol.GetAtomWithIdx(idx)
            oldname = self._get_PDBInfo_atomname(atom, throw=False)
            self.rename_atom(atom, name)  # alters entry
            originals.append(oldname)
        return originals

    def rename_from_str(self, atomnames:str):
        """
        * '1:XX, 3:YY' format: calls ``rename_from_dict``
        * 'XX,YY' format: calls ``rename_from_list``

        :param atomnames: str
        :return:
        """
        if ':' in atomnames:
            self.rename_from_dict({int(k): v for k, v in re.findall(r'(\d+):([\s\w]{1,4})', atomnames)})
        else:
            self.rename_from_list(atomnames.split(','))

    def rename_from_template(self, template: Chem.Mol, overwrite:bool=True):
        """
        Assigns to the atoms in self.mol the names based on the template, which does not need to be a perfect match.
        See ``_fix_atom_names`` for example usage.
        Does not change the Params.

        :param template: mol object with atom names

        :return: None for now.
        """
        AllChem.SanitizeMol(template) #this is where half my issues come from.
        mcs = rdFMCS.FindMCS([self.mol, template],
                             atomCompare=rdFMCS.AtomCompare.CompareElements,
                             bondCompare=rdFMCS.BondCompare.CompareAny,
                             ringMatchesRingOnly=True)
        common = Chem.MolFromSmarts(mcs.smartsString)
        for acceptor, donor in zip(self.mol.GetSubstructMatch(common), template.GetSubstructMatch(common)):
            a_atom = self.mol.GetAtomWithIdx(acceptor)
            d_atom = template.GetAtomWithIdx(donor)
            info = d_atom.GetPDBResidueInfo()
            if info:
                self.rename_atom(a_atom, info.GetName(), overwrite=overwrite)
            else:
                self.log.debug.info(f'No info in template for atom {d_atom.GetSymbol()} #{donor}')

    def rename_from_dict(self, atomnames: Dict[int, str]):
        """
        Renames the ``self.mol`` atom names by a dict that has key atom idx and value name

        :param atomnames: idx -> new name
        :return:
        """
        # sanity check
        assert len(set(atomnames.values())) == len(atomnames), 'Atom Names are repeated.'
        if self.mol.GetNumAtoms() > len(atomnames):
            self.log.info('There are more atoms in mol than were provided.')
        elif self.mol.GetNumAtoms() < len(atomnames):
            #raise ValueError('There are less atoms in mol than were provided.')
            pass # this is fine.
        # operate
        for k, v in atomnames.items():
            self._set_PDBInfo_atomname(self.mol.GetAtomWithIdx(k), v, overwrite=True)
        self.polish_mol()

    def rename_from_list(self, atomnames: List[str]):
        """
        Renames the ``self.mol`` atom names based on the order in list. If None it is skipped.

        :param atomnames: list of unique names.
        :return:
        """
        # sanity check
        assert len(set(atomnames)) == len(atomnames), 'Atom Names are repeated.'
        if self.mol.GetNumAtoms() > len(atomnames):
            self.log.info('There are more atoms in mol than were provided.')
        elif self.mol.GetNumAtoms() < len(atomnames):
            raise ValueError('There are less atoms in mol than were provided.')
        # rename
        r = {i: name for i, name in enumerate(atomnames) if name is not None}
        self.rename_from_dict(r)

    # ============= class method =======================================================================================

    @classmethod
    def add_names(cls, mol: Chem.Mol, atomnames: Union[None, dict, list, str, Chem.Mol], name: Optional[str] = None) -> Chem.Mol:
        """
        Quick way to add atom names to a mol object --adds them the normal way.

        :param mol: Chem.Mol, will actually be edited in place.
        :param names: list/dict/str/chem of unique names.
        :param name: 3letter code for the molecule.
        :return: the mol
        """
        self = cls()

        self.log.debug(f'`add_names` called...')
        if name is not None:
            self.NAME = name
        self.mol = mol
        self.TYPE.append('LIGAND')
        self.fix_mol()
        self.rename(atomnames)
        return self.mol