########################################################################################################################
__doc__ = \
    """
    The main class here is ``Params``. All underscore base classes are not meant to be used standalone.
    ``Entries`` is the class for an list of entries of the same kind.
    """

__author__ = "Matteo Ferla. [Github](https://github.com/matteoferla)"
__email__ = "matteo.ferla@gmail.com"
__date__ = "10 July 2020 A.D."
__license__ = "MIT"
__version__ = "1.1.0"
__citation__ = "None."

########################################################################################################################


from warnings import warn
import os, re, logging
from typing import Union

#################### base classes ######################################################################################

from ._io_mixin import _ParamsIoMixin # in turn inherits _ParamsInitMixin
from .entries import Entries

#################### pyrosetta #########################################################################################
try:
    from ._pyrosetta_mixin import _PoserMixin, pyrosetta
except ImportError:
    warn('PyRosetta is required for the ``test`` method', category=ImportWarning)


    class _PoserMixin:
        pass

#################### rdkit #############################################################################################
try:
    from .rdkitside import _RDKitMixin
    from .constraint import Constraints
    from rdkit import Chem
except ImportError:
    warn('RDkit is required for ``from_mol`` stuff and ``Constraints``', category=ImportWarning)

    class Chem:
        Atom = None

    class _RDKitMixin:
        pass


#################### main class ########################################################################################


class Params(_ParamsIoMixin, _RDKitMixin, _PoserMixin):
    """

    ``Params`` creates and manipulates params files. It can handles several types of params operations,
    such as "atom name surgery" and ``rdkit.Chem.Mol`` to a params file.

    ## Key methods

    * ``Params.load(filename)`` will instantiate from file.
    * ``Params.from_mol(mol)`` will instantiate from ``Chem.Mol``
    * ``p.dump(filename)`` will save a file.
    * ``loads`` and ``dumps``for strings.
    * ``p.fields`` will return all header fields.
    * ``p.test`` tests the params file in PyRosetta.
    * ``p.rename_atom(old, new)`` changes an atom name

    ## Attributes

    The attributes are generally the uppercase line headers, with a few exceptions.

    * `.comments` is for # lines
    * "BOND_TYPE" and "BOND" are merged into ``.BOND``.
    * "UPPER", "LOWER" and "CONNECT" are merged into ``.CONNECT``

    With the exception of ``.NAME`` which depends on ``.IO_STRING`` basically
    all the header type attributes are actually instances of the class `Entries`, which holds a sequence of specific entries.
    see `entry.py` for the properties of each.
    These can be a singleton, such as `.IO_STRING` which when a new line is added it gets overwritten instead, or not like say `.ATOM`.
    That is to say that ``.ATOM[0]`` will give the first atom as expected, but this has to be done for ``.IO_STRING[0]`` too.

    Atomnames...

    * ``p.get_correct_atomname`` will return the 4 letter name of the atom with nice spacing.
    * ``p.rename_atom`` will change one atomname to a new one across all entries.
    * ``BOND``, ``CHI``, ``CUT_BOND`` entries store 4 char atomnames as ``.first``, ``.second``, ``.third``, ``.fourth``.
    * ``ICOOR_INTERNAL`` entries store 5 char atomnames as ``.child``,``.parent``,``.second_parent``,``.third_parent``.
    * ``ATOM_ALIAS``, ``NBR_ATOM``, ``FIRST_SIDECHAIN_ATOM``, ``ADD_RING`` are just ``entries.GenericEntries`` instances, where ``.body`` is a string which will contain the atomname.
    * ``METAL_BINDING_ATOMS``, ``ACT_COORD_ATOMS`` are ``entries.GenericListEntries`` instances where ``.values`` is a list of string with maybe atomnames.

    ## Inheritance

    It inherits several class,
    which are not not mean to be used standalone, except for testing.

    The pyrosetta and rdkit functionality are dependent on these being installed.

    * ``_ParamsIoMixin`` adds read write, and inherits
    * ``_ParamsInitMixin`` which adds the basics.
    * ``_PoserMixin`` is a base that adds pyrosetta functionality if avaliable.
    * ``_RDKitCovertMixin``, which adds rdkit ``from_mol`` conversion functionality, the class is split in two, the other part being
    * ``_RDKitParamsPrepMixin``, which prepares the molecule for `_RDKitCovertMixin.from_mol``.

    """
    log = logging.getLogger(__name__)

    @property
    def NAME(self):
        if len(self.IO_STRING):
            return self.IO_STRING[0].name3
        else:
            self.log.warning('Attempted access to empty IO_STRING/NAME')
            return 'XXX'

    @NAME.setter
    def NAME(self, name):
        if len(self.IO_STRING) == 0:
            self.IO_STRING.append(f'{name} Z')
        else:
            self.IO_STRING[0].name3 = name

    def is_aminoacid(self):
        if len(self.TYPE) == 0:
            self.TYPE.append('LIGAND')
        return self.TYPE[0].body == 'POLYMER'

    def validate(self):
        self.log.critical('This is not finished.')
        if 'CANONICAL_AA' in self.PROPERTIES[0].values:
            assert self.AA != 'UNK', 'CANONICAL_AA property requires a AA type not UNK'
        if 'METALBINDING' in self.PROPERTIES[0].values:
            assert len(self.METAL_BINDING_ATOMS) > 0, 'METALBINDING property requires METAL_BINDING_ATOMS'
        assert os.path.exists(self.PDB_ROTAMERS.strip()), f'PDB_ROTAMERS file {self.PDB_ROTAMERS} does not exist'
        # etc.

    def get_correct_atomname(self, name: str) -> str:
        """
        Given a name, gets the correctly spaced out one as appears in the ATOM entry.
        To pad out a name use pad_name
        This has nothing to do with ``._get_PDBInfo_atomname`` which just returns the atom name from a ``Chem.Atom``.

        :param name: dirty name
        :return: correct name
        """
        name = name.upper()
        # find in atom.
        for atom in self.ATOM:
            if atom.name == name:
                return name
        else:
            for atom in self.ATOM:
                if atom.name.strip() == name.strip():
                    return atom.name
            else:
                raise ValueError(f'{name} is not a valid atom name (does not appear in the entries)')

    def rename_atom(self, atom_or_atomname: Union[str, 'Chem.Atom'], newname: str, overwrite=True) -> Union[str, None]:
        """
        rename an atom by atomname or Chem.Atom (the former just calls ``rename_atom_by_name`` as is just for legacy)

        :param atom_or_atomname:
        :param newname:
        :return:
        """
        # sanity
        if newname is None:
            return None
        try:
            if self.mol:
                atom = self.get_atom_by_name(newname)
                if isinstance(atom_or_atomname, str):
                    raise AssertionError(f'New name {newname} already exists')
                elif isinstance(atom_or_atomname, Chem.Atom) and atom_or_atomname.GetIdx() != atom.GetIdx():
                    raise AssertionError(f'New name {newname} already exists')
                else:
                    pass # already changed.
            if len(self.ATOM) > 0:
                self.get_correct_atomname(newname)
        except ValueError:
            pass # absent
        # change.
        if isinstance(atom_or_atomname, str): #atom name
            oldname = atom_or_atomname
            return self.rename_atom_by_name(oldname, newname)
        elif isinstance(atom_or_atomname, Chem.Atom):
            atom = atom_or_atomname
            oldname = self._get_PDBInfo_atomname(atom, throw=False)
            if oldname:
                return self.rename_atom_by_name(oldname, newname)  # alters entry & rdkit
            else:
                return self._set_PDBInfo_atomname(atom, newname, overwrite=overwrite) # alters rdkit
        else:
            raise TypeError(f'{type(atom_or_atomname)} is not a string or atom')


    def rename_atom_by_name(self, oldname: str, newname: str) -> str:
        """
        Change the atom name from ``oldname`` to ``newname`` and returns the 4 char ``newname``.

        :param oldname: atom name, preferably 4 char long.
        :param newname: atom name, preferably 4 char long.
        :return: 4 char newname
        """
        if newname is None:
            return None
        elif oldname == newname:
            return newname
        else:
            # rdkit mol
            if self.mol:
                atom = self.get_atom_by_name(oldname)
                newname = self.pad_name(newname, atom)
                atom.GetPDBResidueInfo().SetName(newname)
            else:
                newname = self.pad_name(newname)
            # params
            self._rename_atom_in_entries(oldname, newname)
            return newname

    def _rename_atom_in_entries(self, oldname, newname):
        # if params is not filled nothing happens.
        # check if it is a connect atom
        if len(self.ATOM) == 0:
            return None # N/A: unparameterise atm.
        elif oldname.strip() == 'CONN':
            pass  # ...
        elif oldname.strip() in ('CONN1', 'CONN2', 'CONN3', 'LOWER', 'UPPER'):
            for conn in self.CONNECT:
                if conn.connect_name.strip() == oldname.strip():
                    conn.connect_name = newname
                    # TODO fix this properly. LOWER UPPER CONN3 should be the preferred order.
                    break
            else:
                raise ValueError(f'{oldname} does not appear amid the connections')
            for entry in self.ICOOR_INTERNAL:
                for key in 'child', 'parent', 'second_parent', 'third_parent':
                    if getattr(entry, key) == oldname.rjust(5):
                        setattr(entry, key, newname.rjust(5))
        else:
            # fix names
            oldname = self.get_correct_atomname(oldname)
            if len(newname) > 4:
                raise ValueError(f'{newname} is too long.')
            elif newname == 'END':
                self.log.info('I thing END may be an old keyword - What is ``ACT_COORD_ATOMS``?. BEST AVOID IT.')
            newname = self.pad_name(newname).upper()
            # find in atom.
            for atom in self.ATOM:
                if atom.name == newname:
                    raise ValueError(f'{newname} is already taken.')
                elif atom.name == oldname:
                    atom.name = newname
                    break
            # find in entries of the kind with ``first``, ``second``, ``third``, ``fourth`` attributes, which are atom names 4 char
            for attr in ('BOND', 'CHI', 'CUT_BOND'):  # ICOOR_INTERNAL ATOM_ALIAS
                for entry in getattr(self, attr):
                    for key in 'first', 'second', 'third', 'fourth':
                        if hasattr(entry, key) and getattr(entry, key) == oldname:
                            setattr(entry, key, newname)
                            break
            # find ICOOR_INTERNAL entries with ``child``,``parent``,``second_parent``,``third_parent`` attributes,
            # which are atom names 5 char
            for entry in self.ICOOR_INTERNAL:
                for key in 'child', 'parent', 'second_parent', 'third_parent':
                    if getattr(entry, key).strip() == oldname.strip():
                        setattr(entry, key, newname.ljust(5))
            for conn in self.CONNECT:
                if conn.atom_name.strip() == oldname.strip():
                    conn.atom_name = newname
            # find in the Generic entries
            for attr in ('ATOM_ALIAS', 'NBR_ATOM', 'FIRST_SIDECHAIN_ATOM', 'ADD_RING'):
                for entry in getattr(self, attr):
                    if oldname.strip() in entry.body:
                        entry.body = re.sub('(?<!\w)' + oldname.strip() + '(?!\w)', newname, entry.body)
            # find in the Generic list entries
            for attr in ('METAL_BINDING_ATOMS', 'ACT_COORD_ATOMS'):
                for entry in getattr(self, attr):
                    entry.values = [v if v.strip() != oldname.strip() else newname for v in entry.values]

