########################################################################################################################
__doc__ = \
    """
The main class here is `Entries``, which is a fancy list. It gets called for each uppercase attribute 
in the initialisation of ``Params`` (which happens in ``_ParamsInitMixin`` __e.g.__ ``Entries.from_name('IO_STRING')``).
    """
from .version import *

########################################################################################################################

from dataclasses import dataclass
import re, logging
from warnings import warn

from collections import abc
from enum import Enum


# ======================================================================================================================

class Singletony(Enum):
    """
    Is the entry

    1. a singleton (e.g. ``NAME``) which can have only one value
    0. a regular multientry affair (e.g. ``ATOM``)
    2. a singleton (e.g. ``PROPERTY``) which can accept multiple values
    """
    multiton = 0  # non-singleton
    singleton = 1
    list_singleton = 2


# ----------------------------------------------------------------------------------------------------------------------

class Entries(abc.MutableSequence):
    """
    A fancy default list, where the elements are instances of whatver is in ``entry_cls``.
    It can be initialised via the class method ``from_name`` which accepst a string that has to be present in the class attribute ``.choices``.
    The ``.append`` method can work with str, list, dict or instance of the actual class it wants.
    Note that the string for the string way must be without the header to the line.
    The entry classes requires a ``from_str`` classmethod that returns an instance for this.
    They also require __str__ method as this is how the entries are converted into string.

    ``Entries.from_name('BOND')``
    """

    choices = {}  # this gets filled after each class is declared.

    def __init__(self, entry_cls, singleton: Singletony = Singletony.singleton):
        """
        The entries class is a fancy constrained list. The data is actually stored in ``.data``.

        :param entry_cls: what is the allowed class of the entries
        :param singleton: is only one entry allowed?
        """
        self.entry_cls = entry_cls
        self.singleton = Singletony(singleton)
        if self.singleton != Singletony.multiton and isinstance(self.entry_cls, GenericListEntry):
            raise TypeError(f'{type(self.entry_cls)} is incompatible with {self.singleton}')
        self.data = []

    @classmethod
    def from_name(cls, name: str):
        if name in cls.choices:
            cc, singleton = cls.choices[name]
            return cls(entry_cls=cc, singleton=singleton)
        else:
            raise KeyError(f'Name {name} is not one of {", ".join(cls.choices.keys())}')

    def __getitem__(self, index):
        return self.data[index]

    def _assign_value(self, value):
        if isinstance(value, self.entry_cls):
            return value
        elif isinstance(value, str):
            return self.entry_cls.from_str(value)
        elif isinstance(value, list):
            return self.entry_cls(*value)
        elif isinstance(value, dict):
            return self.entry_cls(**value)
        else:
            raise ValueError(f'No idea what to do with {value}')

    def __setitem__(self, index, value):
        if self.singleton == Singletony.multiton:
            self.data[index] = self._assign_value(value)
        elif self.singleton == Singletony.singleton or len(self.data) == 0:
            self.data = [self._assign_value(value)]
        else:  # assumes it's a GenericListEntry and has a value.
            neo = self._assign_value(value)
            self.data[0].values.extend(neo.values)

    def __delitem__(self, index):
        if self.singleton:
            index = 0
        del self.data[index]

    def insert(self, index, value):
        if self.singleton == Singletony.multiton:  # non-singleton
            self.data.insert(index, self._assign_value(value))
        else:
            self[0] = value

    def __len__(self):
        return len(self.data)

    def __str__(self):
        lines = []
        for entry in self.data:
            lines.append(str(entry))
        return '\n'.join(lines)


#########################################################################################################

class GenericEntry:
    """
    This is meant to be inherited. ``header`` is the entry type. body is a string.
    """
    log = logging.getLogger(__name__)

    def __init__(self, header: str, body: str):
        self.header = header.strip().upper()
        self.body = body.rstrip()
        assert self.header, f'Type is empty'
        assert self.body, f'Value is empty'

    def __str__(self) -> str:
        return f'{self.header} {self.body}'

    @classmethod
    def from_str(cls, text):
        return cls(text)


#########################################################################################################

class GenericListEntry:
    """
    This is meant to be inherited. ``header`` is the entry type. ``values`` is a list of strings.
    """

    def __init__(self, header: str, *args: str):
        self.header = header.strip().upper()
        self.values = list(args)

    def __str__(self) -> str:
        v = ' '.join(self.values)
        return f'{self.header} {v}'

    @classmethod
    def from_str(cls, text):
        return cls(*text.split())


#########################################################################################################

class NBR_ATOMEntry(GenericEntry):
    def __init__(self, body: str):
        super().__init__(header='NBR_ATOM', body=body)


Entries.choices['NBR_ATOM'] = (NBR_ATOMEntry, Singletony.singleton)


#########################################################################################################

class NBR_RADIUSEntry(GenericEntry):
    def __init__(self, body: str):
        super().__init__(header='NBR_RADIUS', body=body)


Entries.choices['NBR_RADIUS'] = (NBR_RADIUSEntry, Singletony.singleton)


#########################################################################################################

class MAINCHAIN_ATOMS(GenericListEntry):
    def __init__(self, *args: str):
        super().__init__('MAINCHAIN_ATOMS', *args)


Entries.choices['MAINCHAIN_ATOMS'] = (MAINCHAIN_ATOMS, Singletony.list_singleton)


#########################################################################################################

class CommentEntry(GenericEntry):
    def __init__(self, body: str):
        super().__init__(header='#', body=body)


Entries.choices['#'] = (CommentEntry, Singletony.multiton)
Entries.choices['comment'] = (CommentEntry, Singletony.multiton)


#########################################################################################################


class ATOM_ALIASEntry(GenericEntry):
    def __init__(self, body: str):
        super().__init__(header='ATOM_ALIAS', body=body)


Entries.choices['ATOM_ALIAS'] = (ATOM_ALIASEntry, Singletony.multiton)


#########################################################################################################

@dataclass
class IO_STRINGEntry:
    """
    * ``.name3`` is three letter name. ``Params().NAME`` is actually a dynamic attribute that uses this.
    * ``.name1`` is a one letter name.

    These get checked for length.
    """
    name3: str = 'LIG'
    name1: str = 'Z'

    def __post_init__(self):
        if len(self.name3) == 2:
            self.name3 += ' '
        assert len(self.name3) == 3, f'{self.name3} is not 3 char long'
        assert len(self.name1) == 1, f'{self.name1} is not 1 char long'

    def __str__(self) -> str:
        return f'IO_STRING {self.name3} {self.name1}'

    @classmethod
    def from_str(cls, text):
        name3, name1 = text.strip().split()
        return cls(name3, name1)


Entries.choices['IO_STRING'] = (IO_STRINGEntry, Singletony.singleton)


#########################################################################################################

@dataclass
class CONNECTEntry:
    """
    This is a mess, but it guesses what you mean.
    Deals with UPPER, LOWER and CONNECT.

    """
    atom_name: str
    index: int = 1
    connect_type: str = ''  # | 'CONNECT' | 'UPPER_CONNECT' | 'LOWER_CONNECT'
    connect_name: str = ''  # 'CONN1' | 'UPPER' | 'LOWER'

    def __post_init__(self):
        if self.connect_type and self.connect_name:
            pass
        elif not self.connect_type and not self.connect_name:
            self.connect_type = 'CONNECT'
            self.connect_name = f'CONN{self.index}'
        elif not self.connect_type and 'CONN' not in self.connect_name:
            self.connect_type = f'{self.connect_name}_CONNECT'
        elif not self.connect_type:
            self.connect_type = 'CONNECT'
        elif not self.connect_name and self.connect_type == 'CONNECT':
            self.connect_name = f'CONN{self.index}'
        elif not self.connect_name:
            self.connect_name = self.connect_type.replace('_CONNECT', '')
        else:
            raise ValueError(
                f'I missed this case ({self.connect_name}, {self.connect_type}) in this badly written method')

    def __str__(self) -> str:
        return f'{self.connect_type} {self.atom_name}'

    @classmethod
    def from_str(cls, text):
        return cls(*text.split())


Entries.choices['CONNECT'] = (CONNECTEntry, Singletony.multiton)


#########################################################################################################

@dataclass
class CHIEntry:
    index: int
    first: str
    second: str
    third: str
    fourth: str

    def __post_init__(self):
        self.fourth = self.fourth.ljust(4)

    def __str__(self) -> str:
        return f'CHI {self.index} {self.first} {self.second} {self.third} {self.fourth}'

    @classmethod
    def from_str(cls, text: str):
        # 1  C6   C5   C4   C3
        rex = re.match('(\d+) (.{4}) (.{4}) (.{4}) (.{2,4})', text)
        if rex is None:
            raise ValueError(f'CHI entry "{text}" is not formatted correctly')
        data = dict(zip(('index', 'first', 'second', 'third', 'fourth'), rex.groups()))
        return cls(**data)


Entries.choices['CHI'] = (CHIEntry, Singletony.multiton)


#########################################################################################################


@dataclass
class ICOOR_INTERNALEntry:
    """
    Lines stolen from Rosetta documentation
    >                     Child  Phi Angle    Theta        Distance   Parent  Angle  Torsion
    >   ICOOR_INTERNAL    C14  167.536810   59.880644    1.473042   N2    C11   C12

    * Child atom (A4)
    * phi angle (torsion angle between A1, A2, A3, A4)
    * theta angle (improper angle = (180 - (angle between A4, A3, A2)))
    * distance (between A4 and A3)
    * parent atom (A3)
    * angle atom (A2)
    * torsion atom (A4)
    """
    child: str
    phi: float
    theta: float
    distance: float
    parent: str
    second_parent: str
    third_parent: str

    def __post_init__(self):
        self.third_parent = self.third_parent.ljust(5)

    def __str__(self) -> str:
        return f'ICOOR_INTERNAL  {self.child: <5} {self.phi: >11.6f} {self.theta: >11.6f} {self.distance: >11.6f} ' + \
               f'{self.parent: <5} {self.second_parent: <5} {self.third_parent: <5}'

    @classmethod
    def from_str(cls, text: str):
        # position based.
        rex = re.match(' (.{5}) (.{11}) (.{11}) (.{11}) (.{5}) (.{5}) (.{1,5})', text)
        # space based... bad.
        rex2 = re.search('(\w+)\s+([-\d\.]+)\s+([-\d\.]+)\s+([-\d\.]+)\s+(\w+)\s+(\w+)\s+(\w+)', text)
        if rex:
            data = list(rex.groups())
        elif rex2:
            data = list(rex2.groups())
        else:
            raise ValueError(f'ICOOR_INTERNAL Entry "{text}" is not formatted correctly')
        for i in range(1, 4):
            data[i] = float(data[i].strip())
        return cls(*data)


Entries.choices['ICOOR_INTERNAL'] = (ICOOR_INTERNALEntry, Singletony.multiton)


#########################################################################################################


@dataclass
class BONDEntry:
    """
    dataclass class for both BOND and BOND_ENTRY. The ``__str__`` method will know based on ``.order``.
    The hash is the two atom names sorted. So BOND records with the same names will be equal.
    """
    first: str
    second: str
    order: int = 1  # 2,3, 4|ARO

    def __post_init__(self):
        self.second = self.second.ljust(4)

    def __str__(self) -> str:
        if self.order == 1 or not self.order:
            return f'BOND {self.first: >4} {self.second: >4}'
        else:
            return f'BOND_TYPE {self.first: >4} {self.second: >4} {self.order}'

    def __hash__(self):
        return hash('+'.join(sorted([self.first, self.second])))

    def __eq__(self, other):
        return hash(self) == hash(other)

    @classmethod
    def from_str(cls, text: str):
        rex = re.match('(.{4}) (.{2,4})\s?(.*)', text)
        if rex is None:
            raise ValueError(f'BOND entry "{text}" is not formatted correctly')
        data = dict(zip(('first', 'second', 'order'), rex.groups()))
        data['order'] = data['order'].strip()
        if data['order'] == '':
            data['order'] == 1
        elif data['order'] in ('ARO', '4'):
            data['order'] == 4  # ARO is also acceptable.
        elif isinstance(data['order'], int):
            pass
        else:
            data['order'] = int(data['order'].strip())
        return cls(**data)


Entries.choices['BOND'] = (BONDEntry, Singletony.multiton)


#########################################################################################################


@dataclass
class ATOMEntry:
    # PDB atom name, Rosetta AtomType, MM AtomType, and charge
    name: str
    rtype: str
    mtype: str = 'X'
    partial: float = 0

    def __str__(self) -> str:
        return f'ATOM {self.name: >4} {self.rtype: >4} {self.mtype: >4} {self.partial:.7f}'

    def __eq__(self, other):
        """
        ``atomentry == 'CA'`` will return false because ``atomentry.name.strip() == 'CA'`` will return true.
        """
        if isinstance(other, self.__class__):
            return self.name == other.name
        else:
            return False

    def __hash__(self):
        return hash(self.name)

    @classmethod
    def from_str(cls, text: str):
        rex = re.match('(.{4}) (.{4}) (.{4}) +([-\d\.]+)', text)
        if rex is None:
            raise ValueError(f'ATOM entry "{text}" is not formatted correctly')
        data = dict(zip(('name', 'rtype', 'mtype', 'partial'), rex.groups()))
        data['partial'] = float(data['partial'])
        return cls(**data)


Entries.choices['ATOM'] = (ATOMEntry, Singletony.multiton)


#########################################################################################################

@dataclass
class CUT_BONDEntry:
    """
    No idea what CUT_BOND is for.
    """
    first: str
    second: str

    def __post_init__(self):
        self.second = self.second.ljust(4)

    def __str__(self) -> str:
        return f'CUT_BOND {self.first: >4} {self.second: >4}'

    @classmethod
    def from_str(cls, text: str):
        rex = re.match('(.{4}) (.{2,4})', text)
        if rex is None:
            raise ValueError(f'CUT_BOND entry "{text}" is not formatted correctly')
        data = dict(zip(('first', 'second'), rex.groups()))
        return cls(**data)


Entries.choices['CUT_BOND'] = (CUT_BONDEntry, Singletony.multiton)


#########################################################################################################

@dataclass
class CHARGEEntry:
    """
    No idea if anything respects this.
    """
    atom: str
    charge: int

    def __str__(self) -> str:
        return f'CHARGE {self.atom} FORMAL {self.charge}'

    @classmethod
    def from_str(cls, text: str):
        rex = re.match('(\w+) FORMAL ([-\d])', text)
        if rex is None:
            raise ValueError(f'CHARGE entry "{text}" is not formatted correctly')
        data = dict(zip(('atom', 'charge'), rex.groups()))
        return cls(**data)


Entries.choices['CHARGE'] = (CHARGEEntry, Singletony.multiton)


#########################################################################################################

class PDB_ROTAMERSEntry(GenericEntry):
    """
    This does zero checks for fine existance.
    """

    def __init__(self, body: str):
        super().__init__(header='PDB_ROTAMERS', body=body)


Entries.choices['PDB_ROTAMERS'] = (PDB_ROTAMERSEntry, Singletony.singleton)


#########################################################################################################

class ROTAMER_AAEntry(GenericEntry):
    def __init__(self, body: str):
        super().__init__(header='ROTAMER_AA', body=body)


Entries.choices['ROTAMER_AA'] = (ROTAMER_AAEntry, Singletony.singleton)


#########################################################################################################

class AAEntry(GenericEntry):
    def __init__(self, body: str = 'UNK'):
        if body != 'UNK':
            self.log.info('AA should be UNK... tolerating oddity.')
        super().__init__(header='AA', body=body)


Entries.choices['AA'] = (AAEntry, Singletony.singleton)


#########################################################################################################

class TYPEEntry(GenericEntry):
    """
    LIGAND or POLYMER. No exceptions.
    """

    def __init__(self, body: str = 'LIGAND'):
        assert body in ('POLYMER', 'LIGAND'), f'residue TYPE {body} is neither POLYMER or LIGAND'
        super().__init__(header='TYPE', body=body)


Entries.choices['TYPE'] = (TYPEEntry, Singletony.singleton)


#########################################################################################################

class ADD_RINGEntry(GenericEntry):
    ## To be fixed. Spacing drama...
    def __init__(self, body: str):
        self.log.info('ADD_RING is sloppily coded. the values are stored as an unsplit string!')
        super().__init__(header='ADD_RING', body=body)


Entries.choices['ADD_RING'] = (ADD_RINGEntry, Singletony.multiton)


#########################################################################################################

class PROPERTIESEntry(GenericListEntry):
    # https://graylab.jhu.edu/PyRosetta.documentation/pyrosetta.rosetta.core.chemical.html#pyrosetta.rosetta.core.chemical.ResidueProperty
    def __init__(self, *args: str):
        super().__init__('PROPERTIES', *args)


Entries.choices['PROPERTIES'] = (PROPERTIESEntry, Singletony.list_singleton)


#########################################################################################################

class VARIANTEntry(GenericListEntry):
    # do multiple variant entries get the same line??
    def __init__(self, *args: str):
        super().__init__('VARIANT', *args)


Entries.choices['VARIANT'] = (VARIANTEntry, Singletony.list_singleton)


#########################################################################################################

class FIRST_SIDECHAIN_ATOMEntry(GenericEntry):
    def __init__(self, body: str):
        super().__init__(header='FIRST_SIDECHAIN_ATOM', body=body)


Entries.choices['FIRST_SIDECHAIN_ATOM'] = (FIRST_SIDECHAIN_ATOMEntry, Singletony.singleton)


class BACKBONE_AAEntry(GenericEntry):
    def __init__(self, body: str):
        assert len(body) == 3, f'{body} is not 3 char long'
        super().__init__(header='BACKBONE_AA', body=body.upper())


Entries.choices['BACKBONE_AA'] = (BACKBONE_AAEntry, Singletony.singleton)


#########################################################################################################

class RAMA_PREPRO_FILENAMEEntry(GenericEntry):
    def __init__(self, body: str):
        super().__init__(header='RAMA_PREPRO_FILENAME', body=body)


Entries.choices['RAMA_PREPRO_FILENAME'] = (RAMA_PREPRO_FILENAMEEntry, Singletony.singleton)


#########################################################################################################

class METAL_BINDING_ATOMSEntry(GenericListEntry):
    def __init__(self, *args: str):
        super().__init__('METAL_BINDING_ATOMS', *args)


Entries.choices['METAL_BINDING_ATOMS'] = (METAL_BINDING_ATOMSEntry, Singletony.singleton)


#########################################################################################################

class ACT_COORD_ATOMSEntry(GenericListEntry):
    def __init__(self, *args: str):
        super().__init__('ACT_COORD_ATOMS', *args)


Entries.choices['ACT_COORD_ATOMS'] = (ACT_COORD_ATOMSEntry, Singletony.singleton)

#########################################################################################################
