# rdkit_to_params package

## Submodules

## rdkit_to_params.constraint module

This contains `make_constraint` which is creates a constraint file.
It is completely independent and different in style because it was different.
It is not integral to the conversion, it’s just a utility.


### class rdkit_to_params.constraint.Constraints(smiles, names, ligand_res, target_res)
Bases: `object`


#### \__init__(smiles, names, ligand_res, target_res)
Give smiles of the two sides (with `\*`) and the names, residue numbers convert.
It requires 2 atoms on each side in addition to the attachment point.
Note that the atom names are stored in a non-standard way out of laziness. in the Atom prop ‘_AtomName’.
The instance has the following attributes


* `smiles`: stored input smiles string


* `names`: stored input list of names


* `ligand_res`: stored input liagnd residue


* `target_res`: stored input protein residue


* `cov_template`: Chem.Mol from first smiles.


* `target_template`: Chem.Mol from first smiles


* `combo`: combined templates


* `atom_pair_constraint`: AtomPair


* `angle_constraint`: Angle . NB. this is two lines.


* `dihedral_constaint`: dihedral

Class methods:
\* `assign_names(mol, list)` classmethods that assigns names to a mol in place.
\* `join_by_dummy(molA, molB)` classmethods that returns a joined molecule


* **Parameters**

    
    * **smiles** (`Tuple`[`str`, `str`]) – a tuple/list of two string. The first is the ligand, the second is the peptide.


    * **names** (`List`[`str`]) – a list of atom names. The ‘\*’ will need a name -but will be ignored-, but not the H.


    * **ligand_res** (`Union`[`str`, `int`]) – ligand residue in pose or PDB format (12 vs. 12A)


    * **target_res** (`Union`[`str`, `int`]) – peptide residue in pose or PDB format (12 vs. 12A)



#### classmethod assign_names(mol, names)
Stores names of atoms as given in the list.
totally non-standard way. PDBInfo is correct way. But too much effort.


* **Return type**

    `None`



#### get_atom(mol, name)

* **Return type**

    `Atom`



#### classmethod get_conn(mol)
Get connecting atom of mol.


* **Return type**

    `Atom`



#### classmethod join_by_dummy(a, b)

* **Return type**

    `Mol`


## rdkit_to_params.entries module

The main class here is Entries\`, which is a fancy list. It gets called for each uppercase attribute 
in the initialisation of `Params` (which happens in `_ParamsInitMixin` __e.g.__ `Entries.from_name('IO_STRING')`).


### class rdkit_to_params.entries.AAEntry(body='UNK')
Bases: `rdkit_to_params.entries.GenericEntry`


#### \__init__(body='UNK')
Initialize self.  See help(type(self)) for accurate signature.


### class rdkit_to_params.entries.ACT_COORD_ATOMSEntry(\*args)
Bases: `rdkit_to_params.entries.GenericListEntry`


#### \__init__(\*args)
Initialize self.  See help(type(self)) for accurate signature.


### class rdkit_to_params.entries.ADD_RINGEntry(body)
Bases: `rdkit_to_params.entries.GenericEntry`


#### \__init__(body)
Initialize self.  See help(type(self)) for accurate signature.


### class rdkit_to_params.entries.ATOMEntry(name: str, rtype: str, mtype: str = 'X', partial: float = 0)
Bases: `object`


#### \__init__(name: str, rtype: str, mtype: str = 'X', partial: float = 0)

* **Return type**

    `None`



#### classmethod from_str(text)

#### mtype(: str = 'X')

#### name(: str = None)

#### partial(: float = 0)

#### rtype(: str = None)

### class rdkit_to_params.entries.ATOM_ALIASEntry(body)
Bases: `rdkit_to_params.entries.GenericEntry`


#### \__init__(body)
Initialize self.  See help(type(self)) for accurate signature.


### class rdkit_to_params.entries.BONDEntry(first: str, second: str, order: int = 1)
Bases: `object`

dataclass class for both BOND and BOND_ENTRY. The `__str__` method will know based on `.order`.
The hash is the two atom names sorted. So BOND records with the same names will be equal.


#### \__init__(first: str, second: str, order: int = 1)

* **Return type**

    `None`



#### \__post_init__()

#### first(: str = None)

#### classmethod from_str(text)

#### order(: int = 1)

#### second(: str = None)

### class rdkit_to_params.entries.CHIEntry(index: int, first: str, second: str, third: str, fourth: str)
Bases: `object`


#### \__init__(index: int, first: str, second: str, third: str, fourth: str)

* **Return type**

    `None`



#### \__post_init__()

#### first(: str = None)

#### fourth(: str = None)

#### classmethod from_str(text)

#### index(: int = None)

#### second(: str = None)

#### third(: str = None)

### class rdkit_to_params.entries.CONNECTEntry(atom_name: str, index: int = 1, connect_type: str = '', connect_name: str = '')
Bases: `object`

This is a mess, but it guesses what you mean.
Deals with UPPER, LOWER and CONNECT.


#### \__init__(atom_name: str, index: int = 1, connect_type: str = '', connect_name: str = '')

* **Return type**

    `None`



#### \__post_init__()

#### atom_name(: str = None)

#### connect_name(: str = '')

#### connect_type(: str = '')

#### classmethod from_str(text)

#### index(: int = 1)

### class rdkit_to_params.entries.CUT_BONDEntry(first: str, second: str)
Bases: `object`

No idea what CUT_BOND is for.


#### \__init__(first: str, second: str)

* **Return type**

    `None`



#### \__post_init__()

#### first(: str = None)

#### classmethod from_str(text)

#### second(: str = None)

### class rdkit_to_params.entries.CommentEntry(body)
Bases: `rdkit_to_params.entries.GenericEntry`


#### \__init__(body)
Initialize self.  See help(type(self)) for accurate signature.


### class rdkit_to_params.entries.Entries(entry_cls, singleton=True)
Bases: `collections.abc.MutableSequence`

A fancy default list, where the elements are instances of whatver is in `entry_cls`.
It can be initialised via the class method `from_name` which accepst a string that has to be present in the class attribute `.choices`.
The `.append` method can work with str, list, dict or instance of the actual class it wants.
Note that the string for the string way must be without the header to the line.
The entry classes requires a `from_str` classmethod that returns an instance for this.
They also require __str__ method as this is how the entries are converted into string.

`Entries.from_name('BOND')`


#### \__init__(entry_cls, singleton=True)
The entries class is a fancy constrained list. The data is actually stored in `.data`.


* **Parameters**

    
    * **entry_cls** – what is the allowed class of the entries


    * **singleton** (`bool`) – is only one entry allowed?



#### choices( = {'#': (<class 'rdkit_to_params.entries.CommentEntry'>, False), 'AA': (<class 'rdkit_to_params.entries.AAEntry'>, True), 'ACT_COORD_ATOMS': (<class 'rdkit_to_params.entries.ACT_COORD_ATOMSEntry'>, True), 'ADD_RING': (<class 'rdkit_to_params.entries.ADD_RINGEntry'>, False), 'ATOM': (<class 'rdkit_to_params.entries.ATOMEntry'>, False), 'ATOM_ALIAS': (<class 'rdkit_to_params.entries.ATOM_ALIASEntry'>, False), 'BOND': (<class 'rdkit_to_params.entries.BONDEntry'>, False), 'CHI': (<class 'rdkit_to_params.entries.CHIEntry'>, False), 'CONNECT': (<class 'rdkit_to_params.entries.CONNECTEntry'>, False), 'CUT_BOND': (<class 'rdkit_to_params.entries.CUT_BONDEntry'>, False), 'FIRST_SIDECHAIN_ATOM': (<class 'rdkit_to_params.entries.FIRST_SIDECHAIN_ATOMEntry'>, True), 'ICOOR_INTERNAL': (<class 'rdkit_to_params.entries.ICOOR_INTERNALEntry'>, False), 'IO_STRING': (<class 'rdkit_to_params.entries.IO_STRINGEntry'>, True), 'METAL_BINDING_ATOMS': (<class 'rdkit_to_params.entries.METAL_BINDING_ATOMSEntry'>, True), 'NBR_ATOM': (<class 'rdkit_to_params.entries.NBR_ATOMEntry'>, True), 'NBR_RADIUS': (<class 'rdkit_to_params.entries.NBR_RADIUSEntry'>, True), 'PDB_ROTAMERS': (<class 'rdkit_to_params.entries.PDB_ROTAMERSEntry'>, True), 'PROPERTIES': (<class 'rdkit_to_params.entries.PROPERTIESEntry'>, False), 'RAMA_PREPRO_FILENAME': (<class 'rdkit_to_params.entries.RAMA_PREPRO_FILENAMEEntry'>, True), 'ROTAMER_AA': (<class 'rdkit_to_params.entries.ROTAMER_AAEntry'>, True), 'TYPE': (<class 'rdkit_to_params.entries.TYPEEntry'>, True), 'comment': (<class 'rdkit_to_params.entries.CommentEntry'>, False)})

#### classmethod from_name(name)

#### insert(index, value)
S.insert(index, value) – insert value before index


### class rdkit_to_params.entries.FIRST_SIDECHAIN_ATOMEntry(body)
Bases: `rdkit_to_params.entries.GenericEntry`


#### \__init__(body)
Initialize self.  See help(type(self)) for accurate signature.


### class rdkit_to_params.entries.GenericEntry(header, body)
Bases: `object`

This is meant to be inherited. `header` is the entry type. body is a string.


#### \__init__(header, body)
Initialize self.  See help(type(self)) for accurate signature.


#### classmethod from_str(text)

### class rdkit_to_params.entries.GenericListEntry(header, \*args)
Bases: `object`

This is meant to be inherited. `header` is the entry type. `values` is a list of strings.


#### \__init__(header, \*args)
Initialize self.  See help(type(self)) for accurate signature.


#### classmethod from_str(text)

### class rdkit_to_params.entries.ICOOR_INTERNALEntry(child: str, phi: float, theta: float, distance: float, parent: str, second_parent: str, third_parent: str)
Bases: `object`

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


#### \__init__(child: str, phi: float, theta: float, distance: float, parent: str, second_parent: str, third_parent: str)

* **Return type**

    `None`



#### \__post_init__()

#### child(: str = None)

#### distance(: float = None)

#### classmethod from_str(text)

#### parent(: str = None)

#### phi(: float = None)

#### second_parent(: str = None)

#### theta(: float = None)

#### third_parent(: str = None)

### class rdkit_to_params.entries.IO_STRINGEntry(name3: str = 'LIG', name1: str = 'Z')
Bases: `object`


* `.name3` is three letter name. `Params().NAME` is actually a dynamic attribute that uses this.


* `.name1` is a one letter name.

These get checked for length.


#### \__init__(name3: str = 'LIG', name1: str = 'Z')

* **Return type**

    `None`



#### \__post_init__()

#### classmethod from_str(text)

#### name1(: str = 'Z')

#### name3(: str = 'LIG')

### class rdkit_to_params.entries.METAL_BINDING_ATOMSEntry(\*args)
Bases: `rdkit_to_params.entries.GenericListEntry`


#### \__init__(\*args)
Initialize self.  See help(type(self)) for accurate signature.


### class rdkit_to_params.entries.NBR_ATOMEntry(body)
Bases: `rdkit_to_params.entries.GenericEntry`


#### \__init__(body)
Initialize self.  See help(type(self)) for accurate signature.


### class rdkit_to_params.entries.NBR_RADIUSEntry(body)
Bases: `rdkit_to_params.entries.GenericEntry`


#### \__init__(body)
Initialize self.  See help(type(self)) for accurate signature.


### class rdkit_to_params.entries.PDB_ROTAMERSEntry(body)
Bases: `rdkit_to_params.entries.GenericEntry`

This does zero checks for fine existance.


#### \__init__(body)
Initialize self.  See help(type(self)) for accurate signature.


### class rdkit_to_params.entries.PROPERTIESEntry(\*args)
Bases: `rdkit_to_params.entries.GenericListEntry`


#### \__init__(\*args)
Initialize self.  See help(type(self)) for accurate signature.


### class rdkit_to_params.entries.RAMA_PREPRO_FILENAMEEntry(body)
Bases: `rdkit_to_params.entries.GenericEntry`


#### \__init__(body)
Initialize self.  See help(type(self)) for accurate signature.


### class rdkit_to_params.entries.ROTAMER_AAEntry(body)
Bases: `rdkit_to_params.entries.GenericEntry`


#### \__init__(body)
Initialize self.  See help(type(self)) for accurate signature.


### class rdkit_to_params.entries.TYPEEntry(body='LIGAND')
Bases: `rdkit_to_params.entries.GenericEntry`

LIGAND or POLYMER. No exceptions.


#### \__init__(body='LIGAND')
Initialize self.  See help(type(self)) for accurate signature.

## Module contents

The main class here is `Params`. All underscore base classes are not meant to be used standalone.
`Entries` is the class for an list of entries of the same kind.


### class rdkit_to_params.Params()
Bases: `rdkit_to_params._io_mixin._ParamsIoMixin`, `rdkit_to_params._rdkit_convert._RDKitCovertMixin`, `rdkit_to_params._pyrosetta_mixin._PoserMixin`

`Params` creates and manipulates params files. It can handles several types of params operations,
such as “atom name surgery” and `rdkit.Chem.Mol` to a params file.

## Key methods


* `Params.load(filename)` will instantiate from file.


* `Params.from_mol(mol)` will instantiate from `Chem.Mol`


* `p.dump(filename)` will save a file.


* `loads` and 

```
``
```

dumps\`\`for strings.


* `p.fields` will return all header fields.


* `p.test` tests the params file in PyRosetta.


* `p.rename_atom(old, new)` changes an atom name

## Attributes

The attributes are generally the uppercase line headers, with a few exceptions.


* .comments is for # lines


* “BOND_TYPE” and “BOND” are merged into `.BOND`.


* “UPPER”, “LOWER” and “CONNECT” are merged into `.CONNECT`

With the exception of `.NAME` which depends on `.IO_STRING` basically
all the header type attributes are actually instances of the class Entries, which holds a sequence of specific entries.
see entry.py for the properties of each.
These can be a singleton, such as .IO_STRING which when a new line is added it gets overwritten instead, or not like say .ATOM.
That is to say that `.ATOM[0]` will give the first atom as expected, but this has to be done for `.IO_STRING[0]` too.

Atomnames…


* `p.get_correct_atomname` will return the 4 letter name of the atom with nice spacing as present in the entries


* `p.rename_atom` will change one atomname to a new one across all entries.


* `BOND`, `CHI`, `CUT_BOND` entries store 4 char atomnames as `.first`, `.second`, `.third`, `.fourth`.


* `ICOOR_INTERNAL` entries store 5 char atomnames as `.child`,\`\`.parent\`\`,\`\`.second_parent\`\`,\`\`.third_parent\`\`.


* `ATOM_ALIAS`, `NBR_ATOM`, `FIRST_SIDECHAIN_ATOM`, `ADD_RING` are just `entries.GenericEntries` instances, where `.body` is a string which will contain the atomname.


* `METAL_BINDING_ATOMS`, `ACT_COORD_ATOMS` are `entries.GenericListEntries` instances where `.values` is a list of string with maybe atomnames.

## Inheritance

It inherits several class,
which are not not mean to be used standalone, except for testing.

The pyrosetta and rdkit functionality are dependent on these being installed.


* `_ParamsIoMixin` adds read write, and inherits


* `_ParamsInitMixin` which adds the basics.


* `_PoserMixin` is a base that adds pyrosetta functionality if avaliable.


* `_RDKitCovertMixin`, which adds rdkit `from_mol` conversion functionality, the class is split in two, the other part being


* `_RDKitParamsPrepMixin`, which prepares the molecule for _RDKitCovertMixin.from_mol\`.


#### property NAME()

#### rename_atom(oldname, newname)
Change the atom name from `oldname` to `newname` and returns the 4 char `newname`.


* **Parameters**

    
    * **oldname** (`str`) – atom name, preferably 4 char long.


    * **newname** (`str`) – atom name, preferably 4 char long.



* **Return type**

    `str`



* **Returns**

    4 char newname



#### get_correct_atomname(name)
Given a name, gets the correctly spaced out one.
This has nothing to do with `._get_PDBInfo_atomname` which just returns the atom name from a `Chem.Atom`.


* **Parameters**

    **name** (`str`) – dirty name



* **Return type**

    `str`



* **Returns**

    correct name



#### validate()
