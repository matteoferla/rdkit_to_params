# RDKit to params
Create or modify Rosetta params files (topology files) from scratch, RDKit mols or another params file.

[![Read the Docs](https://img.shields.io/readthedocs/rdkit-to-params)](https://rdkit-to-params.readthedocs.io/en/latest/index.html)
[![https img shields io pypi pyversions rdkit to params logo python](https://img.shields.io/pypi/pyversions/rdkit--to--params?logo=python)](https://pypi.org/project/rdkit-to-params)
[![https img shields io pypi v rdkit to params logo python](https://img.shields.io/pypi/v/rdkit--to--params?logo=python)](https://pypi.org/project/rdkit-to-params)
[![https img shields io pypi dm rdkit to params logo python](https://img.shields.io/pypi/dm/rdkit_to_params?logo=python)](https://pypi.org/project/rdkit-to-params)
[![https img shields io github license matteoferla rdkit_to_params logo github](https://img.shields.io/github/license/matteoferla/rdkit_to_params?logo=github)](https://github.com/matteoferla/rdkit_to_params/raw/master/LICENCE)
[![https img shields io github last commit matteoferla rdkit_to_params logo github](https://img.shields.io/github/last-commit/matteoferla/rdkit_to_params?logo=github)](https://github.com/matteoferla/rdkit_to_params)
[![https img shields io github commit activity m matteoferla rdkit_to_params logo github](https://img.shields.io/github/commit-activity/m/matteoferla/rdkit_to_params?logo=github)](https://github.com/matteoferla/rdkit_to_params)

## Installation

To install from pip type:

    pip install rdkit-to-params

To install the latest version (probably the same) from GitHub

    pip install git+https://github.com/matteoferla/rdkit_to_params
    
RDKit and PyRosetta are optional module, but most of the useful functionality comes from the former.
To install rdkit, `conda install -c conda-forge rdkit` or `apt-get` or `pip install rdkit-pypi`.
To install PyRosetta you need to get a licence (free for academic use) and install it.
One option for installation without visiting the Rosetta Commons site is:

    pip install pyrosetta_help; PYROSETTA_USERNAME='👾👾👾' PYROSETTA_PASSWORD='👾👾👾' install_pyrosetta

## Rationale
This is a fresh rewrite of ``mol_to_params.py``. For three reasons:

* I cannot share my 2to3 port and modd  ed module-version of ``mol_to_params.py`` due to licence.
* I want to modify `params` files and more as opposed to use a standalone script.
* RDKit does not save ``mol2`` files, yet knows about atom names and Gasteiger-Massilli charges and more...

It sounds mad, but did not actually take too long.

## Website

For a web app using this see [https://direvo.mutanalyst.com/params](https://direvo.mutanalyst.com/params).
For the code running the website, see:

* [templates](https://github.com/matteoferla/DirEvo_tools/tree/master/direvo/templates/params)
* [views](https://github.com/matteoferla/DirEvo_tools/blob/master/direvo/views/params.py)


## Roundtrip

Native amino acid params files can be found in the Rosetta folder
`rosetta/main/database/chemical/residue_type_sets/fa_standard/residue_types/l-caa`
Let's do a roundtrip changing an atomname:

    import pyrosetta
    pyrosetta.init(extra_options='-mute all') # required for test
    from rdkit_to_params import Params
    
    p = Params.load('PHE.params')
    p.IO_STRING[0].name3 = 'PHX'
    p.IO_STRING[0].name1 = 'Z'
    p.AA = 'UNK'  #If it's not one of the twenty (plus extras), UNK!
    del p.ROTAMER_AA[0]
    p.rename_atom(' CB ', ' CX ') # this renames
    p.dump('fake.params')
    p.test().dump_pdb('test.pdb')
    
`p.test()` returns a pyrosetta pose. The static method `params_to_pose('something.params', name3)` accepts a params file
    
    import nglview
    pose = Params.params_to_pose('some_topology_I_found.params', name3)
    view = nglview.show_rosetta(pose)
    view

## From mol object
### Requirements
For the sake of sanity, `EmbedMolecule`, `Chem.AddHs(mol)` or any other operation is assumed to have been done beforehand.
And that the user is going to do `Chem.MolToPDBFile(params.mol)` or `Chem.MolToPDBBlock(params.mol)` or use the bound methods of `Params`,
`dump_pdb` and `dump_pdb_conf` (see below).

The molecule should preferably be **not** Kekulised.
3letter name of residue is either from the title row (``_Name``) if a 3letter word or from the PDBInfo or 'LIG'.

Dummy atom (\*/R) is assumed to be a CONNECT —ligand only atm.

Here is a conversion to an amino acid from a SMILES (quickest way):

    import pyrosetta
    pyrosetta.init(extra_options='-mute all')
    from rdkit_to_params import Params
    p = Params.from_smiles('*C(=O)C(Cc1ccccc1)[NH]*', #recognised as amino acid.
            name='PHX', #optional.
            atomnames={3: 'CZ'} #optional, rando atom name as see in previous edit
            )
    print(p.is_aminoacid()) # True
    p.dump('fake.params')
    p.test().dump_pdb('test.pdb')
    Chem.MolToPDBFile(mol, 'ref.pdb')
            
Here is a conversion to a ligand the circuitous way, just for fun:

    import pyrosetta
    pyrosetta.init(extra_options='-mute all')
    # note that pyrosetta needs to be started before rdkit.
    from rdkit_to_params import Params
    # make the molecule in RDKit or chemdraw or download it or whatever.
    mol = Chem.MolFromSmiles('NC(C(=O)O)Cc1ccccc1')
    mol = AllChem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    AllChem.MMFFOptimizeMolecule(mol)
    # add names to the mol beforehand
    Params.add_names(mol, names=['N', 'CA', 'C', 'O', 'OXT', 'CB'], name='PHZ')
    # parameterise
    p = Params.from_mol(mol, name='PHZ')
    p.test().dump_pdb('test.pdb')
    Chem.MolToPDBFile('ref.pdb')
 
The class method `add_names` is based upon atom index
(which is derived from the SMILES or sdf/mol file unless atoms have been replaced).
The instance method `rename_by_substructure` accepts a substructure and a list of atom names in the order they are in the substructure.
    
Note that conformer generation is not fully automatic and is not done by default.

    # make your conformers as you desire
    AllChem.EmbedMultipleConfs(mol, numConfs=10) # or whatever you choose. This is a somewhat important decision.
    AllChem.AlignMolConformers(mol) # I do not know if the conformers need to be aligned for Rosetta
    # params time!
    p = Params.from_mol(mol, name='LIG') 
    p.dump_pdb_conf('LIG_conf.pdb')
    p.PDB_ROTAMERS.append('LIG_conf.pdb')
    p.dump('my_params.params')
    
Note `dump_pdb` and `dump_pdb_conf` will save the molecule(s) without the dummy atoms, to stop this add `stripped=False`.

## From SMILES string
The above is actually a bit convoluted for example purposes as `Params.from_smiles`, accepts a SMILES string.

## From SMILES string and PDB for names    
In some cases one has a PDB file with a ligand that needs parameterising.
Assuming one has also the smiles of the ligand (PubChem has an super easy search), one can do

```python
p = Params.from_smiles_w_pdbfile(pdb_file, smiles, 'XXX') # the name has to match.
```
    
The smiles does not need to match full. It can contain more atoms or even one`*` (CONNECT).
The smiles gets parameterised. So be suse to add correct charges properly —hydrogens are added.
It could be used for scaffold hopping, but if position matters so much,
you may be interested in [Fragmenstein](https://github.com/matteoferla/Fragmenstein).

For more see [autogenerated documentation](sphinx-documentation.md). Sphinx with markdown cannot deal with typehinting,
so checking the code might be clearer.

## Rename

A key part is the atom names ——this can happen at . The following renaming methods are present:

* `p.rename(???)`: "overloaded" method that directs to the others
* `p.rename_from_str('XX,YY,ZZ')` or `p.rename_from_str('0:XX,3:YY')`
* `p.rename_from_list(['XX','YY', 'ZZ'])`
* `p.rename_from_dict({0:'XX',3:'YY'})`
* `p.rename_from_template(Chem.Mol)`
* `p.rename_by_substructure(Chem.Mol, ['XX','YY', 'ZZ'])` where the list is the atom idx in substructure

Note, ``retype_by_name`` does not have all these options (only atomname -> Rosetta atomtype).

The class method ``add_names`` simply uses these, but returns a mol

### DIY

If you have two mol objects from whatever routes, the basic operation is:

```python
p = Params.load_mol(mol, generic=False, name='LIG')
p.rename_from_template(template) # or whatever middle step

p.convert_mol()
```


Note that `convert_mol` should be called once and is already called in the two `from_XXX` classmethods.

```python
p = Params.from_mol(...)
p.convert_mol() # No!!!
p.mol # is the mol...
p2 = Params.load_mol(p.mol)
p2.convert_mol() # Yes
```

## Partial charges

The partial charges in the `rdkit.Chem.Atom property` to use can be set via `pcharge_prop_name` argument in `from_mol` and `load_mol`
among others.
By default it is `'_GasteigerCharge'`, which is the Gasteiger-Marsili charge assigned by RDKit.
Custom partial charges need to be assigned to the atoms beforehand, for example:

```python
import psi4
from rdkit import Chem
from rdkit_to_params import Params

mol: Chem.Mol = ...  # your molecule
pcharge_prop_name: str = 'ESPCharge'  # custom property name for partial charges
mol_string = f"{Chem.GetFormalCharge(mol)} 1\n" + "\n".join(Chem.MolToXYZBlock(mol).split('\n')[2:])  # charge, multiplicity, then coords
psi4_mol = psi4.geometry(mol_string)
psi4.set_options({
    'basis': '6-31G*',
    'scf_type': 'df'
})
energy, wfn = psi4.energy('B3LYP', return_wfn=True)
psi4.set_options({'CUBEPROP_TASKS': ['ESP']})
psi4.cubeprop(wfn)
charges: list[float] = wfn.atomic_point_charges().np.tolist()
for i, atom in enumerate(mol.GetAtoms()):
    atom.SetDoubleProp(pcharge_prop_name, charges[i])
    
p = Params.from_mol(mol, pcharge_prop_name=pcharge_prop_name, name='LIG')
```

## Constraints

The selfstanding class `Constraints` is for generating constraint files, which are a must with covalent attachments
in order to stop janky topologies.
The class is instantiated with a pair of SMILES, each with at least a real atom and with one attachment point,
the first is the ligand and the second is its peptide target. The names of the heavy atoms and the Rosetta residue "numbers".

```python
from rdkit_to_params import Constraints
c = Constraints(smiles=('*C(=N)', '*SC'), names= ['*', 'CX', 'NY', '*', 'SG', 'CB'], ligand_res= '1B', target_res='145A')
c.dump('con.con')
# individual strings can be accessed
c.atom_pair_constraint
c.angle_constraint
c.dihedral_constaint
c.custom_constaint # if you want to add your own before `str`, `.dumps`, `.dump`.
```

Do note that to make covalent links work in Rosetta, NGL and a few other places you need a LINK record, here is a f-string
for it:

    f'LINK        {target_atom: >4} {target_resn: >3} {p_chain[:1]} {target_resi: >3}                '+\
    f'{ligand_atom: >4} {ligand_resn: >3} {ligand_chain[:1]} {ligand_resi: >3}     1555   1555  1.8\n'

This is not to be confused with CCP4 REFMAC's `LINKR`, which are however easy to covert.
Alternatively, you can add it after importing the pose, _cf._ `pose.residue(lig_pos).connect_map`.

## Bond order
It is worth mentioning that the bond order specified in the topology file in the `BOND_ORDER` lines is mostly ignored 
and the bond order is derived from the rosetta types that get assigned. 
To extract and correct a ligand, consider the following

    
```python
# pose to string
buffer = pyrosetta.rosetta.std.stringbuf()
pose.dump_pdb(pyrosetta.rosetta.std.ostream(buffer))
pdbblock = buffer.str()
# get the residue
mol = Chem.MolFromPDBBlock(pdbblock, proximityBonding=False, removeHs=False)
ligand = Chem.SplitMolByPDBResidues(mol, whiteList=[params.NAME])[params.NAME]
# fix bond order
template = AllChem.DeleteSubstructs(params.mol, Chem.MolFromSmiles('*'))
AllChem.AssignBondOrdersFromTemplate(template, ligand)
```

## Amino acids

A `*C(=O)C([*:3])[NH]*` molecule, where R3 is whatever sidechain is automatically converted into an amino acid.
Omitting the hydrogen on the amine is fine (implicit), so `*C(=O)C([*:3])N*` is also automatically accepted.
Likewise, a secondary amine like in proline, `*C(=O)C1CCCN1*`, is automatically determined to be an amino acid.
Omitting the double bond of the carboxyl will result in a hydroxyl backbones amino acid, which will behave like `C=[OH+]`
for properties, but without the partial charge.
The criterion for an amino acid is if the substracture `*NCC(~O)*` is matched (see `_aminoacid_override`).

However, `C(=O)C([*:3])[NH]` will be parsed as `[CH](=O)C([*:3])[NH]`, i.e. with a radical amine and an aldehyde.

Here is an example of making a sequence with a custom residue (without writing to file):
```python
import nglview as nv
from rdkit_to_params import Params
import pyrosetta

# make params
p = Params.from_smiles('CCCCC(N*)C(*)=O', name='NLE')
p.PROPERTIES.append('ALIPHATIC')
p.PROPERTIES.append('HYDROPHOBIC')
# p.test() would test it in isolation.

# add to pose
pose = pyrosetta.Pose()
rst = p.add_residuetype(pose)
pyrosetta.rosetta.core.pose.make_pose_from_sequence(pose, 'AX[NLE]A', rst)

# relax and show
scorefxn = pyrosetta.get_fa_scorefxn()
relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, 15)
relax.apply(pose)
nv.show_rosetta(pose)
```

### Greek
In the amino acid case, the class attribute `greekification` changes the atomnames to CB, CD2 etc.
It is by default `True`. It is called during `fix_mol`, a step in `load_mol`/`load_smiles`,
so should be safe for rename methods.


## Optionals
Installing RDKit with conda is easy (`conda install rdkit`).
With apt-get likewise (`sudo apt-get install python3-rdkit  librdkit1 rdkit-data`).
With brew idem (`brew install rdkit --with-python3`).
But there is always a system where one needs to compile it from source, which is a pain. Hence why it is optional.
For example, I have never tried installing it on a Windows.

Pyrosetta is optional because it has a non-standard installation.

## Terminal caps

To make a cap, there is a quick way:

```python
p = Params.from_smiles('*NCC', name='CAP', atomnames={1: ' N  ', 2: ' CA '})
p.make_C_terminal_cap(mainchain_atoms=['N', 'CA'])

p = Params.from_smiles('*C(=O)C', name='CAP', atomnames={1: ' C  ', 2: ' O  ', 3: ' CA '})
p.make_N_terminal_cap(mainchain_atoms=['C', 'CA'])
```

These methods also accept `connection_idx`, which is the Fortran-style index of the connection that will become a LOWER/UPPER.
i.e. if the cap is further connected but not as a polymer, say `*NCC*`.

Do note:

* `.test()` does not work on a terminal cap and will segfault
* `mainchain_atoms` will change Rosetta atom-types only if the name matches
* The code accepts only cases with 3 or more atoms (so a `*=O` cap is a no go and requires virtual atoms added manually)
    
## Caveat: I do not know many things!

### Chemical
I suspect I am doing stuff weirdly and I am meant to create ligands via ``pyrosetta.rosetta.core.chemical`` and not via params files... 
If this is so let me know. I don't mind knowing I made a mistake!

Edit: it is indeed possible, but it is far from trivial/sane. The class
`pyrosetta.rosetta.core.chemical.MutableResidueType` can be operated upon by various
other classes such as `pyrosetta.rosetta.core.chemical.SetAtomicCharge`.
A mutable residue type can be converted to a regular residue type via `.make()`.

### Generic
I like this generic atom type business, but I am not sure how to use them in RL.
``from_mol(mol, generic=True)`` will make generic atom types.
I made several guesses with the classic atom types and I am sure many things are wrong...

### Rings and cis-trans
* `CUT_BOND` is implemented, but I am not sure it does anything. `CHI` entries cannot work across a cut bond,
even when undeclared, so is likely redudant.
* `ADD_RING` is not implemented in the `from_mol` conversion as I think it's an old command.
* Does a cis-trans tautomer bond (say `C(=O)-C=O`) gets a `CHI` entry? I am assuming no, but not sure.

### Notes

There are some other things to pay attention to:

* Atom names are 4-letters. It is always safer to add the spaces yourself if assigning them.
* CHI struggles with rings, so currently `C1CCCCC1CCC` has only one CHI (C7, C8, C9, H10), even if (C6, C7, C8, C9) most probably counts.

## To Do

The `from_mol` class method recognises `*[NH]CC(~O)*` and assigns it as a backbone properly.
However, `Chem.MolFromSmiles('*[NH]CC(~O)*')` cannot be embedded, so is a bit of a horrible one for users to use.
Throughout the code, dummy atoms (\*/R) are changed to carbons or chlorines and then changed back.
Cystathionine and similar twinned amino acids are the problem as I cannot simply make an amino acid backbone be recognised,
however if protonated as is the case `[NH1]C[CH0](=O)`.
Maybe the `CC(=O)NCC(=O)NC` option may be a better choice after all.

## Footnote

To save a ResidueType in PyRosetta to a params file, the command is:

```python
pyrosetta.rosetta.core.chemical.write_topology_file(residuetype, 'foo.params')  # noqa
```

## Legal Disclaimer
The author, Matteo Ferla, is not affiliated with either Rosetta or RDKit.

[![Matteo Ferla orcid](https://img.shields.io/badge/orcid-0000--0002--5508--4673-a6ce39?logo=orcid)](https://orcid.org/0000--0002--5508--4673)
[![Matteo Ferla googlescholar](https://img.shields.io/badge/google--scholar-gF--bp_cAAAAJ-success?logo=googlescholar)](https://scholar.google.com/citations?user=gF--bp_cAAAAJ&hl=en)
[![Matteo Ferla twitter](https://img.shields.io/twitter/follow/matteoferla?label=Follow&logo=twitter)](https://twitter.com/matteoferla) 
[![Matteo Ferla stackoverflow](https://img.shields.io/stackexchange/stackoverflow/r/4625475?logo=stackoverflow)](https://stackoverflow.com/users/4625475)
[![Matteo Ferla stackexchange](https://img.shields.io/stackexchange/bioinformatics/r/6322?logo=stackexchange)](https://bioinformatics.stackexchange.com/users/6322) 
[![Matteo Ferla googlemail](https://img.shields.io/badge/email-gmail-informational&logo=googlemail)](https://mailhide.io/e/Ey3RNO2G) 
[![Matteo Ferla Oxford](https://img.shields.io/badge/email-Oxford-informational&logo=googlemail)](https://mailhide.io/e/Y1dbgyyE)
