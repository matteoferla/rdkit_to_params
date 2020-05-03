# RDKit to params
Create or modify Rosetta params files (topology files) from scratch, RDKit mols or another params file.

## Rationale
This is a fresh rewrite of ``mol_to_params.py``. For three reasons:

* I cannot share my 2to3 port and modded module-version of ``mol_to_params.py`` due to licence.
* I want to modify `params` files and more as opposed to use a standalone script.
* RDKit does not save ``mol2`` files, yet knows about atom names and Gasteiger-Massilli charges and more...

It sounds mad, but did not actually take too long.

## Caveat: I do not know many things!

### Chemical
I suspect I am doing stuff weirdly and I am meant to create ligands via ``pyrosetta.rosetta.core.chemical`` and not via params files... If this is so let me know. I don't mind knowing I made a mistake!

### Generic
I like this generic atom type business, but I am not sure how to use them in RL.
``from_mol(mol, generic=True)`` will make generic atom types.
I made several guesses with the classic atom types and I am sure many things are wrong...

### Rings and cis-trans
I don't know what `CUT_BOND` does. I think it is to do with rings.
`ADD_RING` is not implemented in the `from_mol` conversion as I think it's an old command.
Does a cis-trans tautomer bond (say `C(=O)-C=O`) gets a `CHI` entry? I am assuming no, but not sure.

## Roundtrip
Native amino acid params files can be found in the Rosetta folder
`rosetta/main/database/chemical/residue_type_sets/fa_standard/residue_types/l-caa`
Let's do a roundtrip changing an atomname:

    from rdkit_to_params import Params, pyrosetta
    pyrosetta.init(extra_options='-mute all')
    
    p = Params.load('PHE.params')
    p.IO_STRING[0].name3 = 'PHX'
    p.IO_STRING[0].name1 = 'Z'
    p.AA = 'UNK'
    del p.ROTAMER_AA[0]
    p.change_atomname(' CB ', ' CX ')
    p.dump('fake.params')
    p.test().dump_pdb('test.pdb')

## From mol object
### Requirements
For the sake of sanity, `EmbedMolecule`, `Chem.AddHs(mol)` or any weird hack is assumed to have been done beforehand.
And that the user is going to do `Chem.MolToPDBFile(mol)` or `Chem.MolToPDBBlock(mol)` _afterwards_.

The molecule should preferably be **not** Kekulised.
3letter name of residue is either from the title row (``_Name``) if a 3letter word or from the PDBInfo or 'LIG'.

Dummy atom (*/R) is assumed to be a CONNECT —ligand only atm.

Here is a conversion:

from rdkit_to_params import Params, pyrosetta
# note that pyrosetta needs to be started before rdkit.
pyrosetta.init(extra_options='-mute all')
from rdkit import Chem
from rdkit.Chem import AllChem

mol = Chem.MolFromSmiles('NCC(C)C(=O)O')
mol = AllChem.AddHs(mol)
display(mol)
AllChem.EmbedMolecule(mol)
AllChem.MMFFOptimizeMolecule(mol)
Params.add_names(mol, names=['N', 'CA', 'CB', 'C', 'O'], name='LIG')

p = Params.from_mol(mol)
p.NBR_ATOM.append('CB')
p.NBR_RADIUS.append('12.3')
p.test().dump_pdb('test.pdb')


## To Do
I have not coded yet, because I forgot:
* an auto-assignment of `NBR_ATOM` and `NBR_RADIUS` for `from_mol`.
* add rotamer line in `from_mol`
* change option to override starting atom.
* tweak the logic of `NAME` after some thinking.
* output constrain file for the CONNECT atom.