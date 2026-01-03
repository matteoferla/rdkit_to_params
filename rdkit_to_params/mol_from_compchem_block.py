"""
Convert PDB Chemical Component Dictionary (CCD) entries to RDKit molecules.

This module provides utilities to parse CCD mmCIF blocks (e.g., from the PDB's
Components.cif or individual ligand files) and convert them to RDKit Mol objects.

Requires optional dependencies: ``pip install rdkit_to_params[extra]``
"""

from warnings import warn

try:
    import gemmi
    import pandas as pd
except ImportError as err:
    warn(
        f"gemmi and pandas are required for chemcomp_to_mol. "
        f"Install with: pip install rdkit_to_params[extra] (ImportError: {err})",
        category=ImportWarning,
    )
    gemmi = None  # type: ignore[assignment]
    pd = None  # type: ignore[assignment]

from rdkit import Chem
from rdkit.Chem import AllChem

BONDTYPES = {
    "SING": Chem.BondType.SINGLE,
    "DOUB": Chem.BondType.DOUBLE,
    "TRIP": Chem.BondType.TRIPLE,
    "QUAD": Chem.BondType.QUADRUPLE,
    "POLY": Chem.BondType.UNSPECIFIED,
    "AROM": Chem.BondType.AROMATIC,
}


def correct_charges(mol: Chem.Mol):
    """
    Fix formal charges for atoms with valence exceptions.

    Handles:
    - Tetravalent nitrogen (N with 4 bonds -> +1 charge, e.g., quaternary amines)
    - Trivalent oxygen (O with 3 bonds -> +1 charge, e.g., oxonium ions)
    """
    mol.UpdatePropertyCache(strict=False)
    ps = Chem.DetectChemistryProblems(mol)
    if not ps:
        Chem.SanitizeMol(mol)
        return mol
    for p in ps:
        if p.GetType() == "AtomValenceException":
            at: Chem.Atom = mol.GetAtomWithIdx(p.GetAtomIdx())
            # Tetravalent nitrogen -> quaternary amine cation
            if (
                at.GetAtomicNum() == 7
                and at.GetFormalCharge() == 0
                and at.GetExplicitValence() == 4
            ):
                at.SetFormalCharge(1)
            # Trivalent oxygen -> oxonium cation
            elif (
                at.GetAtomicNum() == 8
                and at.GetFormalCharge() == 0
                and at.GetExplicitValence() == 3
            ):
                at.SetFormalCharge(1)
    Chem.SanitizeMol(mol)
    return mol


def chemcomp_to_mol(block: gemmi.cif.Block, comp_id: str | None = None) -> Chem.Mol:
    """
    Convert a chemcomp block into a rdkit Mol
    """
    mol = Chem.RWMol()
    if comp_id is None:
        comp_id = block.name
    mol.SetProp("_Name", comp_id)  # title
    atom_df = pd.DataFrame(block.get_mmcif_category("_chem_comp_atom"))
    if "comp_id" in atom_df.columns:  # there could be more than one chem comp in block
        atom_df = atom_df.loc[atom_df.comp_id == comp_id]
    bond_df = pd.DataFrame(block.get_mmcif_category("_chem_comp_bond"))
    if "comp_id" in bond_df.columns:  # there could be more than one chem comp in block
        bond_df = bond_df.loc[bond_df.comp_id == comp_id]
    name2idx = {}
    for _, row in atom_df.iterrows():
        symbol = row.type_symbol.capitalize()
        if symbol == "D":
            symbol = "H"  # deuterium as hydrogen
        atom = Chem.Atom(symbol)
        # https://blog.matteoferla.com/2020/03/atom-names-purely-in-rdkit.html
        info = Chem.AtomPDBResidueInfo(atomName=row.atom_id, residueName=row.comp_id)
        atom.SetMonomerInfo(info)
        atom.SetProp("molFileAlias", row.atom_id)
        atom.SetFormalCharge(int(row.charge))
        name2idx[row.atom_id] = mol.AddAtom(atom)
    for _, row in bond_df.iterrows():
        bondtype = BONDTYPES.get(
            str(row.value_order).upper()[:4],
            Chem.BondType.UNSPECIFIED,
        )
        # `delo`, `pi`, `poly` (polymeric bond) are not encodable
        mol.AddBond(name2idx[row.atom_id_1], name2idx[row.atom_id_2], bondtype)
    correct_charges(mol)
    return mol.GetMol()
