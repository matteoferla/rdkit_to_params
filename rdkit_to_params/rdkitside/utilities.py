__all__ = ["neutralise", "neutralize", "DummyMasker"]

import warnings
from types import TracebackType

from rdkit import Chem
from rdkit.Chem import AllChem


def neutralize(mol: Chem.Mol) -> Chem.Mol:
    """
    Alters the protonation of tne molecule to be that at pH 7.
    Not great, but does the job.

    * Protonates amines, but not aromatic bound amines.
    * Deprotonates carboxylic acid, phosphoric acid and sulfuric acid, without ruining esters.
    """
    mol = Chem.RemoveHs(mol)
    protons_added = 0
    protons_removed = 0
    for indices in mol.GetSubstructMatches(Chem.MolFromSmarts("[N;D1]")):
        atom = mol.GetAtomWithIdx(indices[0])
        if atom.GetNeighbors()[0].GetIsAromatic():
            continue  # aniline
        atom.SetFormalCharge(1)
        atom.SetNumExplicitHs(3)
        protons_added += 1
    for indices in mol.GetSubstructMatches(Chem.MolFromSmarts("C(=O)[O;D1]")):
        atom = mol.GetAtomWithIdx(indices[2])
        # benzoic acid pKa is low.
        atom.SetFormalCharge(-1)
        atom.SetNumExplicitHs(0)
        protons_removed += 1
    for indices in mol.GetSubstructMatches(Chem.MolFromSmarts("P(=O)[Oh1]")):
        atom = mol.GetAtomWithIdx(indices[2])
        # benzoic acid pKa is low.
        atom.SetFormalCharge(-1)
        atom.SetNumExplicitHs(0)
        protons_removed += 1
    for indices in mol.GetSubstructMatches(Chem.MolFromSmarts("S(=O)(=O)[Oh1]")):
        atom = mol.GetAtomWithIdx(indices[3])
        atom.SetNumExplicitHs(0)
        # benzoic acid pKa is low.
        atom.SetFormalCharge(-1)
        protons_removed += 1
    AllChem.ComputeGasteigerCharges(mol)  # type: ignore[attr-defined]
    # return dict(protons_added=protons_added,
    #             protons_removed=protons_removed)
    return mol


def neutralise(*args, **kwargs):
    warnings.warn("The GB spelling has been changed to US", category=DeprecationWarning)
    return neutralize(*args, **kwargs)


class DummyMasker:
    """
    A context manager that allows operations on a mol containing dummy atoms (R/*) that
    otherwise would raise an RDKit error.
    It simply masks and unmasks the dummy atoms.

    >>> mol = Chem.MolFromSmiles('*CCC(C)C')
    >>> with DummyMasker(mol):
    >>>     AllChem.EmbedMolecule(mol)

    The input options for dummy maker are ``mol`` (Chem.Mol),
    ``placekeeper_zahl`` (Z for atomic number),
    and ``blank_Gasteiger`` to make the dummy atom's '_GasteigerCharge' property zero if present.
    The Zahl of the placekeeping element will affect the Gasteiger partial chargers of nearby atoms though.

    A tricker case is when explicit hydrogens are added to the masked molecule,
    in which case the unmasking needs to be done with the static methods
     ``make_masked_with_extra_protons`` and ``unmask_without_proton``.
    """

    def __init__(self, mol: Chem.Mol, placekeeper_zahl: int = 6, blank_Gasteiger: bool = True):
        self.mol = mol
        self.is_masked = False
        self.zahl = int(placekeeper_zahl)
        self.blank_Gasteiger = bool(blank_Gasteiger)
        self.dummies = list(mol.GetAtomsMatchingQuery(Chem.rdqueries.AtomNumEqualsQueryAtom(0)))  # type: ignore[attr-defined]

    def mask(self):
        for dummy in self.dummies:
            dummy.SetAtomicNum(self.zahl)
            dummy.SetBoolProp("dummy", True)
            dummy.SetHybridization(Chem.HybridizationType.SP3)
        self.is_masked = True

    def unmask(self):
        for dummy in self.dummies:
            assert dummy.HasProp("dummy"), (
                "The atoms have changed because of new protons: use `make_masked` and `unmask_without_proton`"
            )
            dummy.SetAtomicNum(0)
            if dummy.HasProp("_GasteigerCharge") and self.blank_Gasteiger:
                dummy.SetDoubleProp("_GasteigerCharge", 0.0)
        self.is_masked = False

    def __enter__(self):
        self.mask()
        return self

    def __exit__(self, exc_type: Exception, exc_value: str, exc_traceback: TracebackType):
        self.unmask()

    @staticmethod
    def make_masked_with_extra_protons(mol: Chem.Mol, placekeeper_zahl: int = 8) -> Chem.Mol:
        # rdkit_to_params.DummyMasker cannot be used due to explicit H
        # mask dummy atoms
        dummies = list(mol.GetAtomsMatchingQuery(Chem.rdqueries.AtomNumEqualsQueryAtom(0)))  # type: ignore[attr-defined]
        for dummy in dummies:
            dummy.SetAtomicNum(placekeeper_zahl)
            dummy.SetBoolProp("dummy", True)
            dummy.SetHybridization(Chem.HybridizationType.SP3)
            dummy.UpdatePropertyCache()
        mol_with_hs: Chem.Mol = AllChem.AddHs(mol, addCoords=True)  # type: ignore[attr-defined]
        return mol_with_hs

    @staticmethod
    def copy_missed_props(source: Chem.Mol, target: Chem.Mol):
        """
        Assumes the index ordering is correct. And one differs by a proton at the end.
        """
        for ori_atom, mod_atom in zip(source.GetAtoms(), target.GetAtoms()):
            for prop, value in ori_atom.GetPropsAsDict().items():
                if prop.startswith("__"):
                    continue
                elif mod_atom.HasProp(prop):
                    continue
                elif isinstance(value, float):
                    mod_atom.SetDoubleProp(prop, value)
                elif isinstance(value, int):
                    mod_atom.SetIntProp(prop, value)
                elif isinstance(value, bool):
                    mod_atom.SetBoolProp(prop, value)
                else:
                    mod_atom.SetProp(prop, str(value))

    @classmethod
    def unmask_without_proton(cls, mol: Chem.Mol):
        """
        The mask has had protons added to it,
        this removes them and unmasks the dummy atoms,
        returns the molecule without them.
        """
        atom: Chem.Atom
        to_remove_bonds = []  # atom idxs of bond
        for idx, atom in enumerate(mol.GetAtoms()):
            if atom.HasProp("dummy") and atom.GetBoolProp("dummy"):
                dummy = atom
                dummy.SetAtomicNum(0)
                hydrogen_idxs = [
                    neigh.GetIdx() for neigh in dummy.GetNeighbors() if neigh.GetSymbol() == "H"
                ]
                for hydrogen_idx in hydrogen_idxs:
                    to_remove_bonds.append((dummy.GetIdx(), hydrogen_idx))

        rwmol = Chem.RWMol(mol)
        rwmol.BeginBatchEdit()
        for idx1, idx2 in to_remove_bonds:
            rwmol.RemoveBond(idx1, idx2)
        rwmol.CommitBatchEdit()
        split = AllChem.GetMolFrags(rwmol.GetMol(), asMols=True, sanitizeFrags=False)[0]  # type: ignore[attr-defined]
        cls.copy_missed_props(mol, split)
        return split
