__all__ = ['neutralise', 'neutralize', 'DummyMasker']

from rdkit import Chem
from rdkit.Chem import AllChem
import warnings


def neutralize(mol: Chem) -> Chem:
    """
    Alters the protonation of tne molecule to be that at pH 7.
    Not great, but does the job.

    * Protonates amines, but not aromatic bound amines.
    * Deprotonates carboxylic acid, phosphoric acid and sulfuric acid, without ruining esters.
    """
    mol = Chem.RemoveHs(mol)
    protons_added = 0
    protons_removed = 0
    for indices in mol.GetSubstructMatches(Chem.MolFromSmarts('[N;D1]')):
        atom = mol.GetAtomWithIdx(indices[0])
        if atom.GetNeighbors()[0].GetIsAromatic():
            continue  # aniline
        atom.SetFormalCharge(1)
        atom.SetNumExplicitHs(3)
        protons_added += 1
    for indices in mol.GetSubstructMatches(Chem.MolFromSmarts('C(=O)[O;D1]')):
        atom = mol.GetAtomWithIdx(indices[2])
        # benzoic acid pKa is low.
        atom.SetFormalCharge(-1)
        atom.SetNumExplicitHs(0)
        protons_removed += 1
    for indices in mol.GetSubstructMatches(Chem.MolFromSmarts('P(=O)[Oh1]')):
        atom = mol.GetAtomWithIdx(indices[2])
        # benzoic acid pKa is low.
        atom.SetFormalCharge(-1)
        atom.SetNumExplicitHs(0)
        protons_removed += 1
    for indices in mol.GetSubstructMatches(Chem.MolFromSmarts('S(=O)(=O)[Oh1]')):
        atom = mol.GetAtomWithIdx(indices[3])
        atom.SetNumExplicitHs(0)
        # benzoic acid pKa is low.
        atom.SetFormalCharge(-1)
        protons_removed += 1
    AllChem.ComputeGasteigerCharges(mol)
    # return dict(protons_added=protons_added,
    #             protons_removed=protons_removed)
    return mol

def neutralise(*args, **kwargs):
    warnings.warn('The GB spelling has been changed to US', category=DeprecationWarning)
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
    """

    def __init__(self,
                 mol: Chem.Mol,
                 placekeeper_zahl:int=6,
                 blank_Gasteiger:bool=True):
        self.mol = mol
        self.is_masked = False
        self.zahl = int(placekeeper_zahl)
        self.blank_Gasteiger = bool(blank_Gasteiger)
        self.dummies = list(  mol.GetAtomsMatchingQuery(Chem.rdqueries.AtomNumEqualsQueryAtom(0))  )

    def mask(self):
        for dummy in self.dummies:
            dummy.SetAtomicNum(self.zahl)
            dummy.SetBoolProp('dummy', True)
            dummy.SetHybridization(Chem.HybridizationType.SP3)
        self.is_masked = True

    def unmask(self):
        for dummy in self.dummies:
            assert dummy.HasProp('dummy'), 'The atoms have changed somehow? (weird cornercase)'
            dummy.SetAtomicNum(0)
            if dummy.HasProp('_GasteigerCharge') and self.blank_Gasteiger:
                dummy.SetDoubleProp('_GasteigerCharge', 0.)
        self.is_masked = False

    def __enter__(self):
        self.mask()
        return self

    def __exit__(self, exc_type: Exception, exc_value: str, exc_traceback: 'bultins.traceback'):
        self.unmask()
