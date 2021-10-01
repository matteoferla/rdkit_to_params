__all__ = ['neutralise']

from rdkit import Chem
from rdkit.Chem import AllChem


def neutralise(mol: Chem) -> Chem:
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

class DummyMasker:
    """
    A context manager that allows operations on a mol containing dummy atoms (R/*) that
    otherwise would raise an RDKit error.
    It simply masks and unmasks the dummy atoms.

    >>> mol = Chem.MolFromSmiles('*CCC(C)C')
    >>> with DummyMasker(mol):
    >>>     AllChem.EmbedMolecule(mol)
    """

    def __init__(self, mol: Chem.Mol):
        self.mol = mol
        self.is_masked = False
        self.dummies = list(  mol.GetAtomsMatchingQuery(Chem.rdqueries.AtomNumEqualsQueryAtom(0))  )

    def mask(self):
        for dummy in self.dummies:
            dummy.SetAtomicNum(6)
        self.is_masked = True

    def unmask(self):
        for dummy in self.dummies:
            dummy.SetAtomicNum(0)
        self.is_masked = False

    def __enter__(self):
        self.mask()
        return self

    def __exit__(self, exc_type: Exception, exc_value: str, exc_traceback: 'bultins.traceback'):
        self.unmask()
