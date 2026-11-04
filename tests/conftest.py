"""Pytest configuration and fixtures for rdkit_to_params tests."""
from pathlib import Path

import _pytest.config
import pyrosetta
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

PytestConfig = _pytest.config.Config


@pytest.fixture(scope="session", autouse=True)
def initialize_pyrosetta():
    """Initialize pyrosetta once for all tests in the session."""
    pyrosetta.init(extra_options='-mute all -load_PDB_components false')
    yield
    # Cleanup code can go here if needed


def path_to_file(pytestconfig: PytestConfig, *other: str) -> Path:
    """Get an absolute path to the given target file."""
    return pytestconfig.rootpath.joinpath(*other)


@pytest.fixture(scope="module")
def assets_dir(pytestconfig: PytestConfig) -> Path:
    """Get an absolute path to the test assets directory."""
    return path_to_file(pytestconfig, "tests", "assets")


@pytest.fixture(scope="module")
def official_phe_params(assets_dir: Path) -> Path:
    """Get an absolute path to the official PHE params file."""
    return assets_dir / "official_PHE.params"


@pytest.fixture
def uncharged_explicit_ethanol() -> Chem.Mol:
    """Build explicit-H ethanol without partial-charge properties."""
    mol = Chem.MolFromSmiles("CCO")
    if mol is None:
        raise ValueError("Could not build the ethanol test molecule")
    return Chem.AddHs(mol)


@pytest.fixture
def explicit_acetate() -> Chem.Mol:
    """Build explicit-H acetate with net formal charge -1."""
    mol = Chem.MolFromSmiles("CC(=O)[O-]")
    if mol is None:
        raise ValueError("Could not build the acetate test molecule")
    return Chem.AddHs(mol)


@pytest.fixture
def explicit_methylammonium() -> Chem.Mol:
    """Build explicit-H methylammonium with net formal charge +1."""
    mol = Chem.MolFromSmiles("C[NH3+]")
    if mol is None:
        raise ValueError("Could not build the methylammonium test molecule")
    return Chem.AddHs(mol)


@pytest.fixture(params=[("CC(=O)[O-].[Na+]", -1), ("CCO.[Na+]", 0)])
def disconnected_charge_case(request: pytest.FixtureRequest) -> tuple[Chem.Mol, int]:
    """Supply real disconnected inputs with an existing conformer and charges."""
    smiles, retained_charge = request.param
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Could not build the disconnected test molecule")
    mol = Chem.AddHs(mol)
    AllChem.Compute2DCoords(mol)
    for atom in mol.GetAtoms():
        atom.SetDoubleProp("_GasteigerCharge", 0.0)
    return mol, retained_charge
