"""Pytest configuration and fixtures for rdkit_to_params tests."""
from pathlib import Path

import _pytest.config
import pyrosetta
import pytest

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

