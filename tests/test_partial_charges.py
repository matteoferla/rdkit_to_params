"""Tests for partial-charge preparation and fallback handling."""

from __future__ import annotations

import math

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from rdkit_to_params import Params


class TestPartialChargeFallback:
    """Check generated charges and pre-existing charge properties."""

    def test_from_mol_preserves_computed_finite_charges(
        self, uncharged_explicit_ethanol: Chem.Mol
    ) -> None:
        baseline = Chem.Mol(uncharged_explicit_ethanol)
        AllChem.ComputeGasteigerCharges(baseline, throwOnParamFailure=True)
        expected = [atom.GetDoubleProp("_GasteigerCharge") for atom in baseline.GetAtoms()]

        params = Params.from_mol(uncharged_explicit_ethanol, name="ETH")
        actual = [atom.partial for atom in params.ATOM]

        assert actual == pytest.approx(expected, abs=1e-12)
        assert any(abs(charge) > 1e-6 for charge in actual)
        assert sum(actual) == pytest.approx(0.0, abs=1e-12)

    def test_from_mol_preserves_preassigned_charges(
        self, uncharged_explicit_ethanol: Chem.Mol
    ) -> None:
        expected = [0.1 - index * 0.01 for index in range(uncharged_explicit_ethanol.GetNumAtoms())]
        for atom, charge in zip(uncharged_explicit_ethanol.GetAtoms(), expected, strict=True):
            atom.SetDoubleProp("_GasteigerCharge", charge)

        params = Params.from_mol(uncharged_explicit_ethanol, name="ETH")

        assert [atom.partial for atom in params.ATOM] == pytest.approx(expected, abs=1e-12)

    def test_from_smiles_keeps_its_precomputed_gasteiger_charges(self) -> None:
        params = Params.from_smiles("CCO", name="ETH")
        charges = [atom.partial for atom in params.ATOM]

        assert any(abs(charge) > 1e-6 for charge in charges)
        assert all(math.isfinite(charge) for charge in charges)
        assert sum(charges) == pytest.approx(0.0, abs=1e-12)

    def test_from_mol_keeps_dummy_atoms_zero_and_preserves_real_atom_charges(self) -> None:
        mol = Chem.MolFromSmiles("*CCO")
        if mol is None:
            raise ValueError("Could not build the dummy-atom ethanol test molecule")
        mol = Chem.AddHs(mol)

        params = Params.from_mol(mol, name="DUM")
        dummy_index = next(
            atom.GetIdx() for atom in params.mol.GetAtoms() if atom.GetAtomicNum() == 0
        )
        charges = [atom.partial for atom in params.ATOM]

        assert params.mol.GetAtomWithIdx(dummy_index).GetDoubleProp(
            "_GasteigerCharge"
        ) == pytest.approx(0.0, abs=1e-12)
        assert any(abs(charge) > 1e-6 for charge in charges)

    @pytest.mark.parametrize("nonfinite_charge", [math.nan, math.inf, -math.inf])
    def test_from_mol_zeros_nonfinite_computed_charges(
        self,
        uncharged_explicit_ethanol: Chem.Mol,
        monkeypatch: pytest.MonkeyPatch,
        nonfinite_charge: float,
    ) -> None:
        def assign_nonfinite_charges(mol: Chem.Mol, throwOnParamFailure: bool = False) -> None:
            for atom in mol.GetAtoms():
                atom.SetDoubleProp("_GasteigerCharge", nonfinite_charge)

        monkeypatch.setattr(AllChem, "ComputeGasteigerCharges", assign_nonfinite_charges)

        params = Params.from_mol(uncharged_explicit_ethanol, name="ETH")

        assert [atom.partial for atom in params.ATOM] == [0.0] * params.mol.GetNumAtoms()
