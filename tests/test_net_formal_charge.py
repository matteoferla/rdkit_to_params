"""Tests for molecule formal-charge emission into Rosetta params."""

from __future__ import annotations

import pytest
from rdkit import Chem

from rdkit_to_params import Params


class TestNetFormalChargeEmission:
    """Check params metadata, serialization and the attached ResidueType."""

    def test_net_charge_describes_retained_fragment(
        self, disconnected_charge_case: tuple[Chem.Mol, int]
    ) -> None:
        mol, expected_charge = disconnected_charge_case
        params = Params.from_mol(mol, name="FRG")

        assert len(Chem.GetMolFrags(params.mol)) == 1
        assert Chem.GetFormalCharge(params.mol) == expected_charge
        assert params.NET_FORMAL_CHARGE[0].charge == expected_charge
        loaded = Params.loads(params.dumps())
        assert loaded.to_pose().residue(1).type().net_formal_charge() == expected_charge

    def test_neutral_from_mol_emits_zero_and_roundtrips(
        self, uncharged_explicit_ethanol: Chem.Mol
    ) -> None:
        params = Params.from_mol(uncharged_explicit_ethanol, name="ETH")

        assert len(params.NET_FORMAL_CHARGE) == 1
        assert params.NET_FORMAL_CHARGE[0].charge == 0
        loaded = Params.loads(params.dumps())
        assert len(loaded.NET_FORMAL_CHARGE) == 1
        assert loaded.NET_FORMAL_CHARGE[0].charge == 0
        assert loaded.to_pose().residue(1).type().net_formal_charge() == 0

    def test_charged_from_mol_emits_formal_charge_and_roundtrips(
        self, explicit_acetate: Chem.Mol
    ) -> None:
        assert Chem.GetFormalCharge(explicit_acetate) == -1
        params = Params.from_mol(explicit_acetate, name="ACT")

        assert len(params.NET_FORMAL_CHARGE) == 1
        assert params.NET_FORMAL_CHARGE[0].charge == -1
        loaded = Params.loads(params.dumps())
        assert len(loaded.NET_FORMAL_CHARGE) == 1
        assert loaded.NET_FORMAL_CHARGE[0].charge == -1
        assert loaded.to_pose().residue(1).type().net_formal_charge() == -1

    def test_cation_from_mol_emits_formal_charge_and_roundtrips(
        self, explicit_methylammonium: Chem.Mol
    ) -> None:
        assert Chem.GetFormalCharge(explicit_methylammonium) == 1
        params = Params.from_mol(explicit_methylammonium, name="MAM")

        assert len(params.NET_FORMAL_CHARGE) == 1
        assert params.NET_FORMAL_CHARGE[0].charge == 1
        loaded = Params.loads(params.dumps())
        assert len(loaded.NET_FORMAL_CHARGE) == 1
        assert loaded.NET_FORMAL_CHARGE[0].charge == 1
        assert loaded.to_pose().residue(1).type().net_formal_charge() == 1

    def test_from_mol_uses_formal_charges_not_partial_charge_sum(
        self, explicit_acetate: Chem.Mol
    ) -> None:
        for atom in explicit_acetate.GetAtoms():
            atom.SetDoubleProp("_GasteigerCharge", 0.0)

        params = Params.from_mol(explicit_acetate, name="ACT")

        assert sum(atom.partial for atom in params.ATOM) == pytest.approx(0.0)
        assert len(params.NET_FORMAL_CHARGE) == 1
        assert params.NET_FORMAL_CHARGE[0].charge == -1

    def test_from_smiles_emits_exactly_one_net_formal_charge(self) -> None:
        params = Params.from_smiles("CC(=O)[O-]", name="ACT")

        assert len(params.NET_FORMAL_CHARGE) == 1
        assert params.NET_FORMAL_CHARGE[0].charge == -1

    def test_neutral_connection_dummy_does_not_change_net_charge(self) -> None:
        params = Params.from_smiles("*CCO", name="DUM")

        assert len(params.NET_FORMAL_CHARGE) == 1
        assert params.NET_FORMAL_CHARGE[0].charge == 0

    @pytest.mark.parametrize("smiles", ["[*+]CCO", "[*-]CCO"])
    def test_charged_connection_dummy_is_rejected(self, smiles: str) -> None:
        with pytest.raises(ValueError, match="charged dummy atoms are unsupported"):
            Params.from_smiles(smiles, name="DUM")

    def test_charged_dummy_is_rejected_by_from_mol(self) -> None:
        mol = Chem.MolFromSmiles("[*+]CCO")
        if mol is None:
            raise ValueError("Could not build the charged-dummy test molecule")
        mol = Chem.AddHs(mol)

        with pytest.raises(ValueError, match="charged dummy atoms are unsupported"):
            Params.from_mol(mol, name="DUM")

    def test_manual_net_formal_charge_override_is_preserved(
        self, uncharged_explicit_ethanol: Chem.Mol
    ) -> None:
        params = Params.load_mol(uncharged_explicit_ethanol, name="ETH")
        params.NET_FORMAL_CHARGE.append(-3)
        params.convert_mol()

        assert len(params.NET_FORMAL_CHARGE) == 1
        assert params.NET_FORMAL_CHARGE[0].charge == -3
        assert params.to_pose().residue(1).type().net_formal_charge() == -3
