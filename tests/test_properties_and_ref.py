"""Tests for NUMERIC_PROPERTY, STRING_PROPERTY entries and ref_nc estimation."""
import pytest

from rdkit_to_params import Params
from rdkit_to_params.entries import (
    Entries,
    NUMERIC_PROPERTYEntry,
    STRING_PROPERTYEntry,
)


class TestNumericPropertyEntry:
    """NUMERIC_PROPERTYEntry parsing and roundtrip."""

    def test_roundtrip(self):
        entry = NUMERIC_PROPERTYEntry.from_str("REFERENCE -1.5")
        assert entry.tag == "REFERENCE"
        assert entry.value == -1.5
        assert "NUMERIC_PROPERTY REFERENCE -1.5" == str(entry)

    def test_positive_value(self):
        entry = NUMERIC_PROPERTYEntry(tag="WEIGHT", value=3.14)
        reparsed = NUMERIC_PROPERTYEntry.from_str(str(entry).replace("NUMERIC_PROPERTY ", ""))
        assert reparsed.tag == "WEIGHT"
        assert abs(reparsed.value - 3.14) < 0.01

    def test_registry(self):
        assert "NUMERIC_PROPERTY" in Entries.choices

    def test_bad_format(self):
        with pytest.raises(ValueError):
            NUMERIC_PROPERTYEntry.from_str("OOPS")


class TestStringPropertyEntry:
    """STRING_PROPERTYEntry parsing and roundtrip."""

    def test_roundtrip(self):
        entry = STRING_PROPERTYEntry.from_str("SOME_TAG some_value")
        assert entry.tag == "SOME_TAG"
        assert entry.value == "some_value"
        assert "STRING_PROPERTY SOME_TAG some_value" == str(entry)

    def test_registry(self):
        assert "STRING_PROPERTY" in Entries.choices


class TestRefEnergy:
    """Test estimate_ref_energy and auto_ref behaviour."""

    def test_estimate_nle(self):
        """NLE (norleucine) should get a finite reference energy."""
        p = Params.from_smiles("CCCCC(N*)C(*)=O", name="NLE")
        ref = p.estimate_ref_energy()
        assert isinstance(ref, float)
        assert -10 < ref < 10  # sanity bounds

    def test_auto_ref_ncaa(self):
        """auto_ref=True should add NUMERIC_PROPERTY REFERENCE for ncAAs."""
        p = Params.from_smiles("CCCCC(N*)C(*)=O", name="NLE")
        refs = [e for e in p.NUMERIC_PROPERTY if e.tag == "REFERENCE"]
        assert len(refs) == 1

    def test_auto_ref_ligand_excluded(self):
        """auto_ref should NOT add REFERENCE for plain ligands."""
        p = Params.from_smiles("CCCCCO", name="PNL")
        refs = [e for e in p.NUMERIC_PROPERTY if e.tag == "REFERENCE"]
        assert len(refs) == 0

    def test_auto_ref_covalent_excluded(self):
        """auto_ref should NOT add REFERENCE for covalent ligands."""
        p = Params.from_smiles("*C(=O)CCCC", name="COV")
        refs = [e for e in p.NUMERIC_PROPERTY if e.tag == "REFERENCE"]
        assert len(refs) == 0

    def test_add_ref_energy_explicit(self):
        """add_ref_energy with explicit value should use that value."""
        p = Params.from_smiles("CCCCC(N*)C(*)=O", name="NLE")
        p.NUMERIC_PROPERTY.clear()
        val = p.add_ref_energy(ref_value=-2.5)
        assert val == -2.5
        refs = [e for e in p.NUMERIC_PROPERTY if e.tag == "REFERENCE"]
        assert refs[0].value == -2.5

    def test_proline_like_correction(self):
        """Proline-like (aliphatic ring) ncAA should get a lower ref energy."""
        # proline: has aliphatic ring
        p_pro = Params.from_smiles("*C(=O)C1CCCN1*", name="PRO")
        # alanine: no ring
        p_ala = Params.from_smiles("*C(=O)C(C)N*", name="ALA")
        ref_pro = p_pro.estimate_ref_energy()
        ref_ala = p_ala.estimate_ref_energy()
        # proline correction is -2.8
        assert ref_pro < ref_ala

    def test_loads_roundtrip(self):
        """NUMERIC_PROPERTY REFERENCE should survive dumps/loads."""
        p = Params.from_smiles("CCCCC(N*)C(*)=O", name="NLE")
        text = p.dumps()
        p2 = Params.loads(text)
        refs_orig = [e for e in p.NUMERIC_PROPERTY if e.tag == "REFERENCE"]
        refs_loaded = [e for e in p2.NUMERIC_PROPERTY if e.tag == "REFERENCE"]
        assert len(refs_orig) == len(refs_loaded)
        assert abs(refs_orig[0].value - refs_loaded[0].value) < 0.001
