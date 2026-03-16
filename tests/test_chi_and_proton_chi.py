"""Tests for CHI dihedral generation and PROTON_CHI support."""
import pytest

from rdkit_to_params import Params
from rdkit_to_params.entries import Entries, PROTON_CHIEntry

# ### PROTON_CHIEntry roundtrip -------------------------------------------------------------------------

class TestProtonChiEntry:
    """Test PROTON_CHIEntry parsing and serialisation."""

    def test_roundtrip(self):
        """from_str ↔ __str__ should be idempotent."""
        original = "3 SAMPLES 18 0 20 40 60 80 100 120 140 160 180 200 220 240 260 280 300 320 340 EXTRA 0"
        entry = PROTON_CHIEntry.from_str(original)
        assert entry.chi_index == 3
        assert len(entry.samples) == 18
        assert entry.extra == [0]
        # roundtrip
        reparsed = PROTON_CHIEntry.from_str(str(entry).replace("PROTON_CHI ", ""))
        assert reparsed.chi_index == entry.chi_index
        assert reparsed.samples == entry.samples
        assert reparsed.extra == entry.extra

    def test_amine_roundtrip(self):
        """Amine-style PROTON_CHI with EXTRA 1 20."""
        original = "2 SAMPLES 2 0 180 EXTRA 1 20"
        entry = PROTON_CHIEntry.from_str(original)
        assert entry.chi_index == 2
        assert entry.samples == [0, 180]
        assert entry.extra == [1, 20]

    def test_bad_format_raises(self):
        with pytest.raises(ValueError):
            PROTON_CHIEntry.from_str("garbage input")

    def test_entries_registry(self):
        """PROTON_CHI should be registered in Entries.choices."""
        assert "PROTON_CHI" in Entries.choices
        entries = Entries.from_name("PROTON_CHI")
        assert entries is not None

    def test_load_from_params_string(self):
        """PROTON_CHI lines in a params string should be parsed correctly."""
        params_text = """\
NAME SER
IO_STRING SER S
TYPE POLYMER
AA UNK
ATOM  N   Nbb  NH1  -0.4700000
ATOM  CA  CAbb CT1   0.0700000
ATOM  C   CObb C     0.5100000
ATOM  O   OCbb O    -0.5100000
ATOM  CB  CH2  CT2  -0.0100000
ATOM  OG  OH   OH1  -0.6600000
ATOM  HG  Hpol H     0.4300000
ATOM 1HB  Hapo HA    0.0900000
ATOM 2HB  Hapo HA    0.0900000
ATOM  H   HNbb H     0.3100000
ATOM  HA  Hapo HB    0.0900000
BOND  N    CA
BOND  N    H
BOND  CA   C
BOND  CA   CB
BOND  CA   HA
BOND  C    O
BOND  CB   OG
BOND  CB  1HB
BOND  CB  2HB
BOND  OG   HG
CHI 1  N    CA   CB   OG
CHI 2  CA   CB   OG   HG
PROTON_CHI 2 SAMPLES 18 0 20 40 60 80 100 120 140 160 180 200 220 240 260 280 300 320 340 EXTRA 0
CONNECT  N
CONNECT  C
NBR_ATOM CB
NBR_RADIUS 3.4473
ICOOR_INTERNAL    N      0.000000    0.000000    0.000000   N     CA    C
ICOOR_INTERNAL    CA     0.000000  180.000000    1.458001   N     CA    C
ICOOR_INTERNAL    C      0.000000   68.800003    1.523258  CA     N     C
ICOOR_INTERNAL    O   -180.000000   59.200001    1.231015   C    CA     N
ICOOR_INTERNAL    CB  -122.800000   69.625412    1.529508  CA     N     C
ICOOR_INTERNAL    OG     0.000000   66.393341    1.433430  CB    CA     N
ICOOR_INTERNAL    HG     0.000000   70.500000    0.960000  OG    CB    CA
ICOOR_INTERNAL   1HB   121.400000   70.500000    1.090000  CB    CA    OG
ICOOR_INTERNAL   2HB   117.600000   70.500000    1.090000  CB    CA   1HB
ICOOR_INTERNAL    H   -180.000000   60.850000    1.010000   N    CA     C
ICOOR_INTERNAL    HA  -119.000000   71.900000    1.090000  CA     N    CB
"""
        p = Params.loads(params_text)
        assert len(p.PROTON_CHI) == 1
        pchi = p.PROTON_CHI[0]
        assert pchi.chi_index == 2
        assert len(pchi.samples) == 18
        assert pchi.extra == [0]


# ### CHI generation from SMILES -----------------------------------------------------------------------

class TestChiGeneration:
    """Test that _parse_rotatables produces correct CHI entries."""

    def test_pentanol_chi_count(self):
        """Pentanol (CCCCCO) should have ≥3 heavy-atom CHIs."""
        p = Params.from_smiles("CCCCCO", name="PNL")
        heavy_chis = [c for c in p.CHI if c.fourth.strip()[0] != "H"]
        assert len(heavy_chis) >= 3, f"Expected ≥3 heavy CHIs, got {len(heavy_chis)}: {[str(c) for c in p.CHI]}"

    def test_pentanol_proton_chi(self):
        """Pentanol should have ≥1 PROTON_CHI with 18 samples (hydroxyl)."""
        p = Params.from_smiles("CCCCCO", name="PNL")
        assert len(p.PROTON_CHI) >= 1, f"Expected ≥1 PROTON_CHI, got {len(p.PROTON_CHI)}"
        has_18 = any(len(pc.samples) == 18 for pc in p.PROTON_CHI)
        assert has_18, "Expected a PROTON_CHI with 18 samples for hydroxyl"

    def test_hexanol_chi_overlap(self):
        """Hexanol (CCCCCCO) should have ≥4 CHIs with overlap between adjacent ones."""
        p = Params.from_smiles("CCCCCCO", name="HXL")
        heavy_chis = [c for c in p.CHI if c.fourth.strip()[0] != "H"]
        assert len(heavy_chis) >= 4, f"Expected ≥4 heavy CHIs, got {len(heavy_chis)}"
        # check overlap: atom 4 of CHI k == atom 2 of CHI k+1
        for i in range(len(heavy_chis) - 1):
            # adjacent CHIs should share atoms (not strictly 4→2, but overlap)
            atoms_k = {heavy_chis[i].first.strip(), heavy_chis[i].second.strip(),
                       heavy_chis[i].third.strip(), heavy_chis[i].fourth.strip()}
            atoms_k1 = {heavy_chis[i + 1].first.strip(), heavy_chis[i + 1].second.strip(),
                        heavy_chis[i + 1].third.strip(), heavy_chis[i + 1].fourth.strip()}
            overlap = atoms_k & atoms_k1
            assert len(overlap) >= 2, f"CHI {i + 1} and {i + 2} should share ≥2 atoms, got {overlap}"

    def test_ethanol_proton_chi(self):
        """Ethanol (CCO) should have ≥1 CHI ending in H and 1 PROTON_CHI."""
        p = Params.from_smiles("CCO", name="ETH")
        h_chis = [c for c in p.CHI if c.fourth.strip().startswith("H")]
        assert len(h_chis) >= 1, f"Expected ≥1 CHI ending in H, got {len(h_chis)}"
        assert len(p.PROTON_CHI) >= 1, f"Expected ≥1 PROTON_CHI, got {len(p.PROTON_CHI)}"

    def test_thiol_proton_chi(self):
        """Ethanethiol (CCS) should have ≥1 PROTON_CHI."""
        p = Params.from_smiles("CCS", name="THI")
        assert len(p.PROTON_CHI) >= 1, f"Expected ≥1 PROTON_CHI for thiol, got {len(p.PROTON_CHI)}"

    def test_branched_chi(self):
        """Branched molecule CC(C)CO should not crash and produce ≥1 CHI."""
        p = Params.from_smiles("CC(C)CO", name="BRA")
        assert len(p.CHI) >= 1, f"Expected ≥1 CHI for branched mol, got {len(p.CHI)}"

    def test_aa_serine_chi(self):
        """Serine-like (*C(=O)C(CO)[NH]*) should have ≥2 CHIs and ≥1 PROTON_CHI."""
        p = Params.from_smiles("*C(=O)C(CO)[NH]*", name="SRX")
        assert len(p.CHI) >= 2, f"Expected ≥2 CHIs for serine, got {len(p.CHI)}"
        assert len(p.PROTON_CHI) >= 1, f"Expected ≥1 PROTON_CHI for serine, got {len(p.PROTON_CHI)}"

    def test_no_crash_on_ring(self):
        """Cyclohexane-methanol should not crash."""
        p = Params.from_smiles("C1CCCCC1CO", name="CHM")
        # should have at least 1 CHI for the -CH2-OH tail
        assert len(p.CHI) >= 1

    def test_amine_proton_chi(self):
        """Ethylamine (CCN) should get PROTON_CHI with 2 samples."""
        p = Params.from_smiles("CCN", name="ETA")
        assert len(p.PROTON_CHI) >= 1, "Expected ≥1 PROTON_CHI for amine"
        # amine with 2H → samples should be [0, 180]
        amine_pchis = [pc for pc in p.PROTON_CHI if len(pc.samples) == 2]
        assert len(amine_pchis) >= 1, "Expected amine PROTON_CHI with 2 samples"

    def test_aib_no_chi(self):
        """AIB (gem-dimethyl) should produce 0 CHI and 0 PROTON_CHI.

        Both methyl groups are terminal (all non-parent neighbours are H),
        and the backbone N-CA-C bond must not leak through either.
        """
        p = Params.from_smiles("CC(C)(N*)C(*)=O", name="AIB")
        assert len(p.CHI) == 0, f"AIB should have 0 CHI, got {len(p.CHI)}: {[str(c) for c in p.CHI]}"
        assert len(p.PROTON_CHI) == 0, f"AIB should have 0 PROTON_CHI, got {len(p.PROTON_CHI)}"

    def test_valine_one_chi(self):
        """Valine-like should produce exactly 1 CHI (N-CA-CB-CG1), no PROTON_CHI.

        The isopropyl sidechain has one heavy-atom rotatable bond (CA-CB).
        The two terminal methyl groups (CG1, CG2) should be skipped.
        """
        p = Params.from_smiles("*[NH]C(C(*)=O)C(C)C", name="VAX")
        assert len(p.CHI) == 1, f"VAL-like should have 1 CHI, got {len(p.CHI)}: {[str(c) for c in p.CHI]}"
        assert len(p.PROTON_CHI) == 0, f"VAL-like should have 0 PROTON_CHI, got {len(p.PROTON_CHI)}"

    def test_no_backbone_chi_for_aminoacids(self):
        """Amino acid backbone bonds (N-CA, CA-C) must never appear as CHI."""
        p = Params.from_smiles("*C(=O)C(CO)[NH]*", name="SRX")
        backbone = {"N", "CA", "C", "O"}
        for chi in p.CHI:
            bc = {chi.second.strip(), chi.third.strip()}
            assert not bc.issubset(backbone), (
                f"CHI {chi} has both central atoms in backbone: {bc}"
            )


# ### dumps() includes PROTON_CHI -----------------------------------------------------------------------

class TestProtonChiInDumps:
    """Verify PROTON_CHI appears in serialised output."""

    def test_proton_chi_in_dumps(self):
        """Pentanol params string should contain PROTON_CHI line."""
        p = Params.from_smiles("CCCCCO", name="PNL")
        text = p.dumps()
        assert "PROTON_CHI" in text, "PROTON_CHI should appear in dumps() output"

    def test_proton_chi_roundtrip_via_dumps(self):
        """Generate params with PROTON_CHI, dumps, loads, check preserved."""
        p = Params.from_smiles("CCO", name="ETH")
        text = p.dumps()
        p2 = Params.loads(text)
        assert len(p2.PROTON_CHI) == len(p.PROTON_CHI)
        for a, b in zip(p.PROTON_CHI, p2.PROTON_CHI):
            assert a.chi_index == b.chi_index
            assert a.samples == b.samples
            assert a.extra == b.extra
