"""Tests for RAMA_PREPRO_FILENAME and RAMA_PREPRO_RESNAME entries."""
import pytest

from rdkit_to_params import Params
from rdkit_to_params.entries import (
    Entries,
    RAMA_PREPRO_FILENAMEEntry,
    RAMA_PREPRO_RESNAMEEntry,
)

# ### Entry classes -------------------------------------------------------------------------

class TestRamaPreProFilenameEntry:
    def test_roundtrip(self):
        entry = RAMA_PREPRO_FILENAMEEntry(
            general="scoring/score_functions/rama/alpha_ncaa/AIB_general.rama",
            prepro="scoring/score_functions/rama/alpha_ncaa/AIB_prepro.rama",
        )
        text = str(entry).replace("RAMA_PREPRO_FILENAME ", "")
        reparsed = RAMA_PREPRO_FILENAMEEntry.from_str(text)
        assert reparsed.general == entry.general
        assert reparsed.prepro == entry.prepro

    def test_from_str_two_paths(self):
        entry = RAMA_PREPRO_FILENAMEEntry.from_str("path/general.rama path/prepro.rama")
        assert entry.general == "path/general.rama"
        assert entry.prepro == "path/prepro.rama"

    def test_from_str_single_path_duplicates(self):
        """Single path is accepted and used for both general and prepro."""
        entry = RAMA_PREPRO_FILENAMEEntry.from_str("path/shared.rama")
        assert entry.general == "path/shared.rama"
        assert entry.prepro == "path/shared.rama"

    def test_from_str_empty_raises(self):
        with pytest.raises(ValueError):
            RAMA_PREPRO_FILENAMEEntry.from_str("")

    def test_registry(self):
        assert "RAMA_PREPRO_FILENAME" in Entries.choices
        entries = Entries.from_name("RAMA_PREPRO_FILENAME")
        assert entries is not None


class TestRamaPreProResnameEntry:
    def test_roundtrip(self):
        entry = RAMA_PREPRO_RESNAMEEntry("GENERIC_ALPHA_AMINOISOBUTYRIC_AA")
        text = str(entry)
        assert "RAMA_PREPRO_RESNAME GENERIC_ALPHA_AMINOISOBUTYRIC_AA" == text

    def test_registry(self):
        assert "RAMA_PREPRO_RESNAME" in Entries.choices


# ### Load from params string ---------------------------------------------------------------

class TestRamaPreProInParams:
    AIB_SNIPPET = """\
NAME AIB
IO_STRING AIB X
TYPE POLYMER
AA UNK
ATOM  N   Nbb  NH1  -0.4700000
ATOM  CA  CAbb CT1  -0.0700000
ATOM  C   CObb C     0.5100000
ATOM  O   OCbb O    -0.5100000
BOND  N    CA
BOND  CA   C
BOND  C    O
PROPERTIES PROTEIN ALPHA_AA ACHIRAL_BACKBONE
RAMA_PREPRO_FILENAME scoring/score_functions/rama/alpha_ncaa/AIB_general.rama scoring/score_functions/rama/alpha_ncaa/AIB_prepro.rama
RAMA_PREPRO_RESNAME GENERIC_ALPHA_AMINOISOBUTYRIC_AA
CONNECT  N
CONNECT  C
NBR_ATOM  N
NBR_RADIUS 6.99
ICOOR_INTERNAL    N      0.000000    0.000000    0.000000   N     CA    C
ICOOR_INTERNAL    CA     0.000000  180.000000    1.458001   N     CA    C
ICOOR_INTERNAL    C      0.000000   68.800003    1.523258  CA     N     C
ICOOR_INTERNAL    O   -180.000000   59.200001    1.231015   C    CA     N
"""

    def test_loads_filename(self):
        p = Params.loads(self.AIB_SNIPPET)
        assert len(p.RAMA_PREPRO_FILENAME) == 1
        entry = p.RAMA_PREPRO_FILENAME[0]
        assert "AIB_general" in entry.general
        assert "AIB_prepro" in entry.prepro

    def test_loads_resname(self):
        p = Params.loads(self.AIB_SNIPPET)
        assert len(p.RAMA_PREPRO_RESNAME) == 1
        assert "GENERIC_ALPHA_AMINOISOBUTYRIC_AA" in str(p.RAMA_PREPRO_RESNAME[0])

    def test_dumps_roundtrip(self):
        p = Params.loads(self.AIB_SNIPPET)
        text = p.dumps()
        assert "RAMA_PREPRO_FILENAME" in text
        assert "RAMA_PREPRO_RESNAME" in text
        assert "AIB_general" in text
        assert "GENERIC_ALPHA_AMINOISOBUTYRIC_AA" in text
        # re-parse
        p2 = Params.loads(text)
        assert p2.RAMA_PREPRO_FILENAME[0].general == p.RAMA_PREPRO_FILENAME[0].general
        assert str(p2.RAMA_PREPRO_RESNAME[0]) == str(p.RAMA_PREPRO_RESNAME[0])

    def test_convert_mol_preserves_rama(self):
        """Setting RAMA_PREPRO before convert_mol should survive the conversion."""
        p = Params.from_smiles("CC(C)(N*)C(*)=O", name="AIB")
        p.RAMA_PREPRO_FILENAME.append(
            "scoring/score_functions/rama/alpha_ncaa/AIB_general.rama "
            "scoring/score_functions/rama/alpha_ncaa/AIB_prepro.rama"
        )
        p.RAMA_PREPRO_RESNAME.append("GENERIC_ALPHA_AMINOISOBUTYRIC_AA")
        # re-convert should not blow away these entries
        p.convert_mol()
        assert len(p.RAMA_PREPRO_FILENAME) == 1
        assert len(p.RAMA_PREPRO_RESNAME) == 1
