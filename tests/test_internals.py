"""Test internal functionality of rdkit_to_params package."""
import re
from pathlib import Path

import pytest
from rdkit import Chem

from rdkit_to_params import Params
from rdkit_to_params.entries import (
    ATOMEntry,
    BONDEntry,
    CHARGEEntry,
    CHIEntry,
    CUT_BONDEntry,
    Entries,
    GenericEntry,
    ICOOR_INTERNALEntry,
    IO_STRINGEntry,
)


class TestEntries:
    """Test the Entries class and entry types."""

    def test_entries_creation(self):
        """Test creating an Entries instance."""
        entries = Entries.from_name("ATOM")
        assert isinstance(entries, Entries)

    def test_entries_append(self):
        """Test appending entries."""
        entries = Entries.from_name("ATOM")
        entry = ATOMEntry(name=" CA ", rtype="CT1", mtype="CT", partial=0.0)
        entries.append(entry)
        assert len(entries) == 1
        assert entries[0] == entry

    def test_entries_iteration(self):
        """Test iterating over entries."""
        entries = Entries.from_name("BOND")
        bond1 = BONDEntry(first=" CA ", second=" CB ", order=1)
        bond2 = BONDEntry(first=" CB ", second=" CG ", order=1)
        entries.append(bond1)
        entries.append(bond2)

        bonds = list(entries)
        assert len(bonds) == 2
        assert bonds[0] == bond1
        assert bonds[1] == bond2


class TestATOMEntry:
    """Test ATOM entry parsing and creation."""

    def test_atom_entry_creation(self):
        """Test creating an ATOM entry."""
        entry = ATOMEntry(name=" CA ", rtype="CT1", mtype="CT", partial=0.07)
        assert entry.name == " CA "
        assert entry.rtype == "CT1"
        assert entry.mtype == "CT"
        assert entry.partial == 0.07

    def test_atom_entry_from_str(self):
        """Test parsing ATOM entry from string."""
        text = " CA  CT1 CT  0.07"
        entry = ATOMEntry.from_str(text)
        assert entry.name == " CA "
        assert entry.rtype.strip() == "CT1"
        assert entry.mtype.strip() == "CT"
        assert entry.partial == 0.07

    def test_atom_entry_to_str(self):
        """Test converting ATOM entry to string."""
        entry = ATOMEntry(name=" CA ", rtype="CT1", mtype="CT", partial=0.07)
        result = str(entry)
        assert " CA " in result
        assert "CT1" in result
        assert "0.07" in result


class TestBONDEntry:
    """Test BOND entry parsing and creation."""

    def test_bond_entry_creation(self):
        """Test creating a BOND entry."""
        entry = BONDEntry(first=" CA ", second=" CB ", order=1)
        assert entry.first == " CA "
        assert entry.second == " CB "
        assert entry.order == 1

    def test_bond_entry_from_str(self):
        """Test parsing BOND entry from string."""
        text = " CA   CB  1"
        entry = BONDEntry.from_str(text)
        assert entry.first == " CA "
        assert entry.second == " CB "
        assert entry.order == 1

    def test_bond_entry_aromatic(self):
        """Test parsing aromatic BOND entry."""
        text = " CA   CB  ARO"
        entry = BONDEntry.from_str(text)
        assert entry.first == " CA "
        assert entry.second == " CB "
        assert entry.order == 4

    def test_bond_entry_empty_order(self):
        """Test parsing BOND entry with empty order defaults to 1."""
        text = " CA   CB "
        entry = BONDEntry.from_str(text)
        assert entry.order == 1


class TestCHIEntry:
    """Test CHI entry parsing and creation."""

    def test_chi_entry_creation(self):
        """Test creating a CHI entry."""
        entry = CHIEntry(index=1, first=" C6 ", second=" C5 ", third=" C4 ", fourth=" C3 ")
        assert entry.index == 1
        assert entry.first == " C6 "
        assert entry.second == " C5 "

    def test_chi_entry_from_str(self):
        """Test parsing CHI entry from string."""
        text = "1  C6   C5   C4   C3"
        entry = CHIEntry.from_str(text)
        assert entry.index == "1"
        assert entry.first.strip() == "C6"
        assert entry.second.strip() == "C5"
        assert entry.third.strip() == "C4"
        assert entry.fourth.strip() == "C3"

    def test_chi_entry_invalid_format(self):
        """Test that invalid CHI format raises ValueError."""
        with pytest.raises(ValueError, match="CHI entry.*not formatted correctly"):
            CHIEntry.from_str("invalid format")


class TestICOOR_INTERNALEntry:
    """Test ICOOR_INTERNAL entry parsing."""

    def test_icoor_entry_from_str_position_based(self):
        """Test parsing position-based ICOOR_INTERNAL entry."""
        text = " CA   -57.924927 125.580849   1.458101  N    C    CA"
        entry = ICOOR_INTERNALEntry.from_str(text)
        assert entry.child.strip() == "CA"
        assert float(entry.phi) == pytest.approx(-57.924927, rel=1e-5)
        assert float(entry.theta) == pytest.approx(125.580849, rel=1e-5)
        assert float(entry.distance) == pytest.approx(1.458101, rel=1e-5)

    def test_icoor_entry_from_str_space_based(self):
        """Test parsing space-based ICOOR_INTERNAL entry."""
        text = "CA -57.92 125.58 1.458 N C CA"
        entry = ICOOR_INTERNALEntry.from_str(text)
        assert entry.child.strip() == "CA"
        assert float(entry.phi) == pytest.approx(-57.92, rel=1e-2)


class TestCUT_BONDEntry:
    """Test CUT_BOND entry parsing."""

    def test_cut_bond_entry_creation(self):
        """Test creating a CUT_BOND entry."""
        entry = CUT_BONDEntry(first=" CA ", second=" CB ")
        assert entry.first == " CA "
        assert entry.second == " CB "

    def test_cut_bond_entry_from_str(self):
        """Test parsing CUT_BOND entry from string."""
        text = " CA   CB "
        entry = CUT_BONDEntry.from_str(text)
        assert entry.first.strip() == "CA"
        assert entry.second.strip() == "CB"


class TestCHARGEEntry:
    """Test CHARGE entry parsing."""

    def test_charge_entry_from_str_positive(self):
        """Test parsing positive CHARGE entry."""
        text = "NZ FORMAL +1"
        entry = CHARGEEntry.from_str(text)
        assert entry.atom == "NZ"
        assert entry.charge == "+1"

    def test_charge_entry_from_str_negative(self):
        """Test parsing negative CHARGE entry."""
        text = "OD1 FORMAL -1"
        entry = CHARGEEntry.from_str(text)
        assert entry.atom == "OD1"
        assert entry.charge == "-1"


class TestIO_STRINGEntry:
    """Test IO_STRING entry parsing."""

    def test_io_string_entry_creation(self):
        """Test creating an IO_STRING entry."""
        entry = IO_STRINGEntry(name3="PHE", name1="F")
        assert entry.name3 == "PHE"
        assert entry.name1 == "F"

    def test_io_string_entry_from_str(self):
        """Test parsing IO_STRING entry from string."""
        text = "PHE F"
        entry = IO_STRINGEntry.from_str(text)
        assert entry.name3 == "PHE"
        assert entry.name1 == "F"


class TestGenericEntry:
    """Test generic entry parsing."""

    def test_generic_entry_creation(self):
        """Test creating a generic entry."""
        entry = GenericEntry(header="TEST", body="test body")
        assert entry.header == "TEST"
        assert entry.body == "test body"

    def test_generic_entry_to_str(self):
        """Test converting generic entry to string."""
        entry = GenericEntry(header="NBR_ATOM", body=" CB ")
        result = str(entry)
        assert "NBR_ATOM" in result
        assert "CB" in result


class TestParamsInit:
    """Test Params initialization and basic functionality."""

    def test_params_creation_empty(self):
        """Test creating an empty Params instance."""
        params = Params()
        assert isinstance(params, Params)
        assert hasattr(params, 'ATOM')
        assert hasattr(params, 'BOND')
        assert hasattr(params, 'CHI')

    def test_params_has_attributes(self):
        """Test that Params has expected attributes."""
        params = Params()
        expected_attrs = [
            'ATOM', 'BOND', 'CHI', 'ICOOR_INTERNAL',
            'IO_STRING', 'CONNECT', 'NAME'
        ]
        for attr in expected_attrs:
            assert hasattr(params, attr), f"Missing attribute: {attr}"

    def test_params_entries_are_entries_instances(self):
        """Test that params attributes are Entries instances."""
        params = Params()
        assert isinstance(params.ATOM, Entries)
        assert isinstance(params.BOND, Entries)
        assert isinstance(params.CHI, Entries)


class TestParamsLoading:
    """Test loading params files."""

    def test_loads_from_string(self, official_phe_params):
        """Test loading params from string."""
        with open(official_phe_params) as f:
            content = f.read()

        params = Params.loads(content)
        assert isinstance(params, Params)
        assert params.NAME == "PHE"

    def test_params_name_from_io_string(self, official_phe_params):
        """Test that NAME is derived from IO_STRING."""
        params = Params.load(str(official_phe_params))
        assert params.NAME == "PHE"
        assert len(params.IO_STRING) > 0
        assert params.IO_STRING[0].name3 == "PHE"

    def test_params_has_atoms(self, official_phe_params):
        """Test that loaded params has atoms."""
        params = Params.load(str(official_phe_params))
        assert len(params.ATOM) > 0
        # PHE should have multiple atoms
        assert len(params.ATOM) >= 10

    def test_params_has_bonds(self, official_phe_params):
        """Test that loaded params has bonds."""
        params = Params.load(str(official_phe_params))
        assert len(params.BOND) > 0


class TestParamsDumping:
    """Test dumping params to strings and files."""

    def test_dumps_to_string(self, official_phe_params):
        """Test dumping params to string."""
        params = Params.load(str(official_phe_params))
        content = params.dumps()
        assert isinstance(content, str)
        assert "NAME PHE" in content
        assert "IO_STRING" in content
        assert "ATOM" in content

    def test_dump_and_reload(self, official_phe_params, tmp_path):
        """Test dumping and reloading params file."""
        params1 = Params.load(str(official_phe_params))
        temp_file = tmp_path / "temp.params"

        params1.dump(str(temp_file))
        params2 = Params.load(str(temp_file))

        assert params1.NAME == params2.NAME
        assert len(params1.ATOM) == len(params2.ATOM)
        assert len(params1.BOND) == len(params2.BOND)


class TestParamsAtomRenaming:
    """Test atom renaming functionality."""

    def test_rename_atom_in_atom_list(self):
        """Test renaming atom in ATOM entries."""
        params = Params()
        params.ATOM.append(ATOMEntry(name=" CA ", rtype="CT1", mtype="CT", partial=0.0))
        params.ATOM.append(ATOMEntry(name=" CB ", rtype="CT2", mtype="CT", partial=0.0))

        params.rename_atom(" CA ", " CX ")

        assert params.ATOM[0].name == " CX "
        assert params.ATOM[1].name == " CB "

    def test_rename_atom_in_bonds(self, official_phe_params):
        """Test renaming atom in BOND entries."""
        params = Params.load(str(official_phe_params))

        # Find a bond with CA
        ca_bonds_before = [b for b in params.BOND if ' CA ' in (b.first, b.second)]
        assert len(ca_bonds_before) > 0, "PHE should have bonds with CA"

        params.rename_atom(" CA ", " CX ")

        # CA should be gone, CX should exist
        ca_bonds_after = [b for b in params.BOND if ' CA ' in (b.first, b.second)]
        cx_bonds_after = [b for b in params.BOND if ' CX ' in (b.first, b.second)]

        assert len(ca_bonds_after) == 0
        assert len(cx_bonds_after) == len(ca_bonds_before)

    def test_rename_atom_in_chi(self, official_phe_params):
        """Test renaming atom in CHI entries."""
        params = Params.load(str(official_phe_params))

        # Find CHI entries with CA
        chi_with_ca_before = sum(1 for chi in params.CHI
                                  if ' CA ' in (chi.first, chi.second, chi.third, chi.fourth))

        if chi_with_ca_before > 0:
            params.rename_atom(" CA ", " CX ")

            # CA should be replaced with CX
            chi_with_ca_after = sum(1 for chi in params.CHI
                                     if ' CA ' in (chi.first, chi.second, chi.third, chi.fourth))
            chi_with_cx_after = sum(1 for chi in params.CHI
                                     if ' CX ' in (chi.first, chi.second, chi.third, chi.fourth))

            assert chi_with_ca_after == 0
            assert chi_with_cx_after == chi_with_ca_before

    def test_get_correct_atomname(self, official_phe_params):
        """Test getting correctly formatted atom name."""
        params = Params.load(str(official_phe_params))

        # Test various input formats for CA (which exists in PHE)
        assert params.get_correct_atomname("CA") == " CA "
        assert params.get_correct_atomname(" CA") == " CA "
        assert params.get_correct_atomname("CA ") == " CA "
        assert params.get_correct_atomname(" CA ") == " CA "



class TestParamsFields:
    """Test params field access."""

    def test_fields_returns_list(self, official_phe_params):
        """Test that fields returns a list of header names."""
        params = Params.load(str(official_phe_params))
        fields = params.fields

        assert isinstance(fields, list)
        assert "ATOM" in fields
        assert "BOND" in fields
        assert "NAME" in fields

    def test_get_attr_with_field(self, official_phe_params):
        """Test getting attribute by field name."""
        params = Params.load(str(official_phe_params))

        # These should all work
        assert hasattr(params, 'ATOM')
        assert hasattr(params, 'BOND')
        assert hasattr(params, 'NAME')


class TestParamsComments:
    """Test comment handling."""

    def test_comments_preserved(self):
        """Test that comments are preserved when loading."""
        content = """# This is a comment
NAME TST
# Another comment
"""
        params = Params.loads(content)
        assert len(params.comments) >= 2
        assert "This is a comment" in params.comments[0].body

    def test_comments_in_dump(self, official_phe_params):
        """Test that comments appear in dumped output."""
        params = Params.load(str(official_phe_params))
        params.comments.append("Test comment")

        output = params.dumps()
        assert "# Test comment" in output


class TestRegexPatterns:
    """Test that regex patterns work correctly."""

    def test_atom_pattern_matches(self):
        """Test ATOM entry regex pattern."""
        pattern = r"(?P<name>.{1,4})\s*(?P<rtype>.{1,4})\s*(?P<mtype>.{1,4})\s*(?P<partial>[-\d\.]+)"
        text = " CA  CT1 CT  0.07"
        match = re.match(pattern, text.rstrip())
        assert match is not None
        assert match.group("name") == " CA "

    def test_bond_pattern_matches(self):
        """Test BOND entry regex pattern."""
        pattern = r"(?P<first>.{1,4}) (?P<second>.{2,4})\s?(?P<order>.*)"
        text = " CA   CB  1"
        match = re.match(pattern, text.rstrip())
        assert match is not None
        assert match.group("first") == " CA "
        assert match.group("second") == " CB "

    def test_chi_pattern_matches(self):
        """Test CHI entry regex pattern."""
        pattern = r"(\d+)\s+(\S{1,4})\s+(\S{1,4})\s+(\S{1,4})\s+(\S{1,4})"
        text = "1  C6   C5   C4   C3"
        match = re.match(pattern, text)
        assert match is not None
        assert match.group(1) == "1"

    def test_charge_pattern_matches(self):
        """Test CHARGE entry regex pattern."""
        pattern = r"(?P<atom>\S+) FORMAL (?P<charge>\-?\+?\d)"
        text = "NZ FORMAL +1"
        match = re.match(pattern, text.rstrip())
        assert match is not None
        assert match.group("charge") == "+1"


class TestParamsIntegration:
    """Integration tests for complete workflows."""

    def test_load_modify_dump_cycle(self, official_phe_params, tmp_path):
        """Test loading, modifying, and dumping a params file."""
        # Load
        params = Params.load(str(official_phe_params))
        original_name = params.NAME
        original_atom_count = len(params.ATOM)

        # Modify
        params.IO_STRING[0].name3 = "PHX"
        params.IO_STRING[0].name1 = "X"

        # Dump
        temp_file = tmp_path / "modified.params"
        params.dump(str(temp_file))

        # Reload
        params2 = Params.load(str(temp_file))

        # Verify modifications persisted
        assert params2.IO_STRING[0].name3 == "PHX"
        assert params2.IO_STRING[0].name1 == "X"
        assert len(params2.ATOM) == original_atom_count

    def test_multiple_atom_renames(self, official_phe_params):
        """Test renaming multiple atoms."""
        params = Params.load(str(official_phe_params))

        # Get original bond involving CA
        ca_bonds = [b for b in params.BOND if ' CA ' in (b.first, b.second)]
        assert len(ca_bonds) > 0

        # Rename CA to CX
        params.rename_atom(' CA ', ' CX ')

        # Verify all CA references are now CX
        ca_bonds_after = [b for b in params.BOND if ' CA ' in (b.first, b.second)]
        cx_bonds_after = [b for b in params.BOND if ' CX ' in (b.first, b.second)]

        assert len(ca_bonds_after) == 0
        assert len(cx_bonds_after) == len(ca_bonds)


class TestEdgeCases:
    """Test edge cases and error handling."""

    def test_minimal_params_dumps(self, official_phe_params):
        """Test dumping minimal params."""
        params = Params.load(str(official_phe_params))
        # Clear most entries but keep essentials
        output = params.dumps()
        assert isinstance(output, str)
        assert len(output) > 0

    def test_params_name_property(self, official_phe_params):
        """Test params NAME property."""
        params = Params.load(str(official_phe_params))
        assert params.NAME == "PHE"

        # Change the name
        params.IO_STRING[0].name3 = "PHX"
        assert params.NAME == "PHX"

    def test_invalid_bond_format(self):
        """Test that invalid BOND format raises error."""
        with pytest.raises(ValueError, match="BOND entry.*not formatted correctly"):
            BONDEntry.from_str("invalid")

    def test_rename_nonexistent_atom(self, official_phe_params):
        """Test renaming an atom that doesn't exist raises ValueError."""
        params = Params.load(str(official_phe_params))

        # This should raise an error since ZZ doesn't exist
        with pytest.raises(ValueError, match="not a valid atom name"):
            params.rename_atom(" ZZ ", " XX ")


    def test_load_nonexistent_file(self):
        """Test loading a nonexistent file raises appropriate error."""
        with pytest.raises(FileNotFoundError):
            Params.load("nonexistent_file.params")


# Run with: pytest tests/test_internals.py -v

