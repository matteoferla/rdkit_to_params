"""Tests for virtual shadow atom generation in ring systems.

Virtual shadow atoms bridge CUT_BOND gaps in rings so that ring torsion (NU)
angles and the ``ring_close`` score term can maintain ring geometry during
relaxation.

Reference: GTS params at ``PROTDES-421/glucosidase/GTS_proto.params``
"""

import pytest

from rdkit_to_params import Params

# ### Fixtures ########################################################################

GLUCOSE_SMILES = "OC1C(O)C(OC(O)C1O)CO"  # ← standard glucose
THF_SMILES = "C1CCOC1"  # ← tetrahydrofuran (5-membered ring)
BENZENE_SMILES = "c1ccccc1"  # ← fully aromatic — should be skipped


@pytest.fixture
def glucose_params() -> Params:
    """Glucose params with virtual shadow atoms."""
    p = Params.from_smiles(GLUCOSE_SMILES, name="GLC")
    return p


@pytest.fixture
def thf_params() -> Params:
    """Tetrahydrofuran params with virtual shadow atoms."""
    p = Params.from_smiles(THF_SMILES, name="THF")
    return p


@pytest.fixture
def benzene_params() -> Params:
    """Benzene params — aromatic ring, should NOT get virtual shadows by default."""
    p = Params.from_smiles(BENZENE_SMILES, name="BNZ")
    return p


# ### TestMakeVirtualName ##############################################################


class TestMakeVirtualName:
    """Test the virtual name generation helper."""

    def test_standard_name(self, glucose_params: Params) -> None:
        """O5 → VO5 (prepend V, strip whitespace, pad to 4 chars)."""
        result = glucose_params._make_virtual_name("O5")
        assert result.strip() == "VO5"

    def test_name_too_long_raises(self, glucose_params: Params) -> None:
        """A 4-char real name would produce a 5-char virtual name → ValueError."""
        with pytest.raises(ValueError, match="too long"):
            glucose_params._make_virtual_name("LONG")


# ### TestGlucoseVirtualShadows ########################################################


class TestGlucoseVirtualShadows:
    """Glucose (6-membered ring) should produce the full virtual shadow infrastructure."""

    def test_has_virtual_shadow_entries(self, glucose_params: Params) -> None:
        """Exactly 2 VIRTUAL_SHADOW entries (one per cut-bond endpoint)."""
        assert len(glucose_params.VIRTUAL_SHADOW) == 2

    def test_has_nu_entries(self, glucose_params: Params) -> None:
        """A 6-membered ring should produce 5 NU angles."""
        assert len(glucose_params.NU) == 5

    def test_has_add_ring_entry(self, glucose_params: Params) -> None:
        """Exactly 1 ADD_RING entry listing all 6 ring atoms."""
        assert len(glucose_params.ADD_RING) == 1
        ring_body = glucose_params.ADD_RING[0].body
        # Ring size indicator "1" plus 6 atom names
        parts = ring_body.split()
        assert len(parts) == 7  # ring_idx + 6 atoms

    def test_virtual_atoms_in_atom_list(self, glucose_params: Params) -> None:
        """Virtual atoms should appear in the ATOM entries with VIRT type."""
        atom_names = [a.name.strip() for a in glucose_params.ATOM]
        vs_names = [vs.virtual_atom.strip() for vs in glucose_params.VIRTUAL_SHADOW]
        for vn in vs_names:
            assert vn in atom_names, f"Virtual atom {vn} not found in ATOM list"

    def test_virtual_atoms_have_virt_rtype(self, glucose_params: Params) -> None:
        """Virtual atoms must have rtype=VIRT and mtype=VIRT."""
        vs_names = {vs.virtual_atom.strip() for vs in glucose_params.VIRTUAL_SHADOW}
        for atom_entry in glucose_params.ATOM:
            if atom_entry.name.strip() in vs_names:
                assert atom_entry.rtype.strip() == "VIRT"

    def test_virtual_bonds_present(self, glucose_params: Params) -> None:
        """Virtual atoms should be bonded to the real atom across the cut."""
        vs_pairs = [
            (vs.virtual_atom.strip(), vs.shadow_atom.strip())
            for vs in glucose_params.VIRTUAL_SHADOW
        ]
        bond_pairs = set()
        for b in glucose_params.BOND:
            bond_pairs.add((b.first.strip(), b.second.strip()))
            bond_pairs.add((b.second.strip(), b.first.strip()))
        # For each virtual/shadow pair, the virtual should be bonded to
        # the atom across the cut (not to its shadow directly)
        for virt, _shadow in vs_pairs:
            # The virtual atom must appear in at least one bond
            bonded = any(virt in pair for pair in bond_pairs)
            assert bonded, f"Virtual atom {virt} has no BOND entry"

    def test_virtual_icoor_present(self, glucose_params: Params) -> None:
        """Each virtual atom must have an ICOOR_INTERNAL entry."""
        vs_names = {vs.virtual_atom.strip() for vs in glucose_params.VIRTUAL_SHADOW}
        icoor_children = {ic.child.strip() for ic in glucose_params.ICOOR_INTERNAL}
        for vn in vs_names:
            assert vn in icoor_children, f"Virtual atom {vn} has no ICOOR_INTERNAL"

    def test_nu_first_starts_with_virtual(self, glucose_params: Params) -> None:
        """The first NU angle should start with the first virtual atom."""
        first_nu = glucose_params.NU[0]
        vs_names = {vs.virtual_atom.strip() for vs in glucose_params.VIRTUAL_SHADOW}
        assert first_nu.first.strip() in vs_names

    def test_nu_last_ends_with_virtual(self, glucose_params: Params) -> None:
        """The last NU angle should end with the second virtual atom."""
        last_nu = glucose_params.NU[-1]
        vs_names = {vs.virtual_atom.strip() for vs in glucose_params.VIRTUAL_SHADOW}
        assert last_nu.fourth.strip() in vs_names

    def test_cyclic_property(self, glucose_params: Params) -> None:
        """PROPERTIES should contain CYCLIC."""
        assert len(glucose_params.PROPERTIES) > 0
        assert "CYCLIC" in glucose_params.PROPERTIES[0].values

    def test_roundtrip_dumps_loads(self, glucose_params: Params) -> None:
        """Dump and re-load should preserve virtual shadow entries."""
        text = glucose_params.dumps()
        reloaded = Params.loads(text)
        assert len(reloaded.VIRTUAL_SHADOW) == len(glucose_params.VIRTUAL_SHADOW)
        assert len(reloaded.NU) == len(glucose_params.NU)
        assert len(reloaded.ADD_RING) == len(glucose_params.ADD_RING)


# ### TestFiveRing #####################################################################


class TestFiveRing:
    """THF (5-membered ring) should produce 2 VIRTUAL_SHADOW and 4 NU entries."""

    def test_has_virtual_shadows(self, thf_params: Params) -> None:
        assert len(thf_params.VIRTUAL_SHADOW) == 2

    def test_has_correct_nu_count(self, thf_params: Params) -> None:
        """5-membered ring → N-1 = 4 NU angles."""
        assert len(thf_params.NU) == 4

    def test_has_add_ring(self, thf_params: Params) -> None:
        assert len(thf_params.ADD_RING) == 1
        parts = thf_params.ADD_RING[0].body.split()
        assert len(parts) == 6  # ring_idx + 5 atoms


# ### TestAromaticSkip #################################################################


class TestAromaticSkip:
    """Fully aromatic rings should be skipped by default."""

    def test_benzene_no_virtual_shadows(self, benzene_params: Params) -> None:
        """Benzene should NOT get VIRTUAL_SHADOW entries by default."""
        assert len(benzene_params.VIRTUAL_SHADOW) == 0

    def test_benzene_no_nu(self, benzene_params: Params) -> None:
        """Benzene should NOT get NU entries by default."""
        assert len(benzene_params.NU) == 0

    def test_benzene_forced_aromatic(self) -> None:
        """With aromatic=True, benzene SHOULD get virtual shadows."""
        p = Params.from_smiles(BENZENE_SMILES, name="BNZ")
        p.add_ring_virtual_shadows(aromatic=True)
        # This will add a second set since convert_mol already ran with aromatic=False
        # So just check they exist
        assert len(p.VIRTUAL_SHADOW) >= 2


# ### TestRenameAtomVirtualShadow ######################################################


class TestRenameAtomVirtualShadow:
    """rename_atom must update NU and VIRTUAL_SHADOW entries."""

    def test_rename_updates_nu_entries(self, glucose_params: Params) -> None:
        """Renaming an atom appearing in NU entries should propagate."""
        # Find an atom that appears in NU entries
        first_nu = glucose_params.NU[0]
        old_second = first_nu.second.strip()
        new_name = "X99"
        glucose_params.rename_atom(old_second, new_name)
        # Verify the NU entry was updated
        updated_nu = glucose_params.NU[0]
        assert updated_nu.second.strip() == new_name.strip()

    def test_rename_updates_virtual_shadow_entries(self, glucose_params: Params) -> None:
        """Renaming a shadow atom should update VIRTUAL_SHADOW."""
        shadow_name = glucose_params.VIRTUAL_SHADOW[0].shadow_atom.strip()
        new_name = "Y99"
        glucose_params.rename_atom(shadow_name, new_name)
        updated_shadow = glucose_params.VIRTUAL_SHADOW[0].shadow_atom.strip()
        assert updated_shadow == new_name.strip()


# ### TestIdempotency ##################################################################


class TestIdempotency:
    """Calling add_ring_virtual_shadows twice should not duplicate entries."""

    def test_idempotent(self) -> None:
        p = Params.from_smiles(GLUCOSE_SMILES, name="GLC")
        n_vs = len(p.VIRTUAL_SHADOW)
        n_nu = len(p.NU)
        # Call again — should be a no-op
        p.add_ring_virtual_shadows()
        assert len(p.VIRTUAL_SHADOW) == n_vs
        assert len(p.NU) == n_nu


# ### TestFromFileReconstruction #######################################################


class TestFromFileReconstruction:
    """Params loaded from file (no mol) can reconstruct rings from CUT_BOND + BOND."""

    def test_file_loaded_params_get_rings(self, glucose_params: Params) -> None:
        """Roundtrip through dumps/loads loses the mol, but rings are still detectable."""
        text = glucose_params.dumps()
        reloaded = Params.loads(text)
        assert reloaded.mol is None
        # The ring info is already present from the original, so check it persisted
        ring_atoms = reloaded._get_ring_atom_names()
        assert len(ring_atoms) >= 1
        assert len(ring_atoms[0]) >= 5  # at least 5 atoms in the ring
