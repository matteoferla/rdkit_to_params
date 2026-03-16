########################################################################################################################
__doc__ = """
    THis the core conversion. the method ``convert_mol`` does something you'd never have guessed.
    There are probably loads of mistakes here ---or more correctly
    in ``_RDKitPrepMixin()._add_rtypes`` or ``_RDKitPrepMixin()._add_genrtypes``.
    """

########################################################################################################################

import random
from collections import defaultdict, deque, namedtuple
from typing import Any

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from rdkit_to_params.entries import Entries
from rdkit_to_params.rdkitside._rdkit_prep import _RDKitPrepMixin
from rdkit_to_params.rdkitside.utilities import DummyMasker


class _RDKitCovertMixin(_RDKitPrepMixin):
    _Measure = namedtuple("_Measure", ["distance", "angle", "torsion"])

    # Type hints for attributes from other mixins
    comments: Entries
    ICOOR_INTERNAL: Entries
    CONNECT: Entries
    CHARGE: Entries
    CHI: Entries
    PROTON_CHI: Entries
    CUT_BOND: Entries
    BOND: Entries
    NBR_ATOM: Entries
    NBR_RADIUS: Entries
    PDB_ROTAMERS: Entries
    ATOM: Entries
    VIRTUAL_SHADOW: Entries
    NU: Entries
    ADD_RING: Entries
    PROPERTIES: Entries

    greekification = True  # controls whether to change the atom names to greek for amino acids for CB and beyond.
    auto_ref = True  # auto-add NUMERIC_PROPERTY REFERENCE for ncAAs

    _CONVERT_SKIP = frozenset({
        "#", "TYPE", "comment", "IO_STRING", "ROTAMER_AA",
        "AA", "PROPERTIES", "VARIANT", "<UNKNOWN>",
        "RAMA_PREPRO_FILENAME", "RAMA_PREPRO_RESNAME",
    })
    # ↓ backbone atoms that should never appear in CHI entries for amino acids
    _BACKBONE_ATOMS = frozenset({"N", "CA", "C", "O", "H", "HA", "UPPER", "LOWER", "OXT"})

    def convert_mol(self, pcharge_prop_name: str = "_GasteigerCharge") -> None:
        """Convert the RDKit mol into Rosetta params entries."""
        self._ensure_conformer()
        title = self._resolve_name()
        self.log.debug(f"{title} is being converted (`.convert_mol`)")
        self._handle_fragments(title)
        self.comments.append(Chem.MolToSmiles(self.mol))
        if len(self.TYPE) == 0:
            self.TYPE.append("LIGAND")
        self._clear_regenerated_entries()
        self._build_topology(pcharge_prop_name)
        self._find_centroid()
        if self.mol.GetRingInfo().NumRings() > 0:
            self.add_ring_virtual_shadows()
        if self.auto_ref and self.is_aminoacid():
            self.add_ref_energy()  # type: ignore[attr-defined]  # ← mixin on Params

    def _ensure_conformer(self) -> None:
        """Embed + optimise if the mol has no conformers."""
        if self.mol.GetConformers():
            return
        self.log.warning("The molecule passed has no conformers... Adding.")
        with DummyMasker(self.mol):
            AllChem.EmbedMolecule(self.mol)  # type: ignore[attr-defined]
            AllChem.MMFFOptimizeMolecule(self.mol)  # type: ignore[attr-defined]

    def _resolve_name(self) -> str:
        """Set self.NAME from mol _Name prop or PDBInfo; return the title."""
        title = self.mol.GetProp("_Name").strip() if self.mol.HasProp("_Name") else "Unnamed molecule"
        if len(title) == 3:
            self.NAME = title
        else:
            self.comments.append(title)
            self.NAME = self._get_resn_from_PDBInfo()
        return title

    def _handle_fragments(self, title: str) -> None:
        """Keep only the largest fragment if the mol is fragmented."""
        if len(Chem.GetMolFrags(self.mol)) > 1:
            self.log.warning(f"{title} is split in {len(Chem.GetMolFrags(self.mol))}")
            self.mol = Chem.GetMolFrags(self.mol, asMols=True)[0]

    def _clear_regenerated_entries(self) -> None:
        """Clear all entries that will be regenerated, preserving metadata entries."""
        for attr_name in Entries.choices:
            if attr_name in self._CONVERT_SKIP or not hasattr(self, attr_name):
                continue
            getattr(self, attr_name).clear()

    def _build_topology(self, pcharge_prop_name: str) -> None:
        """Build ICOOR, ATOM, CONNECT, CHI, and BOND entries."""
        if self.mol.GetNumAtoms() >= 3:
            self._parse_icoors()
            self._parse_atoms(pcharge_prop_name)
            self._parse_rotatables()
        else:
            self._parse_w_virtuals()
            self._parse_atoms(pcharge_prop_name)
        self._parse_bonds()

    ## ============= internal coordinates ==============================================================================

    def _parse_icoors(self) -> None:
        self.log.debug("Filling ICOOR")
        self._undescribed = deque(self.mol.GetAtoms())
        self.ordered_atoms: list[Chem.Atom] = []
        self._rotation_count = 0
        if self.is_aminoacid():
            self._icoor_aminoacid_bootstrap()
        else:
            self._icoor_ligand_bootstrap()
        self.log.debug("Parsing forth line onwards")
        self._icoor_remaining_atoms()

    def _icoor_aminoacid_bootstrap(self) -> None:
        """Emit the first 6 ICOOR lines for an amino acid (N, CA, C, UPPER, O, LOWER)."""
        N = self.get_atom_by_name("N")
        CA = self.get_atom_by_name("CA")
        C = self.get_atom_by_name("C")
        self._add_icoor([N, N, CA, C])
        self._add_icoor([CA, N, CA, C])
        self._add_icoor([C, CA, N, C])
        self._add_icoor([self.get_atom_by_name("UPPER"), C, CA, N])
        self._add_icoor([self.get_atom_by_name("O"), C, CA, self.get_atom_by_name("UPPER")])
        self._add_icoor([self.get_atom_by_name("LOWER"), N, CA, C])

    def _icoor_ligand_bootstrap(self) -> None:
        """Find a valid root triplet and emit the first 3 ICOOR lines for a ligand."""
        while self._undescribed:
            self._prevent_overcycling()
            atom = self._undescribed[0]
            if atom.GetSymbol() == "*" or atom.GetProp("_rType") == "VIRT":
                self.log.debug("Dummy or Virtual not allowed in first three lines...")
                self._undescribed.rotate(-1)
                self._rotation_count += 1
                continue
            triplet = self._find_root_triplet(atom)
            if triplet is None:
                self._undescribed.rotate(-1)
                self._rotation_count += 1
                continue
            root, parent, grandparent = triplet
            self._add_icoor([root, root, parent, grandparent])
            self._add_icoor([parent, root, parent, grandparent])
            self._add_icoor([grandparent, parent, root, grandparent])
            break

    def _find_root_triplet(self, atom: Chem.Atom) -> tuple[Chem.Atom, Chem.Atom, Chem.Atom] | None:
        """Find a root-parent-grandparent triplet for ICOOR bootstrap, or None."""
        for parent_atom in self._get_unseen_neighbors(atom, [], True):
            siblings = self._get_unseen_neighbors(parent_atom, [atom], True)
            if siblings:
                return (atom, parent_atom, siblings[0])
        self.log.debug(f"First row: No siblings for {self._get_PDBInfo_atomname(atom)}!")
        return None

    def _icoor_remaining_atoms(self) -> None:
        """BFS-walk remaining undescribed atoms and emit ICOOR entries."""
        while self._undescribed:
            self._prevent_overcycling()
            atom = self._undescribed[0]
            quadruplet = self._find_icoor_quadruplet(atom)
            if quadruplet is None:
                self._undescribed.rotate(-1)
                self._rotation_count += 1
                continue
            self._add_icoor(quadruplet)
            self._rotation_count -= 1

    def _described_neighbors(self, atom: Chem.Atom, seen: list[Chem.Atom]) -> list[Chem.Atom]:
        """Return neighbours of ``atom`` that are already in ``ordered_atoms``."""
        described_idx = {d.GetIdx() for d in self.ordered_atoms}
        return [n for n in self._get_unseen_neighbors(atom, seen, False)
                if n.GetIdx() in described_idx]

    def _find_icoor_quadruplet(self, atom: Chem.Atom) -> list[Chem.Atom] | None:
        """Find [child, parent, sibling, ref] for an ICOOR entry, or None if not yet possible."""
        parents = self._described_neighbors(atom, [])
        if not parents:
            self.log.debug(f"Other row: No parent for {self._get_PDBInfo_atomname(atom)}!")
            return None
        parent = parents[0]
        siblings = self._described_neighbors(parent, [atom])
        if not siblings:
            self.log.debug(f"Other row: No siblings for  {self._get_PDBInfo_atomname(atom)}!")
            return None
        if len(siblings) >= 2:
            return [atom, parent, siblings[0], siblings[1]]
        # only 1 sibling — need a cousin
        cousins = self._described_neighbors(siblings[0], [parent, atom])
        if not cousins:
            self.log.debug(f"Other row: No cousins for  {self._get_PDBInfo_atomname(atom)}")
            return None
        return [atom, parent, siblings[0], cousins[0]]

    def _add_icoor(self, atoms: list[Chem.Atom]) -> None:
        # due to madness the same atom objects differ.
        self._undescribed.remove(
            [this for this in self._undescribed if atoms[0].GetIdx() == this.GetIdx()][0]
        )
        self.ordered_atoms.append(atoms[0])
        m = self._get_measurements(self.mol.GetConformer(), *atoms)
        self.ICOOR_INTERNAL.append(
            dict(
                child=self._get_PDBInfo_atomname(atoms[0]),
                phi=m.torsion,
                theta=m.angle,
                distance=m.distance,
                parent=self._get_PDBInfo_atomname(atoms[1]),
                second_parent=self._get_PDBInfo_atomname(atoms[2]),
                third_parent=self._get_PDBInfo_atomname(atoms[3]),
            )
        )

    def _prevent_overcycling(self) -> None:
        # prevent overcycling...
        self.log.debug(f"Current deque rotation count {self._rotation_count}")
        if self._rotation_count > 5 * self.mol.GetNumAtoms():
            d = [f"{a.GetIdx()}: {self._get_PDBInfo_atomname(a)}" for a in self.ordered_atoms]
            u = [f"{a.GetIdx()}: {self._get_PDBInfo_atomname(a)}" for a in self._undescribed]
            raise StopIteration(f"Too many cycles. described: {d}, undescribed: {u}")
        elif self._rotation_count % self.mol.GetNumAtoms() == 0:
            # it has done a full rotation.
            if self._rotation_count != 0:
                self.log.debug("Ordering atoms with a shuffle first")
                random.shuffle(self._undescribed)  # mix up equal tier ones
            else:
                self.log.debug("Ordering atoms without a shuffle first")
            if len(self.ordered_atoms) != 0:
                self.log.debug(
                    "Ordering atoms by distance from root, and being not hydrogen or dummy"
                )
                root = self.ordered_atoms[0]

                def distance_penaltifier(atom):
                    return len(Chem.GetShortestPath(self.mol, root.GetIdx(), atom.GetIdx()))
            else:
                self.log.debug("Ordering atoms by being not hydrogen or dummy")

                def distance_penaltifier(atom):
                    return 0

            def hydrogen_penaltifier(atom):
                return 0 if atom.GetSymbol() != "H" else 100

            def dummy_penaltifier(atom):
                return 0 if atom.GetAtomicNum() == 0 else 200

            def scorer(atom):
                return (
                    distance_penaltifier(atom)
                    + hydrogen_penaltifier(atom)
                    + dummy_penaltifier(atom)
                )

            self._undescribed = deque(sorted(self._undescribed, key=scorer))

    def _get_measurements(
        self,
        conf: Chem.Conformer,
        a: Chem.Atom,
        b: Chem.Atom,
        c: Chem.Atom,
        d: Chem.Atom,
    ) -> Any:
        ai = a.GetIdx()
        bi = b.GetIdx()
        ci = c.GetIdx()
        di = d.GetIdx()
        dist = 0
        angle = 0
        tor = 0
        try:
            dist = Chem.rdMolTransforms.GetBondLength(conf, ai, bi)  # type: ignore[attr-defined]
            angle = 180 - Chem.rdMolTransforms.GetAngleDeg(conf, ai, bi, ci)  # type: ignore[attr-defined]
            tor = Chem.rdMolTransforms.GetDihedralDeg(conf, ai, bi, ci, di)  # type: ignore[attr-defined]
        except Exception:
            pass
        if str(tor) == "nan":  # quicker than isnan.
            tor = 0
        return self._Measure(distance=dist, angle=angle, torsion=tor)

    ## ============= Atom entries ======================================================================================

    def _parse_atoms(self, pcharge_prop_name: str = "_GasteigerCharge") -> None:
        self.log.debug("Filling ATOM records...")
        for atom in self.ordered_atoms:
            self._parse_atom(atom, pcharge_prop_name)

    def _parse_atom(self, atom: Chem.Atom, pcharge_prop_name: str = "_GasteigerCharge") -> None:
        self.log.debug(f"Parsing {atom.GetSymbol()} at position {atom.GetIdx()}")
        if atom.GetAtomicNum() == 0:  # * (Chem.Atom) or R (Chem.QueryAtom)
            neighbor = atom.GetNeighbors()[0]
            n_name = self._get_PDBInfo_atomname(neighbor)
            if self.is_aminoacid() and neighbor.GetSymbol() == "N":
                # atom_name, index, connect_type, connect_name
                self.CONNECT.append([n_name, 1, "LOWER_CONNECT", "LOWER"])
            elif self.is_aminoacid() and neighbor.GetSymbol() == "C":
                # atom_name, index, connect_type, connect_name
                self.CONNECT.append([n_name, 2, "UPPER_CONNECT", "UPPER"])
            elif self.is_aminoacid():
                i = max(3, len(self.CONNECT) + 1)
                self.CONNECT.append([n_name, i, "CONNECT"])
            else:
                self.CONNECT.append([n_name, len(self.CONNECT) + 1, "CONNECT"])
        else:
            d = self._get_atom_descriptors(
                atom, pcharge_prop_name
            )  # dict of 'name', 'rtype': 'mtype', 'partial'
            self.ATOM.append(d)
            formal = atom.GetFormalCharge()
            if formal != 0:
                self.CHARGE.append([d["name"], formal])

    ## ============= Chi entries =======================================================================================

    def _parse_rotatables(self) -> None:
        """Bond-centric CHI generation with PROTON_CHI support.

        Iterates mol bonds to find rotatable single bonds between non-ring,
        non-dummy heavy atoms. For each, builds a 4-atom CHI quadruplet.
        Polar H neighbours get an additional CHI + PROTON_CHI entry.
        For CHI a-b-c-d, Rosetta requires atom_base(c)==b in the ICOOR tree.
        """
        self.log.debug("Filling CHI")
        ring_idxs = self._collect_ring_atom_indices()
        icoor_parent = self._build_icoor_parent_map()
        rotatable_bonds = self._find_rotatable_bonds(ring_idxs, icoor_parent)
        for b_atom, c_atom in rotatable_bonds:
            self._emit_chi(b_atom, c_atom, ring_idxs, icoor_parent)

    def _build_icoor_parent_map(self) -> dict[int, int]:
        """Build atom_idx → parent_idx map from ICOOR_INTERNAL entries.

        For CHI a-b-c-d, Rosetta requires atom_base(c)==b, i.e. the ICOOR parent
        of atom c must be atom b. This map lets us orient bonds correctly.
        """
        name_to_idx: dict[str, int] = {}
        for atom in self.mol.GetAtoms():
            info = atom.GetPDBResidueInfo()
            if info is not None:
                name_to_idx[info.GetName().strip()] = atom.GetIdx()
        parent_map: dict[int, int] = {}
        for ic in self.ICOOR_INTERNAL:
            child_name = ic.child.strip()
            parent_name = ic.parent.strip()
            if child_name == parent_name:
                continue  # ↑ skip self-referential root line
            child_idx = name_to_idx.get(child_name)
            parent_idx = name_to_idx.get(parent_name)
            if child_idx is not None and parent_idx is not None:
                parent_map[child_idx] = parent_idx
        return parent_map

    def _collect_ring_atom_indices(self) -> set[int]:
        """Return set of atom indices that belong to any ring."""
        ring_idxs: set[int] = set()
        for ring_set in self.mol.GetRingInfo().AtomRings():
            ring_idxs.update(ring_set)
        return ring_idxs

    def _is_chi_eligible_bond(self, bond: Chem.Bond, ring_idxs: set[int],
                              ordered_set: set[int]) -> bool:
        """Return True if bond is a rotatable single bond between non-ring heavy atoms.

        For amino acids, bonds where both atoms are backbone are excluded —
        backbone torsions (e.g. N-CA-C-O) should never produce CHI entries.
        """
        if bond.GetBondType() != Chem.BondType.SINGLE:
            return False
        ai, bi = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if ai in ring_idxs or bi in ring_idxs:
            return False
        a_atom, b_atom = bond.GetBeginAtom(), bond.GetEndAtom()
        if a_atom.GetAtomicNum() in (0, 1) or b_atom.GetAtomicNum() in (0, 1):
            return False
        if not (ai in ordered_set and bi in ordered_set):
            return False
        # ↓ reject backbone-only bonds for amino acids (e.g. N-CA, CA-C)
        if self.is_aminoacid():
            a_name = self._get_PDBInfo_atomname(a_atom).strip()
            b_name = self._get_PDBInfo_atomname(b_atom).strip()
            if a_name in self._BACKBONE_ATOMS and b_name in self._BACKBONE_ATOMS:
                return False
        return True

    def _orient_bond_by_icoor(self, bond: Chem.Bond,
                              icoor_parent: dict[int, int]) -> tuple[Chem.Atom, Chem.Atom] | None:
        """Orient a bond (b, c) so icoor_parent[c]==b. Returns None if impossible."""
        ai, bi = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if icoor_parent.get(bi) == ai:
            return bond.GetBeginAtom(), bond.GetEndAtom()
        elif icoor_parent.get(ai) == bi:
            return bond.GetEndAtom(), bond.GetBeginAtom()
        return None

    def _find_rotatable_bonds(self, ring_idxs: set[int],
                              icoor_parent: dict[int, int]) -> list[tuple[Chem.Atom, Chem.Atom]]:
        """Return rotatable single bonds as (b, c) pairs sorted by ICOOR tree position.

        Oriented so that icoor_parent[c] == b — satisfying Rosetta's atom_base constraint.
        """
        ordered_indices = [a.GetIdx() for a in self.ordered_atoms]
        ordered_set = set(ordered_indices)
        seen: set[tuple[int, int]] = set()
        results: list[tuple[int, Chem.Atom, Chem.Atom]] = []
        for bond in self.mol.GetBonds():
            if not self._is_chi_eligible_bond(bond, ring_idxs, ordered_set):
                continue
            pair = (min(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()),
                    max(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))
            if pair in seen:
                continue
            seen.add(pair)
            oriented = self._orient_bond_by_icoor(bond, icoor_parent)
            if oriented is None:
                continue
            b_out, c_out = oriented
            results.append((ordered_indices.index(b_out.GetIdx()), b_out, c_out))
        results.sort(key=lambda t: t[0])
        return [(b, c) for _, b, c in results]

    @staticmethod
    def _is_terminal_methyl(c_atom: Chem.Atom, b_idx: int) -> bool:
        """True if c is a carbon whose only non-b neighbours are H (CH3 rotation)."""
        if c_atom.GetAtomicNum() != 6:
            return False
        return all(n.GetAtomicNum() == 1 or n.GetIdx() == b_idx
                   for n in c_atom.GetNeighbors())

    def _emit_chi(self, b_atom: Chem.Atom, c_atom: Chem.Atom,
                  ring_idxs: set[int], icoor_parent: dict[int, int]) -> None:
        """Build CHI quadruplet a-b-c-d and append. Also emits PROTON_CHI if applicable.

        For atom d, prefer a child of c in the ICOOR tree (so atom_base(d)==c).
        PROTON_CHI references the emitted CHI — no extra CHI entry needed.
        Follows the canonical pattern: SER ``CA CB OG HG`` + ``PROTON_CHI 2``,
        LYS ``CG CD CE NZ`` + ``PROTON_CHI 4``.
        """
        # ↓ terminal methyl rotations (C→H only) are never meaningful CHIs
        if self._is_terminal_methyl(c_atom, b_atom.GetIdx()):
            return
        a_atom = self._pick_flanking_atom(
            b_atom, exclude={c_atom.GetIdx()}, ring_idxs=ring_idxs,
            prefer_heavy=True, icoor_parent=icoor_parent, want_child_of=None,
        )
        if a_atom is None:
            return
        # ↓ for d, prefer an ICOOR child of c (atom_base(d)==c) and heavy
        d_atom = self._pick_flanking_atom(
            c_atom, exclude={b_atom.GetIdx()}, ring_idxs=ring_idxs,
            prefer_heavy=True, icoor_parent=icoor_parent, want_child_of=c_atom.GetIdx(),
        )
        if d_atom is None:
            d_atom = self._pick_flanking_atom(
                c_atom, exclude={b_atom.GetIdx()}, ring_idxs=ring_idxs,
                prefer_heavy=False, icoor_parent=icoor_parent, want_child_of=c_atom.GetIdx(),
            )
        if d_atom is None:
            return
        chi_idx = len(self.CHI) + 1
        self.CHI.append(dict(
            index=chi_idx,
            first=self._get_PDBInfo_atomname(a_atom),
            second=self._get_PDBInfo_atomname(b_atom),
            third=self._get_PDBInfo_atomname(c_atom),
            fourth=self._get_PDBInfo_atomname(d_atom),
        ))
        # ↓ PROTON_CHI: reference this CHI if it involves a polar atom with H
        self._maybe_add_proton_chi(chi_idx, c_atom, d_atom, ring_idxs)

    @staticmethod
    def _count_h_neighbors(atom: Chem.Atom, ring_idxs: set[int]) -> int:
        """Count non-ring hydrogen neighbours of an atom."""
        return sum(1 for n in atom.GetNeighbors()
                   if n.GetAtomicNum() == 1 and n.GetIdx() not in ring_idxs)

    _POLAR_ATOMIC_NUMS = frozenset({7, 8, 16})  # N, O, S

    def _maybe_add_proton_chi(self, chi_idx: int, c_atom: Chem.Atom,
                              d_atom: Chem.Atom, ring_idxs: set[int]) -> None:
        """Emit PROTON_CHI if the CHI involves a polar atom bearing hydrogens.

        Two canonical patterns:
        - d is H and c is polar (e.g. SER: CA-CB-OG-HG)
        - d is polar with H children (e.g. LYS: CG-CD-CE-NZ)
        """
        # pattern 1: terminal polar-H
        if d_atom.GetAtomicNum() == 1 and c_atom.GetAtomicNum() in self._POLAR_ATOMIC_NUMS:
            self._add_proton_chi(chi_idx, c_atom, self._count_h_neighbors(c_atom, ring_idxs))
            return
        # pattern 2: CHI ends at polar atom with H
        if d_atom.GetAtomicNum() not in self._POLAR_ATOMIC_NUMS:
            return
        n_h = self._count_h_neighbors(d_atom, ring_idxs)
        if n_h > 0:
            self._add_proton_chi(chi_idx, d_atom, n_h)

    @staticmethod
    def _score_flanking_candidate(atom: Chem.Atom, ring_idxs: set[int],
                                  icoor_parent: dict[int, int] | None,
                                  want_child_of: int | None) -> tuple[int, int]:
        """Scoring key for flanking atom selection: (ICOOR-child-priority, ring-priority)."""
        is_child = (0 if (want_child_of is not None and icoor_parent is not None
                          and icoor_parent.get(atom.GetIdx()) == want_child_of) else 1)
        is_ring = 0 if atom.GetIdx() not in ring_idxs else 1
        return (is_child, is_ring)

    def _pick_flanking_atom(self, center: Chem.Atom, exclude: set[int],
                            ring_idxs: set[int], prefer_heavy: bool,
                            icoor_parent: dict[int, int] | None = None,
                            want_child_of: int | None = None) -> Chem.Atom | None:
        """Choose the best flanking atom for a CHI quadruplet.

        Prefers ICOOR children of ``want_child_of``, then non-ring atoms.
        """
        candidates = [n for n in center.GetNeighbors()
                      if n.GetIdx() not in exclude and n.GetAtomicNum() != 0]
        if prefer_heavy:
            candidates = [c for c in candidates if c.GetAtomicNum() != 1]
        if not candidates:
            return None
        candidates.sort(key=lambda c: self._score_flanking_candidate(
            c, ring_idxs, icoor_parent, want_child_of))
        return candidates[0]  # type: ignore[no-any-return]  # ← RDKit Atom stubs incomplete

    def _add_proton_chi(self, chi_index: int, polar_atom: Chem.Atom, n_hydrogens: int) -> None:
        """Append a PROTON_CHI entry with sampling appropriate to the polar atom type."""
        symbol = polar_atom.GetSymbol()
        if symbol in ("O", "S"):
            # hydroxyl / thiol: 18 samples every 20°
            samples = list(range(0, 360, 20))
            extra = [0]
        elif symbol == "N" and n_hydrogens >= 3:
            # NH3+ (e.g. lysine): 6 samples
            samples = [60, -60, 180, 0, 120, -120]
            extra = [1, 20]
        else:
            # amine with 1-2 H
            samples = [0, 180]
            extra = [1, 20]
        self.PROTON_CHI.append(dict(chi_index=chi_index, samples=samples, extra=extra))

    def _get_unseen_neighbors(
        self, atom: Chem.Atom, seen: list[Chem.Atom], nondummy: bool = True
    ) -> list[Chem.Atom]:
        neighbors_list = list(atom.GetNeighbors())
        if nondummy:
            neighbors_list = [
                neighbor for neighbor in neighbors_list if neighbor.GetSymbol() != "*"
            ]
        seenIdx = {a.GetIdx() for a in seen}
        return [neighbor for neighbor in neighbors_list if neighbor.GetIdx() not in seenIdx]

    ## ============= Bond entries ======================================================================================

    def _parse_bonds(self) -> None:
        self.log.debug("Filling BOND")
        for bond in self.mol.GetBonds():
            self._parse_bond(bond)
        for ring_set in self.mol.GetRingInfo().AtomRings():
            self.log.debug("Mol has a ring, adding a `CUT_BOND` entry (`.convert_mol`)")
            self.CUT_BOND.append(
                {
                    "first": self._get_PDBInfo_atomname(self.mol.GetAtomWithIdx(ring_set[0])),
                    "second": self._get_PDBInfo_atomname(self.mol.GetAtomWithIdx(ring_set[-1])),
                }
            )

    def _parse_bond(self, bond: Chem.Bond) -> None:
        if any([atom.GetAtomicNum() == 0 for atom in (bond.GetBeginAtom(), bond.GetEndAtom())]):
            return None  # CONNECT.
        order = bond.GetBondTypeAsDouble()
        if order == 0:
            order = 1  # it's a virtual atom, but actually a vanadium
        if order == 1.5:
            order = 4
        else:
            order = int(order)
        begin = bond.GetBeginAtom()
        end = bond.GetEndAtom()
        self.BOND.append(
            [self._get_PDBInfo_atomname(begin), self._get_PDBInfo_atomname(end), order]
        )

    def _get_atom_descriptors(
        self, atom: Chem.Atom, pcharge_prop_name: str
    ) -> dict[str, str | float]:
        return {
            "name": self._get_PDBInfo_atomname(atom),
            "rtype": atom.GetProp("_rType"),
            "mtype": " X  ",
            "partial": atom.GetDoubleProp(pcharge_prop_name),
        }

    def _get_nondummy_neighbors(self, atom: Chem.Atom) -> list[str]:
        # returns list of names!
        return [
            self._get_PDBInfo_atomname(neighbor)
            for neighbor in atom.GetNeighbors()
            if neighbor.GetSymbol() != "*"
        ]

    def _find_centroid(self) -> None:
        a = Chem.Get3DDistanceMatrix(self.mol)
        n = np.linalg.norm(a, axis=0)  # Frobenius norm
        i = int(np.argmin(n))
        s = np.max(a, axis=0)[i]
        self.NBR_ATOM.append(self.mol.GetAtomWithIdx(i).GetPDBResidueInfo().GetName())
        self.NBR_RADIUS.append(str(s))

    # ============ Virtual Shadow Atoms for Rings ========================================================================

    def add_ring_virtual_shadows(self, aromatic: bool = False) -> None:
        """Add virtual shadow atoms, NU angles, and ADD_RING entries for ring systems.

        Virtual shadow atoms are VIRT-type atoms bonded across the CUT_BOND gap.
        They track the real atoms via the ``ring_close`` score term (harmonic, sigma=0.1A).
        Without them, both dihedral and Cartesian relaxation produce ring distortions.

        Called automatically from ``convert_mol()`` for non-aromatic rings.

        :param aromatic: if True, also process fully aromatic rings (default: skip them)
        """
        rings = self._get_ring_atom_names()
        for ring_idx, ring_atoms in enumerate(rings):
            if self._should_skip_ring(ring_atoms, aromatic):
                continue
            self._add_virtual_shadow_for_ring(ring_idx + 1, ring_atoms)
        self._ensure_cyclic_property()

    def _should_skip_ring(self, ring_atoms: list[str], aromatic: bool) -> bool:
        """Return True if ring already has virtual shadows or is aromatic and unwanted."""
        existing = {vs.shadow_atom.strip() for vs in self.VIRTUAL_SHADOW}
        if ring_atoms[0] in existing or ring_atoms[-1] in existing:
            return True
        if not aromatic and self.mol is not None and self._ring_is_aromatic(ring_atoms):
            return True
        return False

    def _ensure_cyclic_property(self) -> None:
        """Add CYCLIC to PROPERTIES if virtual shadows were generated."""
        if len(self.VIRTUAL_SHADOW) == 0:
            return
        if len(self.PROPERTIES) == 0:
            self.PROPERTIES.append("CYCLIC")
        elif "CYCLIC" not in self.PROPERTIES[0].values:
            self.PROPERTIES[0].values.append("CYCLIC")

    def _ring_is_aromatic(self, ring_atom_names: list[str]) -> bool:
        """Check whether ALL bonds in a ring are aromatic using the RDKit mol."""
        if self.mol is None:
            return False
        name_to_idx = self._build_name_to_idx_map()
        for i in range(len(ring_atom_names)):
            a_idx = name_to_idx.get(ring_atom_names[i].strip())
            b_idx = name_to_idx.get(ring_atom_names[(i + 1) % len(ring_atom_names)].strip())
            if a_idx is None or b_idx is None:
                return False
            bond = self.mol.GetBondBetweenAtoms(a_idx, b_idx)
            if bond is None or not bond.GetIsAromatic():
                return False
        return True

    def _build_name_to_idx_map(self) -> dict[str, int]:
        """Build a mapping from stripped PDB atom names to RDKit atom indices."""
        result: dict[str, int] = {}
        for atom in self.mol.GetAtoms():
            info = atom.GetPDBResidueInfo()
            if info is not None:
                result[info.GetName().strip()] = atom.GetIdx()
        return result

    def _get_ring_atom_names(self) -> list[list[str]]:
        """Get ring atom names, ordered so that CUT_BOND is at [0]--[-1].

        Uses RDKit ring info when mol is available, otherwise reconstructs
        from BOND + CUT_BOND entries (for file-loaded params without mol).
        """
        if self.mol is not None:
            return self._get_ring_atom_names_from_mol()
        return self._get_ring_atom_names_from_entries()

    def _get_ring_atom_names_from_mol(self) -> list[list[str]]:
        """Extract ring atom names from RDKit mol, ordered per AtomRings()."""
        rings: list[list[str]] = []
        for ring_set in self.mol.GetRingInfo().AtomRings():
            names = [
                self._get_PDBInfo_atomname(self.mol.GetAtomWithIdx(idx)).strip() for idx in ring_set
            ]
            rings.append(names)
        return rings

    def _get_ring_atom_names_from_entries(self) -> list[list[str]]:
        """Reconstruct ring atom names from BOND + CUT_BOND entries (no mol)."""
        # ← build adjacency from BOND entries
        adj: dict[str, set[str]] = defaultdict(set)
        for bond in self.BOND:
            a, b = bond.first.strip(), bond.second.strip()
            adj[a].add(b)
            adj[b].add(a)
        rings: list[list[str]] = []
        for cut in self.CUT_BOND:
            start, end = cut.first.strip(), cut.second.strip()
            path = self._find_ring_path(adj, start, end)
            if path:
                rings.append(path)
        return rings

    def _find_ring_path(self, adj: dict[str, set[str]], start: str, end: str) -> list[str]:
        """BFS to find a path from start to end through BOND adjacency, excluding the direct edge."""
        queue: deque[list[str]] = deque([[start]])
        visited: set[str] = {start}
        while queue:
            path = queue.popleft()
            current = path[-1]
            for neighbour in adj.get(current, set()):
                # ← skip the direct cut-bond edge
                if current == start and neighbour == end:
                    continue
                if current == end and neighbour == start:
                    continue
                if neighbour == end:
                    return path + [end]
                if neighbour not in visited:
                    visited.add(neighbour)
                    queue.append(path + [neighbour])
        return []

    def _add_virtual_shadow_for_ring(self, ring_idx: int, ring_atoms: list[str]) -> None:
        """Assemble all virtual shadow entries for one ring.

        For ring [a0, a1, ..., a(N-1)] with CUT_BOND at a(N-1)--a0:
        - V_last shadows a(N-1), bonded to a0
        - V_first shadows a0, bonded to a(N-1)
        """
        n = len(ring_atoms)
        first_atom = ring_atoms[0]  # ← a0 — one end of the cut
        last_atom = ring_atoms[-1]  # ← a(N-1) — other end of the cut
        v_last = self._make_virtual_name(last_atom)  # shadows last, bonded to first
        v_first = self._make_virtual_name(first_atom)  # shadows first, bonded to last
        # ← handle name collisions between the two virtuals
        if v_last.strip() == v_first.strip():
            v_first = self._make_virtual_name(first_atom, suffix="2")
        # ### ATOM entries — VIRT type, zero charge
        self.ATOM.append(dict(name=v_last, rtype="VIRT", mtype="VIRT", partial=0.0))
        self.ATOM.append(dict(name=v_first, rtype="VIRT", mtype="VIRT", partial=0.0))
        # ### BOND entries — virtual bonded to real atom across the cut
        self.BOND.append([first_atom, v_last, 1])  # a0 -- V_last
        self.BOND.append([last_atom, v_first, 1])  # a(N-1) -- V_first
        # ### VIRTUAL_SHADOW entries
        self.VIRTUAL_SHADOW.append(dict(virtual_atom=v_last, shadow_atom=last_atom))
        self.VIRTUAL_SHADOW.append(dict(virtual_atom=v_first, shadow_atom=first_atom))
        # ### NU angles — N-1 torsions spanning the ring
        # NU 1:   V_last, a0,   a1,   a2
        # NU k:   a(k-2), a(k-1), a(k), a(k+1)   for k=2..N-2
        # NU N-1: a(N-3), a(N-2), a(N-1), V_first
        extended = [v_last] + ring_atoms + [v_first]
        for i in range(n - 1):
            self.NU.append(
                dict(
                    index=i + 1,
                    first=extended[i],
                    second=extended[i + 1],
                    third=extended[i + 2],
                    fourth=extended[i + 3],
                )
            )
        # ### ADD_RING entry
        ring_body = f"{ring_idx} " + " ".join(ring_atoms)
        self.ADD_RING.append(ring_body)
        # ### ICOOR_INTERNAL for virtual atoms
        self._add_virtual_icoor(
            v_last, last_atom, first_atom, ring_atoms[1], ring_atoms[2] if n > 2 else ring_atoms[0]
        )
        self._add_virtual_icoor(
            v_first,
            first_atom,
            last_atom,
            ring_atoms[-2],
            ring_atoms[-3] if n > 2 else ring_atoms[-1],
        )
        self.log.info(
            f"Ring {ring_idx}: added virtual shadows {v_last.strip()}/{v_first.strip()}. "
            "LOWEST_RING_CONFORMER and LOW_RING_CONFORMERS not auto-set "
            "(requires Cremer-Pople analysis) — set manually if needed."
        )

    def _make_virtual_name(self, real_name: str, suffix: str = "") -> str:
        """Generate a virtual atom name from a real atom name.

        Prepends 'V' to the stripped real name and pads to 4 chars.
        Raises ValueError if the result exceeds 4 characters.
        """
        stripped = real_name.strip()
        vname = f"V{stripped}{suffix}"
        if len(vname) > 4:
            raise ValueError(
                f"Virtual name '{vname}' is too long (>4 chars) for real atom '{stripped}'"
            )
        return vname.ljust(4)

    def _add_virtual_icoor(
        self, virt_name: str, shadow_name: str, bonded_to: str, ref2: str, ref3: str
    ) -> None:
        """Add an ICOOR_INTERNAL entry for a virtual atom.

        Uses 3D coordinates from the mol when available, otherwise sensible defaults
        matching typical ring geometry (phi=60, theta=70, distance=1.44).
        """
        phi, theta, distance = self._measure_virtual_icoor(shadow_name, bonded_to, ref2, ref3)
        self.ICOOR_INTERNAL.append(
            dict(
                child=virt_name,
                phi=phi,
                theta=theta,
                distance=distance,
                parent=bonded_to,
                second_parent=ref2,
                third_parent=ref3,
            )
        )

    def _measure_virtual_icoor(
        self, shadow_name: str, bonded_to: str, ref2: str, ref3: str
    ) -> tuple[float, float, float]:
        """Compute (phi, theta, distance) from mol 3D coords, or return defaults."""
        defaults = (60.0, 70.0, 1.44)
        if self.mol is None:
            return defaults
        try:
            name_to_idx = self._build_name_to_idx_map()
            indices = [name_to_idx.get(n.strip()) for n in (shadow_name, bonded_to, ref2, ref3)]
            if any(idx is None for idx in indices):
                return defaults
            si, bi, r2i, r3i = indices  # type: ignore[misc]
            conf = self.mol.GetConformer()
            distance = Chem.rdMolTransforms.GetBondLength(conf, si, bi)  # type: ignore[attr-defined]
            theta = 180 - Chem.rdMolTransforms.GetAngleDeg(conf, si, bi, r2i)  # type: ignore[attr-defined]
            tor = Chem.rdMolTransforms.GetDihedralDeg(conf, si, bi, r2i, r3i)  # type: ignore[attr-defined]
            phi = tor if str(tor) != "nan" else 60.0
            return (phi, theta, distance)
        except Exception:
            return defaults

    def polish_mol(self, resi: int = 1, chain: str = "X") -> None:
        """
        The mol may be inconsistent in its PDBResidueInfo

        :return:
        """
        resn = self.NAME[:3]
        self.log.debug(f"Making all PDBResidueInfo resn={resn} resi={resi} chain={chain}")
        for atom in self.mol.GetAtoms():
            info = atom.GetPDBResidueInfo()
            if info is not None:
                pass
            if atom.HasProp("molFileAlias"):
                info = Chem.AtomPDBResidueInfo()
                info.SetName(atom.GetProp("molFileAlias"))
            else:
                info = Chem.AtomPDBResidueInfo()
                info.SetName("")  # this is the default but illegal as per dejavu
            info.SetResidueName(resn)
            info.SetResidueNumber(resi)
            info.SetChainId(chain)
            # info.SetSegmentNumber(segi) commented out because it is an int not string
            if self.is_aminoacid():
                info.SetIsHeteroAtom(False)
            else:
                info.SetIsHeteroAtom(True)
        # is the name unique?
        self.rename_repeated_atoms()

    def rename_repeated_atoms(self) -> None:
        dejavu: set[str] = {""}  # set of already seen stripped names plus the illegal '' name
        to_be_renamed: list[Chem.Atom] = []
        # partial repetition of ``_fix_atom_names``
        for atom in self.mol.GetAtoms():
            name: str = atom.GetPDBResidueInfo().GetName()
            if not isinstance(name, str) or name.strip() in dejavu:
                to_be_renamed.append(atom)
            dejavu.add(name.strip())
        # rename
        for atom in to_be_renamed:
            info: Chem.AtomPDBResidueInfo = atom.GetPDBResidueInfo()
            element: str = atom.GetSymbol()
            for i in range(1, 100):
                name = f"{element: >2}{i: <2}"
                if name.strip() in dejavu:
                    continue
                print(name)
                info.SetName(name)
                break
            else:
                raise RuntimeError(f"Could not find a unique name for {atom}")

    # ============ VIRT ================================================================================================
    def _parse_w_virtuals(self) -> None:
        """
        Add 1 or 2 vanadium (virtual) atoms, 0.1A away
        :return:
        """
        target_atom_count = 3  # this is fixed due to icoor. But the code works for more if add_icoor part is corrected.
        nonvirtuals = self.mol.GetNumAtoms()
        virtuals = target_atom_count - nonvirtuals
        if virtuals <= 0:
            raise ValueError("Human called this when there are 3+ atoms.")
        elif virtuals == target_atom_count:
            raise ValueError("There are no atoms")
        # 1 or more virtuals
        anchor = sorted(self.mol.GetAtoms(), key=lambda atom: atom.GetAtomicNum(), reverse=True)[0]
        anchor_idx = anchor.GetIdx()  # either 0 or 1
        mol = Chem.RWMol(self.mol)
        for i in range(virtuals):
            virtual = Chem.Atom("V")
            virtual.SetDoubleProp("_GasteigerCharge", 0.0)
            virtual.SetProp("_rType", "VIRT")
            mol.AddAtom(virtual)
            virtual_idx = mol.GetNumAtoms() - 1
            virtual = mol.GetAtomWithIdx(virtual_idx)  # the PDBResidue info does not get set?
            mol.AddBond(anchor_idx, virtual_idx, Chem.BondType.ZERO)
            virtual.SetMonomerInfo(
                Chem.AtomPDBResidueInfo(
                    atomName=f"V{i + 1}",
                    serialNumber=virtual_idx,
                    residueName=self.NAME,
                    isHeteroAtom=not self.is_aminoacid(),
                )
            )
        # get the coords off the original (without vanadium "virtual" atoms)
        coordMap = {i: mol.GetConformer().GetAtomPosition(i) for i in range(self.mol.GetNumAtoms())}
        AllChem.EmbedMolecule(mol, coordMap=coordMap)  # type: ignore[attr-defined]
        for i in range(virtuals):
            AllChem.SetBondLength(mol.GetConformer(), anchor_idx, i + nonvirtuals, 0.1)  # type: ignore[attr-defined]
        self.mol = mol.GetMol()
        self._undescribed = deque(self.mol.GetAtoms())
        self.ordered_atoms = []
        atoms = mol.GetAtoms()
        self._add_icoor([atoms[0], atoms[0], atoms[1], atoms[2]])
        self._add_icoor([atoms[1], atoms[0], atoms[1], atoms[2]])
        self._add_icoor([atoms[2], atoms[1], atoms[0], atoms[2]])
        self.ordered_atoms = self.mol.GetAtoms()
        # no idea why but the atoms in ordered_atoms get altered and this segfaults
        # calling any getter results in a segfault

    def write_conformers(self, filename: str) -> None:
        """
        Adds all conformers to a PDB file and
        adds `PDB_ROTAMERS`.

        NB. Whereas the conformers for covalents need to be generated with DummyMasker
        this file needs to be outside of that context manager to preserve the virtual atoms.
        """
        assert self.mol.GetNumConformers() > 1, (
            "Molecule needs multiple conformers for `PDB_ROTAMERS"
        )
        with open(filename, "w") as fh:
            for c in range(self.mol.GetNumConformers()):
                fh.write(Chem.MolToPDBBlock(self.mol, confId=c, flavor=1 + 2 + 8))
        self.PDB_ROTAMERS.append(filename)  # noqa
