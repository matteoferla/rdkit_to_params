########################################################################################################################
__doc__ = """
    THis the core conversion. the method ``convert_mol`` does something you'd never have guessed.
    There are probably loads of mistakes here ---or more correctly
    in ``_RDKitPrepMixin()._add_rtypes`` or ``_RDKitPrepMixin()._add_genrtypes``.
    """

########################################################################################################################

import random
from collections import deque, namedtuple
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
    CUT_BOND: Entries
    BOND: Entries
    NBR_ATOM: Entries
    NBR_RADIUS: Entries
    PDB_ROTAMERS: Entries

    greekification = True  # controls whether to change the atom names to greek for amino acids for CB and beyond.

    def convert_mol(self, pcharge_prop_name: str = "_GasteigerCharge") -> None:
        """
        This method does the actual conversion to params entries

        :return:
        """
        if not self.mol.GetConformers():
            self.log.warning("The molecule passed has no conformers... Adding.")
            # add_conformer does also the gasteiger charges :/
            with DummyMasker(self.mol):
                AllChem.EmbedMolecule(self.mol)  # type: ignore[attr-defined]
                AllChem.MMFFOptimizeMolecule(self.mol)  # type: ignore[attr-defined]
        # NAME
        if self.mol.HasProp("_Name"):
            title = self.mol.GetProp("_Name").strip()
        else:
            title = "Unnamed molecule"
        if len(title) == 3:
            self.NAME = title
        else:
            self.comments.append(title)
            self.NAME = self._get_resn_from_PDBInfo()
        self.log.debug(f"{title} is being converted (`.convert_mol`)")
        if len(Chem.GetMolFrags(self.mol)) > 1:
            self.log.warning(f"{title} is split in {len(Chem.GetMolFrags(self.mol))}")
            self.mol = Chem.GetMolFrags(self.mol, asMols=True)[0]
        # SMILES
        self.comments.append(Chem.MolToSmiles(self.mol))
        if len(self.TYPE) == 0:
            self.TYPE.append("LIGAND")
        # remove all previous entries.
        for attr_name in Entries.choices:
            # do not blank these:
            if attr_name in [
                "#",
                "TYPE",
                "comment",
                "IO_STRING",
                "ROTAMER_AA",
                "AA",
                "PROPERTIES",
                "VARIANT",
                "<UNKNOWN>",
            ]:
                continue
            elif not hasattr(self, attr_name):
                continue
            else:
                getattr(self, attr_name).clear()
        # correct for ions
        if self.mol.GetNumAtoms() >= 3:
            # ICOOR
            self._parse_icoors()  # fills self.ordered_atoms
            # ATOM & CONNECT
            self._parse_atoms(pcharge_prop_name)
            # CHI
            self._parse_rotatables()
        else:
            self._parse_w_virtuals()  # adds virt atom to the mol
            self._parse_atoms(pcharge_prop_name)
        # BOND
        self._parse_bonds()
        # NBR
        self._find_centroid()

    ## ============= internal coordinates ==============================================================================

    def _parse_icoors(self) -> None:
        self.log.debug("Filling ICOOR")
        self._undescribed = deque(self.mol.GetAtoms())
        self.ordered_atoms: list[Chem.Atom] = []
        self._rotation_count = 0
        if self.is_aminoacid():
            # you have to start with N.
            N_atom: Chem.Atom = self.get_atom_by_name("N")
            CA_atom: Chem.Atom = self.get_atom_by_name("CA")
            C_atom: Chem.Atom = self.get_atom_by_name("C")
            O_atom: Chem.Atom = self.get_atom_by_name("O")
            UPPER_atom: Chem.Atom = self.get_atom_by_name("UPPER")
            LOWER_atom: Chem.Atom = self.get_atom_by_name("LOWER")
            self._add_icoor([N_atom, N_atom, CA_atom, C_atom])
            self._add_icoor([CA_atom, N_atom, CA_atom, C_atom])
            self._add_icoor([C_atom, CA_atom, N_atom, C_atom])
            self._add_icoor([UPPER_atom, C_atom, CA_atom, N_atom])
            self._add_icoor([O_atom, C_atom, CA_atom, UPPER_atom])
            self._add_icoor([LOWER_atom, N_atom, CA_atom, C_atom])
        else:
            # cycle through until a good root is found!
            while self._undescribed:
                self._prevent_overcycling()
                atom: Chem.Atom = self._undescribed[0]
                if atom.GetSymbol() == "*" or atom.GetProp("_rType") == "VIRT":
                    self.log.debug("Dummy or Virtual not allowed in first three lines...")
                    self._undescribed.rotate(-1)
                    self._rotation_count += 1
                    continue  # This cannot be in the first three lines, whereas the R group is first in a canonical SMILES.
                else:
                    # I assume in some cases Virtual atoms are called for... but for now, no!
                    for parent_atom in self._get_unseen_neighbors(atom, [], True):
                        siblings_atoms = self._get_unseen_neighbors(parent_atom, [atom], True)
                        if len(siblings_atoms) != 0:
                            break
                    else:
                        self.log.debug(
                            f"First row: No siblings for {self._get_PDBInfo_atomname(atom)}!"
                        )
                        self._undescribed.rotate(-1)
                        self._rotation_count += 1
                        continue
                    grandparent_atom: Chem.Atom = siblings_atoms[0]
                    # first line
                    self._add_icoor([atom, atom, parent_atom, grandparent_atom])
                    # Second line
                    self._add_icoor([parent_atom, atom, parent_atom, grandparent_atom])
                    # Third line
                    self._add_icoor([grandparent_atom, parent_atom, atom, grandparent_atom])
                    break
        self.log.debug("Parsing forth line onwards")
        # not first three lines
        while self._undescribed:
            self._prevent_overcycling()
            atom = self._undescribed[0]
            d_idx = [d.GetIdx() for d in self.ordered_atoms]
            d_neighs = [
                n for n in self._get_unseen_neighbors(atom, [], False) if n.GetIdx() in d_idx
            ]
            if len(d_neighs) == 0:
                self.log.debug(f"Other row: No parent for {self._get_PDBInfo_atomname(atom)}!")
                self._undescribed.rotate(-1)
                self._rotation_count += 1
                continue
            parent_atom = d_neighs[0]
            d_neighs = [
                n for n in self._get_unseen_neighbors(parent_atom, [atom], False) if n.GetIdx() in d_idx
            ]
            if len(d_neighs) == 0:
                self._undescribed.rotate(-1)
                self.log.debug(f"Other row: No siblings for  {self._get_PDBInfo_atomname(atom)}!")
                self._rotation_count += 1
                continue
            elif len(d_neighs) == 1:
                sibling_atom: Chem.Atom = d_neighs[0]
                d_neighs = [
                    n
                    for n in self._get_unseen_neighbors(sibling_atom, [parent_atom, atom], False)
                    if n.GetIdx() in d_idx
                ]
                if len(d_neighs) == 0:
                    self._undescribed.rotate(-1)
                    self.log.debug(f"Other row: No cousins for  {self._get_PDBInfo_atomname(atom)}")
                    self._rotation_count += 1
                    continue
                else:
                    atoms = [atom, parent_atom, sibling_atom, d_neighs[0]]
            else:
                atoms = [atom, parent_atom, d_neighs[0], d_neighs[1]]
            self._add_icoor(atoms)
            self._rotation_count -= 1

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
                    return len(
                        Chem.GetShortestPath(self.mol, root.GetIdx(), atom.GetIdx())
                    )
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
        """
        I am not sure if this is right...
        I do not know if a C(=O)-C=O is rotatable (cis-trans).

        :return:
        """
        self.log.debug("Filling CHI")

        def is_single(a, b):
            bond = self.mol.GetBondBetweenAtoms(a.GetIdx(), b.GetIdx())
            return bond.GetBondType() == Chem.BondType.SINGLE

        cuts = [
            tuple(sorted([ring_set[0], ring_set[-1]]))
            for ring_set in self.mol.GetRingInfo().AtomRings()
        ]

        def is_cut(a, b):
            pair = tuple(sorted([a.GetIdx(), b.GetIdx()]))
            return pair in cuts

        ordered_indices = [a.GetIdx() for a in self.ordered_atoms]

        def get_place(atom):
            return ordered_indices.index(atom.GetIdx())

        # mark ring atoms.
        for ring_set in self.mol.GetRingInfo().AtomRings():
            for i in ring_set:
                self.mol.GetAtomWithIdx(i).SetBoolProp("_isRing", True)
        for atom in self.ordered_atoms:
            if atom.HasProp("_isRing"):
                continue
            elif atom.GetSymbol() == "H":
                continue  # PROTON_CHI 3 SAMPLES 2 0 180 EXTRA 1 20 thing...
            elif atom.GetSymbol() == "*":
                continue
            # should there be a `is_single` case?
            for neighbor in self._get_unseen_neighbors(atom, []):
                if neighbor.HasProp("_isRing"):
                    continue
                elif get_place(atom) > get_place(neighbor):
                    continue
                elif is_cut(atom, neighbor):  # redundant as _isRing in effect.
                    continue  # dont cross a cut!
                elif is_single(atom, neighbor):
                    for grandneighbor in self._get_unseen_neighbors(neighbor, [atom]):
                        if grandneighbor.HasProp("_isRing"):
                            continue
                        elif is_cut(neighbor, grandneighbor):
                            continue
                        elif is_single(grandneighbor, neighbor):
                            ggneighbors = self._get_unseen_neighbors(
                                grandneighbor, [atom, neighbor]
                            )
                            if ggneighbors:  # 3-4 bond does not need to be single, nor does it matter which is the 4th atom
                                self.CHI.append(
                                    dict(
                                        index=len(self.CHI) + 1,
                                        first=self._get_PDBInfo_atomname(atom),
                                        second=self._get_PDBInfo_atomname(neighbor),
                                        third=self._get_PDBInfo_atomname(grandneighbor),
                                        fourth=self._get_PDBInfo_atomname(ggneighbors[0]),
                                    )
                                )

    def _get_unseen_neighbors(
        self, atom: Chem.Atom, seen: list[Chem.Atom], nondummy: bool = True
    ) -> list[Chem.Atom]:
        neighbors_list = list(atom.GetNeighbors())
        if nondummy:
            neighbors_list = [neighbor for neighbor in neighbors_list if neighbor.GetSymbol() != "*"]
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
