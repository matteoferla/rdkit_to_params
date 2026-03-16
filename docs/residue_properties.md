# Rosetta Residue Properties Reference

Rosetta `.params` files contain a `PROPERTIES` line that declares one or more
residue properties. These properties control how the scoring, packing, and
protocol machinery treats a residue type.

In `rdkit_to_params` you can set properties directly:

```python
params = Params.from_smiles('CCO', name='ETH')
params.PROPERTIES.append('POLAR')
# or set several at once:
params.PROPERTIES[0].values = ['LIGAND', 'POLAR']
```

The canonical list lives in `source/src/core/chemical/residue_properties/general_properties.list`
in the Rosetta source tree. Every property below is drawn from that file.

> **Auto-set column key:**
> **Yes** = `rdkit_to_params` sets this automatically under certain conditions.
> Empty = you must add it yourself if needed.

---

## General: Residue Type

| Property | Description | Auto-set? |
|----------|-------------|-----------|
| `POLYMER` | Polymeric residue (connected via UPPER/LOWER backbone connections) | |
| `LIGAND` | Non-polymeric small-molecule residue (default for `rdkit_to_params` output) | Yes (default TYPE) |

---

## Residue Family

| Property | Description | Auto-set? |
|----------|-------------|-----------|
| `PROTEIN` | Amino acid residue (implies POLYMER). Checked by most scoring terms, rotamer libraries, and the Ramachandran potential | Yes (amino acid detection) |
| `CANONICAL_AA` | One of the 20 standard amino acids. Enables canonical rotamer libraries and sequence-profile scoring | |
| `CANONICAL_NUCLEIC` | One of the standard DNA/RNA nucleotides | |
| `DNA` | Deoxyribonucleic acid residue | |
| `RNA` | Ribonucleic acid residue | |
| `TNA` | Threose nucleic acid residue | |
| `PNA` | Peptide nucleic acid residue | |
| `PEPTOID` | N-substituted glycine (peptoid) residue | |
| `OLIGOUREA` | Oligourea residue (urea linkage instead of peptide bond; resembles beta-amino acids) | |
| `ARAMID` | Aromatic polyamide residue (generic aramid) | |
| `ORTHO_ARAMID` | Ortho-linked aramid | |
| `META_ARAMID` | Meta-linked aramid | |
| `PARA_ARAMID` | Para-linked aramid | |
| `PRE_METHYLENE_ORTHO_ARAMID` | Ortho-aramid with pre-methylene spacer | |
| `PRE_METHYLENE_META_ARAMID` | Meta-aramid with pre-methylene spacer | |
| `PRE_METHYLENE_PARA_ARAMID` | Para-aramid with pre-methylene spacer | |
| `POST_METHYLENE_ORTHO_ARAMID` | Ortho-aramid with post-methylene spacer | |
| `POST_METHYLENE_META_ARAMID` | Meta-aramid with post-methylene spacer | |
| `POST_METHYLENE_PARA_ARAMID` | Para-aramid with post-methylene spacer | |
| `PRE_METHYLENE_POST_METHYLENE_ORTHO_ARAMID` | Ortho-aramid with both pre- and post-methylene spacers | |
| `PRE_METHYLENE_POST_METHYLENE_META_ARAMID` | Meta-aramid with both pre- and post-methylene spacers | |
| `PRE_METHYLENE_POST_METHYLENE_PARA_ARAMID` | Para-aramid with both pre- and post-methylene spacers | |
| `CARBOHYDRATE` | Saccharide residue. Activates sugar-specific scoring and ring conformer handling | |
| `LIPID` | Lipid residue | |
| `TERPENE` | Terpene residue | |
| `NUCLEOTIDE_DIPHOSPHATE` | Nucleotide diphosphate (e.g. GDP, ATP) | |
| `SRI` | Artificial heteropolymer class (squaramide-based) | |
| `TRIAZOLE_LINKER` | Triazole-linked artificial heteropolymer | |

---

## Special Residue Types

| Property | Description | Auto-set? |
|----------|-------------|-----------|
| `METAL` | Metal ion residue. Used by the metal-binding constraint and scoring machinery | |
| `SURFACE` | Surface residue (for surface simulations) | |
| `WATER` | Water molecule | |
| `TP3` | TIP3P water model specifically. Used by packer residue-type setup to avoid string parsing | |
| `VIRTUALIZABLE_BY_PACKER` | During packing, this residue can be interconverted with a virtual (absent) form | |
| `SOLVENT` | Generic solvent molecule | |
| `VIRTUAL_RESIDUE` | Virtual (non-physical) residue used for coordinate frames, fold-tree jumps, etc. | |
| `VRT1` | A particular virtual residue type, used by packer setup machinery | |
| `INVERTED_VIRTUAL_RESIDUE` | Used by the symmetry machinery for mirror symmetry operations | |

---

## Variants and Terminus

| Property | Description | Auto-set? |
|----------|-------------|-----------|
| `LOWER_TERMINUS` | Polymeric residue that cannot be connected through its lower (N-terminal) connect | |
| `UPPER_TERMINUS` | Polymeric residue that cannot be connected through its upper (C-terminal) connect | |
| `BRANCH_POINT` | Polymeric residue with a non-polymeric connection, or non-polymeric residue with 3+ connections | |
| `TERMINUS` | Contains a single non-polymeric connection | Yes (`_prep_for_terminal()`) |
| `LOWERTERM_TRUNC` | Truncated lower terminus (virtual lower-terminal residue) | |
| `UPPERTERM_TRUNC` | Truncated upper terminus (virtual upper-terminal residue) | |
| `COARSE` | Coarse-grained representation (currently used for coarse RNA only) | |
| `ADDUCT` | Residue is an adduct (small-molecule covalently attached to a residue) | |
| `SC_ORBITALS` | Residue carries side-chain orbital pseudo-atoms for orbital scoring | |
| `FRAGMENT` | Fragment residue (used by GrowLigand in fragment-based ligand design) | |
| `UPPERTERM_CAP` | Upper-terminal capping residue | |
| `LOWERTERM_CAP` | Lower-terminal capping residue | |

---

## Qualitative Properties

| Property | Description | Auto-set? |
|----------|-------------|-----------|
| `POLAR` | Clearly polar side-chain (D, E, H, K, N, Q, R, S, T among canonicals). Excludes G, A, C, P and amphipathic residues like W/Y. Used by layer-design selectors and solvation scoring | |
| `HYDROPHOBIC` | Clearly hydrophobic side-chain (F, M, I, L, Y, V, W among canonicals). Also excludes G, A, C, P. The union of POLAR and HYDROPHOBIC does *not* cover all amino acids | |
| `CHARGED` | Residue carries a net formal charge | |
| `NEGATIVE_CHARGE` | Has a negative charge (includes zwitterionic side-chains like phosphoserine) | |
| `POSITIVE_CHARGE` | Has a positive charge (includes zwitterionic side-chains) | |
| `AROMATIC` | Contains an aromatic ring in the side-chain. Used by aromatic-specific scoring terms and packing heuristics | |
| `ALIPHATIC` | Side-chain contains only C and H with no aromatic rings. Met is included as an honorary aliphatic. Set: A, V, I, L, P, M | |
| `CYCLIC` | Residue contains a ring (e.g. proline, or any ligand with ring virtual-shadow atoms). Used by ring-conformer sampling and `ADD_RING` handling | Yes (ring virtual shadow generation) |
| `BETA_BRANCHED_SIDECHAIN` | Side-chain is beta-branched (e.g. Val, Ile). Affects rotamer-library selection | |
| `METALBINDING` | Residue can coordinate a metal ion. Used by the metal-site constraint generator | |
| `SIDECHAIN_THIOL` | Residue can form disulphide bonds (e.g. cysteine) | |
| `DISULFIDE_BONDED` | Residue is currently participating in a disulphide bond or non-canonical connection | |
| `ELECTROPHILE` | Side-chain can be conjugated to a nucleophile (e.g. contains an alpha-beta unsaturated carbonyl) | |
| `SIDECHAIN_AMINE` | Residue has a primary amine in the side-chain (e.g. lysine) | |
| `N_METHYLATED` | Amide proton replaced by a methyl group (mid-chain, not N-terminal) | |

---

## Membrane

| Property | Description | Auto-set? |
|----------|-------------|-----------|
| `MEMBRANE` | Membrane residue. Used by the membrane scoring framework and implicit-membrane potential | |

---

## Modifications

| Property | Description | Auto-set? |
|----------|-------------|-----------|
| `PHOSPHONATE` | Amino phosphonic acid (phosphonate replaces carboxylate) | |
| `PHOSPHONATE_UPPER` | Upper-connected phosphonate variant | |
| `ACETYLATED_NTERMINUS` | N-terminus is acetylated (ACE capping) | |
| `METHYLATED_CTERMINUS` | C-terminus is methylated (NME capping) | |
| `DIMETHYLATED_CTERMINUS` | C-terminus is dimethylated | |

---

## Amino-Acid-Specific

| Property | Description | Auto-set? |
|----------|-------------|-----------|
| `ALPHA_AA` | Alpha-amino acid (standard backbone: N-CA-C) | Yes (amino acid detection) |
| `BETA_AA` | Beta-amino acid (extra backbone carbon: N-CA-CB-C) | |
| `GAMMA_AA` | Gamma-amino acid (two extra backbone carbons) | |
| `L_AA` | L-chirality at the first chiral mainchain atom closest to LOWER_CONNECT | Yes (amino acid detection) |
| `D_AA` | D-chirality at the first chiral mainchain atom closest to LOWER_CONNECT | |
| `ACHIRAL_BACKBONE` | Backbone centre is achiral (e.g. glycine) | |
| `ACHIRAL_SIDECHAIN` | Side-chain is achiral | |
| `R_PEPTOID` | R-chirality peptoid side-chain (auto-generated from S-chirality counterparts) | |
| `S_PEPTOID` | S-chirality peptoid side-chain (only S-chirality peptoids are in the Rosetta database) | |
| `TAUTOMER` | Residue is a tautomeric form | |

---

## Nucleic-Acid-Specific

| Property | Description | Auto-set? |
|----------|-------------|-----------|
| `PURINE` | Purine base (A, G) | |
| `PYRIMIDINE` | Pyrimidine base (C, T, U) | |
| `L_RNA` | L-ribose RNA | |
| `D_RNA` | D-ribose RNA (natural chirality) | |
| `METHYLATED_NA` | Methylated nucleic acid | |

---

## Carbohydrate-Specific

These properties are set in sugar `.params` files and are generally not relevant
to small-molecule ligands. Included for completeness.

### Ring size and stereochemistry

| Property | Description |
|----------|-------------|
| `TRIOSE` | 3-carbon sugar |
| `TETROSE` | 4-carbon sugar |
| `PENTOSE` | 5-carbon sugar |
| `HEXOSE` | 6-carbon sugar |
| `HEPTOSE` | 7-carbon sugar |
| `OCTOSE` | 8-carbon sugar |
| `NONOSE` | 9-carbon sugar |
| `ALDOSE` | Aldehyde oxidation form |
| `KETOSE` | Ketone oxidation form |
| `L_SUGAR` | L-stereochemistry |
| `D_SUGAR` | D-stereochemistry |
| `OXIROSE` | 3-membered ring (oxirane) |
| `OXETOSE` | 4-membered ring (oxetane) |
| `FURANOSE` | 5-membered ring |
| `PYRANOSE` | 6-membered ring |
| `SEPTANOSE` | 7-membered ring |
| `ALPHA_SUGAR` | Alpha anomeric form |
| `BETA_SUGAR` | Beta anomeric form |
| `SIALIC_ACID` | Sialic acid residue |
| `GLYCOSIDE` | Glycosidic linkage |

### Sugar modifications

| Property | Description |
|----------|-------------|
| `ALDONIC_ACID` | Aldonic acid modification |
| `URONIC_ACID` | Uronic acid modification |
| `DEOXY_SUGAR` | Deoxygenated sugar |
| `AMINO_SUGAR` | Amino sugar (e.g. glucosamine) |
| `ACETYLAMINO_SUGAR` | N-acetylated amino sugar (e.g. GlcNAc) |
| `GLYCOLYLAMINO_SUGAR` | N-glycolylated amino sugar |
| `ACETYL_SUGAR` | O-acetylated sugar |
| `BUTYRYL_SUGAR` | O6-butyrylated sugar |
| `LACTYL_SUGAR` | O-lactylated sugar |
| `R3PRIMEHYDROXYBUTYRYLAMINO_SUGAR` | R-3'-hydroxybutyrylamino sugar |
| `PHOSPHORYLATED_SUGAR` | Phosphorylated sugar |
| `SULFATED_SUGAR` | Sulphated sugar |
| `SULFOAMINO_SUGAR` | Sulphoamino sugar |
| `THIO_SUGAR` | Thio sugar |
| `C_METHYLATED_SUGAR` | C-methylated sugar |
| `PROPARGYL_SUGAR` | Propargyl sugar |
| `FLUORO_SUGAR` | Fluorinated sugar |
| `METHYL_SUGAR` | Methylated sugar |
| `C1_MODIFIED` | Saccharide modified at position 1 |
| `C2_MODIFIED` | Saccharide modified at position 2 |
| `C3_MODIFIED` | Saccharide modified at position 3 |
| `C4_MODIFIED` | Saccharide modified at position 4 |
| `C5_MODIFIED` | Saccharide modified at position 5 |
| `C6_MODIFIED` | Saccharide modified at position 6 |
| `C7_MODIFIED` | Saccharide modified at position 7 |
| `C8_MODIFIED` | Saccharide modified at position 8 |
| `C9_MODIFIED` | Saccharide modified at position 9 |

---

## Miscellaneous

| Property | Description | Auto-set? |
|----------|-------------|-----------|
| `PHOSPHATE` | Phosphate group | |

---

## NMR Spin-Label

| Property | Description | Auto-set? |
|----------|-------------|-----------|
| `SPINLABEL` | NMR spin-label residue | |
| `DIAMAGNETIC` | Diamagnetic spin-label variant | |
| `PARAMAGNETIC` | Paramagnetic spin-label variant | |

---

## Auto-set summary

`rdkit_to_params` automatically assigns these properties under certain conditions:

| Property | Condition | Code location |
|----------|-----------|---------------|
| `LIGAND` | Default `TYPE` for all generated params | `entries.py` `TYPEEntry` default |
| `PROTEIN` | Amino acid backbone detected (N-CA-C-O pattern) | `_rdkit_prep.py:657` |
| `ALPHA_AA` | Set together with PROTEIN when amino acid detected | `_rdkit_prep.py:658` |
| `L_AA` | Set together with PROTEIN when amino acid detected | `_rdkit_prep.py:659` |
| `CYCLIC` | Ring virtual-shadow atoms generated (`ADD_RING`) | `_rdkit_convert.py:660` |
| `TERMINUS` | `make_C_terminal_cap()` or `make_N_terminal_cap()` called | `__init__.py:484` |

---

*Generated from Rosetta `general_properties.list`. Last updated: 2026-03-16.*
