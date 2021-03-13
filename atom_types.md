## Atom Types

In addition to the partial charges, Rosetta AtomTypes are defined.
These basically control hydrogen-bonding and atom size.
In some cases the assignment may have gone wrong and it is worth checking.

These are added to the rdkit.Chem.Mol Property via (`atom.SetProp('_rType', 'xxx')`)
and to the `params.ATOM[n].rtype` string.

Example:

```jupyterpython
from rdkit import Chem
from rdkit_to_params import Params

for smiles in ('CC[NH3+]','CC[NH2]','C=C[NH2]','CC[OH]','CC[OH2+]','CC[O-]','C=CO','CC=O'):
    p = Params.from_smiles(smiles)
    i = 2
    print(smiles, p.mol.GetAtomWithIdx(i).GetHybridization().name, p.ATOM[i].rtype)
```
| SMILES | Hybridisation | AtomType |
| ---- | ---- | ---- |
| CC[NH3+] | SP3 | Nlys |
| CC[NH2] | SP3 | Npro |
| C=C[NH2] | SP2 | NH2O |
| CC[OH] | SP3 | OH |
| CC[OH2+] | SP3 | Oet3 |
| CC[O-] | SP3 | Oet3 |
| C=CO | SP2 | OH |
| CC=O | SP2 | OOC |

To see what the different atom type properties one could inspect a residue in pyrosetta

```jupyterpython
at = pose.residue(1).atom_type(1)
at.get_all_properties()
```

Or look at the file `rosetta/main/database/chemical/atom_type_sets/fa_standard/atom_properties.txt`
but briefly these are the accepted atom types.

## AtomType options

|    | AtomType   | Element   | Properties                             |
|---:|:-----------|:----------|:---------------------------------------|
|  0 | aroC       | C         | AROMATIC ORBITALS                      |
|  1 | Ntrp       | N         | DONOR AROMATIC ORBITALS                |
|  2 | Nhis       | N         | ACCEPTOR AROMATIC RING_HYBRID ORBITALS |
|  3 | NtrR       | N         | DONOR AROMATIC ORBITALS                |
|  4 | NH2O       | N         | DONOR                                  |
|  5 | Narg       | N         | DONOR ORBITALS                         |
|  6 | OH         | O         | ACCEPTOR SP3_HYBRID DONOR ORBITALS     |
|  7 | OW         | O         | ACCEPTOR SP3_HYBRID DONOR ORBITALS     |
|  8 | ONH2       | O         | ACCEPTOR SP2_HYBRID ORBITALS           |
|  9 | OOC        | O         | ACCEPTOR SP2_HYBRID ORBITALS           |
| 10 | Oaro       | O         | ACCEPTOR AROMATIC RING_HYBRID ORBITALS |
| 11 | Oet2       | O         | ACCEPTOR SP2_HYBRID ORBITALS           |
| 12 | Oet3       | O         | ACCEPTOR SP3_HYBRID DONOR ORBITALS     |
| 13 | Nbb        | N         | DONOR                                  |
| 14 | OCbb       | O         | ACCEPTOR SP2_HYBRID ORBITALS           |
| 15 | Hpol       | H         | POLAR_HYDROGEN                         |
| 16 | HS         | H         | POLAR_HYDROGEN                         |
| 17 | HNbb       | H         | POLAR_HYDROGEN                         |
| 18 | Hwat       | H         | POLAR_HYDROGEN                         |
| 19 | Owat       | O         | ACCEPTOR DONOR SP3_HYBRID              |
| 20 | HOH        | O         | ACCEPTOR DONOR SP3_HYBRID              |
| 21 | F          | F         | SP3_HYBRID                             |
| 22 | Cl         | CL        | SP3_HYBRID                             |
| 23 | Br         | BR        | SP3_HYBRID                             |
| 24 | I          | I         | SP3_HYBRID                             |
| 25 | #Zn2p      | ZN        | SP3_HYBRID                             |
| 26 | Zn2p       | ZN        | SP3_HYBRID                             |
| 27 | Co2p       | CO        | SP3_HYBRID                             |
| 28 | Cu2p       | CU        | SP3_HYBRID                             |
| 29 | Fe2p       | FE        | SP3_HYBRID                             |
| 30 | Fe3p       | FE        | SP3_HYBRID                             |
| 31 | Mg2p       | MG        | SP3_HYBRID                             |
| 32 | Ca2p       | CA        | SP3_HYBRID                             |
| 33 | Pha        | P         | SP3_HYBRID                             |
| 34 | OPha       | O         | ACCEPTOR SP3_HYBRID                    |
| 35 | OHha       | O         | ACCEPTOR DONOR SP3_HYBRID              |
| 36 | Hha        | H         | POLAR_HYDROGEN                         |
| 37 | CO3        | C         | SP2_HYBRID                             |
| 38 | OC3        | O         | ACCEPTOR SP2_HYBRID                    |
| 39 | Si         | Si        | SP3_HYBRID                             |
| 40 | OSi        | O         | ACCEPTOR SP3_HYBRID                    |
| 41 | Oice       | O         | ACCEPTOR SP3_HYBRID                    |
| 42 | Hice       | H         | POLAR_HYDROGEN                         |
| 43 | Na1p       | NA        | SP3_HYBRID                             |
| 44 | K1p        | K         | SP3_HYBRID                             |
| 45 | REPL       | Z         | REPULSIVE                              |
| 46 | REPLS      | Z         | REPULSIVE                              |
| 47 | HREPS      | Z         | REPULSIVE                              |
| 48 | VIRT       | X         | VIRTUAL                                |

