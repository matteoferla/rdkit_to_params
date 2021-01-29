#!/usr/bin/env python3

__version__ = '0.1'
__doc__ = """
Command line interface to RDKit to params.
NB. For most applications, the website is better.
""".strip()

from typing import *
import argparse

# ======================================================================================================================

from .. import Params


def smiles():
    parser = argparse.ArgumentParser(description=self.__doc__)
    parser.add_argument('smiles', type=str, help='SMILES string to convert. If unfamiliar, see Wikipedia.')
    parser.add_argument('-n', '--name', type=str, help='Three letter code.')
    parser.add_argument('-o', '--output', type=str, help='Output params file')
    parser.add_argument('-o', '--pdb', type=str, help='Output pdb file')
    args = parser.parse_args()
    p = Params.from_smiles(args.smiles,
                           name=args.name if args.name else 'LIG'
                           )
    # params
    if args.output:
        p.dump(args.output)
    else:
        print(p.dumps())
    # pdb
    if args.pdb:
        p.dump_pdb(args.pdb)


if __name__ == '__main__':
    smiles()
