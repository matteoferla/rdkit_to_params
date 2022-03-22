########################################################################################################################
__doc__ = \
    """
The main class here is ``_ParamsIoMixin``, which adds load and dump to ``_ParamsInitMixin``.
    """

import warnings

from .version import *

########################################################################################################################


import re
from .entries import Entries
from ._init_mixin import _ParamsInitMixin


class _ParamsIoMixin(_ParamsInitMixin):

    @classmethod
    def load(cls, filename: str, skip_unknown: bool=True):
        if hasattr(cls, '__name__'):
            self = cls()
        else:
            self = cls
        with open(filename, 'r') as w:
            for line in w:
                self._parse_line(line, skip_unknown)
        return self

    @classmethod
    def loads(cls, text: str, skip_unknown: bool=True):
        if hasattr(cls, '__name__'):
            self = cls()
        else:
            self = cls
        for line in text.split('\n'):
            self._parse_line(line, skip_unknown)
        return self

    def _parse_line(self, line, skip_unknown: bool=True):
        sline = line.strip()
        if not sline:
            return
        elif re.match('#', sline):
            header, body = re.match('(#+)\s?(.*)', sline).groups()
            self.comments.append(body)
        elif re.match('[A-Z]', sline):
            header, body = re.match('([_\w]+) (.*)', sline).groups()
            if '#' in body:
                body, comment = re.match('(.*?) ?#(.*)', body).groups()
                self.comments.append(comment)
            if 'CONNECT' in header:
                self.CONNECT.append(dict(atom_name=body,
                                         connect_type=header,
                                         index=len(self.CONNECT) + 1))
            elif header == 'NAME':
                self.NAME = body
            elif header in ('BOND', 'BOND_TYPE'):
                self.BOND.append(body)
            elif header in Entries.choices:
                getattr(self, header).append(body)
            else:
                msg = f'Entry {header} is for a not-yet-encoded entry type.'
                if not skip_unknown:
                    raise ValueError(msg +'\nAdditing to `.OTHERS`')
                warnings.warn(msg)
                self.OTHERS.append(self.OTHERS.entry_cls(header, body))
        else:
            self.comments.append(sline)

    def dumps(self, html=False) -> str:
        assert self.NBR_ATOM and self.NBR_ATOM[0], 'Undeclared NBR_ATOM entry'
        assert self.NBR_RADIUS and self.NBR_RADIUS[0], 'Undeclared NBR_RADIUS entry'
        if not html:
            lines = [f'NAME {self.NAME}', str(self.comments)]
        else:
            lines = [f'NAME {self.NAME}', self.comments._repr_html_()]
        for entries in (self.IO_STRING, self.TYPE, self.AA, self.ROTAMER_AA, self.ATOM, self.ATOM_ALIAS, self.BOND,
                        self.CUT_BOND, self.CHARGE, self.ADD_RING, self.PROPERTIES, self.VARIANT, self.METAL_BINDING_ATOMS,
                        self.FIRST_SIDECHAIN_ATOM, self.BACKBONE_AA, self.RAMA_PREPRO_FILENAME,
                        self.ACT_COORD_ATOMS, self.CHI,
                        self.CONNECT, self.NBR_ATOM, self.NBR_RADIUS, self.MAINCHAIN_ATOMS,
                        self.ICOOR_INTERNAL, self.PDB_ROTAMERS,
                        self.OTHERS):
            if not entries:
                continue
            elif not html:
                lines.append(str(entries))
            else:
                lines.append(entries._repr_html_())
        if not html:
            return '\n'.join(lines)
        else:
            return '<br>'.join(lines)

    def dump(self, filename: str) -> None:
        with open(filename, 'w') as w:
            w.write(self.dumps())

    def _repr_html_(self):
        return f'<p>{self.dumps(html=True)}</p>'

