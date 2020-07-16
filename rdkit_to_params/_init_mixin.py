########################################################################################################################
__doc__ = \
    """
The main class here is ``_ParamsInitMixin``, which simply adds the __init__ method...
    """

__author__ = "Matteo Ferla. [Github](https://github.com/matteoferla)"
__email__ = "matteo.ferla@gmail.com"
__date__ = "10 July 2020 A.D."
__license__ = "MIT"
__version__ = "1.1.0"
__citation__ = "None."

########################################################################################################################


from .entries import Entries

class _ParamsInitMixin:

    # This is solely for my sanity/notebook and does not and should not be used.
    ordering = ['NAME', 'IO_STRING', 'TYPE', 'AA', 'ROTAMER_AA', '#' 'ATOM', 'ATOM_ALIAS', 'BOND', 'CUT_BOND',
                'CHI', 'CONNECT', 'ADD_RING', 'PROPERTIES', 'METAL_BINDING_ATOMS', 'FIRST_SIDECHAIN_ATOM',
                'RAMA_PREPRO_FILENAME', 'ACT_COORD_ATOMS',
                'NBR_ATOM', 'NBR_RADIUS', 'ICOOR_INTERNAL', 'PDB_ROTAMERS']

    def __init__(self):
        self.AA = Entries.from_name('AA')  # for ligands
        self.TYPE = Entries.from_name('TYPE') #POLYMER or LIGAND
        self.ROTAMER_AA = Entries.from_name('ROTAMER_AA')
        self.IO_STRING = Entries.from_name('IO_STRING')
        self.comments = Entries.from_name('#')
        self.ATOM = Entries.from_name('ATOM')
        self.ATOM_ALIAS = Entries.from_name('ATOM_ALIAS')
        self.BOND = Entries.from_name('BOND')
        self.CUT_BOND = Entries.from_name('CUT_BOND')
        self.CONNECT = Entries.from_name('CONNECT')
        self.CHI = Entries.from_name('CHI')
        self.NBR_ATOM = Entries.from_name('NBR_ATOM')
        self.NBR_RADIUS = Entries.from_name('NBR_RADIUS')
        self.ICOOR_INTERNAL = Entries.from_name('ICOOR_INTERNAL')
        self.ADD_RING = Entries.from_name('ADD_RING')
        self.PROPERTIES = Entries.from_name('PROPERTIES')
        self.FIRST_SIDECHAIN_ATOM = Entries.from_name('FIRST_SIDECHAIN_ATOM')
        self.RAMA_PREPRO_FILENAME = Entries.from_name('RAMA_PREPRO_FILENAME')
        self.METAL_BINDING_ATOMS = Entries.from_name('METAL_BINDING_ATOMS')
        self.ACT_COORD_ATOMS = Entries.from_name('ACT_COORD_ATOMS')
        self.PDB_ROTAMERS = Entries.from_name('PDB_ROTAMERS') #TODO it would be nice to check this file exists.
        ## RDKit route specific
        self.mol = None
        self.generic = False
        self._rtype = []

    @property
    def fields(self):
        """
        This operates under the assumption that the user may have added extra uppercase entries!

        :return:
        """
        # NAME is a dynamic property.
        return ['NAME'] + [k for k in self.__dict__ if k.upper() == k]