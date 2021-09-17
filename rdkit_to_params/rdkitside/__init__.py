from ._rdkit_inits import _RDKitInitMixin
from ._rdkit_rename import _RDKitRenameMixin
from ._rdkit_convert import _RDKitCovertMixin
from ._rdkit_prep import _RDKitPrepMixin
from .utilities import *

########################################################################################################################
__doc__ = \
    """
    The functionality is too big for a single class so was split into functional units.
    
    * **init** contains the various entry points
    * **rename** deals with atom renaming
    * **prep** labels atoms etc. for conversion
    * **convert** converts the mol into params entries
    """

from ..version import *

########################################################################################################################

class _RDKitMixin(_RDKitInitMixin):
    # this is so that parts can change happily
    pass

