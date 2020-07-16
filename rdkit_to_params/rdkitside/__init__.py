from ._rdkit_inits import _RDKitInitMixin
from ._rdkit_rename import _RDKitRenameMixin
from ._rdkit_convert import _RDKitCovertMixin
from ._rdkit_prep import _RDKitPrepMixin

########################################################################################################################
__doc__ = \
    """
    The functionality is too big for a single class so was split into functional units.
    
    * **init** contains the various entry points
    * **rename** deals with atom renaming
    * **prep** labels atoms etc. for conversion
    * **convert** converts the mol into params entries
    """

__author__ = "Matteo Ferla. [Github](https://github.com/matteoferla)"
__email__ = "matteo.ferla@gmail.com"
__date__ = "10 July 2020 A.D."
__license__ = "MIT"
__version__ = "1.1.0"
__citation__ = "None."

########################################################################################################################

class _RDKitMixin(_RDKitInitMixin):
    # this is so that parts can change happily
    pass