from rdkit_to_params.rdkitside._rdkit_convert import _RDKitCovertMixin as _RDKitCovertMixin
from rdkit_to_params.rdkitside._rdkit_inits import _RDKitInitMixin
from rdkit_to_params.rdkitside._rdkit_prep import _RDKitPrepMixin as _RDKitPrepMixin
from rdkit_to_params.rdkitside._rdkit_rename import _RDKitRenameMixin as _RDKitRenameMixin
from rdkit_to_params.rdkitside.utilities import DummyMasker, neutralize
from rdkit_to_params.version import (
    __author__,
    __citation__,
    __date__,
    __email__,
    __license__,
    __version__,
)

########################################################################################################################
__doc__ = """
    The functionality is too big for a single class so was split into functional units.

    * **init** contains the various entry points
    * **rename** deals with atom renaming
    * **prep** labels atoms etc. for conversion
    * **convert** converts the mol into params entries
    """


########################################################################################################################

__all__ = [
    "__author__",
    "__email__",
    "__date__",
    "__license__",
    "__version__",
    "__citation__",
    "_RDKitMixin",
    "neutralize",
    "DummyMasker",
]


class _RDKitMixin(_RDKitInitMixin):
    # this is so that parts can change happily
    pass
