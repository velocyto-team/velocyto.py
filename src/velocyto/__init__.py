# src/pyplier/__init__.py
from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version(__name__)
except PackageNotFoundError:  # pragma: no cover
    __version__ = "unknown"

from math import isclose

from numpy import arange

from numba.core.errors import NumbaDeprecationWarning, NumbaPendingDeprecationWarning
import warnings

warnings.simplefilter("ignore", category=NumbaDeprecationWarning)
warnings.simplefilter("ignore", category=NumbaPendingDeprecationWarning)

# from .analysis import VelocytoLoom, ixs_thatsort_a2b, load_velocyto_hdf5, scatter_viz
# from .constants import *
# from .counter import ExInCounter
# from .estimation import _fit1_slope, clusters_stats, fit_slope
# from .feature import Feature
# from .gene_info import GeneInfo
# from .indexes import FeatureIndex, TransciptsIndex
# from .logic import *
# from .metadata import Metadata, MetadataCollection
# from .molitem import Molitem
# from .neighbors import BalancedKNN, convolve_by_sparse_weights
# from .read import Read
# from .segment_match import SegmentMatch
# from .serialization import dump_hdf5, load_hdf5
# from .transcript_model import TranscriptModel

# Protect users from a nasty bug in Anaconda
# See https://github.com/velocyto-team/velocyto.py/issues/104
# and https://github.com/ContinuumIO/anaconda-issues/issues/10089


MKL_BUG_ERROR_MSG = """
Your current Python installation is affected by a critical bug in numpy and
MKL, and is going to return wrong results in velocyto and potentially other
scientific packages.

Please try updating your `numpy` version.

For more information, see
https://github.com/velocyto-team/velocyto.py/issues/104
and
https://github.com/ContinuumIO/anaconda-issues/issues/10089
"""

std_check = arange(1000000).std()
expected = 288675.1345946685

if not isclose(std_check, expected):
    raise RuntimeError(MKL_BUG_ERROR_MSG)
