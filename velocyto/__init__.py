from .constants import *
from .gene_info import GeneInfo
from .read import Read
from .feature import Feature
from .transcript_model import TranscriptModel
from .segment_match import SegmentMatch
from .indexes import FeatureIndex, TransciptsIndex
from .molitem import Molitem
from .logic import *
from .counter import ExInCounter
from .metadata import MetadataCollection, Metadata
from .neighbors import BalancedKNN, convolve_by_sparse_weights
from .estimation import fit_slope, _fit1_slope, clusters_stats
from .serialization import dump_hdf5, load_hdf5
from .analysis import VelocytoLoom, scatter_viz, ixs_thatsort_a2b, load_velocyto_hdf5
from ._version import __version__


# Protect users from a nasty bug in Anaconda
# See https://github.com/velocyto-team/velocyto.py/issues/104
# and https://github.com/ContinuumIO/anaconda-issues/issues/10089

import math
import numpy as np

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

std_check = np.arange(1000000).std()
expected = 288675.1345946685

if not math.isclose(std_check, expected):
    raise RuntimeError(MKL_BUG_ERROR_MSG)
