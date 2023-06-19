"""eFISHent."""

__version__ = "0.0.5"

import logging
import warnings

# Silence some troublemakers
warnings.filterwarnings("ignore")
logging.getLogger("luigi-interface").setLevel(level=logging.CRITICAL)
logging.getLogger("matplotlib.font_manager").setLevel(level=logging.CRITICAL)

from . import alignment
from . import analyze
from . import basic_filtering
from . import cleanup
from . import cli
from . import config
from . import constants
from . import generate_probes
from . import kmers
from . import optimization
from . import prepare_sequence
from . import secondary_structure
from . import util
