__name__ = 'CFPy'
__author__ = 'Torsten Noffz, Max G. Rudolph'

try:
    from importlib.metadata import version, PackageNotFoundError
    __version__ = version("CFPy_TUD")
except PackageNotFoundError:
    __version__ = "unknown"

from .cfp import *
from .utils import *