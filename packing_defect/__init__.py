"""
packing_defect
membrane packing defect analysis
"""

from importlib.metadata import version

# Expose submodules so autodoc / autosummary can import them
from . import run_radius
from . import run_defect
from . import utils
from . import data

__all__ = [
    "run_radius",
    "run_defect",
    "utils",
    "data",
]

try:
    __version__ = version("packing_defect")
except Exception:
    __version__ = "0.0.0"
