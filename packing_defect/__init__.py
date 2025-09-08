"""packing_defect: membrane packing defect analysis."""

from importlib.metadata import version

# Keep package imports lightweight; avoid importing CLI modules here.
# Expose version for consumers without side effects.
try:  # pragma: no cover - falls back when not installed
    __version__ = version("packing_defect")
except Exception:
    __version__ = "0.0.0"

__all__ = ["__version__"]
