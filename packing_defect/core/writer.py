"""GRO writers for defect coordinates.

``write_combined_gro`` (the low-level implementation) lives in
:mod:`packing_defect.utils`; it is re-exported here so tests and callers can
monkeypatch ``packing_defect.core.writer.write_combined_gro``.
"""

import os

from packing_defect.utils import write_combined_gro  # noqa: F401


def write_defect_coordinates(protein, defect_atoms, dims, path):
    """Write protein + defect atoms to a ``.gro`` file, creating parent dirs."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    write_combined_gro(protein, defect_atoms, dims, path)


__all__ = ["write_defect_coordinates", "write_combined_gro"]
