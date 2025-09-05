"""Deprecated shim for writer functions.

Re-exports write_defect_coordinates from packing_defect.io.writer for backward compatibility,
and exposes write_combined_gro from utils for existing monkeypatch targets in tests/code.
"""

from packing_defect.utils import write_combined_gro  # noqa: F401


def write_defect_coordinates(protein, defect_atoms, dims, path):
    import os
    os.makedirs(os.path.dirname(path), exist_ok=True)
    write_combined_gro(protein, defect_atoms, dims, path)


__all__ = ["write_defect_coordinates", "write_combined_gro"]

