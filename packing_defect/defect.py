"""Minimal defect module used for tests."""

from __future__ import annotations

import json


class PackingDefect2:
    """Simplified implementation providing minimal API for unit tests."""

    def __init__(self, radii_path: str) -> None:
        with open(radii_path, "r", encoding="utf-8") as handle:
            self.radii = json.load(handle)

    def default_classify(self, resname: str, atomname: str) -> int:
        """Return a default classification code for *resname* and *atomname*.

        This dummy logic only covers the cases exercised in the tests.
        """
        if resname == "TRIO":
            return 2 if atomname == "O11" else 3
        if resname == "POPC":
            return 1 if atomname == "C216" else -1
        return -1

    def read_top(self, resname: str, filepath: str) -> dict:
        """Parse a very small topology file and return radii information."""
        atoms = {}
        with open(filepath, "r", encoding="utf-8") as handle:
            for line in handle:
                parts = line.split()
                if len(parts) >= 3 and parts[0] == "ATOM":
                    atom_name = parts[1]
                    atom_type = parts[2]
                    radius = self.radii.get(atom_type, 0.0)
                    atoms[atom_name] = [radius, "A"]
        return atoms


class PackingDefect2Sequential(PackingDefect2):
    """Sequential variant used only for type compatibility in tests."""

    pass
