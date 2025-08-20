"""
classification.py

Implements the Strategy pattern for atom classification, with a default
that matches your original logic and an optional JSON-driven rule set.
"""

from abc import ABC, abstractmethod
from typing import Dict, Tuple
import json


class ClassificationStrategy(ABC):
    """
    Abstract interface for atom classification strategies.
    """

    @abstractmethod
    def classify(self, resname: str, atom_name: str) -> int:
        """
        Given a residue name and atom name, return an integer code.
        """


class DefaultClassification(ClassificationStrategy):
    """
    Default logic (mirrors your original `default_classify`):
    - For non-TRIO residues: tail atoms → 1, else → -1
    - For TRIO residues: glycerol atoms → 2, else → 3
    """

    def __init__(self):
        self.tails = [f"C2{i}" for i in range(2, 23)] + \
                     [f"C3{i}" for i in range(2, 23)] + \
                     [f"H{i}{s}" for i in range(2, 23) for s in ['R', 'S', 'X', 'Y']] + \
                     ['H16Z', 'H18T', 'H91', 'H101', 'H18Z', 'H20T']
        self.TGglyc = ['O11', 'O21', 'O31', 'O12', 'O22', 'O32',
                       'C1', 'C2', 'C3', 'C11', 'C21', 'C31',
                       'HA', 'HB', 'HS', 'HX', 'HY']
        self.PL_resnames = ('POPC', 'DOPE', 'SAPI')

    def classify(self, resname: str, atom_name: str) -> int:
        if resname in self.PL_resnames:
            return 1 if atom_name in self.tails else -1
        if resname == 'TRIO':
            return 2 if atom_name in self.TGglyc else 3
        return -1


class UserDictClassification(ClassificationStrategy):
    """
    Load classification codes from a JSON file of the form::

        {
          "RES1": {"ATOM1": "heads", "ATOM2": "tails", ...},
          "RES2": { ... }
        }

    Labels ("heads", "tails") are automatically mapped to integer codes.
    """

    def __init__(self, rules: Dict[Tuple[str, str], int], label_to_code: Dict[str, int]):
        self.rules = rules
        self.label_to_code = label_to_code

    @classmethod
    def from_json(cls, json_file: str) -> 'UserDictClassification':
        with open(json_file, 'r') as f:
            data: Dict[str, Dict[str, str]] = json.load(f)

        label_to_code: Dict[str, int] = {}
        next_code = 1
        rules: Dict[Tuple[str, str], int] = {}

        for resname, atom_map in data.items():
            for atom_name, label in atom_map.items():
                # assign a new integer code if label not seen
                if label not in label_to_code:
                    label_to_code[label] = next_code
                    next_code += 1
                rules[(resname, atom_name)] = label_to_code[label]

        return cls(rules, label_to_code)

    def classify(self, resname: str, atom_name: str) -> int:
        return self.rules.get((resname, atom_name), -1)
