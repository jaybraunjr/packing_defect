# packing_defect/core/writer.py

import os
import json
from packing_defect.utils import write_combined_gro

def write_defect_coordinates(protein, defect_atoms, dims, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    write_combined_gro(protein, defect_atoms, dims, path)

