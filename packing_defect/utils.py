
import numpy as np
import MDAnalysis as mda
from MDAnalysis import Universe
import numpy as np

def apply_pbc(positions, box):
    """Apply periodic boundary conditions to positions."""
    box_xy = np.array([box[0], box[1], 0])
    box_xyz = np.array([box[0], box[1], box[2]])
    return positions - box_xy * np.floor(positions / box_xyz)




def compute_pairwise_distances(positions1, positions2):
    """Compute pairwise distances between two sets of positions."""
    diff = positions1[:, np.newaxis, :] - positions2
    return np.sqrt(np.sum(diff**2, axis=2))




def validate_defect_thresholds(defect_types, defect_thresholds):
    for dt in defect_types:
        if dt not in defect_thresholds:
            raise ValueError(f"Missing threshold for defect type: {dt}")




def write_combined_gro(protein_atoms, defect_atoms, dimensions, filepath):
    combined = mda.Merge(protein_atoms, defect_atoms)
    combined.atoms.positions[:len(protein_atoms)] = protein_atoms.positions
    combined.atoms.positions[len(protein_atoms):] = defect_atoms.positions
    combined.trajectory.ts.dimensions = dimensions
    combined.atoms.write(filepath)



def initialize_empty_defect_universe(n_atoms, nframes, dims, dt):
    fac = np.zeros((nframes, n_atoms, 3))

    df = Universe.empty(
        n_atoms=n_atoms,
        n_residues=n_atoms,
        atom_resindex=np.arange(n_atoms),
        residue_segindex=[0] * n_atoms,
        trajectory=True,
    )
    df.add_TopologyAttr('resname', ['O'] * n_atoms)
    df.add_TopologyAttr('name', ['O'] * n_atoms)
    df.add_TopologyAttr('resid', np.arange(n_atoms) + 1)
    df.load_new(fac, order='fac')
    df.trajectory[0].dt = dt

    for i, ts in enumerate(df.trajectory):
        df.trajectory[i].dimensions = dims[i]

    return df