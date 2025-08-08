import numpy as np
import MDAnalysis as mda
from MDAnalysis import Universe
from packing_defect.utils import (
    apply_pbc,
    compute_pairwise_distances,
    validate_defect_thresholds,
    initialize_empty_defect_universe,
    write_combined_gro
)

def test_apply_pbc_wraps_into_box():
    box = np.array([10.0, 10.0, 10.0])
    pos = np.array([[11.0, -1.0, 5.0]])
    wrapped = apply_pbc(pos, box)
    assert np.all(wrapped >= 0) and np.all(wrapped < box)

def test_pairwise_distances_basic():
    a = np.array([[0,0,0]])
    b = np.array([[3,4,0]])
    d = compute_pairwise_distances(a, b)
    assert d.shape == (1,1) and d[0,0] == 5.0

def test_validate_defect_thresholds_ok():
    validate_defect_thresholds(["A","B"], {"A":1, "B":2})

def test_validate_defect_thresholds_missing():
    try:
        validate_defect_thresholds(["A","B"], {"A":1})
        raised = False
    except ValueError:
        raised = True
    assert raised

def test_initialize_empty_defect_universe_and_write(tmp_path):
    n_atoms, nframes = 4, 3
    dims = [np.array([20,20,20,90,90,90]) for _ in range(nframes)]
    u = initialize_empty_defect_universe(n_atoms, nframes, dims, dt=1.0)
    assert len(u.atoms) == n_atoms
    assert u.trajectory.n_frames == nframes
    # Create two Universes to merge/write
    prot = mda.Universe.empty(2, trajectory=False)
    prot.add_TopologyAttr("name", ["P1","P2"])
    prot.atoms.positions = np.array([[1,1,1],[2,2,2]])
    dfx = mda.Universe.empty(2, trajectory=False)
    dfx.add_TopologyAttr("name", ["D1","D2"])
    dfx.atoms.positions = np.array([[3,3,3],[4,4,4]])
    out = tmp_path / "combo.gro"
    write_combined_gro(prot.atoms, dfx.atoms, np.array([20,20,20,90,90,90]), str(out))
    assert out.exists()

from MDAnalysis.coordinates.memory import MemoryReader

def test_initialize_empty_defect_universe_and_write(tmp_path):
    n_atoms, nframes = 4, 3
    dims = [np.array([20,20,20,90,90,90]) for _ in range(nframes)]
    u = initialize_empty_defect_universe(n_atoms, nframes, dims, dt=1.0)
    assert len(u.atoms) == n_atoms
    assert u.trajectory.n_frames == nframes

    # Create two Universes with 1-frame memory trajectories
    prot = mda.Universe.empty(2, n_residues=2, n_segments=1, trajectory=True,
                              atom_resindex=np.arange(2), residue_segindex=np.zeros(2, int))
    prot.add_TopologyAttr("name", ["P1", "P2"])
    prot.trajectory = MemoryReader(np.zeros((1, 2, 3)))
    prot.atoms.positions = np.array([[1, 1, 1], [2, 2, 2]])

    dfx = mda.Universe.empty(2, n_residues=2, n_segments=1, trajectory=True,
                             atom_resindex=np.arange(2), residue_segindex=np.zeros(2, int))
    dfx.add_TopologyAttr("name", ["D1", "D2"])
    dfx.trajectory = MemoryReader(np.zeros((1, 2, 3)))
    dfx.atoms.positions = np.array([[3, 3, 3], [4, 4, 4]])

    out = tmp_path / "combo.gro"
    write_combined_gro(prot.atoms, dfx.atoms, np.array([20,20,20,90,90,90]), str(out))
    assert out.exists()
