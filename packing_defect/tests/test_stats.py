import numpy as np
from packing_defect.core import stats

def test_cluster_positions_and_sizes_basic():
    pts = {(0,0), (0,1), (5,5)}
    clusters = stats.cluster_positions(pts)
    sizes = stats.cluster_sizes(pts)
    assert sorted(sizes) == sorted([2, 1])
    assert sum(len(c) for c in clusters) == 3

import numpy as np
import MDAnalysis as mda
from MDAnalysis.coordinates.memory import MemoryReader
from packing_defect.core import stats

def test_defect_distro_and_global_stats(tmp_path):
    def write_gro(coords, path):
        n = len(coords)
        # Build a tiny Universe with one frame and valid topology
        u = mda.Universe.empty(
            n_atoms=n, n_residues=n, n_segments=1, trajectory=True,
            atom_resindex=np.arange(n), residue_segindex=np.zeros(n, dtype=int)
        )
        u.add_TopologyAttr("name", [f"A{i}" for i in range(n)])
        u.add_TopologyAttr("resname", [f"R{i}" for i in range(n)])
        u.add_TopologyAttr("resid", np.arange(1, n + 1))
        u.trajectory = MemoryReader(coords[np.newaxis, :, :])
        u.trajectory.ts.dimensions = np.array([20., 20., 20., 90., 90., 90.])
        with mda.Writer(str(path), n_atoms=n) as W:
            W.write(u.atoms)

    frames = [
        np.array([[0., 0., 0.], [0., 1., 0.], [10., 10., 0.]]),
        np.array([[10., 10., 0.], [11., 10., 0.]])
    ]
    for i, xyz in enumerate(frames, 1):
        write_gro(xyz, tmp_path / f"frame_{i}.gro")

    distros = stats.defect_distribution(str(tmp_path))
    assert set(distros) == {"frame_1", "frame_2"}
    assert distros["frame_1"][2] == 1 and distros["frame_1"][1] == 1

    per_size = stats.defect_distro_stats(distros)
    assert 1 in per_size and 2 in per_size

    gmean, gstd = stats.global_avg_std_defect_size(distros)
    assert gmean > 0
