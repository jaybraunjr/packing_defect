import numpy as np
import pytest
from MDAnalysis import Universe
from packing_defect.utils import (
    apply_pbc,
    compute_pairwise_distances,
    validate_defect_thresholds,
    write_combined_gro,
    initialize_empty_defect_universe,
)
from packing_defect.core.grid import _dfs, _make_graph, DefectGrid


def test_apply_pbc_wraps_positions():
    box = np.array([10.0, 10.0, 10.0])
    positions = np.array([[11.0, -1.0, 5.0]])
    result = apply_pbc(positions, box)
    assert np.all(result >= 0)
    assert np.all(result < box)


def test_pairwise_distance_basic():
    a = np.array([[0, 0, 0]])
    b = np.array([[3, 4, 0]])
    result = compute_pairwise_distances(a, b)
    assert np.isclose(result[0, 0], 5.0)


def test_validate_defect_thresholds():
    defect_types = ['a', 'b']
    with pytest.raises(ValueError):
        validate_defect_thresholds(defect_types, {'a': 1})


def test_write_combined_gro(tmp_path):
    protein = Universe.empty(2, trajectory=True)
    defect = Universe.empty(1, trajectory=True)
    protein.atoms.positions = np.array([[0, 0, 0], [1, 0, 0]])
    defect.atoms.positions = np.array([[0.5, 0.5, 0]])
    out = tmp_path / 'out.gro'
    dims = np.array([10, 10, 10, 90, 90, 90])
    write_combined_gro(protein.atoms, defect.atoms, dims, str(out))
    assert out.exists()
    u = Universe(str(out))
    assert len(u.atoms) == 3


def test_initialize_empty_defect_universe():
    dims = np.array([[10, 10, 10, 90, 90, 90],
                     [10, 10, 10, 90, 90, 90]])
    u = initialize_empty_defect_universe(5, 2, dims, dt=1.0)
    assert len(u.atoms) == 5
    assert u.trajectory.n_frames == 2
    assert np.allclose(u.trajectory[0].dimensions, dims[0])


def test_defectgrid_cluster_sizes():
    grid = DefectGrid((5, 5), dx=1, dy=1, hz=0)
    grid.update(0, 0, 1, 0, 1, 'up')
    grid.update(1, 0, 1, 0, 1, 'up')
    grid.update(3, 4, 1, 0, 1, 'up')
    sizes = sorted(grid.cluster_sizes('up'))
    assert sizes == [1, 2]


def test_make_graph_and_dfs():
    matrix = np.array([[1, 0], [1, 1]])
    graph = _make_graph(matrix)
    comp = _dfs(graph, 0)
    assert comp == {0, 2, 3}
