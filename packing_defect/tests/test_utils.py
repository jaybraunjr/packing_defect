import numpy as np
import os
import pytest
from MDAnalysis import Universe

from packing_defect.utils import (
    apply_pbc,
    calculate_bounding_box,
    initialize_grid,
    compute_pairwise_distances,
    _dfs,
    _make_graph,
    filter_defects_by_distance,
    get_defect_coordinates,
    write_combined_gro
)

def test_apply_pbc_wraps_positions():
    box = np.array([10.0, 10.0, 10.0])
    positions = np.array([[11.0, -1.0, 5.0]])
    result = apply_pbc(positions, box)
    assert np.all(result >= 0)
    assert np.all(result < box)

def test_calculate_bounding_box():
    u = Universe.empty(3,trajectory=True)
    u.atoms.positions = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    bbox = calculate_bounding_box(u.atoms, padding=1.0)

    assert bbox['x_min'] < 1
    assert bbox['x_max'] > 7
    assert bbox['y_min'] < 2
    assert bbox['y_max'] > 8

def test_initialize_grid_structure():
    box = [5.0, 5.0, 10.0]
    grid = initialize_grid(box, dx=1, dy=1, hz=10.0)

    assert 'xx' in grid and 'yy' in grid
    assert grid['xx'].shape == grid['yy'].shape
    assert np.all(grid['z_up'] == 10.0)

def test_pairwise_distance_basic():
    a = np.array([[0, 0, 0]])
    b = np.array([[3, 4, 0]])
    result = compute_pairwise_distances(a, b)
    assert np.isclose(result[0, 0], 5.0)

def test_dfs_finds_connected():
    graph = {0: {1}, 1: {0, 2}, 2: {1}, 3: set()}
    visited = _dfs(graph, 0)
    assert visited == {0, 1, 2}

def test_make_graph_from_matrix():
    matrix = np.array([[0, 1], [1, 1]])
    graph = _make_graph(matrix)
    assert len(graph) > 0
    assert all(isinstance(v, set) for v in graph.values())

def test_filter_defects_by_distance_removes_close_points():
    defects = np.array([[0, 0, 0], [10, 10, 10]])
    positions = np.array([[0, 0, 1]])
    filtered = filter_defects_by_distance(defects, positions, threshold=2.0)

    assert len(filtered) == 1
    assert np.allclose(filtered[0], [10, 10, 10])

def test_get_defect_coordinates_works():
    grid = np.array([[0, 1], [1, 0]])
    xs, ys = get_defect_coordinates(grid, 1)
    assert set(zip(xs, ys)) == {(0, 1), (1, 0)}


