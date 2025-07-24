
import numpy as np
import os
import pytest
from MDAnalysis import Universe
from packing_defect.defect import PackingDefect2, PackingDefect2Sequential
from packing_defect.utils import apply_pbc, initialize_grid, compute_pairwise_distances




def test_apply_pbc():
    box = np.array([10.0, 10.0, 10.0])
    pos = np.array([[11.0, -1.0, 5.0]])
    corrected = apply_pbc(pos, box)
    assert np.all(corrected >= 0) and np.all(corrected <= box), "PBC application failed"

def test_initialize_grid():
    box = [10.0, 10.0, 10.0]
    grid = initialize_grid(box, dx=1, dy=1)
    assert grid['xx'].shape == grid['yy'].shape, "Grid initialization failed"

def test_compute_pairwise_distances():
    p1 = np.array([[0, 0, 0]])
    p2 = np.array([[3, 4, 0]])
    dist = compute_pairwise_distances(p1, p2)
    assert np.isclose(dist[0][0], 5.0), "Distance calculation is incorrect"

def test_classification_default():
    radii_path = os.path.join(os.path.dirname(__file__), '..', 'data', 'radii.json')
    radii_path = os.path.abspath(radii_path)
    pd = PackingDefect2(radii_path)  # Replace with mock or test file
    assert pd.default_classify('TRIO', 'O11') == 2
    assert pd.default_classify('TRIO', 'XYZ') == 3
    assert pd.default_classify('POPC', 'C216') == 1
    assert pd.default_classify('POPC', 'ZZZ') == -1

def test_read_top_parses_correctly():
    from packing_defect.defect import PackingDefect2
    radii_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'data', 'radii.json'))
    test_top = os.path.abspath(os.path.join(os.path.dirname(__file__), 'test.top'))

    pd = PackingDefect2(radii_path)
    result = pd.read_top('POPC', test_top)

    assert 'C216' in result
    assert isinstance(result['C216'], list)
    assert len(result['C216']) == 2  # [radius, acyl_class]
