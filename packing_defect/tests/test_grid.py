import sys, os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
import numpy as np

from packing_defect.core.grid import DefectGrid


def test_grid_update():
    # 10x10 box with 1x1 bins → 10x10 grid
    g = DefectGrid(box_xy=(10, 10), dx=1.0, dy=1.0, hz=5.0)

    # Place two points in 'up' leaflet with different z; latest with higher z wins
    g.update(x=1.0, y=1.0, z=6.0, r=0.4, code=2, leaflet="up")
    g.update(x=1.0, y=1.0, z=7.0, r=0.4, code=3, leaflet="up")
    assert g.zdepth['up'].max() >= 7.0
    # Cell should hold the later/higher-z code
    assert (g.grid['up'] > 0).any()

    # Down leaflet prefers smaller z
    g.update(x=2.0, y=2.0, z=1.0, r=0.4, code=4, leaflet="dw")
    g.update(x=2.0, y=2.0, z=2.0, r=0.4, code=5, leaflet="dw")
    # Still should keep the first update (z=1.0)
    mask = g.get_binary_mask('dw', 4)
    assert mask.sum() >= 1

    # Coordinates for a given code return matching centers
    x_coords, y_coords = g.get_coordinates('up', 3)
    assert x_coords.shape == y_coords.shape


def test_grid_cluster_sizes():
    g = DefectGrid(box_xy=(5, 5), dx=1.0, dy=1.0, hz=2.5)
    # Paint adjacent cells with nonzero codes so they form a cluster of size >= 2
    g.update(1.0, 1.0, 6.0, 0.6, 1, 'up')
    g.update(1.9, 1.0, 6.1, 0.6, 2, 'up')
    sizes = g.cluster_sizes('up')
    assert any(s >= 2 for s in sizes)
