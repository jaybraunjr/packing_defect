import numpy as np
import pytest
from packing_defect.core.grid import DefectGrid

def test_defectgrid_init_shapes():
    g = DefectGrid(box_xy=(10.0, 8.0), dx=1.0, dy=2.0)
    assert g.xx.shape == g.yy.shape
    assert g.grid["up"].shape == (g.nx, g.ny)
    assert g.grid["dw"].shape == (g.nx, g.ny)
    assert g.zdepth["up"].shape == (g.nx, g.ny)
    assert g.zdepth["dw"].shape == (g.nx, g.ny)

def test_update_sets_code_and_depth_up():
    g = DefectGrid(box_xy=(5, 5), dx=1.0, dy=1.0)
    g.update(x=2.0, y=2.0, z=3.5, r=0.2, code=9, leaflet="up")
    i = int(round(2.0 / g.dx))
    j = int(round(2.0 / g.dy))
    assert g.grid["up"][i, j] == 9
    assert g.zdepth["up"][i, j] == pytest.approx(3.5)

def test_update_respects_leaflet_depth_rule():
    g = DefectGrid(box_xy=(5, 5), dx=1, dy=1)
    i = int(round(2.0 / g.dx))
    j = int(round(2.0 / g.dy))
    g.update(2.0, 2.0, z=1.0, r=0.2, code=1, leaflet="up")
    g.update(2.0, 2.0, z=0.5, r=0.2, code=2, leaflet="up")
    assert g.grid["up"][i, j] == 1  # higher z wins for up

    g = DefectGrid(box_xy=(5, 5), dx=1, dy=1)
    g.update(2.0, 2.0, z=-1.0, r=0.2, code=3, leaflet="dw")
    g.update(2.0, 2.0, z=-0.5, r=0.2, code=4, leaflet="dw")
    assert g.grid["dw"][i, j] == 3  # lower z wins for dw

def test_get_binary_mask_and_coordinates():
    g = DefectGrid(box_xy=(5, 5), dx=1, dy=1)
    g.update(1.0, 1.0, z=1.0, r=0.49, code=7, leaflet="up")
    mask = g.get_binary_mask("up", 7)
    assert mask.dtype == int
    xs, ys = g.get_coordinates("up", 7)
    assert len(xs) == len(ys)
    assert len(xs) >= 1

def test_cluster_sizes_simple_patch():
    g = DefectGrid(box_xy=(5, 5), dx=1, dy=1)
    # fill a 2x2 patch with code 1
    for x in (1.0, 2.0):
        for y in (1.0, 2.0):
            g.update(x, y, z=1.0, r=0.49, code=1, leaflet="up")
    sizes = g.cluster_sizes("up")
    assert sum(sizes) >= 4
