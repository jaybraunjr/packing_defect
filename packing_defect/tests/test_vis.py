# tests/test_leaflet_visualizer.py
import matplotlib
matplotlib.use("Agg")  
import os
import numpy as np
import pytest
import MDAnalysis as mda
from MDAnalysis import Universe
from packing_defect.core.visualization import LeafletVisualizer




def make_universe(n_atoms=10, box=(100, 100, 100)):
    """Create a dummy universe with atoms spread in z-dimension."""
    u = Universe.empty(n_atoms,
                       n_residues=n_atoms,
                       atom_resindex=np.arange(n_atoms),
                       residue_segindex=[0]*n_atoms,
                       trajectory=True)
    u.add_TopologyAttr("name", ["C"]*n_atoms)
    u.add_TopologyAttr("resname", ["LIP"]*n_atoms)
    u.add_TopologyAttr("resid", np.arange(n_atoms)+1)
    # segid not needed

    # random coordinates along z
    coords = np.zeros((n_atoms, 3))
    coords[:, 0] = np.linspace(0, box[0]-1, n_atoms)
    coords[:, 1] = np.linspace(0, box[1]-1, n_atoms)
    coords[:, 2] = np.linspace(0, box[2]-1, n_atoms)
    u.atoms.positions = coords
    u.trajectory.ts.dimensions = np.array([*box, 90, 90, 90])
    return u



@pytest.fixture
def viz(tmp_path):
    return LeafletVisualizer(base_dir=str(tmp_path),
                             lipid="TEST",
                             leaflet="up",
                             defect_color="cyan",
                             protein_color="black")


def test_z_mid_box_method(viz):
    u = make_universe()
    assert np.isclose(viz._z_mid(u), 50.0)


def test_z_mid_data_method(viz):
    u = make_universe()
    viz.method = "data"
    mid = viz._z_mid(u)
    assert 0 < mid < 100


def test_leaflet_mask_up_and_dw(viz):
    u = make_universe()
    mask_up = viz._leaflet_mask(u)
    assert mask_up.sum() > 0
    viz.leaflet = "dw"
    mask_dw = viz._leaflet_mask(u)
    assert mask_dw.sum() > 0
    assert (mask_up & mask_dw).sum() == 0


def test_defect_xy_returns_coordinates(viz):
    u = make_universe()
    x, y = viz._defect_xy(u)
    assert len(x) == len(y)
    assert len(x) > 0


def test_protein_xy_empty_if_no_protein(viz):
    u = make_universe()
    x, y = viz._protein_xy(u)
    assert x.size == 0
    assert y.size == 0


def test_plot_and_render(tmp_path, viz):
    u = make_universe()
    fig, ax = viz._plot(u, frame=0, show=False)
    assert fig is not None and ax is not None

    # simulate .gro file to test render_frame
    gro_path = tmp_path / "TEST" / "TEST_frame_0.gro"
    os.makedirs(gro_path.parent, exist_ok=True)
    u.atoms.write(str(gro_path))

    out_file = viz.render_frame(0)
    assert out_file and os.path.exists(out_file)
