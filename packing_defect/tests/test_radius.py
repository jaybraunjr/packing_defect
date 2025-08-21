import sys, os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
import shutil
import os
import tempfile
import MDAnalysis as mda
import numpy as np
import pytest

from packing_defect.core.analyzers.radius import RadiusDefectAnalyzer
from packing_defect.run_utils import run_analysis


# Use a non-interactive backend for matplotlib inside tests
os.environ.setdefault("MPLBACKEND", "Agg")


def _write_tiny_gro(path, nprot=2, ndef=4, box=(60.0, 60.0, 60.0), zshift=0.0):
    """
    Minimal .gro writer:
      - Header
      - Atom count
      - nprot protein atoms, then ndef defect atoms
      - Box line
    Coordinates are simple and parseable by MDAnalysis.
    """
    with open(path, "w", encoding="utf-8") as f:
        f.write("tiny gro\n")
        f.write(f"{nprot + ndef}\n")
        # Protein atoms at small z
        for i in range(nprot):
            # Columns: resid, resname, atom, atomid, x y z
            f.write(f"{1:5d}{'PROT':>5s}{('C'+str(i+1)):>5s}{(i+1):5d}"
                    f"{(10.0+i):8.3f}{10.0:8.3f}{(10.0+zshift):8.3f}\n")
        # Defect atoms spread in z so we get some up/down relative to hz
        for j in range(ndef):
            z = 20.0 + j + zshift
            f.write(f"{2:5d}{'DEF':>5s}{('D'+str(j+1)):>5s}{(nprot+j+1):5d}"
                    f"{(20.0+j):8.3f}{20.0:8.3f}{z:8.3f}\n")
        f.write(f"{box[0]:10.5f}{box[1]:10.5f}{box[2]:10.5f}\n")


@pytest.fixture
def input_dirs(tmp_path):
    """
    Create synthetic input directory structure with two frames per lipid type.
    Layout:
      <tmp>/input/
        PLacyl/PLacyl_frame_0.gro, PLacyl_frame_1.gro
        TGglyc/TGglyc_frame_0.gro, TGglyc_frame_1.gro
        TGacyl/TGacyl_frame_0.gro, TGacyl_frame_1.gro
    """
    base = tmp_path / "input"
    base.mkdir(parents=True, exist_ok=True)
    lipids = ["PLacyl", "TGglyc", "TGacyl"]
    for lipid in lipids:
        d = base / lipid
        d.mkdir(exist_ok=True)
        _write_tiny_gro(d / f"{lipid}_frame_0.gro", nprot=2, ndef=5, zshift=0.0)
        _write_tiny_gro(d / f"{lipid}_frame_1.gro", nprot=2, ndef=6, zshift=0.5)
    return str(base), lipids


@pytest.fixture
def dummy_universe():
    """
    The RadiusDefectAnalyzer doesn't actually use the passed universe for its
    per-file analysis, but BaseDefectAnalyzer stores it. Provide a tiny one.
    """
    # 1-atom empty universe with no trajectory
    u = mda.Universe.empty(1, atom_resindex=[0], residue_segindex=[0], trajectory=False)
    u.add_TopologyAttr("name", ["X"])
    u.add_TopologyAttr("resname", ["DUM"])
    u.add_TopologyAttr("resid", [1])
    return u


def test_radius_pipeline_end_to_end(tmp_path, input_dirs, dummy_universe, monkeypatch):
    base_dir, lipids = input_dirs
    out_dir = tmp_path / "out"
    out_dir.mkdir(exist_ok=True)

    def _fake_renumber(self, inp, out):
        shutil.copy(inp, out)

    monkeypatch.setattr(RadiusDefectAnalyzer, "_renumber_gro", _fake_renumber, raising=True)

    analyzer = RadiusDefectAnalyzer(
        universe=dummy_universe,
        base_directory=base_dir,
        output_dir=str(out_dir),
        lipid_types=lipids,
        frame_start=0,
        frame_end=1,
        protein_atom_count=2,
        apply_protein_cutoff=False,  # keep it simple for the test
        cutoff_distance=1.5,
    )

    analyzer.run()
    analyzer.plot()
    out_txt = analyzer.save_results("summary.txt", analyzer.results)

    # Assertions
    # 1) Results populated for all lipids with keys
    assert set(analyzer.results.keys()) == set(lipids)
    for lipid in lipids:
        assert set(analyzer.results[lipid].keys()) == {"up", "down", "combined"}
        # There should be at least one measured defect size overall
        assert isinstance(analyzer.results[lipid]["combined"], list)

    # 2) Output artifacts exist
    assert os.path.exists(out_txt)
    assert os.path.exists(out_dir / "top_bottom_defects.png")
    assert os.path.exists(out_dir / "combined_defects.png")

    # 3) Renumbered files were produced for each lipid and frame
    for lipid in lipids:
        ldir = out_dir / lipid
        assert ldir.exists()
        # we expect "renumbered_<lipid>_corrected_frame_X.gro"
        for i in (0, 1):
            rn = ldir / f"renumbered_{lipid}_corrected_frame_{i}.gro"
            assert rn.exists()
            # MDAnalysis can read them
            mda.Universe(str(rn))