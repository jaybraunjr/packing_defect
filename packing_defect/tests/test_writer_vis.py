import sys, os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")

import packing_defect.core.vis as vis_mod
from packing_defect.core.writer import write_defect_coordinates


def test_plot_png(tmp_path: Path):
    # Create three tiny .dat files with two columns
    paths = []
    for name in ("PLacyl.dat", "TGglyc.dat", "TGacyl.dat"):
        p = tmp_path / name
        np.savetxt(p, np.column_stack((np.arange(1, 4), np.array([0.1, 0.01, 0.001]))))
        paths.append(p)

    outpng = tmp_path / "plot.png"
    vis_mod.plot_defect_data([str(paths[0]), str(paths[1]), str(paths[2])],
                             labels=["PL", "TGg", "TGa"],
                             colors=["b", "r", "g"],
                             markers=["o", "s", "^"],
                             title="Test",
                             output_path=str(outpng))
    assert outpng.exists() and outpng.stat().st_size > 0


def test_writer_calls(tmp_path, monkeypatch):
    calls = {}

    def fake_write_combined_gro(protein, defect_atoms, dimensions, filepath):
        calls['args'] = (protein, defect_atoms, dimensions, filepath)
        Path(filepath).write_text("stub")

    # Patch the utility used by writer
    monkeypatch.setattr("packing_defect.core.writer.write_combined_gro", fake_write_combined_gro)

    out = tmp_path / "out" / "frame_0.gro"
    write_defect_coordinates(protein="P", defect_atoms="D", dims=[1, 2, 3, 90, 90, 90], path=str(out))
    assert out.exists()
    assert calls['args'][3].endswith("frame_0.gro")
