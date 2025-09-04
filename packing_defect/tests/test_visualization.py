import sys, os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from pathlib import Path

import packing_defect.core.visualization as visz_mod


def test_lv_render_frame(tmp_path, monkeypatch):
    # Create expected input file path
    base = tmp_path / "base"
    lipid = "PLacyl"
    (base / lipid).mkdir(parents=True)
    gro = base / lipid / f"{lipid}_frame_3.gro"
    gro.write_text("x")

    # Monkeypatch Universe to avoid parsing .gro; and _plot to just return the savepath
    monkeypatch.setattr(visz_mod.mda, "Universe", lambda path: object())
    def fake_plot(self, u, frame, show=True, savepath=None):
        Path(savepath).parent.mkdir(parents=True, exist_ok=True)
        Path(savepath).write_text("png")
        return savepath

    monkeypatch.setattr(visz_mod.LeafletVisualizer, "_plot", fake_plot)

    lv = visz_mod.LeafletVisualizer(base_dir=str(base), lipid=lipid, leaflet="up")
    out = lv.render_frame(3)
    assert out and Path(out).exists()


def test_lv_render_frames(tmp_path, monkeypatch):
    base = tmp_path / "base2"
    lipid = "TGglyc"
    (base / lipid).mkdir(parents=True)
    for i in range(2):
        (base / lipid / f"{lipid}_frame_{i}.gro").write_text("x")

    monkeypatch.setattr(visz_mod.mda, "Universe", lambda path: object())
    def fake_plot(self, u, frame, show=True, savepath=None):
        Path(savepath).parent.mkdir(parents=True, exist_ok=True)
        Path(savepath).write_text("png")
        return savepath

    monkeypatch.setattr(visz_mod.LeafletVisualizer, "_plot", fake_plot)

    lv = visz_mod.LeafletVisualizer(base_dir=str(base), lipid=lipid, leaflet="dw")
    outs = lv.render_frames(0, 2)
    assert len(outs) == 2


def test_lv_plot(tmp_path, monkeypatch):
    base = tmp_path / "base3"
    lipid = "X"
    (base / lipid).mkdir(parents=True)
    (base / lipid / f"{lipid}_frame_1.gro").write_text("x")

    monkeypatch.setattr(visz_mod.mda, "Universe", lambda path: object())
    def fake_plot(self, u, frame, show=True, savepath=None):
        return "FIG", "AX"

    monkeypatch.setattr(visz_mod.LeafletVisualizer, "_plot", fake_plot)

    lv = visz_mod.LeafletVisualizer(base_dir=str(base), lipid=lipid)
    fig, ax = lv.plot_frame(1, show=False)
    assert fig == "FIG" and ax == "AX"
