import sys, os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
import csv
from pathlib import Path

import packing_defect.core.interactions as inter_mod


def test_interact_process(tmp_path, monkeypatch):
    # Create directory structure with dummy files
    lipids = ["PLacyl", "TGglyc"]
    base = tmp_path / "in"
    base.mkdir()
    for lipid in lipids:
        d = base / lipid
        d.mkdir()
        # Filenames drive frame index parsing
        (d / f"{lipid}_frame_0.gro").write_text("x")
        (d / f"{lipid}_frame_1.gro").write_text("x")

    # Monkeypatch MDAnalysis.Universe used in this module to a no-op
    class DummyU:
        pass

    monkeypatch.setattr(inter_mod.mda, "Universe", lambda path: DummyU())

    # Monkeypatch tracker to record deterministic results without MDAnalysis
    def fake_track_frame(self, universe, frame_index):
        # One hit per frame with resid=frame_index and resname depending on lipid
        # We cannot see lipid here, so just produce one generic hit
        self.results.setdefault(frame_index, []).append((int(frame_index), "ALA"))

    monkeypatch.setattr(inter_mod.ProteinDefectInteractions, "track_frame", fake_track_frame)

    outcsv = tmp_path / "combined.csv"
    inter_mod.process_multiple([str(base / l) for l in lipids], str(outcsv), cutoff=1.0)

    # Validate CSV structure
    rows = list(csv.DictReader(open(outcsv, newline='')))
    assert rows and set(rows[0].keys()) == {"type", "frame", "resid", "resname"}
    types = {r["type"] for r in rows}
    assert types == set(lipids)
