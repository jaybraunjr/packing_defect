import sys, os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from collections import Counter

import pytest

import packing_defect.core.stats as stats_mod


def test_stats_counts(tmp_path, monkeypatch):
    # Create a fake directory with pretend .gro filenames
    d = tmp_path / "gro"
    d.mkdir()
    (d / "f0.gro").write_text("x")
    (d / "f1.gro").write_text("x")

    # Predefine positions per file
    positions_map = {
        str(d / "f0.gro"): {(0, 0), (0, 1), (2, 2)},   # cluster sizes: 2 and 1
        str(d / "f1.gro"): {(5, 5), (6, 6), (7, 7)},   # diagonal connects with 8-neighbors → cluster of 3
    }

    def fake_open_gro(path):
        return positions_map[str(path)]

    monkeypatch.setattr(stats_mod, "open_gro", fake_open_gro, raising=True)

    # count_defect_types with cutoff 2
    counts = stats_mod.count_defect_types(str(d), cutoff=2)
    assert counts == [("f0", 1), ("f1", 1)]

    summary = stats_mod.return_stats(str(d), cutoff=2)
    assert pytest.approx(summary['mean']) == 1.0
    assert summary['median'] == 1

    distro = stats_mod.defect_distribution(str(d))
    assert distro["f0"] == Counter({2: 1, 1: 1})
    assert distro["f1"] == Counter({3: 1})

    per_size = stats_mod.defect_distro_stats(distro)
    assert per_size[1]['mean'] == 1
    assert per_size[2]['mean'] == 1
    assert per_size[3]['mean'] == 1

    frame_stats = stats_mod.average_std_per_frame(distro)
    assert set(frame_stats.keys()) == {"f0", "f1"}

    gmean, gstd = stats_mod.global_avg_std_defect_size(distro)
    assert gmean > 0
