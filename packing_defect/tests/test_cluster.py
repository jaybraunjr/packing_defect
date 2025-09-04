import sys, os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from pathlib import Path

import numpy as np

from packing_defect.core.cluster import DefectClustering


def test_cluster_sizes():
    # 5x5 mask with a cluster of 3 and an isolated 1, away from edges
    mask = np.zeros((5, 5), dtype=int)
    mask[2, 2] = 1
    mask[2, 3] = 1
    mask[3, 3] = 1  # diagonal connectivity
    mask[0, 2] = 1  # isolated elsewhere

    sizes = sorted(DefectClustering.cluster_sizes_from_mask(mask))
    assert sizes == [1, 3]


def test_hist_file(tmp_path: Path):
    # Construct two masks: one cluster of size 2 and one of size 1
    m1 = np.zeros((4, 4), dtype=int)
    m1[1, 1] = 1
    m1[1, 2] = 1
    m2 = np.zeros((4, 4), dtype=int)
    m2[3, 0] = 1  # isolated

    out = tmp_path / "hist.dat"
    DefectClustering.defect_size([m1, m2], nbins=5, bin_max=5, fname=str(out), prob=True)
    assert out.exists()
    txt = out.read_text().strip().splitlines()
    # Expect nbins-1 lines with two columns each (histogram bins)
    assert len(txt) == 4
    cols = txt[0].split()
    assert len(cols) == 2


def test_hist_empty(tmp_path: Path):
    out = tmp_path / "empty.dat"
    DefectClustering.defect_size([np.zeros((3, 3), dtype=int)], nbins=4, bin_max=4, fname=str(out), prob=True)
    # No file should be written when there are zero counts
    assert not out.exists()
