
import sys, os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))

# tests/test_packing_analyzer.py

import os
import numpy as np
import MDAnalysis as mda
import pytest

from packing_defect.core.analyzers.packing import PackingDefectAnalyzer
from packing_defect.utils import initialize_empty_defect_universe


@pytest.fixture
def tiny_universe_one_frame():
    """
    Build a minimal universe with:
      - 2 protein atoms near z ~ 25
      - POPC:
          P atoms split across leaflets to set hz and COMs
          one tail atom "C22" in each leaflet
      - TRIO:
          one glycerol-like atom "O11" in the up leaflet
          one acyl-like atom "C31" in the down leaflet
    Box ~ 60 x 60 x 50, one frame.
    """
    n_prot = 2
    n_popc = 5   # P(up), P(dw), C22(up), C22(dw), dummy
    n_trio = 2   # O11(up), C31(dw)
    natoms = n_prot + n_popc + n_trio

    u = mda.Universe.empty(
        n_atoms=natoms,
        n_residues=natoms,
        atom_resindex=np.arange(natoms),
        residue_segindex=[0] * natoms,
        trajectory=True,
    )
    # Topology attributes
    names = []
    resnames = []
    # Protein
    for _ in range(n_prot):
        names.append("CA")
        resnames.append("ALA")
    # POPC block
    names += ["P", "P", "C22", "C22", "C22"]
    resnames += ["POPC"] * 5
    # TRIO block
    names += ["O11", "C31"]
    resnames += ["TRIO", "TRIO"]

    u.add_TopologyAttr("name", names)
    u.add_TopologyAttr("resname", resnames)
    u.add_TopologyAttr("resid", np.arange(natoms) + 1)

    # Coordinates (Å): spread in x,y; choose z to separate leaflets
    xyz = np.zeros((1, natoms, 3), dtype=float)
    # protein near midplane
    xyz[0, 0] = [10, 10, 25]
    xyz[0, 1] = [12, 10, 25.5]
    # POPC P up, P down, tails up/down, extra
    xyz[0, 2] = [20, 20, 35]  # P up
    xyz[0, 3] = [25, 25, 10]  # P down
    xyz[0, 4] = [22, 21, 34]  # C22 up
    xyz[0, 5] = [27, 26, 11]  # C22 down
    xyz[0, 6] = [21, 24, 33]  # dummy POPC tail (still code 1)
    # TRIO glycerol up, acyl down
    xyz[0, 7] = [30, 30, 36]  # O11 up
    xyz[0, 8] = [35, 31, 12]  # C31 down

    u.load_new(xyz, order="fac")
    u.trajectory[0].dt = 2.0
    u.trajectory[0].dimensions = np.array([60.0, 60.0, 50.0, 90.0, 90.0, 90.0])

    return u


@pytest.fixture
def radii_map_minimal():
    """
    Provide radii with embedded classification codes:
      PL tails → code 1
      TRIO glycerol → code 2
      TRIO acyl → code 3
    The PackingDefectAnalyzer reads radius, code from radii[resname][atomname].
    """
    return {
        "POPC": {
            "P": (1.5, -1),    # not a defect code; used to compute hz/COM
            "C22": (1.0, 1),   # PLacyl
        },
        "TRIO": {
            "O11": (1.0, 2),   # TGglyc
            "C31": (1.0, 3),   # TGacyl
        },
        # DOPE present in real runs but not required for this synthetic test
    }


def _build_analyzer(u, radii, tmp_path, leaflet="both", start=None, stop=None, stride=1):
    memb = u.select_atoms("resname POPC TRIO")
    return PackingDefectAnalyzer(
        universe=u,
        atomgroups=[memb],
        radii=radii,
        output_dir=str(tmp_path),
        leaflet=leaflet,
        defect_types=["PLacyl", "TGglyc", "TGacyl"],
        defect_thresholds={"PLacyl": 1, "TGglyc": 2, "TGacyl": 3},
        start=start,
        stop=stop,
        stride=stride,
    )


def test_defaults_and_validation(tiny_universe_one_frame, radii_map_minimal, tmp_path):
    u = tiny_universe_one_frame
    # Use default thresholds via None and ensure mapping 1,2,3 assigned in that order
    analyzer = PackingDefectAnalyzer(
        universe=u,
        atomgroups=[u.select_atoms("resname POPC TRIO")],
        radii=radii_map_minimal,
        output_dir=str(tmp_path),
        leaflet="both",
        defect_types=["PLacyl", "TGglyc", "TGacyl"],
        defect_thresholds=None,  # trigger default {PLacyl:1,TGglyc:2,TGacyl:3}
    )
    assert analyzer.defect_thresholds["PLacyl"] == 1
    assert analyzer.defect_thresholds["TGglyc"] == 2
    assert analyzer.defect_thresholds["TGacyl"] == 3


def test_run_creates_dat_and_frames(tiny_universe_one_frame, radii_map_minimal, tmp_path):
    u = tiny_universe_one_frame
    analyzer = _build_analyzer(u, radii_map_minimal, tmp_path, leaflet="both")

    analyzer.run()

    # .dat files present and non-empty
    for dname in ("PLacyl", "TGglyc", "TGacyl"):
        dat = os.path.join(tmp_path, f"{dname}.dat")
        assert os.path.exists(dat), f"Missing {dat}"
        assert os.path.getsize(dat) > 0, f"{dat} is empty"

    # Per-type GRO frames exist
    for dname in ("PLacyl", "TGglyc", "TGacyl"):
        ddir = os.path.join(tmp_path, dname)
        # at least frame 0
        gro0 = os.path.join(ddir, f"{dname}_frame_0.gro")
        assert os.path.exists(gro0), f"Missing {gro0}"
        # sanity check it can be opened
        mda.Universe(gro0)  # should not throw


def test_leaflet_up_only_reduces_down_contribution(tiny_universe_one_frame, radii_map_minimal, tmp_path):
    u = tiny_universe_one_frame
    # Run "both" and "up" so we can compare masks produced in memory
    a_both = _build_analyzer(u, radii_map_minimal, tmp_path / "both", leaflet="both")
    a_both.run()
    both_masks = {k: [m.copy() for m in v] for k, v in a_both.defect_cluster_masks.items()}

    a_up = _build_analyzer(u, radii_map_minimal, tmp_path / "up", leaflet="up")
    a_up.run()
    up_masks = a_up.defect_cluster_masks

    # For each defect type, "up" run should have zero or equal counts compared to "both" on down-leaflet masks.
    # In this implementation masks are appended per-frame as [up_mask, down_mask].
    for dname in ("PLacyl", "TGglyc", "TGacyl"):
        assert len(both_masks[dname]) == len(up_masks[dname]) == 2
        up_mask_up_both, down_mask_both = both_masks[dname]
        up_mask_up_only, down_mask_up_only = up_masks[dname]
        # Up leaflet presence should be preserved or equal
        assert up_mask_up_only.sum() <= up_mask_up_both.sum()
        # Down leaflet should be zero when leaflet="up"
        assert down_mask_up_only.sum() == 0


def test_multiple_frames_stride(tmp_path):
    """
    Smoke test that multiple frames and stride do not crash.
    We reuse initialize_empty_defect_universe to fabricate a 2-frame dummy
    'defect' universe, then embed into a parent universe with the same dims/dt.
    """
    # Make a parent 'real' universe with 1 protein atom and 2 POPC P atoms
    base = mda.Universe.empty(
        n_atoms=4,
        n_residues=4,
        atom_resindex=np.arange(4),
        residue_segindex=[0] * 4,
        trajectory=True,
    )
    base.add_TopologyAttr("name", ["CA", "P", "P", "C22"])
    base.add_TopologyAttr("resname", ["ALA", "POPC", "POPC", "POPC"])
    base.add_TopologyAttr("resid", [1, 2, 3, 4])

    coords = np.zeros((2, 4, 3))
    # frame 0 up/down at z 35 and 10
    coords[0] = [
        [10, 10, 25],
        [20, 20, 35],
        [25, 20, 10],
        [22, 21, 34],
    ]
    # frame 1 shift slightly
    coords[1] = [
        [10, 10, 25.2],
        [21, 20, 35.1],
        [26, 20, 10.2],
        [23, 22, 34.2],
    ]
    base.load_new(coords, order="fac")
    for i in range(2):
        base.trajectory[i].dt = 2.0
        base.trajectory[i].dimensions = np.array([60.0, 60.0, 50.0, 90.0, 90.0, 90.0])

    radii = {"POPC": {"P": (1.5, -1), "C22": (1.0, 1)}}

    memb = base.select_atoms("resname POPC")
    analyzer = PackingDefectAnalyzer(
        universe=base,
        atomgroups=[memb],
        radii=radii,
        output_dir=str(tmp_path / "mf"),
        leaflet="both",
        defect_types=["PLacyl"],
        defect_thresholds={"PLacyl": 1},
        start=0,
        stop=2,
        stride=1,
    )
    analyzer.run()

    # Two frames worth of GROs
    assert os.path.exists(tmp_path / "mf" / "PLacyl" / "PLacyl_frame_0.gro")
    assert os.path.exists(tmp_path / "mf" / "PLacyl" / "PLacyl_frame_1.gro")
    # Dat present
    assert os.path.getsize(tmp_path / "mf" / "PLacyl.dat") > 0
