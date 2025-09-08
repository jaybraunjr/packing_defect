"""CLI entrypoint for packing-defect analysis.

Provides a convenience builder for a PackingDefectAnalyzer and a
command-line interface that mirrors the common options used in notebooks.

Examples
--------
Run against a topology and trajectory and write outputs to a folder::

    python -m packing_defect.run_defect \
      --top system.gro \
      --traj traj.xtc \
      --out results \
      --leaflet both --start 0 --stop 100 --stride 1

"""

import argparse
import os
import json
from typing import Optional
import MDAnalysis as mda

from packing_defect.core.classification import DefaultClassification, UserDictClassification
from packing_defect.core.topology import TopologyReader
from packing_defect.core.analyzers.packing import PackingDefectAnalyzer
from packing_defect.run_utils import run_analysis


def build_packing_analyzer(
    top_file: str,
    traj_file: str,
    out_dir: str,
    class_json: Optional[str] = None,
    leaflet: str = "both",
    start: Optional[int] = None,
    stop: Optional[int] = None,
    stride: int = 1,
):
    """Construct a PackingDefectAnalyzer configured for a typical workflow.

    Parameters
    ----------
    top_file : str
        Topology file readable by MDAnalysis (e.g., GRO, PSF, PDB).
    traj_file : str
        Trajectory file (e.g., XTC, DCD) aligned with ``top_file``.
    out_dir : str
        Directory where results and artifacts will be written.
    class_json : str, optional
        Path to a user JSON classification mapping. If omitted, uses
        ``DefaultClassification``.
    leaflet : {"both", "up", "dw"}, optional
        Which leaflet(s) to analyze.
    start, stop, stride : int, optional
        Frame slicing parameters (as in ``trajectory[start:stop:stride]``).

    Returns
    -------
    PackingDefectAnalyzer
        Configured analyzer instance ready to ``run()``.
    """
    # load radii
    radii_file = os.path.join(os.path.dirname(__file__), "data", "radii.json")
    with open(radii_file, encoding="utf-8") as f:
        types_radii = json.load(f)

    classifier = (
        UserDictClassification.from_json(class_json)
        if class_json
        else DefaultClassification()
    )
    topo_reader = TopologyReader(types_radii, classifier.classify)

    radii = {
        resname: topo_reader.read(resname, os.path.join(os.path.dirname(__file__), "data", "top", topo))
        for resname, topo in {
            "POPC": "top_all36_lipid.rtf",
            "DOPE": "top_all36_lipid.rtf",
            "TRIO": "TRIO.rtf",
        }.items()
    }

    u = mda.Universe(top_file, traj_file)
    memb = u.select_atoms("resname " + " ".join(radii.keys()))

    return PackingDefectAnalyzer(
        universe=u,
        atomgroups=[memb],
        radii=radii,
        output_dir=out_dir,
        leaflet=leaflet,
        defect_types=["PLacyl", "TGglyc", "TGacyl"],
        defect_thresholds={"PLacyl": 1, "TGglyc": 2, "TGacyl": 3},
        start=start,
        stop=stop,
        stride=stride,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run packing defect analysis (radius stamping + clustering)"
    )
    parser.add_argument("--top", required=True)
    parser.add_argument("--traj", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--class-json", help="Optional classification rules JSON")
    parser.add_argument("--leaflet", choices=["both", "up", "dw"], default="both")
    parser.add_argument("--start", type=int, default=None, help="First frame")
    parser.add_argument("--stop", type=int, default=None, help="Last frame (exclusive)")
    parser.add_argument("--stride", type=int, default=1, help="Stride between frames")

    args = parser.parse_args()

    analyzer = build_packing_analyzer(
        args.top, args.traj, args.out,
        args.class_json, args.leaflet,
        args.start, args.stop, args.stride
    )
    run_analysis(analyzer)
