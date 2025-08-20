import argparse
import os
import json
import MDAnalysis as mda

from packing_defect.core.classification import DefaultClassification, UserDictClassification
from packing_defect.core.topology import TopologyReader
from packing_defect.core.analyzers.packing import PackingDefectAnalyzer
from packing_defect.run_utils import run_analysis


def build_packing_analyzer(top_file, traj_file, out_dir, class_json=None, leaflet="both",
                           start=None, stop=None, stride=1):
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
    parser = argparse.ArgumentParser(description="Run packing defect analysis")
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
