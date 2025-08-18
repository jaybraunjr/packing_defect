# pylint: disable=C0301,C0304,R0913,R0914,R0912,R0915,R1702,W1514,W0707,W0718
import os
import argparse
import json
from pathlib import Path
from pprint import pprint
import MDAnalysis as mda

from packing_defect.core.classification import DefaultClassification, UserDictClassification
from packing_defect.core.topology import TopologyReader
from packing_defect.core.analyzer_defect import PackingDefectAnalyzer
from packing_defect.core.cluster import DefectClustering


def run_defect(
    top_file: str,
    traj_file: str,
    output_dir: str,
    class_json: str = None,
    json_only: bool = False,
    leaflet: str = 'both',
    start: int = None,
    stop: int = None,
    stride: int = 1
):
    radii_file = os.path.join(os.path.dirname(__file__), 'data', 'radii.json')
    with open(radii_file, encoding='utf-8') as f:
        types_radii = json.load(f)

    if class_json:
        classifier = UserDictClassification.from_json(class_json)
        classify_fn = classifier.classify
    else:
        classifier = DefaultClassification()
        classify_fn = classifier.classify

    radii = {}

    if json_only:
        print("\u2705 JSON-only mode: using JSON for classification; skipping topology parsing")
        if not isinstance(classifier, UserDictClassification):
            raise ValueError("--json-only requires a JSON classifier")

        for resname, atom_map in classifier.rules.items():
            radii[resname] = {}
            for atom_name, code in atom_map.items():
                try:
                    radius = types_radii[atom_name]
                except KeyError as exc:
                    raise ValueError(
                        f"No radius found for atom type '{atom_name}' in radii.json"
                    ) from exc
                radii[resname][atom_name] = (radius, code)
    else:
        base_top = os.path.join(os.path.dirname(__file__), 'data', 'top')
        topo_map = {
            'POPC': 'top_all36_lipid.rtf',
            'DOPE': 'top_all36_lipid.rtf',
            'SAPI': 'toppar_all36_lipid_inositol.str',
            'TRIO': 'TRIO.rtf',
            'CHYO': 'CHYO.rtf',
        }

        topo_reader = TopologyReader(types_radii, classify_fn)

        for resname, topo_filename in topo_map.items():
            topopath = os.path.join(base_top, topo_filename)
            try:
                if os.path.exists(topopath):
                    radii[resname] = topo_reader.read(resname, topopath)
                    print(f"\u2705 Loaded topology for {resname}")
                elif isinstance(classifier, UserDictClassification) and resname in classifier.rules:
                    print(f"\u26A0\uFE0F  No topology for {resname}, using JSON fallback")
                    radii[resname] = {}
                    for atom_name, code in classifier.rules[resname].items():
                        try:
                            radius = types_radii[atom_name]
                        except KeyError as exc:
                            raise ValueError(
                                f"No radius found for atom type '{atom_name}' in radii.json"
                            ) from exc
                        radii[resname][atom_name] = (radius, code)
                else:
                    print(f"\u26A0\uFE0F  Skipped {resname}: no topology and no JSON entry")
            except Exception:
                print(f"\u274C Failed to handle {resname}")

    if not radii:
        raise ValueError(
            "No residue radii definitions found. Check topology parsing or JSON mapping."
        )

    u = mda.Universe(top_file)
    u.load_new(traj_file)

    resnames = ' '.join(radii.keys())
    memb = u.select_atoms(f'resname {resnames}')

    defect_types = ['PLacyl', 'TGglyc', 'TGacyl']
    defect_thresholds = {'PLacyl': 1, 'TGglyc': 2, 'TGacyl': 3}

    os.makedirs(output_dir, exist_ok=True)
    analyzer = PackingDefectAnalyzer(
        atomgroups=[memb],
        radii=radii,
        output_prefix=output_dir,
        leaflet=leaflet,
        defect_types=defect_types,
        defect_thresholds=defect_thresholds,
        start=start,
        stop=stop,
        stride=stride
    )
    analyzer.run()

    for defect_type in defect_types:
        masks = analyzer.defect_cluster_masks[defect_type]
        dat_path = os.path.join(output_dir, f"{defect_type}.dat")
        DefectClustering.defect_size(
            masks,
            nbins=600,
            bin_max=150,
            fname=dat_path,
            prob=True
        )




def build_parser():
    parser = argparse.ArgumentParser(description="Run packing defect analysis.")
    parser.add_argument('--top', required=True, help='Topology file (.gro, .psf)')
    parser.add_argument('--traj', required=True, help='Trajectory file (.xtc, .dcd)')
    parser.add_argument('--out', required=True, help='Output directory')
    parser.add_argument('--class', dest='class_json', help='Optional JSON with classification rules')
    parser.add_argument('--json-only', action='store_true', help='Use JSON only; skip topology files')
    parser.add_argument('--leaflet', choices=['both', 'up', 'dw'], default='both', help='choose leaflet?')
    parser.add_argument("--start", type=int, default=None, help="first frame to include")
    parser.add_argument("--stop",  type=int, default=None, help="stop before this frame")
    parser.add_argument("--stride", type=int, default=1, help="take one every N frames")
    return parser

if __name__ == '__main__':
    parser = build_parser()
    args = parser.parse_args()

    run_defect(
        args.top,
        args.traj,
        args.out,
        args.class_json,
        args.json_only,
        args.leaflet,
        start=args.start,
        stop=args.stop,
        stride=args.stride
    )
