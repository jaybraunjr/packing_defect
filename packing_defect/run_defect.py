import os
import argparse
import json
import logging
import MDAnalysis as mda

from packing_defect.core.classification import DefaultClassification, UserDictClassification
from packing_defect.core.topology import TopologyReader
from packing_defect.core.analyzer_defect import PackingDefectAnalyzer
from packing_defect.core.cluster import DefectClustering


# We need to choose the topology files and defect types
# Json files can allow us to use different atoms for different defects

def run_defect(top_file: str,
               traj_file: str,
               output_dir: str,
               class_json: str = None,
               json_only: bool = False,
               leaflet: str = 'both'):
    radii_file = os.path.join(os.path.dirname(__file__), 'data', 'radii.json')
    with open(radii_file) as f:
        types_radii = json.load(f)

    # Classification logic. UserDictClassification uses logic from the input JSON files.
    # DefaultClassification uses the default lipid droplet classification.
    # Classification gives the groups different numbers, depending on type of defect
    if class_json:
        classifier = UserDictClassification.from_json(class_json)
        classify_fn = classifier.classify
    else:
        classifier = DefaultClassification()
        classify_fn = classifier.classify

    # Radii mapping
    radii = {}

    # this needs to be adressed, currently not working:
    if json_only:
        logging.info("\u2705 JSON-only mode: using JSON for classification; skipping topology parsing")
        if not isinstance(classifier, UserDictClassification):
            raise ValueError("--json-only requires a JSON classifier")

        for resname, atom_map in classifier.rules.items():
            radii[resname] = {}
            for atom_name, code in atom_map.items():
                try:
                    radius = types_radii[atom_name]
                except KeyError:
                    raise ValueError(f"No radius found for atom type '{atom_name}' in radii.json")
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
                    logging.info(f"\u2705 Loaded topology for {resname}")
                elif isinstance(classifier, UserDictClassification) and resname in classifier.rules:
                    logging.warning(f"\u26A0\uFE0F  No topology for {resname}, using JSON fallback")
                    radii[resname] = {}
                    for atom_name, code in classifier.rules[resname].items():
                        try:
                            radius = types_radii[atom_name]
                        except KeyError:
                            raise ValueError(f"\u274C No radius found for atom type '{atom_name}' in radii.json")
                        radii[resname][atom_name] = (radius, code)
                else:
                    logging.warning(f"\u26A0\uFE0F  Skipped {resname}: no topology and no JSON entry")
            except Exception as e:
                logging.error(f"\u274C Failed to handle {resname}: {e}")

    if not radii:
        raise ValueError("\u274C No residue radii definitions found. Check topology parsing or JSON mapping.")




    # Load trajectory
    u = mda.Universe(top_file)
    u.load_new(traj_file)

    resnames = ' '.join(radii.keys())
    MEMB = u.select_atoms(f'resname {resnames}')

    defect_types = ['PLacyl', 'TGglyc', 'TGacyl']
    defect_thresholds = {'PLacyl': 1, 'TGglyc': 2, 'TGacyl': 3}

    os.makedirs(output_dir, exist_ok=True)
    analyzer = PackingDefectAnalyzer(
        atomgroups=[MEMB],
        radii=radii,
        output_prefix=output_dir,
        leaflet=leaflet,
        defect_types=defect_types,
        defect_thresholds=defect_thresholds
    )
    analyzer.run()



    for d in defect_types:
        defect_matrix_list = analyzer.defect_cluster_masks[d]
        dat_path = os.path.join(output_dir, f"{d}.dat")
        DefectClustering.defect_size(defect_matrix_list, nbins=600, bin_max=150, fname=dat_path, prob=True)





if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    parser = argparse.ArgumentParser(description="Run packing defect analysis.")
    parser.add_argument('--top', required=True, help='Topology file (.gro, .psf)')
    parser.add_argument('--traj', required=True, help='Trajectory file (.xtc, .dcd)')
    parser.add_argument('--out', required=True, help='Output directory')
    parser.add_argument('--class', dest='class_json', help='Optional JSON with classification rules')
    parser.add_argument('--json-only', action='store_true', help='Use JSON only; skip topology files')
    parser.add_argument('--leaflet',
                        choices=['both','up','dw'],
                        default='both',
                        help = 'choose leaflet?')
    args = parser.parse_args()
    run_defect(args.top, args.traj, args.out, args.class_json, args.json_only, args.leaflet)