import os
import csv
import numpy as np
import MDAnalysis as mda
import argparse

class ProteinDefectInteractions:
    def __init__(self, cutoff=0.1):
        self.cutoff = cutoff
        self.results = {}

    def track_frame(self, universe, frame_index):
        protein = universe.select_atoms("protein")
        defect_atoms = universe.atoms.difference(protein)
        defect_positions = defect_atoms.positions

        hits = []
        for residue in protein.residues:
            com = residue.atoms.center_of_mass()
            if defect_positions.size == 0:
                continue
            dists = np.linalg.norm(defect_positions - com, axis=1)
            if np.any(dists < self.cutoff):
                hits.append((residue.resid, residue.resname))
        self.results[frame_index] = hits


def process_multiple(input_dirs, output_csv, cutoff):
    combined = []
    # record lipid types in given order
    lipid_order = [os.path.basename(d) for d in input_dirs]

    for input_dir in input_dirs:
        lipid = os.path.basename(input_dir)
        tracker = ProteinDefectInteractions(cutoff=cutoff)
        for fname in sorted(os.listdir(input_dir)):
            if not (fname.endswith(".gro") or fname.endswith(".pdb")):
                continue
            path = os.path.join(input_dir, fname)
            print(f"Processing {lipid} frame {fname}")
            u = mda.Universe(path)
            # extract numeric frame index
            try:
                frame_index = int(os.path.splitext(fname)[0].rsplit('_', 1)[-1])
            except ValueError:
                frame_index = os.path.splitext(fname)[0]
            tracker.track_frame(u, frame_index)
        # collect hits with lipid tag and numeric frame
        for frame_idx, hits in tracker.results.items():
            for resid, resname in hits:
                combined.append((lipid, frame_idx, resid, resname))

    # sort combined by lipid order, then numeric frame, then resid
    combined.sort(key=lambda x: (lipid_order.index(x[0]), x[1], x[2]))

    # write out
    os.makedirs(os.path.dirname(output_csv) or ".", exist_ok=True)
    with open(output_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["type", "frame", "resid", "resname"])
        for lipid, frame, resid, resname in combined:
            writer.writerow([lipid, frame, resid, resname])
    print(f"Combined CSV written to {output_csv}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Track protein-defect interactions for multiple lipid types"
    )
    parser.add_argument(
        "-i", "--input-dirs",
        required=True,
        nargs='+',
        help="List of directories containing .gro/.pdb frames (one per lipid)"
    )
    parser.add_argument(
        "-o", "--output-csv",
        default=None,
        help="Path for the combined CSV (default: ./combined_interactions.csv)"
    )
    parser.add_argument(
        "-c", "--cutoff",
        type=float,
        default=1.0,
        help="Distance cutoff in Å for interaction"
    )
    args = parser.parse_args()

    output_csv = args.output_csv or "combined_interactions.csv"
    process_multiple(args.input_dirs, output_csv, args.cutoff)
