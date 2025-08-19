"""
interactions.py

Tools for detecting protein–defect interactions in lipid membrane simulations.
This module can process multiple lipid types, identify protein residues that 
interact with packing defects within a distance cutoff, and write results to a CSV.
"""


import os
import csv
import numpy as np
import MDAnalysis as mda
import argparse


class ProteinDefectInteractions:
    """
    Tracks protein–defect interactions in a single trajectory.

    Parameters
    ----------
    cutoff : float, optional
        Minimum distance in Å between a protein residue center of mass
        and any defect atom for the residue to be counted as interacting.
        Default is 0.1 Å.

    Attributes
    ----------
    cutoff : float
        Interaction distance threshold in Å.
    results : dict
        Mapping from frame index → list of (resid, resname) tuples for interacting residues.
    """

    def __init__(self, cutoff=0.1):
        self.cutoff = cutoff
        self.results = {}

    def track_frame(self, universe, frame_index):
        """
        Analyze a single frame and record protein residues near defects.

        Parameters
        ----------
        universe : MDAnalysis.Universe
            Universe containing the current frame.
        frame_index : int or str
            Index or identifier for the current frame.

        Notes
        -----
        - All protein residues whose center of mass lies within `cutoff` Å
          of any defect atom are recorded.
        - Defect atoms are defined as all atoms not belonging to the protein.
        """
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
    """
    Process multiple lipid-type directories and record protein–defect interactions.

    Parameters
    ----------
    input_dirs : list of str
        Paths to directories containing `.gro` or `.pdb` files for each lipid type.
        Each directory is treated as a separate lipid type, and its basename
        is used as the lipid label in the output.
    output_csv : str
        Path to the CSV file to be written. Parent directories will be created if necessary.
    cutoff : float
        Distance cutoff in Å for protein–defect interactions.

    Output CSV Format
    -----------------
    type, frame, resid, resname
        - type : lipid type (directory name)
        - frame : frame index (int or str extracted from filename)
        - resid : residue ID of interacting protein residue
        - resname : residue name of interacting protein residue

    Notes
    -----
    - Frame index is parsed from the filename by extracting the last underscore-delimited
      component before the extension. If parsing fails, the filename stem is used as-is.
    - The CSV rows are sorted by lipid type (order given in `input_dirs`), then frame index,
      then residue ID.
    """
    combined = []
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
            try:
                frame_index = int(os.path.splitext(fname)[0].rsplit('_', 1)[-1])
            except ValueError:
                frame_index = os.path.splitext(fname)[0]
            tracker.track_frame(u, frame_index)
        # Collect hits with lipid tag
        for frame_idx, hits in tracker.results.items():
            for resid, resname in hits:
                combined.append((lipid, frame_idx, resid, resname))

    combined.sort(key=lambda x: (lipid_order.index(x[0]), x[1], x[2]))

    os.makedirs(os.path.dirname(output_csv) or ".", exist_ok=True)
    with open(output_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["type", "frame", "resid", "resname"])
        for lipid, frame, resid, resname in combined:
            writer.writerow([lipid, frame, resid, resname])
    print(f"Combined CSV written to {output_csv}")


def run(input_dirs, cutoff=1.0, output_csv="combined_interactions.csv", as_dataframe=True):
    """
    Compute protein–defect interactions, write the CSV, and return results.

    Parameters
    ----------
    input_dirs : list of str
        Directories containing .gro/.pdb frames, one per lipid type.
    cutoff : float, optional
        Distance cutoff in Å. Default 1.0.
    output_csv : str, optional
        Output CSV path. Default "combined_interactions.csv".
    as_dataframe : bool, optional
        If True and pandas is available, return a pandas DataFrame.
        Otherwise return a list of dicts.

    Returns
    -------
    pandas.DataFrame or list[dict]
        Rows with columns: type, frame, resid, resname.
    """
    # Write CSV using the existing logic
    os.makedirs(os.path.dirname(output_csv) or ".", exist_ok=True)
    process_multiple(input_dirs, output_csv, cutoff)

    if as_dataframe:
        try:
            import pandas as pd
            return pd.read_csv(output_csv)
        except Exception:
            pass

    rows = []
    with open(output_csv, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                row["frame"] = int(row["frame"])
            except Exception:
                pass
            try:
                row["resid"] = int(row["resid"])
            except Exception:
                pass
            rows.append(row)
    return rows




if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Track protein–defect interactions for multiple lipid types"
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
