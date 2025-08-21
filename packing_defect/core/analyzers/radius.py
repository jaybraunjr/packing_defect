# packing_defect/core/analyzers/radius.py

import os
import subprocess
import numpy as np
import MDAnalysis as mda

from packing_defect.core.analyzers.base import BaseDefectAnalyzer
from packing_defect.core.grid import DefectGrid
from packing_defect.utils import apply_pbc, compute_pairwise_distances


class RadiusDefectAnalyzer(BaseDefectAnalyzer):
    """
    Defect analysis based on radius-cutoff GRO filtering.
    """

    def __init__(
        self,
        base_directory,
        output_dir,
        lipid_types,
        frame_start,
        frame_end,
        protein_atom_count,
        universe=None,
        apply_protein_cutoff=True,
        cutoff_distance=1.5,
    ):
        if universe is None:
            universe = mda.Universe.empty(0)
        super().__init__(universe, output_dir, lipid_types)

        self.base_directory = base_directory
        self.frame_start = frame_start
        self.frame_end = frame_end
        self.protein_atom_count = protein_atom_count
        self.apply_protein_cutoff = apply_protein_cutoff
        self.cutoff_distance = cutoff_distance

        self.defects_up = {}
        self.defects_down = {}
        self.defects_combined = {}

    # ------------------- Internal helpers -------------------

    def _write_filtered_gro(self, input_file, output_file):
        """Filter defect atoms by protein cutoff and write new GRO."""
        with open(input_file, "r", encoding="utf-8") as f:
            lines = f.readlines()
        header, footer = lines[0], lines[-1]
        prot = lines[2 : 2 + self.protein_atom_count]
        defect = lines[2 + self.protein_atom_count : -1]

        if self.apply_protein_cutoff:
            ppos = np.array([list(map(float, l[20:44].split())) for l in prot])
            dpos = np.array([list(map(float, l[20:44].split())) for l in defect])
            mins = np.min(compute_pairwise_distances(dpos, ppos), axis=1)
            defect = [d for i, d in enumerate(defect) if mins[i] > self.cutoff_distance]

        with open(output_file, "w", encoding="utf-8") as f:
            f.write(header)
            f.write(f"{len(prot)+len(defect)}\n")
            f.writelines(prot + defect)
            f.write(footer)

    def _renumber_gro(self, input_file, output_file):
        """Renumber atoms using GROMACS genconf."""
        try:
            subprocess.run(
                ["gmx", "genconf", "-f", input_file, "-o", output_file, "-renumber"],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
        except subprocess.CalledProcessError:
            print(f" Failed renumbering {input_file}")


    def _calculate_defects(self, universe, combine=False):
        """Calculate defect cluster sizes for one Universe."""
        up_sizes, dw_sizes = [], []
        for ts in universe.trajectory:
            box = ts.dimensions[:3]
            apply_pbc(universe.atoms.positions, box)
            half = universe.select_atoms("prop z > 0")
            hz = np.mean(half.positions[:, 2])

            grid = DefectGrid(box_xy=(box[0], box[1]), dx=1.0, dy=1.0, hz=hz)
            for atom in universe.select_atoms(f"prop z > {hz}"):
                grid.update(*atom.position, r=0.0, code=1, leaflet="up")
            for atom in universe.select_atoms(f"prop z < {hz}"):
                grid.update(*atom.position, r=0.0, code=1, leaflet="dw")

            up_sizes.extend(grid.cluster_sizes("up"))
            dw_sizes.extend(grid.cluster_sizes("dw"))
        return (up_sizes + dw_sizes) if combine else (up_sizes, dw_sizes)

    # ------------------- Required overrides -------------------

    def run(self):
        """Main loop: filter, renumber, analyze defects for each lipid type."""
        for lipid_type in self.lipid_types:
            print(f"\n▶ Processing {lipid_type}...")

            input_dir = os.path.join(self.base_directory, lipid_type)
            output_dir = os.path.join(self.output_dir, lipid_type)
            os.makedirs(output_dir, exist_ok=True)

            # Step 1: filter frames
            filtered_files = []
            for i in range(self.frame_start, self.frame_end + 1):
                inp = os.path.join(input_dir, f"{lipid_type}_frame_{i}.gro")
                out = os.path.join(output_dir, f"{lipid_type}_corrected_frame_{i}.gro")
                if os.path.exists(inp):
                    self._write_filtered_gro(inp, out)
                    filtered_files.append(out)

            # Step 2: renumber
            renumbered_files = []
            for f in filtered_files:
                out = os.path.join(output_dir, "renumbered_" + os.path.basename(f))
                self._renumber_gro(f, out)
                renumbered_files.append(out if os.path.exists(out) else f)

            # Step 3: calculate defects
            up, down, combined = [], [], []
            for f in renumbered_files:
                u = mda.Universe(f)
                up_frame, down_frame = self._calculate_defects(u, combine=False)
                up.extend(up_frame)
                down.extend(down_frame)
                combined.extend(up_frame + down_frame)

            # store
            self.results[lipid_type] = {"up": up, "down": down, "combined": combined}
            self.defects_up[lipid_type] = up
            self.defects_down[lipid_type] = down
            self.defects_combined[lipid_type] = combined

            print(f"Done: {lipid_type}, Frames processed: {len(renumbered_files)}")

    def plot(self):
        """Plot histograms of defect sizes."""
        import matplotlib.pyplot as plt

        def plot_hist(defects, label, color):
            h, bins = np.histogram(defects, bins=np.linspace(0, 150, 600))
            h[0] = 0
            total = h.sum()
            if total == 0:
                return
            centers = 0.5 * (bins[1:] + bins[:-1])
            plt.scatter(centers, h / total, label=label, color=color)
            plt.yscale("log")

        # Up/Down
        plt.figure(figsize=(8, 6))
        for lipid, color in zip(self.lipid_types, ["red", "blue", "green"]):
            plot_hist(self.defects_up[lipid], f"{lipid} Up", color)
            plot_hist(self.defects_down[lipid], f"{lipid} Down", color)
        plt.legend()
        plt.title("Top and Bottom Defects")
        plt.xlabel("Defect Size")
        plt.ylabel("Probability (log)")
        plt.savefig(os.path.join(self.output_dir, "top_bottom_defects.png"), dpi=300)
        plt.close()

        # Combined
        plt.figure(figsize=(8, 6))
        for lipid, color in zip(self.lipid_types, ["red", "blue", "green"]):
            plot_hist(self.defects_combined[lipid], f"{lipid} Combined", color)
        plt.legend()
        plt.title("Combined Defects")
        plt.xlabel("Defect Size")
        plt.ylabel("Probability (log)")
        plt.savefig(os.path.join(self.output_dir, "combined_defects.png"), dpi=300)
        plt.close()
