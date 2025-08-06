# pylint: disable=C0301,C0413,R0913,R0914,W0612
import os
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import MDAnalysis as mda
from packing_defect.core.analyzer_radius import (
    process_frames,
    renumber_all_gro_files,
    calculate_defects,
)


def plot_histogram(defects, label, color=None):
    h, _ = np.histogram(defects, bins=np.linspace(0, 150, 600))
    h[0] = 0
    total = h.sum()
    if total == 0:
        print(f"Skipping plot for {label} (no defect data)")
        return
    binp = 0.5 * (_[1:] + _[:-1])
    plt.scatter(binp, h / total, label=label, color=color)
    plt.yscale('log')


def run_radius(
    base_directory: str,
    output_base_dir: str,
    lipid_types: list,
    frame_start: int,
    frame_end: int,
    protein_atom_count: int,
    apply_protein_cutoff: bool = True,
    cutoff_distance: float = 1.5,
):
    processed_defects_up = {}
    processed_defects_down = {}
    processed_defects_combined = {}

    for lipid_type in lipid_types:
        directory_prefix = os.path.join(base_directory, lipid_type)
        output_dir = os.path.join(output_base_dir, lipid_type)

        print(f"Processing {lipid_type}...")
        output_files = process_frames(
            frame_start, frame_end, protein_atom_count,
            directory_prefix, lipid_type, output_dir,
            min_cutoff_distance=cutoff_distance, use_cutoff=apply_protein_cutoff
        )
        print(f"Renumbering {lipid_type} files...")
        _renumbered_files = renumber_all_gro_files(output_files)
        print(f"Completed processing for {lipid_type}.")

        defects_up_all, defects_down_all, defects_combined_all = [], [], []
        for frame_idx in range(frame_start, frame_end + 1):
            processed_file_name = (
                f"renumbered_{lipid_type}_corrected_frame_{frame_idx}.gro"
            )
            processed_file_path = os.path.join(
                output_dir, processed_file_name
            )

            if os.path.exists(processed_file_path):
                u = mda.Universe(processed_file_path)
                defects_up, defects_down = calculate_defects(
                    u, combine=False
                )
                defects_up_all.extend(defects_up)
                defects_down_all.extend(defects_down)
                defects_combined_all.extend(defects_up + defects_down)

        processed_defects_up[lipid_type] = defects_up_all
        processed_defects_down[lipid_type] = defects_down_all
        processed_defects_combined[lipid_type] = defects_combined_all

    plt.figure(figsize=(8, 6))
    for lipid_type, color in zip(lipid_types, ['red', 'blue', 'green']):
        plot_histogram(
            processed_defects_up[lipid_type],
            label=f'{lipid_type} Up',
            color=color,
        )
        plot_histogram(
            processed_defects_down[lipid_type],
            label=f'{lipid_type} Down',
            color=color,
        )
    plt.legend()
    plt.title('Top and Bottom Defects')
    plt.xlabel('Defect Size')
    plt.ylabel('Frequency (log scale)')
    plt.savefig(
        'top_bottom_defect_histogram2.png',
        dpi=300,
        bbox_inches='tight',
    )
    plt.close()

    plt.figure(figsize=(8, 6))
    for lipid_type, color in zip(lipid_types, ['red', 'blue', 'green']):
        plot_histogram(
            processed_defects_combined[lipid_type],
            label=f'{lipid_type} Combined',
            color=color,
        )
    plt.legend()
    plt.title('Combined Defects')
    plt.xlabel('Defect Size')
    plt.ylabel('Frequency (log scale)')
    plt.savefig(
        'combined_defect_histogram2.png',
        dpi=300,
        bbox_inches='tight',
    )
    plt.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run radius-based defect analysis."
    )
    parser.add_argument(
        '--input', required=True,
        help='Input directory with *_frame_X.gro files'
    )
    parser.add_argument(
        '--output', required=True,
        help='Output directory'
    )
    parser.add_argument(
        '--lipids', nargs='+', required=True,
        help='List of lipid prefixes (e.g. PLacyl TGacyl)'
    )
    parser.add_argument(
        '--start', type=int, required=True,
        help='Starting frame index'
    )
    parser.add_argument(
        '--end', type=int, required=True,
        help='Ending frame index'
    )
    parser.add_argument(
        '--protein-count', type=int, required=True,
        help='Number of protein atoms in each .gro file'
    )
    parser.add_argument(
        '--cutoff', type=float, default=1.5,
        help='Minimum cutoff distance to protein (default: 1.5 Å)'
    )
    parser.add_argument(
        '--no-cutoff', action='store_true',
        help='Disable distance cutoff filter'
    )
    args = parser.parse_args()

    run_radius(
        base_directory=args.input,
        output_base_dir=args.output,
        lipid_types=args.lipids,
        frame_start=args.start,
        frame_end=args.end,
        protein_atom_count=args.protein_count,
        apply_protein_cutoff=not args.no_cutoff,
        cutoff_distance=args.cutoff,
    )
