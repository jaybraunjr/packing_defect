# scripts/run_radius.py
import os
import numpy as np
import MDAnalysis as mda
import matplotlib.pyplot as plt
from defect import radius2

def plot_histogram(defects, label, color=None):
    h, _ = np.histogram(defects, bins=np.linspace(0, 150, 600))
    h[0] = 0
    binp = 0.5 * (_[1:] + _[:-1])
    plt.scatter(binp, h / np.sum(h), label=label, color=color)
    plt.yscale('log')

def run_radius(base_directory, output_base_dir, lipid_types, frame_start, frame_end,
               protein_atom_count, apply_protein_cutoff, cutoff_distance=1.5):

    processed_defects_up = {}
    processed_defects_down = {}
    processed_defects_combined = {}

    for lipid_type in lipid_types:
        directory_prefix = os.path.join(base_directory, lipid_type)
        output_dir = os.path.join(output_base_dir, lipid_type)

        print(f"Processing {lipid_type}...")
        output_files = radius2.process_frames(
            frame_start, frame_end, protein_atom_count,
            directory_prefix, lipid_type, output_dir,
            min_cutoff_distance=cutoff_distance, use_cutoff=apply_protein_cutoff
        )
        print(f"Renumbering {lipid_type} files...")
        renumbered_files = radius2.renumber_all_gro_files(output_files)
        print(f"Completed processing for {lipid_type}.")

        defects_up_all, defects_down_all, defects_combined_all = [], [], []
        for frame_idx in range(frame_start, frame_end + 1):
            processed_file_name = f"renumbered_{lipid_type}_corrected_frame_{frame_idx}.gro"
            processed_file_path = os.path.join(output_dir, processed_file_name)

            if os.path.exists(processed_file_path):
                u = mda.Universe(processed_file_path)
                defects_up, defects_down = radius2.calculate_defects(u, combine=False)
                defects_up_all.extend(defects_up)
                defects_down_all.extend(defects_down)
                defects_combined_all.extend(defects_up + defects_down)

        processed_defects_up[lipid_type] = defects_up_all
        processed_defects_down[lipid_type] = defects_down_all
        processed_defects_combined[lipid_type] = defects_combined_all

    plt.figure(figsize=(8, 6))
    for lipid_type, color in zip(lipid_types, ['red', 'blue', 'green']):
        plot_histogram(processed_defects_up[lipid_type], label=f'{lipid_type} Up', color=color)
        plot_histogram(processed_defects_down[lipid_type], label=f'{lipid_type} protein', color=color)
    plt.legend()
    plt.title('Top and Bottom Defects')
    plt.xlabel('Defect Size')
    plt.ylabel('Frequency (log scale)')
    plt.savefig('top_bottom_defect_histogram2.png', dpi=300, bbox_inches='tight')
    plt.close()

    plt.figure(figsize=(8, 6))
    for lipid_type, color in zip(lipid_types, ['red', 'blue', 'green']):
        plot_histogram(processed_defects_combined[lipid_type], label=f'{lipid_type} Combined', color=color)
    plt.legend()
    plt.title('Combined Defects')
    plt.xlabel('Defect Size')
    plt.ylabel('Frequency (log scale)')
    plt.savefig('combined_defect_histogram2.png', dpi=300, bbox_inches='tight')
    plt.close()