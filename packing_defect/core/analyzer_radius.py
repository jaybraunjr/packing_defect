# packing_defect/core/analyzer_radius.py

import numpy as np
import MDAnalysis as mda
import os
import subprocess
from MDAnalysis.lib import distances
from .grid import DefectGrid
from ..utils import (
    apply_pbc,
    compute_pairwise_distances,
)


def calculate_defects(u, combine=False):
    """
    Scan each frame of the Universe `u`, stamp a binary defect (code=1) wherever
    there’s an atom in each leaflet, then return the connected-component sizes.

    Parameters
    ----------
    u : MDAnalysis.Universe
        Universe containing only the atoms to analyze.
    combine : bool, optional
        If True, returns a single combined list of sizes. Default is False.

    Returns
    -------
    up_sizes, dw_sizes : list of int
        Size of each defect cluster in upper and lower leaflets (if combine=False).
    combined : list of int
        Single list of all defect sizes (if combine=True).
    """
    up_sizes = []
    dw_sizes = []

    for ts in u.trajectory:
        box = ts.dimensions[:3]
        apply_pbc(u.atoms.positions, box)
        half = u.select_atoms("prop z > 0")
        hz = np.mean(half.positions[:, 2])

        Lx, Ly = box[0], box[1]
        grid = DefectGrid(box_xy=(Lx, Ly), dx=1.0, dy=1.0, hz=hz)

        # Stamp defects: code=1, radius=0 (binary) for each leaflet
        for atom in u.select_atoms(f"prop z > {hz}"):
            x, y, z = atom.position
            grid.update(x, y, z, r=0.0, code=1, leaflet="up")
        for atom in u.select_atoms(f"prop z < {hz}"):
            x, y, z = atom.position
            grid.update(x, y, z, r=0.0, code=1, leaflet="dw")

        up_sizes.extend(grid.cluster_sizes("up"))
        dw_sizes.extend(grid.cluster_sizes("dw"))

    if combine:
        return up_sizes + dw_sizes
    return up_sizes, dw_sizes



def calculate_defects_from_gro(directory_prefix, file_prefix, start_frame=0):
    defects_up = []
    defects_down = []
    print(f'Working in directory: {directory_prefix}')

    frame_idx = start_frame
    while True:
        gro_file_path = os.path.join(directory_prefix, f"{file_prefix}_frame_{frame_idx}.gro")
        print(f'Checking for file: {gro_file_path}')

        if not os.path.exists(gro_file_path):
            print(f'File not found: {gro_file_path}')
            break

        u = mda.Universe(gro_file_path)
        up, down = calculate_defects(u)
        defects_up.extend(up)
        defects_down.extend(down)

        frame_idx += 1

    return defects_up, defects_down



def process_directories(base_directory):
    suffix_pairs = {
        'resultsPLacyl': 'PLacyl',
        'resultsTGacyl': 'TGacyl',
        'resultsTGglyc': 'TGglyc'
    }

    all_defects_up = {}
    all_defects_down = {}

    for dir_suffix, file_suffix in suffix_pairs.items():
        directory_prefix = os.path.join(base_directory, dir_suffix)
        defects_up, defects_down = calculate_defects_from_gro(directory_prefix, file_suffix)
        all_defects_up[dir_suffix] = defects_up
        all_defects_down[dir_suffix] = defects_down

    return all_defects_up, all_defects_down



def compute_distances(positions1, positions2):
    return compute_pairwise_distances(positions1, positions2)



def write_filtered_gro_by_atom_count(input_file, output_file, cutoff_distance=1.5, protein_atom_count=627, use_cutoff=True):
    with open(input_file, 'r') as f:
        lines = f.readlines()

    header = lines[0].strip()
    footer = lines[-1].strip()

    protein_atoms = lines[2:2 + protein_atom_count]
    defect_atoms = lines[2 + protein_atom_count:-1]

    if use_cutoff:
        protein_positions = np.array([[float(line[20:28].strip()), float(line[28:36].strip()), float(line[36:44].strip())] for line in protein_atoms])
        defect_positions = np.array([[float(line[20:28].strip()), float(line[28:36].strip()), float(line[36:44].strip())] for line in defect_atoms])
        min_distances = np.min(compute_distances(defect_positions, protein_positions), axis=1)
        filtered_defect_atoms = [atom for i, atom in enumerate(defect_atoms) if min_distances[i] > cutoff_distance]
    else:
        filtered_defect_atoms = defect_atoms

    with open(output_file, 'w') as f:
        f.write(header + '\n')
        f.write(f"{len(protein_atoms) + len(filtered_defect_atoms)}\n")
        f.writelines(protein_atoms)
        f.writelines(filtered_defect_atoms)
        f.write(footer + '\n')



def process_frames(frame_start, frame_end, protein_atom_count, directory_prefix, lipid_type, output_dir, min_cutoff_distance=1.0, use_cutoff=True):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    output_files = []
    for frame_idx in range(frame_start, frame_end + 1):
        input_file_path = f"{directory_prefix}/{lipid_type}_frame_{frame_idx}.gro"
        output_file_path = f"{output_dir}/{lipid_type}_corrected_frame_{frame_idx}.gro"
        write_filtered_gro_by_atom_count(input_file_path, output_file_path, min_cutoff_distance, protein_atom_count, use_cutoff)
        output_files.append(output_file_path)
    return output_files



def renumber_gro(input_file, output_file):
    try:
        subprocess.run(['gmx', 'genconf', '-f', input_file, '-o', output_file, '-renumber'], check=True)
        print(f'Renumbered gro saved to {output_file}')
    except subprocess.CalledProcessError as e:
        print(f'Error in renumbering: {e}')



def renumber_all_gro_files(input_files):
    renumbered_files = []
    for input_file in input_files:
        output_file = os.path.join(os.path.dirname(input_file), "renumbered_" + os.path.basename(input_file))
        renumber_gro(input_file, output_file)
        renumbered_files.append(output_file)
        try:
            os.remove(input_file)
            print(f"Deleted file: {input_file}")
        except OSError as e:
            print(f"Error deleting file {input_file}: {e}")
    return renumbered_files
