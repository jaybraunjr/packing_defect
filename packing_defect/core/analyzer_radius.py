# pylint: disable=C0301,E0402,W1514,R0913,R0914,R0915,C0411,C0413
import os
import subprocess
import numpy as np
import MDAnalysis as mda
from packing_defect.core.grid import DefectGrid
from packing_defect.utils import apply_pbc, compute_pairwise_distances


def calculate_defects(u, combine=False):
    """Calculate defect sizes for upper and lower leaflets."""
    up_sizes = []
    dw_sizes = []

    for ts in u.trajectory:
        box = ts.dimensions[:3]
        apply_pbc(u.atoms.positions, box)
        half = u.select_atoms('prop z > 0')
        hz = np.mean(half.positions[:, 2])

        grid = DefectGrid(box_xy=(box[0], box[1]), dx=1.0, dy=1.0, hz=hz)

        # Stamp defects for each leaflet (binary)
        for atom in u.select_atoms(f"prop z > {hz}"):
            x, y, z = atom.position
            grid.update(x, y, z, r=0.0, code=1, leaflet='up')
        for atom in u.select_atoms(f"prop z < {hz}"):
            x, y, z = atom.position
            grid.update(x, y, z, r=0.0, code=1, leaflet='dw')

        up_sizes.extend(grid.cluster_sizes('up'))
        dw_sizes.extend(grid.cluster_sizes('dw'))

    return (up_sizes + dw_sizes) if combine else (up_sizes, dw_sizes)


def calculate_defects_from_gro(directory_prefix, file_prefix, start_frame=0):
    """Process .gro frames in a directory to accumulate defect sizes."""
    defects_up = []
    defects_down = []
    frame = start_frame
    while True:
        gro_file = os.path.join(directory_prefix, f"{file_prefix}_frame_{frame}.gro")
        if not os.path.exists(gro_file):
            break
        u = mda.Universe(gro_file)
        up, down = calculate_defects(u)
        defects_up.extend(up)
        defects_down.extend(down)
        frame += 1
    return defects_up, defects_down


def process_directories(base_directory):
    """Batch process PLacyl, TGacyl, TGglyc directories."""
    suffix = {
        'resultsPLacyl': 'PLacyl',
        'resultsTGacyl': 'TGacyl',
        'resultsTGglyc': 'TGglyc',
    }
    all_up, all_down = {}, {}
    for dir_key, lip in suffix.items():
        prefix = os.path.join(base_directory, dir_key)
        up, down = calculate_defects_from_gro(prefix, lip)
        all_up[dir_key] = up
        all_down[dir_key] = down
    return all_up, all_down


def write_filtered_gro_by_atom_count(
    input_file, output_file, cutoff_distance=1.5,
    protein_atom_count=627, use_cutoff=True
):
    """Filter defect atoms by distance to protein and write a new GRO."""
    with open(input_file, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    header, footer = lines[0], lines[-1]
    prot = lines[2:2 + protein_atom_count]
    defect = lines[2 + protein_atom_count:-1]
    if use_cutoff:
        ppos = np.array([list(map(float, l[20:44].split())) for l in prot])
        dpos = np.array([list(map(float, l[20:44].split())) for l in defect])
        mins = np.min(compute_pairwise_distances(dpos, ppos), axis=1)
        defect = [d for i, d in enumerate(defect) if mins[i] > cutoff_distance]
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(header)
        f.write(f"{len(prot)+len(defect)}\n")
        f.writelines(prot + defect)
        f.write(footer)


def renumber_gro(input_file, output_file):
    """Use GROMACS genconf to renumber atoms in a GRO file."""
    try:
        subprocess.run(
            ['gmx', 'genconf', '-f', input_file, '-o', output_file, '-renumber'],
            check=True
        )
    except subprocess.CalledProcessError:
        print(f"Failed renumbering {input_file}")


def renumber_all_gro_files(input_files):
    """Renumber multiple GRO files and return list of new filenames."""
    new_files = []
    for inp in input_files:
        out = os.path.join(os.path.dirname(inp), 'renumbered_' + os.path.basename(inp))
        renumber_gro(inp, out)
        new_files.append(out)
    return new_files


def process_frames(
    frame_start: int,
    frame_end: int,
    protein_atom_count: int,
    directory_prefix: str,
    lipid_type: str,
    output_dir: str,
    min_cutoff_distance: float = 1.0,
    use_cutoff: bool = True,
) -> list[str]:
    """
    For each frame X in [frame_start..frame_end], reads
    {directory_prefix}/{lipid_type}_frame_X.gro, filters its defect atoms
    by distance to protein, writes to
    {output_dir}/{lipid_type}_corrected_frame_X.gro, and returns the list
    of all written file paths.
    """
    import os
    from packing_defect.core.analyzer_radius import write_filtered_gro_by_atom_count

    os.makedirs(output_dir, exist_ok=True)
    output_files: list[str] = []

    for frame_idx in range(frame_start, frame_end + 1):
        inp = f"{directory_prefix}/{lipid_type}_frame_{frame_idx}.gro"
        out = os.path.join(output_dir, f"{lipid_type}_corrected_frame_{frame_idx}.gro")
        write_filtered_gro_by_atom_count(
            inp, out, min_cutoff_distance, protein_atom_count, use_cutoff
        )
        output_files.append(out)

    return output_files

