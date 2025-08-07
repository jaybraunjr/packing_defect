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


class RadiusDefectAnalyzer:
    def __init__(self, base_directory, output_base_dir, lipid_types,
                 frame_start, frame_end, protein_atom_count,
                 apply_protein_cutoff=True, cutoff_distance=1.5):
        """
        Set up the analysis configuration.
        """
        self.base_directory = base_directory
        self.output_base_dir = output_base_dir
        self.lipid_types = lipid_types
        self.frame_start = frame_start
        self.frame_end = frame_end
        self.protein_atom_count = protein_atom_count
        self.apply_protein_cutoff = apply_protein_cutoff
        self.cutoff_distance = cutoff_distance
        self.defects_up = {}
        self.defects_down = {}
        self.defects_combined = {}

    def run(self):
        """
        Process files, renumber, and extract defect sizes.
        """
        for lipid_type in self.lipid_types:
            print(f"\n▶ Processing {lipid_type}...")

            input_dir = os.path.join(self.base_directory, lipid_type)
            output_dir = os.path.join(self.output_base_dir, lipid_type)

            from packing_defect.core.analyzer_radius import process_frames
            files = process_frames(
                self.frame_start, self.frame_end, self.protein_atom_count,
                input_dir, lipid_type, output_dir,
                min_cutoff_distance=self.cutoff_distance,
                use_cutoff=self.apply_protein_cutoff
            )

            from packing_defect.core.analyzer_radius import renumber_all_gro_files
            renumbered = renumber_all_gro_files(files)
            from packing_defect.core.analyzer_radius import calculate_defects
            import MDAnalysis as mda

            up, down, combined = [], [], []
            for i in range(self.frame_start, self.frame_end + 1):
                fname = f"renumbered_{lipid_type}_corrected_frame_{i}.gro"
                path = os.path.join(output_dir, fname)

                if not os.path.exists(path):
                    print(f"Missing file: {path}")
                    continue

                u = mda.Universe(path)
                up_frame, down_frame = calculate_defects(u, combine=False)
                up.extend(up_frame)
                down.extend(down_frame)
                combined.extend(up_frame + down_frame)

            # 5. Store results
            self.defects_up[lipid_type] = up
            self.defects_down[lipid_type] = down
            self.defects_combined[lipid_type] = combined

            print(f"Done: {lipid_type}, Frames: {self.frame_end - self.frame_start + 1}")



    def plot(self):
        """
        Create defect histograms for each lipid.
        """
        import matplotlib.pyplot as plt
        import numpy as np

        def plot_hist(defects, label, color):
            h, bins = np.histogram(defects, bins=np.linspace(0, 150, 600))
            h[0] = 0
            total = h.sum()
            if total == 0:
                return
            centers = 0.5 * (bins[1:] + bins[:-1])
            plt.scatter(centers, h / total, label=label, color=color)
            plt.yscale('log')

        # Plot up + down
        plt.figure(figsize=(8, 6))
        for lipid, color in zip(self.lipid_types, ['red', 'blue', 'green']):
            plot_hist(self.defects_up[lipid], f'{lipid} Up', color)
            plot_hist(self.defects_down[lipid], f'{lipid} Down', color)
        plt.legend()
        plt.title('Top and Bottom Defects')
        plt.xlabel('Defect Size')
        plt.ylabel('Probability (log)')
        plt.savefig('top_bottom_defects.png', dpi=300)
        plt.close()

        # Plot combined
        plt.figure(figsize=(8, 6))
        for lipid, color in zip(self.lipid_types, ['red', 'blue', 'green']):
            plot_hist(self.defects_combined[lipid], f'{lipid} Combined', color)
        plt.legend()
        plt.title('Combined Defects')
        plt.xlabel('Defect Size')
        plt.ylabel('Probability (log)')
        plt.savefig('combined_defects.png', dpi=300)
        plt.close()


