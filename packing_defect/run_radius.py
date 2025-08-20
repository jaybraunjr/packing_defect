import argparse
import MDAnalysis as mda

from packing_defect.core.analyzers.radius import RadiusDefectAnalyzer
from packing_defect.run_utils import run_analysis


def build_radius_analyzer(input_dir, output_dir, lipids, start, end, protein_count, cutoff, no_cutoff):
    # dummy universe just for consistency (BaseDefectAnalyzer requires one)
    u = mda.Universe()  # creates empty placeholder
    return RadiusDefectAnalyzer(
        universe=u,
        base_directory=input_dir,
        output_dir=output_dir,
        lipid_types=lipids,
        frame_start=start,
        frame_end=end,
        protein_atom_count=protein_count,
        apply_protein_cutoff=not no_cutoff,
        cutoff_distance=cutoff,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run radius defect analysis")
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--lipids", nargs="+", required=True)
    parser.add_argument("--start", type=int, required=True)
    parser.add_argument("--end", type=int, required=True)
    parser.add_argument("--protein-count", type=int, required=True)
    parser.add_argument("--cutoff", type=float, default=1.5)
    parser.add_argument("--no-cutoff", action="store_true")
    args = parser.parse_args()

    analyzer = build_radius_analyzer(
        args.input,
        args.output,
        args.lipids,
        args.start,
        args.end,
        args.protein_count,
        args.cutoff,
        args.no_cutoff,
    )
    run_analysis(analyzer)
