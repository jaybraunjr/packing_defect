"""CLI entrypoint for radius-based GRO filtering and defect sizing.

This utility filters GRO frames by a protein-distance cutoff, optionally
renumbers atoms, and measures defect cluster sizes per leaflet.
"""

import argparse
import MDAnalysis as mda

from packing_defect.core.analyzers.radius import RadiusDefectAnalyzer
from packing_defect.run_utils import run_analysis


def build_radius_analyzer(
    input_dir: str,
    output_dir: str,
    lipids: list[str],
    start: int,
    end: int,
    protein_count: int,
    cutoff: float,
    no_cutoff: bool,
):
    """Construct a RadiusDefectAnalyzer for filtered GRO inputs.

    Parameters
    ----------
    input_dir : str
        Directory containing per-lipid subfolders with ``*_frame_N.gro``.
    output_dir : str
        Directory to write corrected/renumbered outputs and plots.
    lipids : list[str]
        Lipid residue-name prefixes for subfolders.
    start, end : int
        Inclusive range of frame indices to process.
    protein_count : int
        Number of protein atoms at the top of each GRO file.
    cutoff : float
        Distance cutoff (in same units as GRO) to exclude near-protein atoms.
    no_cutoff : bool
        If True, skip the protein-distance filtering step.

    Returns
    -------
    RadiusDefectAnalyzer
        Configured analyzer instance ready to ``run()``.
    """
    # Dummy universe to satisfy BaseDefectAnalyzer interface
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
