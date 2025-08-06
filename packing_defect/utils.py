"""Utility functions for packing defect calculations."""

from __future__ import annotations

import numpy as np
import MDAnalysis as mda
from MDAnalysis import Universe


def apply_pbc(positions: np.ndarray, box: np.ndarray) -> np.ndarray:
    """Apply periodic boundary conditions to *positions* within *box* dimensions."""
    box_xy = np.array([box[0], box[1], 0])
    box_xyz = np.array([box[0], box[1], box[2]])
    return positions - box_xy * np.floor(positions / box_xyz)


def calculate_bounding_box(atoms, padding: float = 0.0) -> dict:
    """Return an axis-aligned bounding box around *atoms* with optional *padding*."""
    coords = atoms.positions
    return {
        "x_min": coords[:, 0].min() - padding,
        "x_max": coords[:, 0].max() + padding,
        "y_min": coords[:, 1].min() - padding,
        "y_max": coords[:, 1].max() + padding,
        "z_min": coords[:, 2].min() - padding,
        "z_max": coords[:, 2].max() + padding,
    }


def initialize_grid(box, dx: float, dy: float, hz: float = 0.0) -> dict:
    """Create a regularly spaced grid covering *box* with spacing *dx* and *dy*."""
    x = np.arange(0, box[0] + dx, dx)
    y = np.arange(0, box[1] + dy, dy)
    xx, yy = np.meshgrid(x, y)
    z_up = np.full_like(xx, hz, dtype=float)
    return {"xx": xx, "yy": yy, "z_up": z_up}


def compute_pairwise_distances(positions1: np.ndarray, positions2: np.ndarray) -> np.ndarray:
    """Compute pairwise Euclidean distances between two point clouds."""
    diff = positions1[:, np.newaxis, :] - positions2
    return np.sqrt(np.sum(diff ** 2, axis=2))


def _dfs(graph: dict, start: int) -> set:
    """Depth first search of *graph* starting at *start* returning visited nodes."""
    visited = set()
    stack = [start]
    while stack:
        node = stack.pop()
        if node not in visited:
            visited.add(node)
            stack.extend(graph.get(node, []))
    return visited


def _make_graph(matrix: np.ndarray) -> dict:
    """Convert adjacency *matrix* into a graph represented as an adjacency list."""
    graph = {}
    for i, row in enumerate(matrix):
        graph[i] = {j for j, val in enumerate(row) if val}
    return graph


def filter_defects_by_distance(defects: np.ndarray, positions: np.ndarray, threshold: float) -> np.ndarray:
    """Remove defect points closer than *threshold* to any given *positions*."""
    if len(defects) == 0 or len(positions) == 0:
        return defects
    dists = compute_pairwise_distances(defects, positions)
    mask = np.all(dists >= threshold, axis=1)
    return defects[mask]


def get_defect_coordinates(grid: np.ndarray, value: int):
    """Return coordinate arrays for cells in *grid* equal to *value*."""
    coords = np.argwhere(grid == value)
    return coords[:, 0], coords[:, 1]


def validate_defect_thresholds(defect_types, defect_thresholds):
    """Ensure each defect type has a corresponding threshold defined."""
    for defect_type in defect_types:
        if defect_type not in defect_thresholds:
            raise ValueError(f"Missing threshold for defect type: {defect_type}")


def write_combined_gro(protein_atoms, defect_atoms, dimensions, filepath):
    """Write a combined GRO file from *protein_atoms* and *defect_atoms*."""
    combined = mda.Merge(protein_atoms, defect_atoms)
    combined.atoms.positions[: len(protein_atoms)] = protein_atoms.positions
    combined.atoms.positions[len(protein_atoms) :] = defect_atoms.positions
    combined.trajectory.ts.dimensions = dimensions
    combined.atoms.write(filepath)


def initialize_empty_defect_universe(n_atoms, nframes, dims, dt):
    """Create an empty :class:`MDAnalysis.Universe` for defect coordinates."""
    fac = np.zeros((nframes, n_atoms, 3))
    df = Universe.empty(
        n_atoms=n_atoms,
        n_residues=n_atoms,
        atom_resindex=np.arange(n_atoms),
        residue_segindex=[0] * n_atoms,
        trajectory=True,
    )
    df.add_TopologyAttr("resname", ["O"] * n_atoms)
    df.add_TopologyAttr("name", ["O"] * n_atoms)
    df.add_TopologyAttr("resid", np.arange(n_atoms) + 1)
    df.load_new(fac, order="fac")
    df.trajectory[0].dt = dt
    for i, _ in enumerate(df.trajectory):
        df.trajectory[i].dimensions = dims[i]
    return df
