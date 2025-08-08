# pylint: disable=C0301,C0304,C0411,C0413,R0902,R0913,R0914,R0915,W0201,W0612
import numpy as np
import os
from MDAnalysis import AtomGroup
from packing_defect.core.grid import DefectGrid
from packing_defect.core.cluster import DefectClustering
from packing_defect.core.writer import write_defect_coordinates
from packing_defect.utils import (
    apply_pbc,
    initialize_empty_defect_universe,
    validate_defect_thresholds,
)


class PackingDefectAnalyzer:
    def __init__(
        self,
        atomgroups,
        radii,
        output_prefix: str = './',
        leaflet: str = 'both',
        defect_types=None,
        defect_thresholds=None,
    ):
        self.N = 50000
        self.universe = atomgroups[0].universe
        self.dt = self.universe.trajectory[0].dt
        self.atomgroups = atomgroups
        self.dx = 1.0
        self.dy = 1.0
        self.radii = radii
        self.output_prefix = output_prefix
        self.leaflet = leaflet
        self.protein_atoms = self.universe.select_atoms('protein', updating=True)
        self._results = []

        self.defect_types = defect_types or ['PLacyl', 'TGglyc', 'TGacyl']
        self.defect_thresholds = (
            defect_thresholds if defect_thresholds
            else {t: i + 1 for i, t in enumerate(self.defect_types)}
        )
        validate_defect_thresholds(self.defect_types, self.defect_thresholds)

    def run(self):
        for ts in self.universe.trajectory:
            print(
                f"Processing frame {ts.frame}, time: {ts.time:.3f}, "
                f"pbc: {ts.dimensions[:3]}"
            )
            result = self._analyze_frame(ts)
            if result:
                self._results.append(result)
        if self._results:
            self._finalize()
        else:
            print('No frames processed.')

    def _analyze_frame(self, ts):
        ag = self.atomgroups[0]
        dim = ts.dimensions.copy()
        pbc = dim[:3]
        ag.universe.atoms.positions = apply_pbc(
            ag.universe.atoms.positions, pbc
        )

        hz = np.average(ag.select_atoms('name P').positions[:, 2])
        grid = DefectGrid(
            box_xy=(pbc[0], pbc[1]), dx=self.dx, dy=self.dy, hz=hz
        )
        zlim, PL = self._classify_leaflets(ag, grid)
        return grid, PL['up'] + 2, PL['dw'] - 2, dim

    def _classify_leaflets(
        self, ag: AtomGroup, grid: DefectGrid
    ):
        hz = grid.hz
        PL = {
            'up': ag.select_atoms(f'name P and prop z > {hz}').center_of_mass()[2],
            'dw': ag.select_atoms(f'name P and prop z < {hz}').center_of_mass()[2],
        }

        atoms = {}
        if self.leaflet in ('both', 'up'):
            atoms['up'] = ag.select_atoms(
                f'prop z > {PL["up"] - 10}'
            )
        if self.leaflet in ('both', 'dw'):
            atoms['dw'] = ag.select_atoms(
                f'prop z < {PL["dw"] + 10}'
            )

        for leaflet, group in atoms.items():
            for atom in group:
                try:
                    radius, code = self.radii[atom.resname][atom.name]
                except KeyError:
                    continue
                x, y, z = atom.position
                grid.update(x, y, z, radius, code, leaflet)

        return {
            'up': np.max(ag.positions[:, 2]),
            'dw': np.min(ag.positions[:, 2]),
        }, PL



    def _finalize(self):
        grids, zlimup, zlimdw, dims = zip(*self._results)
        df = initialize_empty_defect_universe(
            self.N, len(dims), dims, self.dt
        )
        defect_uni = {d: df.copy() for d in self.defect_types}
        defect_clu = {d: [] for d in self.defect_types}

        for d in self.defect_types:
            threshold = self.defect_thresholds[d]
            for i, ts in enumerate(defect_uni[d].trajectory):
                num = 0
                mask_up = grids[i].get_binary_mask('up', threshold)
                defect_clu[d].append(mask_up.astype(int))
                num = self._populate_defect_coords(
                    threshold, grids[i], zlimup[i], defect_uni[d], num, 'up'
                )
                mask_dw = grids[i].get_binary_mask('dw', threshold)
                defect_clu[d].append(mask_dw.astype(int))
                self._populate_defect_coords(
                    threshold, grids[i], zlimdw[i], defect_uni[d], num, 'dw'
                )

        self._write_outputs(defect_uni)
        self.defect_cluster_masks = defect_clu

        for d in self.defect_types:
            dat_path = os.path.join(
                self.output_prefix, f'{d}.dat'
            )
            DefectClustering.defect_size(
                defect_clu[d], nbins=600, bin_max=150,
                fname=dat_path, prob=True
            )

    def _populate_defect_coords(
        self, threshold, grid: DefectGrid, zlim, universe, num, leaflet
    ):
        xs, ys = grid.get_coordinates(leaflet, threshold)
        for x1, y1 in zip(xs, ys):
            if num >= self.N:
                break
            universe.atoms[num].position = np.array([x1, y1, zlim])
            num += 1
        return num

    def _write_outputs(self, defect_uni):
        for d, u in defect_uni.items():
            outdir = os.path.join(self.output_prefix, d)
            os.makedirs(outdir, exist_ok=True)
            for i, ts in enumerate(u.trajectory):
                _ = self.protein_atoms.universe.trajectory[i]
                path = os.path.join(
                    outdir, f'{d}_frame_{i}.gro'
                )
                write_defect_coordinates(
                    self.protein_atoms, u.atoms, ts.dimensions, path
                )
