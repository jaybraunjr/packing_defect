import numpy as np
import math
from packing_defect.core.cluster import DefectClustering


class DefectGrid:
    def __init__(self, box_xy, dx=1.0, dy=1.0, hz=None):
        self.dx = dx
        self.dy = dy
        self.Lx, self.Ly = box_xy
        self.hz = hz

        self.nx = int(np.ceil(self.Lx / self.dx))
        self.ny = int(np.ceil(self.Ly / self.dy))
        self.xbins = self.nx
        self.ybins = self.ny

        self.xx, self.yy = np.meshgrid(
            np.linspace(0, self.Lx, self.nx),
            np.linspace(0, self.Ly, self.ny),
            indexing='ij'
        )

        self.grid = {
            'up': np.zeros((self.nx, self.ny), dtype=int),
            'dw': np.zeros((self.nx, self.ny), dtype=int)
        }

        self.zdepth = {
            'up': np.full((self.nx, self.ny), -np.inf),
            'dw': np.full((self.nx, self.ny), np.inf)
        }

    def cluster_sizes(self, leaflet):
        # ensure a binary mask; codes may be >0
        mask = (self._binary_mask(leaflet) != 0).astype(int)
        return DefectClustering.cluster_sizes_from_mask(mask)


    def get_binary_mask(self, leaflet: str, threshold: int) -> np.ndarray:
        """
        Return a binary mask where the defect grid matches the given threshold.

        """
        return (self.grid[leaflet] == threshold).astype(int)



    def update(self, x, y, z, r, code, leaflet):
        # only up or dw matter
        if leaflet not in ('up','dw'):
            return

        r_eff = r + math.sqrt(self.dx**2 + self.dy**2)/2.0
        i0 = int(round(x / self.dx))
        j0 = int(round(y / self.dy))
        max_bin = int(math.ceil(r_eff / self.dx))

        for di in range(-max_bin, max_bin+1):
            ii = i0 + di
            if not (0 <= ii < self.xbins):
                continue
            x_c = (ii + 0.5) * self.dx

            for dj in range(-max_bin, max_bin+1):
                jj = j0 + dj
                if not (0 <= jj < self.ybins):
                    continue
                y_c = (jj + 0.5) * self.dy

                if (x_c - x)**2 + (y_c - y)**2 > r_eff**2:
                    continue

                if leaflet == 'up':
                    if z > self.zdepth['up'][ii, jj]:
                        self.zdepth['up'][ii, jj] = z
                        self.grid   ['up'][ii, jj] = code
                else:  # leaflet == 'dw'
                    if z < self.zdepth['dw'][ii, jj]:
                        self.zdepth['dw'][ii, jj] = z
                        self.grid   ['dw'][ii, jj] = code




    def get_coordinates(self, leaflet, code):
        mask = self.grid[leaflet] == code
        x_coords = self.xx[mask]
        y_coords = self.yy[mask]
        return x_coords, y_coords



    def _binary_mask(self, leaflet):
        return self.grid[leaflet]


