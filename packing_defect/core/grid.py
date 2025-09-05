import numpy as np
import math
from packing_defect.core.cluster import DefectClustering


class DefectGrid:
    """
    2D grid accumulator for membrane packing defects.

    The grid is defined over the XY plane and tracks two leaflets
    ("up" and "dw"). Each grid cell stores an integer code identifying
    the type of defect stamped there and a per-leaflet depth value used
    to keep only the most exposed atom per cell.

    Parameters
    ----------
    box_xy : tuple[float, float]
        Box lengths along X and Y in the same units as positions.
    dx, dy : float, optional
        Grid spacing along X and Y.
    hz : float, optional
        Mid-plane Z used to separate leaflets; when omitted, callers
        should ensure consistent leaflet assignment when calling ``update``.

    Attributes
    ----------
    grid : dict[str, np.ndarray]
        Integer codes for each leaflet, shape ``(nx, ny)``.
    zdepth : dict[str, np.ndarray]
        Depth tracker per leaflet used to keep the outermost atom per cell.
    nx, ny : int
        Number of grid bins along X and Y.
    dx, dy : float
        Grid spacings along X and Y.
    """

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
        """Return connected-component sizes for a leaflet.

        Parameters
        ----------
        leaflet : {"up", "dw"}
            Leaflet to analyze.

        Returns
        -------
        list[int]
            Sizes (in grid cells) of each connected component.
        """
        mask = (self._binary_mask(leaflet) != 0).astype(int)
        return DefectClustering.cluster_sizes_from_mask(mask)


    def get_binary_mask(self, leaflet: str, threshold: int) -> np.ndarray:
        """Return a binary mask where grid values equal ``threshold``.

        Parameters
        ----------
        leaflet : {"up", "dw"}
            Leaflet to extract.
        threshold : int
            Code value considered as a defect in the mask.

        Returns
        -------
        np.ndarray
            Binary array of shape ``(nx, ny)`` where matches are 1 else 0.
        """
        return (self.grid[leaflet] == threshold).astype(int)



    def update(self, x, y, z, r, code, leaflet):
        """Stamp a circular defect into the grid.

        For the selected ``leaflet``, all grid cells whose centers fall
        within the effective radius are assigned ``code``. The most exposed
        atom along Z wins per cell: for "up" we keep the maximum Z, for
        "dw" the minimum Z.

        Parameters
        ----------
        x, y, z : float
            Atom coordinates.
        r : float
            Stamp radius.
        code : int
            Integer code to write into covered cells.
        leaflet : {"up", "dw"}
            Target leaflet.
        """
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
        """Coordinates of grid cell centers matching ``code``.

        Parameters
        ----------
        leaflet : {"up", "dw"}
            Leaflet to query.
        code : int
            Integer code to match.

        Returns
        -------
        tuple[np.ndarray, np.ndarray]
            Arrays of X and Y coordinates of matching grid centers.
        """
        mask = self.grid[leaflet] == code
        x_coords = self.xx[mask]
        y_coords = self.yy[mask]
        return x_coords, y_coords



    def _binary_mask(self, leaflet):
        """Return the raw grid values for ``leaflet``.

        Notes
        -----
        Values may be greater than zero; callers often convert to 0/1.
        """
        return self.grid[leaflet]


