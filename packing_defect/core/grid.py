

import numpy as np


def _make_graph(matrix):
    """Convert a binary mask into an adjacency map of occupied nodes."""
    graph = {}
    nx, ny = matrix.shape
    for (i, j), val in np.ndenumerate(matrix):
        if val == 0:
            continue
        idx = i * ny + j
        neighbors = []
        for di in (-1, 0, 1):
            for dj in (-1, 0, 1):
                ni, nj = (i + di) % nx, (j + dj) % ny
                if (di or dj) and matrix[ni, nj]:
                    neighbors.append(ni * ny + nj)
        graph[idx] = set(neighbors)
    return graph


def _dfs(graph, start):
    """Depth-first search returning one connected component."""
    visited, stack = set(), [start]
    while stack:
        v = stack.pop()
        if v not in visited:
            visited.add(v)
            stack.extend(graph[v] - visited)
    return visited


class DefectGrid:
    def __init__(self, box_xy, dx=1.0, dy=1.0, hz=None):
        self.dx = dx
        self.dy = dy
        self.Lx, self.Ly = box_xy
        self.hz = hz

        self.xbins = int(np.ceil(self.Lx / self.dx))
        self.ybins = int(np.ceil(self.Ly / self.dy))

        self.xx, self.yy = np.meshgrid(
            np.linspace(0, self.Lx, self.xbins),
            np.linspace(0, self.Ly, self.ybins),
            indexing='ij'
        )

        self.up = np.full((self.xbins, self.ybins), 0, dtype=int)
        self.dw = np.full((self.xbins, self.ybins), 0, dtype=int)
        self.z_up = np.full((self.xbins, self.ybins), -np.inf)
        self.z_dw = np.full((self.xbins, self.ybins), np.inf)

    def update(self, x, y, z, r, code, leaflet):
        if leaflet not in ['up', 'dw']:
            return

        i = int(round(x / self.dx))
        j = int(round(y / self.dy))
        if 0 <= i < self.xbins and 0 <= j < self.ybins:
            if leaflet == 'up':
                if z > self.z_up[i, j]:
                    self.up[i, j] = code
                    self.z_up[i, j] = z
            elif leaflet == 'dw':
                if z < self.z_dw[i, j]:
                    self.dw[i, j] = code
                    self.z_dw[i, j] = z

    def get_coordinates(self, leaflet, code):
        if leaflet == 'up':
            mask = self.up == code
        elif leaflet == 'dw':
            mask = self.dw == code
        else:
            return np.empty(0), np.empty(0)

        x_coords = self.xx[mask]
        y_coords = self.yy[mask]
        return x_coords, y_coords

    def _binary_mask(self, leaflet):
        return self.up if leaflet == 'up' else self.dw

    def detect_clusters(self, leaflet):
        mask = self._binary_mask(leaflet)
        graph = _make_graph(mask)
        visited = set()
        clusters = []
        for node in graph:
            if node not in visited:
                component = _dfs(graph, node)
                clusters.append(component)
                visited.update(component)
        return clusters

    def cluster_sizes(self, leaflet):
        return [len(c) for c in self.detect_clusters(leaflet)]