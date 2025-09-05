# packing_defect/core/cluster.py

import numpy as np


class DefectClustering:
    """Connected-component analysis for defect masks.

    Methods operate on 2D integer masks, typically produced from
    ``DefectGrid``. Non-zero entries are treated as part of a defect.
    """

    @staticmethod
    def defect_size(matrices, nbins, bin_max, fname, prob=True):
        """Compute and write a histogram of defect cluster sizes.

        Parameters
        ----------
        matrices : sequence[np.ndarray]
            Iterable of 2D binary masks (1 for defect cell, 0 otherwise).
        nbins : int
            Number of histogram bins.
        bin_max : float
            Upper edge of the histogram range.
        fname : str
            Output path to write two-column text: center, value.
        prob : bool, optional
            If True, normalize counts to probabilities.
        """
        bins = np.linspace(0, bin_max, nbins)
        defects = []
        for matrix in matrices:
            graph = DefectClustering._make_graph(matrix)
            visited = set()
            for node in graph:
                if node not in visited:
                    cluster = DefectClustering._dfs(graph, node)
                    visited |= cluster
                    defects.append(len(cluster))

        hist, _ = np.histogram(defects, bins)
        hist = hist.astype(np.float64)
        binp = 0.5 * (_[1:] + _[:-1])
        if np.sum(hist) == 0:
            return

        if prob:
            hist /= np.sum(hist)

        np.savetxt(fname, np.column_stack((binp, hist)), fmt="%8.5f")


    @staticmethod
    def cluster_sizes_from_mask(matrix: np.ndarray) -> list[int]:
        """Return sizes of all connected components in a binary mask.

        Parameters
        ----------
        matrix : np.ndarray
            2D array with 1 for defect cells and 0 otherwise.

        Returns
        -------
        list[int]
            Sizes (in cells) of each connected component.
        """
        graph = DefectClustering._make_graph(matrix)
        visited, sizes = set(), []
        for node in graph:
            if node not in visited:
                comp = DefectClustering._dfs(graph, node)
                visited |= comp
                sizes.append(len(comp))
        return sizes


    @staticmethod
    def _dfs(graph, start):
        """Depth-first search on an adjacency graph from ``start``."""
        visited, stack = set(), [start]
        while stack:
            node = stack.pop()
            if node not in visited:
                visited.add(node)
                stack.extend(graph[node] - visited)
        return visited


    @staticmethod
    def _make_graph(matrix):
        """Build an 8-connected adjacency graph from a binary mask.

        Nodes are linear indices; neighbors wrap around the domain (PBC).
        """
        graph = {}
        rows, cols = matrix.shape
        for (x, y), val in np.ndenumerate(matrix):
            if val == 0:
                continue
            idx = x * cols + y
            neighbors = set()
            for dx in [-1, 0, 1]:
                for dy in [-1, 0, 1]:
                    nx = (x + dx) % rows
                    ny = (y + dy) % cols
                    if matrix[nx, ny] == 1:
                        nidx = nx * cols + ny
                        neighbors.add(nidx)
            neighbors.discard(idx)
            graph[idx] = neighbors
        return graph
