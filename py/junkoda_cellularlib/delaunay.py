"""
delaunay
"""
import numpy as np
import scipy.spatial
from matplotlib import collections as mc


class Delaunay:
    """
    delaunay[i]: [x1 y1 x2 x2] of ith edge

    Methods
      edges (np.ndarray): integer array of vertex indices n_edges x 2
      plot()
    """
    def __init__(self, v):
        """
        Args:
          v (np.array): vertices v[:, 0] = x, v[:, 1] = y
        """
        self._v = v
        self._d = scipy.spatial.Delaunay(v)
        self.edges = self._get_edges(self._d.simplices)
        self.edge_coords = self._get_edge_coords(self._v, self.edges)

    def __repr__(self):
        triangles = self._d.simplices
        return 'Delaunay(%d triangles, %d edges)' % (len(triangles),
                                                     len(self.edges))

    def __getitem__(self, idx):
        return self.edge_coords[idx, ]

    @staticmethod
    def _get_edges(triangles):
        """
        Get edges from triangles

        Returns:
          edges (np.array int): n_edges x 2

        Note:
          A triangle is a triplet of vertex index
          An edge is a pair of vertex index
        """
        s = set()

        for t in triangles:
            i, j, k = t
            if not ((i, j) in s or (j, i) in s):
                s.add((i, j))
            if not ((j, k) in s or (k, j) in s):
                s.add((j, k))
            if not ((k, i) in s or (i, j) in s):
                s.add((k, i))

        return np.array(list(s))

    @staticmethod
    def _get_edge_coords(v, edge_indices):
        n = len(edge_indices)
        a = np.empty((n, 4))
        a[:, :2] = v[edge_indices[:, 0]]
        a[:, 2:] = v[edge_indices[:, 1]]
        return a

    def plot(self, **kwargs):
        import matplotlib.pyplot as plt

        v = self._v

        lines = []
        for e in self.edges:
            i, j = e
            lines.append(((v[i, 0], v[i, 1]), (v[j, 0], v[j, 1])))

        lc = mc.LineCollection(lines, **kwargs)
        ax = plt.gca()
        ax.add_collection(lc)

        if ax.get_xlim() == (0, 1):
            ax.autoscale()

        return self
