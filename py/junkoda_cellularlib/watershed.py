"""
Class for watershed analysis
"""

import numpy as np
import pandas as pd
from . import _cellularlib as c
from .clusters import Clusters
from .graph import Graph


class Watershed:
    """
    Watershed(img=None, pixel_theshold=0.0, merge_threshold=-1)

    Args:
      img (array):             2D array of float64
      pixel_threshold (float): construct graph for pixels above
      merge_threshold (int):   do not merge two large clusters above this size
                               no such threshold if -1
      seed_random_direction (int): use random first edge while graph
                               contruction with this seed; no randomness if 0.

    Methods:
      edges

      plot.edges(idx=None, *, color='black', cmap=None, vmin, vmax)
      obtain_graph()
      obtain_clusters(pixel_threshold=0.0,
                      edge_threshold=None,
                      size_threshold=0)
    """
    def __init__(self, img=None, pixel_threshold=0.0, *,
                 merge_threshold=-1,
                 seed_random_direction=0):
        self._watershed = c._watershed_alloc()
        self.img = None
        self.graph = None

        if img is not None:
            self.construct(img, pixel_threshold, merge_threshold,
                           seed_random_direction)

    def __repr__(self):
        s = 'Watershed'
        if self.img is not None:
            s += '(%d %d)' % (self.img.shape[0], self.img.shape[1])

        s += (', pixel_threshold=%.3f, merge_threshold=%d'
              % (self.pixel_threshold, self.merge_threshold))
        return s

    def construct(self, img, pixel_threshold, merge_threshold,
                  seed_random_direction):
        """
        Construct watershed graph

        Args:
          img (array): 2D array of float64
          pixel_threshold

        """
        self.img = img
        self.pixel_threshold = float(pixel_threshold)

        if merge_threshold < 0:
            self.merge_threshold = img.shape[0] * img.shape[1] + 1
        else:
            self.merge_threshold = int(merge_threshold)

        self.seed_random_direction = int(seed_random_direction)

        img1D = img.flatten()
        arg_sort = np.argsort(img1D)

        c._watershed_construct(self._watershed, img, arg_sort,
                               self.pixel_threshold,
                               self.merge_threshold,
                               self.seed_random_direction)

        return self

    @property
    def edges(self):
        """
        Returns: edges (pd.DataFrame)
        """
        ei = self.edge_indices
        return pd.DataFrame({'index1': ei[:, 0],
                             'index2': ei[:, 1],
                             'value': self.edge_values})

    def obtain_graph(self):
        """
        Returns: graph (cellularlib.Graph)
        """
        nx = self.img.shape[0]
        ny = self.img.shape[1]

        return Graph(self.edge_indices, self.edge_values, nx, ny)

    @property
    def edge_indices(self):
        """
        Returns:
          edge_indices (int array): n_edges x 2
            e[i, 0]: index of lower pixel
            e[i, 1]: index of higher pixel
        """
        if self.img is None:
            raise RuntimeError('Graph is not constructed yet')
        return c._watershed_get_edges(self._watershed)

    @property
    def edge_values(self):
        """
        Returns:
          edge values (float array): length n_edges
            f[i]: value of lower pixel
        """
        if self.img is None:
            raise RuntimeError('Graph is not constructed yet')
        return c._watershed_get_edge_values(self._watershed)

    def cluster_sizes(self, *, pixel_threshold=0.0, size_threshold=0):
        """
        Compute the number of clusters with size >= size_threshold
        pixels below pixel_threshold are neglected

        Args:
          pixel_threshold (float)
          size_threshold (int)

        Returns:
          sizes (array of int): sizes of clusters
        """

        if self.img is None:
            raise RuntimeError('Graph is not constructed yet')
        return c._watershed_obtain_cluster_sizes(self._watershed,
                                                 float(pixel_threshold),
                                                 int(size_threshold))

    def obtain_clusters(self, *,
                        pixel_threshold=0.0,
                        edge_threshold=None,
                        size_threshold=0):
        """
        Find connected components in the watershed grapch

        Args:
          pixel_threshold (float): pixel.value >= is added to vertex
          edge_threshold (float): edge.value >= is used
          size_threshold (int): cluster size >= is added to clusters

        Returns:
          clusters (Clusters)

        Note:
          edge_thresholds affect how pixels are connected, while
          pixel_threshold only affects if that pixel is included to the
          cluster, not affecting the graph connectivity.
        """
        if self.img is None:
            raise RuntimeError('Graph is not constructed yet')

        if edge_threshold is None:
            edge_threshold = pixel_threshold

        clusters = Clusters()
        c._watershed_obtain_clusters(self._watershed,
                                     float(pixel_threshold),
                                     float(edge_threshold),
                                     int(size_threshold),
                                     clusters._clusters)

        return clusters

    def plot_edges(self, **kwargs):
        if self.graph is None:
            self.graph = self.obtain_graph()

        self.graph.plot_edges(**kwargs)
