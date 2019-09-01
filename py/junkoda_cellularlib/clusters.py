import numpy as np
import junkoda_cellularlib._cellularlib as c  # library in C++

from .graph import Graph


class Cluster:
    def __init__(self, _cluster, nx, ny):
        self._cluster = _cluster
        self._graph = None
        self.nx = nx
        self.ny = ny

    def __repr__(self):
        return 'Cluster (%d vertices, %d edges)' % (len(self), self.n_edges)

    def __len__(self):
        """
        Number of vertices in the cluster
        """
        return self.n_vertices

    def __lt__(self, other):
        return len(self) > len(other)

    @property
    def n_edges(self):
        return c._clusters_cluster_nedges(self._cluster)

    @property
    def n_vertices(self):
        return c._clusters_cluster_nvertices(self._cluster)

    @property
    def edge_indices(self):
        """
        Returns:
          np.array (int): n_edges x 2
                          [pixel index1, pixel index2] of the endpoints
        """
        return c._clusters_cluster_get_edges(self._cluster)

    @property
    def edge_values(self):
        """
        Returns:
          np.array (float [n_edge,]): the value of lower pixel
        """
        return c._clusters_cluster_get_edge_values(self._cluster)

    def obtain_graph(self):
        """
        Returns: graph (cellularlib.Graph)
        """
        return Graph(self.edge_indices, self.edge_values, self.nx, self.ny)

    def plot_edges(self, colour=None, *, cmap='OrRd', vmin=None, vmax=None,
                   **kwargs):
        """
        Plot the graph of the cluster
        """
        if self._graph is None:
            self._graph = self.obtain_graph()

        self._graph.plot_edges(colour, cmap=cmap, vmin=vmin, vmax=vmax,
                               **kwargs)


class Clusters:
    """
    Clusters(img, pixel_threshold, *, size_threshold=0)

    len(clusters): number of clusters
    clusters[i]: ith cluster

    Methods:
      plot_edges
    """
    def __init__(self, img, pixel_threshold, *, size_threshold=0):
        self._clusters = c._clusters_alloc()

        if img is not None:
            self.obtain(img, pixel_threshold, size_threshold)

    def __len__(self):
        """
        Number of cluster in the clusters
        """
        return c._clusters_len(self._clusters)

    def __repr__(self):
        return 'Clusters %d' % len(self)

    def __getitem__(self, i):
        """
        Returns cellularlib.Cluster

        Args:
          i (int): index of the cluster; can be negative
                   slice not supported

        Exception: raise IndexError when i is out of range
        """
        _cluster, nx, ny = c._clusters_get_cluster(self._clusters, i)
        return Cluster(_cluster, nx, ny)

    def obtain(self, img, pixel_threshold, size_threshold=0):
        if img.ndim != 2:
            raise TypeError('Expeceted a 2-dimensional array for img: '
                            '%d' % img.ndim)

        c._clusters_obtain(self._clusters, img,
                           pixel_threshold, size_threshold)
        return self

    def plot_edges(self, colour=None, *, cmap='OrRd', vmin=None, vmax=None,
                   **kwargs):
        """
        Plot the graph of all clusters
        """
        for cl in self:
            cl.plot_edges(colour, cmap=cmap, vmin=vmin, vmax=vmax, **kwargs)

    @property
    def sizes(self):
        out = np.empty(len(self), dtype=int)
        c._clusters_get_sizes(self._clusters, out)
        return out


def obtain(img, pixel_threshold, size_threshold=0):
    """
    Obtain clusters from image

    cluster is defined as a connected component of pixels >= pixel_threshold

    Args:
      img (array): 2D array of float64
      pixel_threshold (array): 1D array of thresholds (float64)
      size_threshold (int): neglect clusters smaller than this

    Retuns: Clusters

    Exception:
      TypeError
    """

    clusters = Clusters()

    clusters.obtain(img, pixel_threshold, size_threshold)

    return clusters
