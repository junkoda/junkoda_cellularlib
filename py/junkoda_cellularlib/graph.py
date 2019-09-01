import numpy as np
import pandas as pd


class Graph:
    """
    Graph struture in a 2D image

    Methods:
      vetices
      edges
      plot_vertices()
      plot_edges()
    """
    def __init__(self, ei, ev, nx, ny):
        """
        ei: edge indices (index = ix*ny + iy)
        ev: edge values

        nx, ny (int)                : image size nx, ny
        """

        # Convert indices to corrdinate (x, y)

        # Vertices
        v = np.array(list(set(ei.flatten())))
        v.sort()

        n_vertices = len(v)
        v_coords = np.empty((n_vertices, 2), dtype=float)
        v_coords[:, 0] = v // ny  # x
        v_coords[:, 1] = v % ny   # y

        # Edges
        n_edges = len(ei)
        edge_coords = np.empty((n_edges, 5), dtype=float)
        edge_coords[:, 0] = ei[:, 0] // ny  # x1
        edge_coords[:, 1] = ei[:, 0] % ny   # y1
        edge_coords[:, 2] = ei[:, 1] // ny  # x2
        edge_coords[:, 3] = ei[:, 1] % ny   # x2
        edge_coords[:, 4] = ev

        self.v = v_coords
        self.e = edge_coords
        self.nx = nx
        self.ny = ny

    def __repr__(self):
        return 'Graph (%d vertices, %d edges)' % (len(self.v), len(self.e))

    def plot_vertices(self, marker='o', **kwargs):
        """
        Args:
          marker (str): 3rd argmument in plt.plot
          **kwargs: any keyword argments are passed to plt.plot
        """
        import matplotlib.pyplot as plt
        plt.xlim(0, self.nx)
        plt.ylim(0, self.ny)
        plt.gca().invert_yaxis()

        plt.plot(self.v[:, 0], self.v[:, 1], marker, **kwargs)

    def plot_edges(self, idx=None, *,
                   color='black',
                   cmap=None, vmin=None, vmax=None,
                   **kwargs):
        """
        Args:
          idx: edge indices, e.g., plot_edges(range(100))
          color (str): color, e.g., str 'black', rgba (0, 0, 0, 0)
          cmap : matplotlib map name (str) or matplotlib.colors.Colormap
          vmin:    minimum value for colour map
          vmax:    maximum value for colour map
          **kwargs: any keyword argmuments are padded to plt.plot

        Note:
          x is ploted on vertical axis and y on horizontal axis to match
          plt.imshow
        """
        import matplotlib.pyplot as plt
        import matplotlib
        from matplotlib import collections as mc

        ax = plt.gca()

        # Change xlim, ylim if they not set yet
        if ax.get_xlim() == (0, 1):
            plt.xlim(0, self.nx)
        if ax.get_ylim() == (0, 1):
            plt.ylim(0, self.ny)

        # Change xlabel, ylabel if they are not set yet
        if ax.xaxis.get_label().get_text() == '':
            plt.xlabel('$x$')
        if ax.yaxis.get_label().get_text() == '':
            plt.ylabel('$y$')

        lines = []

        # norm convers scaler value to colour
        if vmin is None:
            vmin = np.min(self.e[:, 4])
        if vmax is None:
            vmin = np.max(self.e[:, 4])
        norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)

        if idx is None:
            edges = self.e
        else:
            edges = self.e[idx, ]

        for e in edges:
            lines.append([(e[0], e[1]), (e[2], e[3])])

        if 'linewidths' not in kwargs:
            kwargs['linewidths'] = 2

        if cmap is not None:
            # use cmap
            if isinstance(cmap, str):
                # color name
                cmap = matplotlib.cm.get_cmap(cmap)

            colours = []
            for e in edges:
                colours.append(cmap(norm(e[4])))

            lc = mc.LineCollection(lines, colors=colours,
                                   **kwargs)
        else:
            # use colour
            lc = mc.LineCollection(lines, colors=color,
                                   **kwargs)
        plt.gca().add_collection(lc)

    @property
    def edges(self):
        """
        Returns: edges (pd.DataFrame)
                 x1, y1, x2, y2, value
        """
        ec = self.edge_coords
        return pd.DataFrame({'x1': ec[:, 0],
                             'y1': ec[:, 1],
                             'x2': ec[:, 2],
                             'y2': ec[:, 3],
                             'value': ec[:, 4]})

    @property
    def vetices(self):
        """
        Returns: vetices (pd.DataFrame)
                 x, y
        """
        return pd.DataFrame({'x': self.v[:, 0],
                             'y': self.v[:, 1]})
