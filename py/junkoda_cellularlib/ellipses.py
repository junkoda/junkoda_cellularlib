"""
Cluster finding and ellipse computation
"""

import numpy as np
from . import _cellularlib as c  # library in C++


def obtain(img, pixel_threshold, size_threshold=0):
    """
    Obtain ellipse parameters for clusters

    cluster is defined as a connected component of pixels >= pixel_threshold

    Args:
      img (array): 2D array of float64
      pixel_threshold (array): 1D array of thresholds (float64)
      size_threshold (int): neglect clusters smaller than this

    Retuns: a (np.array)
      a[:, 0]  size: number of pixels in the cluster
      a[:, 1]  x: centre of mass of the cluster
      a[:, 2]  y
      a[:, 3]  a: semi-major axis; a^2 is the eigen value of cov
      a[:, 4]  b: semi-minor axis; b^2 is the eigen value of cov (a >= b)
      a[:, 5]  theta: angle between major axis and x axis in radians [-pi, pi]

    Exception:
      TypeError
    """

    # Checks
    if img.ndim != 2:
        raise TypeError('Expeceted a 2-dimensional array for img: '
                        '%d' % img.ndim)

    es = c._ellipses_obtain(img, float(pixel_threshold), int(size_threshold))
    assert(len(es) % 6 == 0)

    return es.reshape(-1, 6)


def plot(es, img=None, *, axis_factor=2.0, **kwargs):
    """
    ellipses (array): output of ellipses.obtain()
    axis_factor (float): arbitrary factor of expansion
    """

    import matplotlib
    import matplotlib.pyplot as plt

    # default ellipse style
    if 'color' not in kwargs:
        kwargs['color'] = 'red'
    if 'lw' not in kwargs and 'linewidth' not in kwargs:
        kwargs['lw'] = 2

    ax = plt.gca()
    plt.xlabel('$x$')
    plt.ylabel('$y$')

    if img is not None:
        plt.xlim(0, img.shape[1])
        plt.ylim(0, img.shape[0])
        plt.imshow(img.T, origin='lower')

    for e in es:
        ellipse = matplotlib.patches.Ellipse(
            (e[1], e[2]),
            axis_factor * e[3], axis_factor * e[4],
            np.degrees(e[5]),
            **kwargs, fill=False)
        ax.add_patch(ellipse)

    if img is None and ax.get_xlim() == (0, 1):
        ax.autoscale()
