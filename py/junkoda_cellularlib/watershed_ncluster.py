"""
Class for measuring number of clusters for an array of thresholds
"""

import numpy as np
import numbers
from . import _cellularlib as c


def compute_nclusters(img, thresholds=None, *,
                      size_threshold=0, seed_random_direction=0):
    """
    Compute the number of clusters for given array of thresholds

    Args:
      img (array): 2D array of float64
      thresholds (array): 1D array of thresholds (float64)
      size_threshold (int): count clusters larger or equal than this number
      seed_random_direction (int): introduce randomness in neighbour
                                   selection; no randomness with 0

    Retuns: d (dict)
      d['thresholds']: array of thresholds (sorted)
      d['nclusters']:  number of clusters for the threshold at same index

    Exception:
      TypeError
    """

    # Checks
    if img.ndim != 2:
        raise TypeError('Expeceted a 2-dimensional array for img: '
                        '%d' % img.ndim)

    # convert thresholds to np.array if necessary
    # `threshold` can be a real number or list/tuple of numbers
    if thresholds is None:
        thresholds = (0.5 + np.arange(255)) / 256
    elif isinstance(thresholds, numbers.Real):
        thresholds = np.array([float(thresholds), ])  # one number
    elif not isinstance(thresholds, np.ndarray):
        thresholds = np.array(thresholds)             # e.g., list, tuple
    else:
        thresholds = thresholds.copy()

    thresholds.sort()
    thresholds = thresholds[::-1]

    if thresholds.ndim != 1:
        raise TypeError('Expected an 1-dimensional array for tresholds: '
                        '%d' % thresholds.ndim)

    # TODO: merge threshold can be an option
    # if merge_threshold < 0:
    #    merge_threshold = img.shape[0] * img.shape[1] + 1
    # else:
    #    merge_threshold = int(merge_threshold)
    #

    img1D = img.flatten()
    arg_sort = np.argsort(img1D)

    # results
    nclusters = np.zeros(len(thresholds), dtype=int)

    c._watershed_ncluster_compute(img, arg_sort,
                                  thresholds, nclusters,
                                  size_threshold, seed_random_direction)

    return thresholds, nclusters
