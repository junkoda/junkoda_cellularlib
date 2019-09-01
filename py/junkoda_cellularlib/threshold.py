"""
Threshold

Functions:
  mean2
"""

import numpy as np
import numbers
from . import _cellularlib as c  # library in C++
from .watershed_ncluster import compute_nclusters


def median_quarter_maximum(img):
    """
    Args:
      img (np.array): 2D array of image, assuming value in [0, 1]

    Median quatre maximum threshold
      median(threshold) for threshold > 0.25*max(ncluster),
    where,
      ncluster(threshold): number of clusters as a function of threshold

    Exception:
      RuntimeError  -- when no cluster exists
    """

    if img.ndim != 2:
        raise TypeError('Expected a 2-dimensional image')

    thresholds, nclusters = compute_nclusters(img, size_threshold=5)

    quarter_maximum = 0.25 * np.max(nclusters)
    idx = nclusters > quarter_maximum

    if np.any(idx):
        return np.median(thresholds[idx])

    raise RuntimeError('No cluster found')


def mean2(img, *, iter=5):
    if img.ndim != 2:
        raise TypeError('Expected a 2-dimensional image')

    thresholds, nclusters = compute_nclusters(img, size_threshold=5)

    m = np.average(thresholds, weights=nclusters)

    for i in range(1, iter):
        idx = thresholds < m
        m1 = np.average(thresholds[idx], weights=nclusters[idx])

        idx = np.logical_not(idx)
        m2 = np.average(thresholds[idx], weights=nclusters[idx])

        m = 0.5 * (m1 + m2)

    return m


def obtain_nuclei_pixels(img, size_min, size_max, *, thresholds=None):
    assert(img.ndim == 2)

    img1D = img.flatten()
    arg_sort = np.argsort(img1D)

    # Prepare thresholds
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

    # Ouput array
    n = len(img1D)
    nuclei = np.zeros(n, dtype=bool)

    c._watershed_nuclei_obtain(img, arg_sort, thresholds,
                               size_min, size_max, nuclei)

    return nuclei.reshape(img.shape[0], img.shape[1])
