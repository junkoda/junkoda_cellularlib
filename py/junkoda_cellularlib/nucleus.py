"""
Module nucleus locate clusters of neculei in channel 1
"""

import numpy as np
import math

from scipy import ndimage
from .watershed_ncluster import compute_nclusters
from .ellipses import obtain


def median_quarter_maximum_threshold(img):
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

    thresholds = (0.5 + np.arange(255)) / 256
    thresholds, nclusters = compute_nclusters(img,
                                              thresholds, size_threshold=5)

    quarter_maximum = 0.25 * np.max(nclusters)
    idx = nclusters > quarter_maximum

    if np.any(idx):
        return np.median(thresholds[idx])

    raise RuntimeError('No cluster found')


def obtain_clips(img, n_clips, *, clip_size=128, threshold=None):
    """
    Return n_clips image clips centred on a cluster of nuclei.
    The major axis of the cluster is aligned with the x axis

    Args:
      img (np.array): img[ix, iy, ichannel]
      n_clips (int):  number of maximum random clips in the output

    Returns:
      clips[iclip, ix, iy, ichannel], ellipses[iclip, 3]

      ellipse in original image
        ellipses[:, 0]: x
        ellipses[:, 1]: y
        ellipses[:, 2]: theta, angle between major axis and x axis (degree)

    Excption:
      RuntimeError: when no ellipse is found

    Note:
      The number of clips may be smaller than n_clips
      The clips may overlap, but the centers is different
    """

    if threshold is None:
        threshold = median_quarter_maximum_threshold(img[:, :, 0])

    ellipses = obtain(img[:, :, 0], threshold, size_threshold=5)
    n_ellipses = len(ellipses)
    if n_ellipses == 0:
        raise RuntimeError()

    nx = img.shape[0]
    ny = img.shape[1]
    half_clipsize = clip_size // 2
    half_clipsize2 = math.ceil(1.415 * half_clipsize)  # sqrt(2)*half_clipsize

    # outputs
    img_clips = np.zeros((n_clips, clip_size, clip_size, 6))
    meta_data = np.zeros((n_clips, 3))

    ii = np.arange(n_ellipses)
    np.random.shuffle(ii)

    n = 0
    for i in ii:
        e = ellipses[i, :]
        x = int(e[1])
        y = int(e[2])

        theta = e[5] / math.pi * 180.0

        # The clip is within the image for any rotation
        if ((0 <= x - half_clipsize2 and x + half_clipsize2 < nx - 1) and
            (0 <= y - half_clipsize2 and y + half_clipsize2 < ny - 1)):

            img_clip1 = img[(x - half_clipsize2):(x + half_clipsize2),
                            (y - half_clipsize2):(y + half_clipsize2), :]

            img_rot = ndimage.rotate(img_clip1, -theta, reshape=False)
            mar = half_clipsize2 - half_clipsize
            img_clip2 = img_rot[mar:(mar + clip_size),
                                mar:(mar + clip_size), :]

            img_clips[n, :, :, :] = img_clip2
            meta_data[n, 0] = x
            meta_data[n, 1] = y
            meta_data[n, 2] = theta

            n += 1

        if n == n_clips:
            break

    if n == 0:
        raise RuntimeError()

    if n < n_clips:
        img_clips = img_clips[:n, :, :, :]
        meta_data = meta_data[:n, :]

    return img_clips, meta_data
