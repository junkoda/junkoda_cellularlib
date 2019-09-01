import numpy as np
import pandas as pd
import matplotlib.image as mpimg
from typing import Tuple

from .config import _dir


_data_dir = '%s/input' % _dir
cell_types = ['HEPG2', 'HUVEC', 'RPE', 'U2OS']
positive_control = 1108
negative_control = 1138
nsirna = 1108  # excluding 30 positive_control + 1 negative_control
plate_shape = (14, 22)


def set_data_dir(d):
    global _data_dir
    _data_dir = d


def load_header(set_type):
    """
    Args:
      set_type (str): train or test

    id_code        experiment plate well sirna well_type        cell_type
    HEPG2-01_1_B03 HEPG2-01   1     B03  513   posotive_control HEPG2
    """

    df = pd.read_csv('%s/%s.csv' % (_data_dir, set_type))
    df_controls = pd.read_csv('%s/%s_controls.csv' % (_data_dir, set_type))

    if set_type == 'train':
        df.insert(5, 'well_type', 'train')
    else:
        df.insert(4, 'sirna', -1)
        df.insert(5, 'well_type', 'test')

    df_all = pd.concat([df, df_controls], sort=False)

    # Add cell_type; HPEG2-01 -> HPEG2
    df_all['cell_type'] = df_all['experiment'].str.replace('-.*', '',
                                                           regex=True)

    # Add cell_type: HPEG2-01 -> 1
    ex_no = df_all['experiment'].str.replace('^.*-', '', regex=True)

    df_all['experiment_no'] = ex_no.astype(int)

    return df_all


def load(set_type, id_code, site):
    """
    Load image specified by the id_code

    Args:
      set_type (str): train or test
      id_code (str): U2OS-03_4_O19  (or, can be a pd.Series with id_code)
      site (int): 1 or 2

    Returns:
      img (np.array): 6 x 512 x 512

    """

    if not isinstance(id_code, str):
        if isinstance(id_code, pd.Series):
            id_code = id_code.id_code
        else:
            raise TypeError('id_code is not str', id_code)

    if not (site == 1 or site == 2):
        raise ValueError('site must be 1 or 2')

    # Example:
    #  id_code:   U2OS-03_4_O19
    #  filename:  U2OS-03/Plate4/O19_s<site>_w<channel>.png
    v = id_code.split('_')
    batch = v[0]  # U2OS-03
    plate = v[1]  # 4
    well = v[2]   # O19

    nc = 512

    X = np.empty((6, nc, nc))
    for ichannel in range(6):
        filename = ('%s/%s/%s/Plate%s/%s_s%d_w%d.png' % (_data_dir,
                    set_type, batch, plate, well, site, ichannel + 1))

        img = mpimg.imread(filename)
        X[ichannel, :, :] = img[:, :]

    return X


def well_coordinate(well: str) -> Tuple[int]:
    """
    Row: B  -> 0
    Col: 02 -> 0

    Args:
      well (str): e.g. B02

    Returns
      (row, col)
    """

    assert(len(well) == 3)

    row = ord(well[0].lower()) - 96
    col = int(well[1:])

    return (row, col)


def well_index(well: str) -> int:
    row, col = well_coordinate(well)
    return row * 22 + col
