import os
from .cellularroot import _dir

if 'CELLULAR_ROOT' in os.environ:
    _dir = os.environ['CELLULAR_ROOT']
