"""
Load appropriate tqdm Jupyter Notebook or terminal
"""

import sys

if 'ipykernel' in sys.modules:
    # Jupyter notebook
    from tqdm import tqdm_notebook as tqdm
else:
    # Terminal
    from tqdm import tqdm

__all__ = ['tqdm']
