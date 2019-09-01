#
# $ make
#

#from distutils.core import setup, Extension
from setuptools import setup, Extension
import numpy as np
import os

# Set default cellularroot dirctory
path, _ = os.path.split(os.path.abspath(__file__))
libroot = os.path.abspath(os.path.join(path, '..', '..'))

with open('junkoda_cellularlib/cellularroot.py', 'w') as f:
    f.write("_dir = '%s'\n" % libroot)

# read patch version
with open('ver', 'r') as f:
    ver = int(f.readline().rstrip())

    

# Remove the "-Wstrict-prototypes" compiler option, which isn't valid for C++.
# https://stackoverflow.com/questions/8106258
import distutils.sysconfig
cfg_vars = distutils.sysconfig.get_config_vars()
for key, value in cfg_vars.items():
    if type(value) == str:
        cfg_vars[key] = value.replace("-Wstrict-prototypes", "")

setup(name='junkoda_cellularlib',
      version='0.0.%d' % ver,
      author='Jun Koda',
      py_modules=['.cellularroot',
                  '.clusters',
                  '.data',
                  '.delaunay',
                  '.ellipses',
                  '.graph',
                  '.watershed',
      ],
      ext_modules=[
          Extension('junkoda_cellularlib._cellularlib',
                    ['py_package.cpp',                     
                     'py_clusters.cpp',
                     'ellipses.cpp',
                     'graph.cpp',
                     'np_array.cpp',
                     'py_watershed.cpp',
                     'watershed_ncluster.cpp',
                     'watershed_nuclei.cpp',
                    ],
                    depends = ['np_array.h',
                               'buffer.h',
                               'py_clusters.h',
                               'ellipses.h',
                               'error.h',
                               'graph.h',
                               'grid.h',
                               'py_util.h',
                               'py_watershed.h',
                               'watershed_ncluster.h',
                               'watershed_nuclei.h',
                    ],
                    extra_compile_args = ['-std=c++11'],
                    include_dirs = [np.get_include(), ],
                    # libraries = ['gsl', 'gslcblas'],
                    undef_macros = ['NDEBUG'],
          )
      ],
      packages=['junkoda_cellularlib'],
      url='https://github.com/junkoda/junkoda_cellularlib',
      license = "MIT",
)
