#ifndef PY_WATERSHED_H
#define PY_WATERSHED_H 1

#include "Python.h"

PyObject* py_watershed_alloc(PyObject* self, PyObject* args);
PyObject* py_watershed_construct(PyObject* self, PyObject* args);
PyObject* py_watershed_get_edges(PyObject* self, PyObject* args);
PyObject* py_watershed_get_edge_values(PyObject* self, PyObject* args);
PyObject* py_watershed_obtain_cluster_sizes(PyObject* self, PyObject* args);
PyObject* py_watershed_obtain_clusters(PyObject* self, PyObject* args);
#endif
