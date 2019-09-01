#include "Python.h"
#include "np_array.h"
#include "py_clusters.h"
#include "ellipses.h"
#include "py_watershed.h"
#include "watershed_ncluster.h"
#include "watershed_nuclei.h"

//
// List of all functions callable from Python
//

static PyMethodDef methods[] = {
  {"_watershed_alloc", py_watershed_alloc, METH_VARARGS,
   "_watershed_alloc()"},
  {"_watershed_construct",  py_watershed_construct, METH_VARARGS,
   "_watershed_construct(_watershed, img, argsort, threshold)"},
  {"_watershed_get_edges", py_watershed_get_edges, METH_VARARGS,
   "_watershed_get_edges(_watershed)"},
  {"_watershed_get_edge_values", py_watershed_get_edge_values, METH_VARARGS,
   "_watershed_get_edge_values(_watershed)"},
  {"_watershed_obtain_cluster_sizes", py_watershed_obtain_cluster_sizes,
   METH_VARARGS, "_watershed_obtain_cluster_size(_watershed,"
   "pixel_threshold, size_threshold)"},
  {"_watershed_obtain_clusters", py_watershed_obtain_clusters, METH_VARARGS,
   "_watershed_obtain_clusters(_watershed, pixel_threshold, size_threshold)"},

  {"_clusters_alloc", py_clusters_alloc, METH_VARARGS,
   "_clusters_alloc()"},
  {"_clusters_len", py_clusters_len, METH_VARARGS,
   "_clusters_len(_clusters)"},
  {"_clusters_get_cluster", py_clusters_get_cluster, METH_VARARGS,
   "_clusters_get_cluster(_clusters, i)"},
  {"_clusters_obtain", py_clusters_obtain, METH_VARARGS,
   "_clusters_obtain(_clusters, img, pixel_threshold, size_threshold)"},
  {"_clusters_get_sizes", py_clusters_get_sizes, METH_VARARGS,
   "_clusters_get_sizes(_clusters, sizes)"},
  {"_clusters_cluster_nvertices", py_clusters_cluster_nvertices, METH_VARARGS,
   "_clusters_cluster_nvertices(_cluster)"},
  {"_clusters_cluster_nedges", py_clusters_cluster_nedges, METH_VARARGS,
   "_clusters_cluster_nedges(_cluster)"},
  {"_clusters_cluster_get_edges", py_clusters_cluster_get_edges, METH_VARARGS,
   "_clusters_cluster_get_edges(_cluster)"},
  {"_clusters_cluster_get_edge_values", py_clusters_cluster_get_edge_values,
   METH_VARARGS, "_clusters_cluster_get_edge_values(_cluster)"},

  {"_watershed_ncluster_compute", watershed_ncluster::py_compute, METH_VARARGS,
   "_watershed_ncluster_compute()"},
  {"_watershed_nuclei_obtain", watershed_nuclei::obtain, METH_VARARGS,
   "_watershed_nuclei_obtain()"},
  
  {"_ellipses_obtain", ellipses::obtain, METH_VARARGS,
   "_ellipses_obtain(img, pixel_threshold, size_threshold)"},
  
  {NULL, NULL, 0, NULL}
};


static struct PyModuleDef module = {
  PyModuleDef_HEAD_INIT,
  "_junkoda_cellularlib", // name of this module
  "A module for kaggle cellular", // Doc String
  -1,
  methods
};

PyMODINIT_FUNC
PyInit__cellularlib(void) {
  np_array::module_init();
  
  return PyModule_Create(&module);
}
