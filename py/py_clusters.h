#ifndef PY_CLUSTER_H
#define PY_CLUSTER_H

#include <vector>
#include <memory>

#include "Python.h"
#include "graph.h"


struct Cluster {
  std::vector<int> pixels;
  std::vector<Edge> edges;
  int centre[2];
  bool empty() const noexcept {
    return pixels.empty() && edges.empty();
  }
};


class Clusters : public std::vector<Cluster> {
 public:
  Clusters();
  Clusters(Clusters const&) = delete;
  Clusters& operator=(Clusters const&) = delete;

  void construct(PyObject* const py_img,
                 const double pixel_threshold,
                 const int size_threshold);
  
  
  int _nx, _ny;
 private:

};


PyObject* py_clusters_alloc(PyObject* self, PyObject* args);
PyObject* py_clusters_len(PyObject* self, PyObject* args);
PyObject* py_clusters_obtain(PyObject* self, PyObject* args);
PyObject* py_clusters_get_cluster(PyObject* self, PyObject* args);
PyObject* py_clusters_get_sizes(PyObject* self, PyObject* args);
 
PyObject* py_clusters_cluster_nvertices(PyObject* self, PyObject* args);
PyObject* py_clusters_cluster_nedges(PyObject* self, PyObject* args);
PyObject* py_clusters_cluster_get_edges(PyObject* self, PyObject* args);
PyObject* py_clusters_cluster_get_edge_values(PyObject* self, PyObject* args);
#endif
