#include <vector>
#include <queue>
#include <algorithm> 

#include "np_array.h"
#include "py_clusters.h"

//using namespace std;
using std::vector;
using std::queue;
using std::min;

//
// static functions
//
static void py_clusters_free(PyObject *obj);

//
// Clusters
//
Clusters::Clusters() :
  _nx(0), _ny(0)
{

}


//
// C++ code
//
void Clusters::construct(PyObject* const py_img,
                         const double pixel_threshold,
                         const int size_threshold)
{
  /*
   * Args:
   *   py_img (2D array float64): 2D image array
   *   pixel_threshold: pixel value < are neglected
   *   size_threshold: cluster size < are neglected
   *   
   * Exceptions:
   *   TypeError
   */

  // Remove existing cluster in this clusters
  clear();

  // Buffer may throw Type error
  Buffer<double> buf_img(py_img, "py_img");  // image/2D pixels;
  
  // image size
  const int nx = _nx = static_cast<int>(buf_img.shape[0]);
  const int ny = _ny = static_cast<int>(buf_img.shape[1]);

  // number of pixels
  const int n = nx*ny;

  // Define neighbour pixel
  const int dx_list[] = {1, -1, 0, 0}; // right, left, up, down
  const int dy_list[] = {0,  0, -1, 1};

  // remember visited pixels
  vector<bool> visited(n, false);

  // queue of pixel indices in the same cluster
  queue<int> q;

  // Travers all pixels
  for(int index0=0; index0<n; ++index0) {
    int ix0 = index0 / ny;
    int iy0 = index0 % ny;

    if(visited[index0] || buf_img(ix0, iy0) < pixel_threshold)
      continue;

    visited[index0] = true;
    assert(q.empty());

    // First pixel in a new cluster
    q.push(index0);
    emplace_back();
    Cluster& c = back();

    int sum = 0;
    while(!q.empty()) {
      int index1 = q.front(); q.pop();
      assert(0 <= index1 && index1 < n);

      c.pixels.push_back(index1);
      int ix1 = index1 / ny;
      int iy1 = index1 % ny;
      double f1 = buf_img(ix1, iy1);

      ++sum;
      
      // Loop over 4 neighbors
      for(int j=0; j<4; ++j) {
        int ix2 = ix1 + dx_list[j];
        int iy2 = iy1 + dy_list[j];
        double f2 = buf_img(ix2, iy2);
        
        if(!(0 <= ix2 && ix2 < nx && 0 <= iy2 && iy2 < ny))
          continue;  // Outside the image

        int index2 = ix2*ny + iy2;
        if(visited[index2] || f2 < pixel_threshold)
          continue;

        // Add a connected pixel to the queue
        visited[index2] = true;
        q.push(index2);
        c.edges.emplace_back(index1, index2, min(f1, f2));
      }
    }

    if(sum < size_threshold) {
      pop_back();
      continue;
    }
  } // goto to next pixel for a new cluster
}


//
// Python interface
//
PyObject* py_clusters_alloc(PyObject* self, PyObject* args)
{
  Clusters* const c = new Clusters();

  return PyCapsule_New(c, "_Clusters", py_clusters_free);
}

void py_clusters_free(PyObject *obj)
{
  // Delete object, called automatically by Python
  Clusters* const c = (Clusters*) PyCapsule_GetPointer(obj, "_Clusters");
  assert(c);

  delete c;
}

PyObject* py_clusters_len(PyObject* self, PyObject* args)
{
  PyObject *py_clusters;
  if(!PyArg_ParseTuple(args, "O", &py_clusters)) {
    return NULL;
  }
  
  Clusters const * const c =
    (Clusters const *) PyCapsule_GetPointer(py_clusters, "_Clusters");
  assert(c);
    
  return Py_BuildValue("k", static_cast<unsigned long>(c->size()));
}

PyObject* py_clusters_get_cluster(PyObject* self, PyObject* args)
{
  PyObject *py_clusters;
  int i;
  if(!PyArg_ParseTuple(args, "Oi", &py_clusters, &i)) {
    return NULL;
  }
  
  Clusters * const clusters =
    (Clusters*) PyCapsule_GetPointer(py_clusters, "_Clusters");
  assert(clusters);

  if(i < 0) {
    i = static_cast<int>(clusters->size()) + i;
  }

  
  try {
    Cluster& cluster = clusters->at(i);  // may throw std::out_of_range


    return Py_BuildValue("Oii", PyCapsule_New(&cluster, "_Cluster", NULL),
                       clusters->_nx, clusters->_ny);
  }
  catch(std::out_of_range& err) {
    PyErr_SetNone(PyExc_IndexError);
  }

  return NULL;
}


PyObject* py_clusters_obtain(PyObject* self, PyObject* args)
{
  // _clusters_obtain(img, pixel_threshold, size_threshold)
  PyObject *py_clusters;
  PyObject *py_img;
  double pixel_threshold;
  int size_threshold;
  if(!PyArg_ParseTuple(args, "OOdi", &py_clusters, &py_img,
                       &pixel_threshold, &size_threshold)) {
    return NULL;
  }
  
  Clusters* const c =
    (Clusters*) PyCapsule_GetPointer(py_clusters, "_Clusters");
  assert(c);

  try {
    c->construct(py_img, pixel_threshold, size_threshold);
  }
  catch (TypeError e) {
    return NULL;
  }

  Py_RETURN_NONE;
}


PyObject* py_clusters_get_sizes(PyObject* self, PyObject* args)
{
  // Get sizes of the clusters
  // _clusters_get_sizes(_clusters, cluster_sizes)
  //
  PyObject *py_clusters, *py_sizes;
  if(!PyArg_ParseTuple(args, "OO", &py_clusters, &py_sizes)) {
    return NULL;
  }
  
  Clusters const * const clusters =
    (Clusters const*) PyCapsule_GetPointer(py_clusters, "_Clusters");
  assert(clusters);
  const int n_clusters = static_cast<int>(clusters->size());
  
  Buffer<long> buf_sizes(py_sizes, "py_sizes");  // cluster sizes (output)
  assert(buf_sizes.ndim == 1);
  assert(buf_sizes.shape[0] == clusters->size());
  
  for(int i=0; i<n_clusters; ++i) {
    const Cluster& c = (*clusters)[i];
    buf_sizes[i] = static_cast<long>(c.pixels.size());
  }

  Py_RETURN_NONE;
}





//
// Cluster
//
PyObject* py_clusters_cluster_nvertices(PyObject* self, PyObject* args)
{
  PyObject *py_cluster;
  if(!PyArg_ParseTuple(args, "O", &py_cluster)) {
    return NULL;
  }
  
  Cluster const * const c =
    (Cluster const*) PyCapsule_GetPointer(py_cluster, "_Cluster");
  assert(c);

  return Py_BuildValue("k", static_cast<unsigned long>(c->pixels.size()));  
}


PyObject* py_clusters_cluster_nedges(PyObject* self, PyObject* args)
{
  PyObject *py_cluster;
  if(!PyArg_ParseTuple(args, "O", &py_cluster)) {
    return NULL;
  }
  
  Cluster const * const c =
    (Cluster const*) PyCapsule_GetPointer(py_cluster, "_Cluster");
  assert(c);

  return Py_BuildValue("k", static_cast<unsigned long>(c->pixels.size()));  
}


PyObject* py_clusters_cluster_get_edges(PyObject* self, PyObject* args)
{
  PyObject *py_cluster;
  if(!PyArg_ParseTuple(args, "O", &py_cluster)) {
    return NULL;
  }
  
  Cluster* const c =
    (Cluster *) PyCapsule_GetPointer(py_cluster, "_Cluster");
  assert(c);

  return np_array::view_from_vector_struct(c->edges.front().index,
                                           c->edges.size(), 2,sizeof(Edge));
}


PyObject* py_clusters_cluster_get_edge_values(PyObject* self, PyObject* args)
{
  PyObject *py_cluster;
  if(!PyArg_ParseTuple(args, "O", &py_cluster)) {
    return NULL;
  }
  
  Cluster* const c =
    (Cluster*) PyCapsule_GetPointer(py_cluster, "_Cluster");
  assert(c);
  
  return np_array::view_from_vector_struct(&(c->edges.front().value),
                                           c->edges.size(), 1, sizeof(Edge));
}

