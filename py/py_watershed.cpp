//
// Watershed algorithm of cluster identification
//
#include <iostream>
#include <vector>
#include <queue>
#include <random>
#include <cmath>
#include <cassert>

#include "buffer.h"
#include "np_array.h"
#include "graph.h"
#include "py_clusters.h"
#include "py_watershed.h"

//using namespace std;
using std::vector;

//
// C++ structure/class
//

class Watershed {
public:
  Watershed();
  Watershed(Watershed const&) = delete;
  Watershed& operator=(Watershed const&) = delete;

  void construct_graph(PyObject * const py_img,
                       PyObject * const py_argsort,
                       const double pixel_threshold,
                       const int merge_threshold,
		       const int seed_random_direction);
 
  vector<int>& obtain_cluster_sizes(const double pixel_threshold,
                                    const int size_threshold);

  vector<int> v_sizes;
  int _nx, _ny;

  std::shared_ptr<vector<Vertex>> ptr_pixels;
  std::shared_ptr<vector<Edge>> ptr_edges;
};


//
// static functions
//

static void obtain_clusters(const vector<Vertex>& v_pixel,
                            const vector<Edge>& v_edge,
                            const double pixel_threshold,
                            const double edge_threshold,
                            const size_t size_threshold,
                            Clusters& clusters);


// Python deconstructor
static void py_watershed_free(PyObject *obj);




//
// Watershed members
//

Watershed::Watershed() :
  _nx(0), _ny(0),
  ptr_pixels(new vector<Vertex>()),
  ptr_edges(new vector<Edge>())
{

}


void Watershed::construct_graph(PyObject * const py_img,
                                PyObject * const py_argsort,
                                const double pixel_threshold,
                                const int merge_threshold,
				const int seed_random_direction)
{
  // Args:
  //   py_img (2D array float64): 2D image array
  //   py_argsort (1D array int):  argsort indices
  //   pixel_threshold: pixel value < are neglected
  //   merge_threshold: if two clusters have sizes >= merge_threshold,
  //                    they are not merged to one cluster
  //   seed_first_direction: if > 0, select first neighbor randomly
  //   
  // Exceptions:
  //   TypeError

  vector<Vertex>& v = *ptr_pixels;
  vector<Edge>& v_edge = *ptr_edges;
  
  v.clear();
  v_edge.clear();


  // Buffer may throw Type error
  Buffer<double> buf_img(py_img, "py_img");         // image/2D pixels;
  Buffer<long>   buf_arg(py_argsort, "py_argsort"); // sorted order

  // image size
  const int nx = _nx = static_cast<int>(buf_img.shape[0]);
  const int ny = _ny = static_cast<int>(buf_img.shape[1]);
  const int n = static_cast<int>(buf_arg.shape[0]);
  assert(nx*ny == n);


  // Define neighbour pixel
  const int dx_list[] = {0, 1, 0, -1}; // up, right, down, left
  const int dy_list[] = {1, 0, -1, 0};

  // Copy img to vector<Vertex>
  v = graph::obtain_vertices(buf_img);

  // randomly select first neighbour if seed_random_direction > 0
  std::mt19937 mt(seed_random_direction);
  std::uniform_int_distribution<int> rand4(0, 3);  // generates 0, 1, 2, 3

  int n_edges = 0;
  
  // Loop over all pixel from that with largest value to lower
  // The `water level` is going down
  for(int i=n-1; i>=0; --i) {
    // <1> is the lowest land above the water level now
    int index1 = buf_arg[i];
    int ix1 = index1 / ny;
    int iy1 = index1 % ny;
    double f1 = buf_img(ix1, iy1);

    assert(0 <= index1 && index1 < n); // DEBUG!

    if(f1 < pixel_threshold)
      break;

    
    // If this pixel does not link to neighbour pixels,
    // the link points to itself
    v[index1].next = index1;
    v[index1].size = 1;

    // First neibour direction (Always 0 if seed == 0)
    int random_direction = seed_random_direction == 0 ?  0 : rand4(mt);

    int another_top = -1;

    // 4 directions  // up, right, down, left from <1>
    for(int j1=0; j1<4; ++j1) {
      // index of the neighbor cell
      int inbr = (random_direction + j1) % 4;

      // 2D index of neighbour cell
      int ix2 = ix1 + dx_list[inbr];
      int iy2 = iy1 + dy_list[inbr];
      
      if(!(0 <= ix2 && ix2 < nx && 0 <= iy2 && iy2 < ny))
        continue;  // Outside the image
      
      int index2 = ix2*ny + iy2;
        
      if(v[index2].next < 0)
        continue;  // Not obove waterlevel yet.
      
      // <2> is a neighbour above water level, higher than <1>
      // by construction

      // Find top pixel of this neighbor
      int top = graph::get_top(index2, v);
      int j2 = (j1 + 2) % 4;  // direction viewed from <2>

      if(another_top == -1) {
        // This pixel joins this first cluster
        v[index1].next = top;
        v[top].size += 1;
        another_top = top;
        
        assert(v[index1].edge[j1] == -1); // DEBUG!
        assert(v[index2].edge[j2] == -1);
            
        v[index1].edge[j1] = n_edges;
        v[index2].edge[j2] = n_edges;

        // add edge
        v_edge.push_back(Edge(index1, index2, f1));
        n_edges++;
      }
      else if(another_top >= 0 && another_top != top) {
        // New cluster is connected to the `another` existing cluster

        // Do not merge two large clusters
        if(v[top].size >= merge_threshold &&
           v[another_top].size >= merge_threshold)
          continue;

        if(v[top].value > v[another_top].value) {
          // top > another_top
          //   another_top is connected under top,
          //   and top become the new another_top 
          v[another_top].next = top;
          v[top].size += v[another_top].size;
          another_top = top;
        }
        else {
          // another_top >= top
          //   top is connected under another_top
          v[top].next = another_top;
          v[another_top].size += v[top].size;
        }

        // vertex -> edge information
        assert(v[index1].edge[j1] == -1);
        assert(v[index2].edge[j2] == -1);
        
        v[index1].edge[j1] = n_edges;
        v[index2].edge[j2] = n_edges;

        // add edge
        v_edge.push_back(Edge(index1, index2, f1));
        n_edges++;
      }
    }
  }      
}


vector<int>& Watershed::obtain_cluster_sizes(
             const double pixel_threshold, const int size_threshold)
{
  v_sizes.clear(); // or create a new vector??
  const int n = _nx*_ny;

  vector<Vertex>& v = *ptr_pixels;
  
  for(int i=0; i<n; ++i) {
    const Vertex& p = v[i];
    if(p.next == i && p.value >= pixel_threshold && p.size >= size_threshold)
      v_sizes.push_back(p.size);
  }
  
  return v_sizes;
}

//
// Clusters
//
// DEBUG!!!
void print_edge(char msg[], const Edge& ee)
{
  int x0 = ee.index[0] / 64;
  int y0 = ee.index[0] % 64;
  int x1 = ee.index[1] / 64;
  int y1 = ee.index[1] % 64;
  
  printf("%-16s  %2d %2d - %2d %2d\n", msg, x0, y0, x1, y1);
}


// Find clusters in a graph
//   v_pixel: array of verticis
//   v_edge:  array of edges
void obtain_clusters(const vector<Vertex>& v_pixel,
                     const vector<Edge>& v_edge,
                     const double pixel_threshold,
                     const double edge_threshold,
                     const size_t size_threshold,
                     Clusters& clusters)
{
  // Thresholds
  //   pixels < pixel_threshold are neglected
  //   edges < edge_threshold are negelected
  //   clusters with sizez >= size_threshold are in result
  if(v_edge.size() == 0)
    return;
  
  const int n_edges = v_edge.size();
  const int img_size = v_pixel.size();

  // indices of exlpred edges
  vector<bool> edge_explored(n_edges, false);
  vector<bool> pixel_explored(img_size, false);

  // edges to be explored
  std::queue<int> q;

  // empty cluster
  Cluster c_init;
  
  // traverse all edges
  for(int i_new_edge=0; i_new_edge<n_edges; ++i_new_edge) {
    // skip if the edge is alreay explored or
    // both endpoint pixels are already explored
    if(edge_explored[i_new_edge] ||
       (pixel_explored[v_edge[i_new_edge].index[0]] &&
        pixel_explored[v_edge[i_new_edge].index[1]]))
      continue;

    // First edge in this new cluster
    assert(q.empty());
    q.push(i_new_edge);

    //print_edge("New edge", v_edge[i_new_edge]);

    // A new cluster
    clusters.push_back(c_init);
    Cluster& c = clusters.back();

    while(!q.empty()) {
      // Pick up an unexplored edge in this cluster
      const int j_edge = q.front(); q.pop();
      const Edge edge = v_edge[j_edge];
      assert(edge_explored[j_edge] == false); // ERROR!!! DEBUG!!!
      edge_explored[j_edge] = true;

      if(edge.value < edge_threshold)
        continue;
      
      c.edges.push_back(edge);

      //print_edge("  pop edge", v_edge[j_edge]);

      for(int k=0; k<2; ++k) {
        // for each end-point pixel
        int index = edge.index[k];
        if(!pixel_explored[index]) {
          // Add a new pixel to the cluster
          if(v_pixel[index].value >= pixel_threshold)
            c.pixels.push_back(index);
          pixel_explored[index] = true;
          
          // Add the adjacent edges to the queue
          const Vertex& p = v_pixel[index];
          for(int j=0; j<4; ++j) { // loop over 4 neighbours
            int adj_edge = p.edge[j];
            if(adj_edge >= 0 && (!edge_explored[adj_edge])) {
              q.push(adj_edge);
            }
          }
        }
      }
    } // all connected vertices are added to the cluster

    // Only keep cluster with size >= size_threshold
    const size_t cluster_size = c.pixels.size();
    if(0 == cluster_size || cluster_size < size_threshold) {
      clusters.pop_back();
    }
  } // all edges explored
}


//
// Python interface
//


//
// Watershed Python interface
//
PyObject* py_watershed_alloc(PyObject* self, PyObject* args)
{
  // _watershed_alloc()
  // Create a new watershed object
  Watershed* const w = new Watershed();

  return PyCapsule_New(w, "_Watershed", py_watershed_free);
}

void py_watershed_free(PyObject *obj)
{
  // Delete Watershed object, called automatically by Python
  Watershed* const w = (Watershed*) PyCapsule_GetPointer(obj, "_Watershed");
  assert(w);

  delete w;
}


PyObject* py_watershed_construct(PyObject* self, PyObject* args)
{
  // _watershed_get_edges(_watershed)
  PyObject *py_watershed, *py_img, *py_argsort;
  double pixel_threshold;
  int merge_threshold;
  int seed_random_direction;
  if(!PyArg_ParseTuple(args, "OOOdii", &py_watershed, &py_img, &py_argsort,
                       &pixel_threshold, &merge_threshold,
		       &seed_random_direction)) {
    return NULL;
  }

  Watershed* const w =
    (Watershed*) PyCapsule_GetPointer(py_watershed, "_Watershed");
  assert(w);

  try {
    w->construct_graph(py_img, py_argsort,
                       pixel_threshold, merge_threshold,
		       seed_random_direction);
  }
  catch (TypeError e) {
    return NULL;
  }

  Py_RETURN_NONE;
}

PyObject* py_watershed_get_edges(PyObject* self, PyObject* args)
{
  // _watershed_get_edges(_watershed)
  PyObject *py_watershed;
  if(!PyArg_ParseTuple(args, "O", &py_watershed)) {
    return NULL;
  }

  Watershed* const w =
    (Watershed*) PyCapsule_GetPointer(py_watershed, "_Watershed");
  assert(w);

  return np_array::view_from_vector_struct(w->ptr_edges->front().index,
                                           w->ptr_edges->size(), 2,
                                           sizeof(Edge));
}


PyObject* py_watershed_get_edge_values(PyObject* self, PyObject* args)
{
  // _watershed_get_edge_values(_watershed)
  PyObject *py_watershed;
  if(!PyArg_ParseTuple(args, "O", &py_watershed)) {
    return NULL;
  }

  Watershed* const w =
    (Watershed*) PyCapsule_GetPointer(py_watershed, "_Watershed");
  assert(w);
  
  return np_array::view_from_vector_struct(&(w->ptr_edges->front().value),
                                           w->ptr_edges->size(), 1,
                                           sizeof(Edge));
}

PyObject* py_watershed_obtain_cluster_sizes(PyObject* self, PyObject* args)
{
  // _watershed_get_edge_values(_watershed)
  PyObject *py_watershed;
  double pixel_threshold;
  int size_threshold;
  if(!PyArg_ParseTuple(args, "Odi",
                       &py_watershed, &pixel_threshold, &size_threshold)) {
    return NULL;
  }

  Watershed* const w =
    (Watershed*) PyCapsule_GetPointer(py_watershed, "_Watershed");
  assert(w);
  
  return np_array::view_from_vector(
           w->obtain_cluster_sizes(pixel_threshold, size_threshold));  
}

PyObject* py_watershed_obtain_clusters(PyObject* self, PyObject* args)
{
  // _watershed_get_edge_values(_watershed)
  PyObject *py_watershed, *py_clusters;
  double pixel_threshold, edge_threshold;
  int size_threshold;
  if(!PyArg_ParseTuple(args, "OddiO", &py_watershed,
                       &pixel_threshold, &edge_threshold,
                       &size_threshold, &py_clusters)) {
    return NULL;
  }

  Watershed const * const w =
    (Watershed const *) PyCapsule_GetPointer(py_watershed, "_Watershed");
  assert(w);

  Clusters* const clusters =
    (Clusters*) PyCapsule_GetPointer(py_clusters, "_Clusters");
  assert(clusters);

  clusters->_nx = w->_nx;
  clusters->_ny = w->_ny;

  obtain_clusters(*w->ptr_pixels, *w->ptr_edges,
                  pixel_threshold, edge_threshold, size_threshold,
                  *clusters);

  //clusters->ptr_pixels = w->ptr_pixels;

  Py_RETURN_NONE;
}

