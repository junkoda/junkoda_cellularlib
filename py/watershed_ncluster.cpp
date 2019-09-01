/*
 * Computes the number of clusters for an array of thresholds
 * Designed to compute nclusters for multiple threshols efficiently.
*/

#include <vector>
#include <random>

#include "buffer.h"
#include "watershed_ncluster.h"

using namespace std;

//
// static functions
//

static inline int get_top(int i, const vector<int>& v_next)
{
  assert(0 <= i && i < static_cast<int>(v_next.size())); // DEBUG!

  while(i != v_next[i]) {
    i = v_next[i];
    assert(0 <= i && i < static_cast<int>(v_next.size())); // DEBUG!
  }

  return i;
}


//
// Main data analysis
//

static void compute_nclusters(PyObject * const py_img,
                              PyObject * const py_argsort,
                              PyObject * const py_thresholds,
                              PyObject * const py_nclusters,
                              const int size_threshold,
                              const int seed_random_direction)
{
  /*
   * Args:
   *   py_img (2D array float64): 2D image array
   *   py_argsort (1D array int):  argsort indices
   *   pixel_threshold: pixel value < are neglected
   *   size_threshold: cluster size < are neglected
   *   seed_first_direction: if > 0, select first neighbor randomly
   *   
   * Exceptions:
   *   TypeError
   */

  // Buffer may throw TypeError
  Buffer<double> buf_img(py_img, "py_img");         // image/2D pixels;
  Buffer<long>   buf_arg(py_argsort, "py_argsort"); // sorted order
  Buffer<double> buf_thresholds(py_thresholds, "py_thresholds");
  Buffer<long>   buf_nclusters(py_nclusters, "py_nclusters"); // result

  assert(buf_img.ndim == 2);
  assert(buf_arg.ndim == 1);
  assert(buf_thresholds.ndim == 1);
  assert(buf_nclusters.ndim == 1);

  // image size
  const int nx = static_cast<int>(buf_img.shape[0]);
  const int ny = static_cast<int>(buf_img.shape[1]);

  // number of pixels
  const int n = static_cast<int>(buf_arg.shape[0]);
  assert(nx*ny == n);

  // thresholds

  const int n_thresholds = static_cast<int>(buf_thresholds.shape[0]);
  int i_threshold = 0; // current threshold in consideration
  
  // Define neighbour direction
  const int dx_list[] = {0, 1, 0, -1}; // up, right, down, left
  const int dy_list[] = {1, 0, -1, 0};

  // Copy img to C++ vector<Vertex>
  //   initial value: vertex.next = -1 and edge[k] = -1
  vector<int> v_next(n, -1);  // link list pointing `next` pixel
  vector<int> v_size(n, 0);   // size of the cluster if this pixel
                              // is a `top` pixel

  // randomly select first neighbour
  // not used if seed_random_direction = 0
  std::mt19937 mt(seed_random_direction);
  std::uniform_int_distribution<int> rand4(0, 3);  // generates 0, 1, 2, 3


  // The result of this function; the number of clusters larger or equal
  // to size_treshold
  int n_clusters = 0;
  
  // Loop over all pixel from the largest value to lower
  // The `water level` is going down
  for(int i=n-1; i>=0; --i) {
    // <1> is the lowest land above the water level now
    int index1 = buf_arg[i];
    int ix1 = index1 / ny;
    int iy1 = index1 % ny;
    double f1 = buf_img(ix1, iy1);

    assert(0 <= index1 && index1 < n); // DEBUG!

    // Set i_threshold
    for(; i_threshold < n_thresholds; ++i_threshold) {
      if(buf_thresholds(i_threshold) <= f1)
        break;
    }
    if(i_threshold == n_thresholds)
      break;

    double pixel_threshold = buf_thresholds(i_threshold);
    assert(pixel_threshold <= f1);

    
    // If this pixel does not link to neighbour pixels,
    // the link points to itself
    v_next[index1] = index1;
    v_size[index1] = 1;

    // First neibour direction (Always 0 if seed == 0)
    int random_direction = seed_random_direction == 0 ?  0 : rand4(mt);
    int the_cluster = -1;  // the cluster this pixel belongs to

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
        
      //if(v[index2].next < 0)
      if(v_next[index2] < 0)
        continue;  // This neighbour is not obove waterlevel yet.
      
      // <2> is a neighbour above water level, higher than <1>
      // by construction

      // The cluster that neighbor <2> belogs to.
      // A cluster is a connected component above the waterlevel
      int nbr_cluster = get_top(index2, v_next);
      //int j2 = (j1 + 2) % 4;  // direction viewed from <2> (180 deg opposite)

      if(the_cluster == -1) {
        // This is the first cluster that this pixel meets
        // This pixel joins this cluster
                
        the_cluster = nbr_cluster;
        v_next[index1] = the_cluster;
        v_size[the_cluster] += 1;
          
        // Just crossed the size threshold
        if(v_size[the_cluster] == size_threshold) {
          ++n_clusters;
        }
      }
      else if(the_cluster >= 0 && the_cluster != nbr_cluster) {
        // This pixel is a bridge between the_cluster and the other nbr_cluster
        // The nbr_cluster will be combined with the_cluster
        v_next[nbr_cluster] = the_cluster;

        // sizes of two clusters
        int s1 = v_size[the_cluster];
        int s2 = v_size[nbr_cluster];

        if(s1 < size_threshold && s2 < size_threshold &&
           s1 + s2 >= size_threshold) {
          // The cluster just crossed the size_threshold
          ++n_clusters;
        }
        else if(s1 >= size_threshold && s2 >= size_threshold) {
          // Two clusters merged to one
          --n_clusters;
        }

        v_size[the_cluster] += s2;      
      }
    } // loop for 4 neighbor pixels

    
    if(v_next[index1] == index1) {
      // This is a new isolated pixel with size == 1
      assert(v_size[index1] == 1);
      if(1 >= size_threshold)
         n_clusters++;
    }

    assert(0 <= i_threshold && i_threshold < n_thresholds);
    buf_nclusters(i_threshold) = n_clusters;
  } // end of loop over all pixels
}


//
// Python interface
//

namespace watershed_ncluster {

PyObject* py_compute(PyObject* self, PyObject* args)
{
  // _watershed_get_edges(img, argsort, nclusters,
  //                      size_threshold, seed_romdom_direction)
  // Exception
  //   TypeError
  PyObject *py_img, *py_argsort, *py_thresholds, *py_ncluster;
  int size_threshold, seed_random_direction;
  if(!PyArg_ParseTuple(args, "OOOOii",
                       &py_img, &py_argsort,
                       &py_thresholds, &py_ncluster,
                       &size_threshold, &seed_random_direction)) {
    return NULL;
  }

  try {
    compute_nclusters(py_img, py_argsort, py_thresholds, py_ncluster,
                      size_threshold, seed_random_direction);
  }
  catch (TypeError e) {
    return NULL;
  }

  Py_RETURN_NONE;
}

}
