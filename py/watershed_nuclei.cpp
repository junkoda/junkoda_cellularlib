/*
 * Nuclei detection according to cluster size
*/

#include <vector>
#include <deque>
#include <set>
#include <chrono>
#include "buffer.h"
#include "watershed_nuclei.h"

using std::vector;
using std::deque;


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


static void merge_pixels(vector<deque<int>>& v_pixels,
                         const int index1, const int index2)
{
  // merge 2 -> 1
  deque<int>& q1 = v_pixels[index1];
  deque<int>& q2 = v_pixels[index2];

  q1.insert(q1.begin(), q2.begin(), q2.end());
  q2.clear();
}


static void mark_pixels(const deque<int>& q,
                        Buffer<bool>& buf_nuclei)
{
  for(int i : q) {
    buf_nuclei(i) = true;
  }
}


//
// Main data analysis
//
static double mark_nuclei(PyObject * const py_img,
                          PyObject * const py_argsort,
                          PyObject * const py_thresholds,
                          const size_t size_min,
                          const size_t size_max,
                          PyObject * const py_nuclei)
{
  /*
   * Args:
   *   py_img (2D array float64): 2D image array
   *   py_argsort (1D array int): argsort indices
   *   py_thresholds (1D array double0: array of thredholds
   *             cluster sizes are evaluated for each threshold
   *   size_min, size_max (int); size range of nuclei
   *   py_nuclei (2D array int):  [output] pixel is in nuclei or not (output)
   *   
   * Exceptions:
   *   TypeError
   */

  auto ts = std::chrono::high_resolution_clock::now();
  
  // Buffer may throw TypeError
  Buffer<double> buf_img(py_img, "py_img");         // image/2D pixels;
  Buffer<long>   buf_arg(py_argsort, "py_argsort"); // sorted order
  Buffer<double> buf_thresholds(py_thresholds, "py_thresholds");
  Buffer<bool>   buf_nuclei(py_nuclei, "py_nuclei");

  assert(buf_img.ndim == 2);
  assert(buf_arg.ndim == 1);
  assert(buf_thresholds.ndim == 1);
  assert(buf_nuclei.ndim == 1);
  assert(size_min <= size_max);

  // image size
  const int nx = static_cast<int>(buf_img.shape[0]);
  const int ny = static_cast<int>(buf_img.shape[1]);

  // number of pixels
  const int n = static_cast<int>(buf_arg.shape[0]);
  assert(nx*ny == n);
  assert(buf_nuclei.shape[0] == buf_arg.shape[0]);

  // Define neighbour direction
  const int dx_list[] = {-1, 1,  0, 0}; // left, right, top, down
  const int dy_list[] = { 0, 0, -1, 1};

  vector<int> v_next(n, -1);  // link list pointing `next` pixel
  vector<deque<int>> v_pixels(n);
  std::set<int> updated_clusters;

  const int n_thresholds = static_cast<int>(buf_thresholds.shape[0]);
  
  // Loop over all pixels from the largest pixel birghtness to lower
  // The `water level` is going down
  int i = n-1;
  for(int i_threshold=0; i_threshold < n_thresholds; ++i_threshold) {
    double pixel_threshold = buf_thresholds(i_threshold);
    updated_clusters.clear();

    // Find clusters for pixels with value >= pixel_threshold
    while(i >= 0) {
      // <1> is the lowest land above the water level now
      int index1 = buf_arg[i];  assert(0 <= index1 && index1 < n);
      int ix1 = index1 / ny;
      int iy1 = index1 % ny;
      double f1 = buf_img(ix1, iy1);

      if(f1 < pixel_threshold)
        break;

      --i;

      // If this pixel does not link to neighbour pixels,
      // the link points to itself
      v_next[index1] = index1;

      int the_cluster = -1;  // the cluster this pixel belongs to

      // 4 directions  // left, right, up, down
      for(int j1=0; j1<4; ++j1) {
        // 2D index of neighbour cell
        int ix2 = ix1 + dx_list[j1];
        int iy2 = iy1 + dy_list[j1];
        
        if(!(0 <= ix2 && ix2 < nx && 0 <= iy2 && iy2 < ny))
          continue;  // Outside the image
        
        int index2 = ix2*ny + iy2;
        
        if(v_next[index2] < 0)
          continue;  // This neighbour is not obove waterlevel yet.
        
        // <2> is a neighbour above water level, higher than <1>
        // by construction
        
        // The cluster that neighbor <2> belogs to.
        // A cluster is a connected component above the waterlevel
        int nbr_cluster = get_top(index2, v_next);
        assert(nbr_cluster >= 0);
        
        if(the_cluster == -1) {
          // This is the first cluster that this pixel meets
          // This pixel joins this cluster
          the_cluster = nbr_cluster;
          v_next[index1] = the_cluster;
          v_pixels[the_cluster].push_back(index1);
          size_t s1 = v_pixels[the_cluster].size();

          if(size_min <= s1 && s1 <= size_max)
            updated_clusters.insert(the_cluster);
        }
        else if(the_cluster >= 0 && the_cluster != nbr_cluster) {
          // This pixel is a bridge between the_cluster and the other nbr_cluster
          // This pixel is already a member of the_cluster
          // The nbr_cluster will be merged with the_cluster
          v_next[nbr_cluster] = the_cluster;
          
          // sizes of two clusters
          size_t s1 = v_pixels[the_cluster].size();
          size_t s2 = v_pixels[nbr_cluster].size();
          size_t s = s1 + s2;
          
          merge_pixels(v_pixels, the_cluster, nbr_cluster);          
          assert(v_pixels[nbr_cluster].empty());

          if(size_min <= s && s < size_max)
            updated_clusters.insert(the_cluster);
        }
      } // loop for 4 neighbor pixels
      
      
      if(v_next[index1] == index1) {
        // This is a new cluster with this pixel only
        v_pixels[index1].push_back(index1);
      }
    } // endo of loop over pixels >= pixel_threshod
    
    // update nuclei mask
    for(int c : updated_clusters) {
      size_t s = v_pixels[c].size();
      if(size_min <= s && s <= size_max)
        mark_pixels(v_pixels[c], buf_nuclei);
    }    
  } // end of loop over all thresholds

  auto te = std::chrono::high_resolution_clock::now();
  return std::chrono::duration<double>(te - ts).count();
}


//
// Python interface
//

namespace watershed_nuclei {

PyObject* obtain(PyObject* self, PyObject* args)
{
  // Returns:
  //   t (double): computation time [sec]
  // Exception
  //   TypeError
  PyObject *py_img, *py_argsort, *py_thresholds, *py_out;
  int size_min, size_max;
  if(!PyArg_ParseTuple(args, "OOOiiO",
                       &py_img, &py_argsort, &py_thresholds,
                       &size_min, &size_max, &py_out)) {
    return NULL;
  }


  try {
    double t = mark_nuclei(py_img, py_argsort, py_thresholds,
                           size_min, size_max, py_out);
    return Py_BuildValue("d", t);
  }
  catch (TypeError e) {
    return NULL;
  }

  Py_RETURN_NONE;
}

}
