/*
Compute ellipses of clusters

A cluster is a connected component of pixels >= pixel_threshold
Pixels are connected with 4 adjacent neibours
*/

#include <iostream>  // DEBUG!!!
#include <vector>
#include <queue>
#include <random>
#include <utility>  // swap

#include <eigen3/Eigen/Dense>

#include "np_array.h"
#include "buffer.h"
#include "ellipses.h"

using std::vector;

//
// C++ definitions
//

// 5.991
// radius sqrt(5.901)*sigma_x

//
// Static functions
//
//
// C++ implementations
//

static PyObject* obtain_ellipses(PyObject* const py_img,
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

  // Buffer may throw Type error
  Buffer<double> buf_img(py_img, "py_img");  // image/2D pixels;

  // image size
  const int nx = static_cast<int>(buf_img.shape[0]);
  const int ny = static_cast<int>(buf_img.shape[1]);

  // number of pixels
  const int n = nx*ny;

  // Define neighbour pixel
  const int dx_list[] = {1, -1, 0, 0}; // right, left, up, down
  const int dy_list[] = {0,  0, -1, 1};

  // remember visited pixels
  vector<bool> visited(n, false);

  // queue of pixel indices in the same cluster
  std::queue<int> q;

  // Eigen-value solver
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> e;

  // Return value
  vector<double> ellipses;

  constexpr double ellipse_factor = 5.991;
  // a, b = sqrt(5.991*eigen_value)
  // This is the 95% contour for Gaussian


    
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

    // Mean and covariance of the pixels in the cluster
    int sum = 0;
    Eigen::Vector2d mu(0.0, 0.0);
    Eigen::Matrix2d cov; cov << 0.0, 0.0, 0.0, 0.0;

    const Eigen::DiagonalMatrix<double, 2> diag(1.0/12.0, 1.0/12.0);

    while(!q.empty()) {
      int index1 = q.front(); q.pop();
      assert(0 <= index1 && index1 < n);
      
      int ix1 = index1 / ny;
      int iy1 = index1 % ny;

      Eigen::Vector2d x(ix1, iy1);;
      mu += x;
      cov += x*x.transpose();
      ++sum;
      // Note: use value?
      
      // Loop over 4 neighbors
      for(int j=0; j<4; ++j) {
        int ix2 = ix1 + dx_list[j];
        int iy2 = iy1 + dy_list[j];
        
        if(!(0 <= ix2 && ix2 < nx && 0 <= iy2 && iy2 < ny))
          continue;  // Outside the image

        int index2 = ix2*ny + iy2;
        if(visited[index2] || buf_img(ix2, iy2) < pixel_threshold)
          continue;

        // Add a connected pixel to the queue
        visited[index2] = true;
        q.push(index2);
      }      
    }

    if(sum < size_threshold)
      continue;

    
    // Normalise mean and covariance
    mu /= sum;
    cov /= sum;
    cov -= mu*mu.transpose();
    cov += diag; // Add 1/12 to diagonal, variance of the square pixel

    
    // Solve for eigen vectors and eigen values 
    e.compute(cov);
    assert(e.info() == Eigen::Success);

    double a = e.eigenvalues()[0];
    double b = e.eigenvalues()[1];

    // (ex, ey) is the eigen vector along major axis
    double ex, ey;

    if(a >= b) {
      ex = e.eigenvectors().col(0)[0];
      ey = e.eigenvectors().col(0)[1];
    }
    else {
      std::swap(a, b);
      ex = e.eigenvectors().col(1)[0];
      ey = e.eigenvectors().col(1)[1];
    }

    assert(a >= b);

    // angle between x axis and the major axis in radians
    double theta = ex >=0 ? asin(ey) : M_PI - asin(ey);

    a = sqrt(ellipse_factor*a); // semi-major axis
    b = sqrt(ellipse_factor*b); // semi-minor axis;  a >= b

    ellipses.push_back(sum);    
    ellipses.push_back(mu[0]);
    ellipses.push_back(mu[1]);
    ellipses.push_back(a);
    ellipses.push_back(b);
    ellipses.push_back(theta);
  } // goto to next pixel for a new cluster

  return np_array::copy_from_vector(ellipses);
}

//
// Python interface
//

namespace ellipses {

PyObject* obtain(PyObject* self, PyObject* args)
{
  // _ellipses_obtain()
  PyObject *py_img;
  double pixel_threshold;
  int size_threshold;
  if(!PyArg_ParseTuple(args, "Odi", &py_img,
                       &pixel_threshold, &size_threshold)) {
    return NULL;
  }

  try {
    return obtain_ellipses(py_img, pixel_threshold, size_threshold);
  }
  catch (TypeError e) {
    return NULL;
  }

  Py_RETURN_NONE;
}

}
