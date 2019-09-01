#include "graph.h"

using namespace std;

namespace graph {

// Copy img to vector<Vertex>
vector<Vertex> obtain_vertices(const Buffer<double>& buf_img,
                              const int size_init)
{
  const int nx = static_cast<int>(buf_img.shape[0]);
  const int ny = static_cast<int>(buf_img.shape[1]);

  const int n = nx*ny;

  vector<Vertex> v;
  v.reserve(n);

  Vertex p;
  p.next = -1; // The pixel is 'under the water'
  p.edge[0] = p.edge[1] = p.edge[2] = p.edge[3] = -1;
  p.size = size_init;

  int index=0;
  for(int ix=0; ix<nx; ++ix) {
    for(int iy=0; iy<ny; ++iy) {
      //      p.index = index;
      p.value = buf_img(ix, iy);
      v.push_back(p);
      
      index++;
    }
  }

  return v; // C++11 move
}

} // namespace graph
