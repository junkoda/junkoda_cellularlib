#ifndef GRAPH_H
#define GRAPH_H 1

#include <vector>
#include "buffer.h"


// Vertex
struct Vertex {
  double value;   // pixel value
  int next;       // pointing 'next' pixel with the same group
  int size;       // size of the cluster
  int edge[4];    // edge indices, -1 for no edge
};

struct Edge {
  Edge() : index{-1, -1}, value(0.0) {}
  Edge(const int i1, const int i2, const double val) :
    index{i1, i2}, value(val) {
  }
  int index[2];
  double value;
};

namespace graph {

//
// Functions
//
std::vector<Vertex> obtain_vertices(const Buffer<double>& buf_img,
                                   const int size_init=0);


//
// Inline functions
//

  
// Find the `top` of pixes
static inline int get_top(int i, const std::vector<Vertex>& v)
{
  assert(0 <= i && i < static_cast<int>(v.size())); // DEBUG!

  while(i != v[i].next) {
    i = v[i].next;
    assert(0 <= i && i < static_cast<int>(v.size())); // DEBUG!
  }

  return i;
}

} // namespace graph

#endif
