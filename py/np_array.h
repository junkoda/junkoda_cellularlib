#ifndef ARRAY_H
#define ARRAY_H 1

#include <string>
#include <vector>
#include "Python.h"

namespace np_array {
  
PyMODINIT_FUNC
module_init();

int dtype(const std::string name);

//This is not working
//extern template PyObject* array_from_vector<int>(const std::vector<int>&);

PyObject* copy_from_vector(const std::vector<int>& v);
PyObject* copy_from_vector(const std::vector<long>& v);
PyObject* copy_from_vector(const std::vector<float>& v);
PyObject* copy_from_vector(const std::vector<double>& v);

PyObject* view_from_vector(std::vector<int>& v);
PyObject* view_from_vector(std::vector<long>& v);
PyObject* view_from_vector(std::vector<double>& v);

PyObject* view_from_vector_struct(int* const p,
                      const int n, const int ncol, const size_t sizeof_struct);
PyObject* view_from_vector_struct(long* const p,
                      const int n, const int ncol, const size_t sizeof_struct);
PyObject* view_from_vector_struct(float* const p,
                      const int n, const int ncol, const size_t sizeof_struct);
PyObject* view_from_vector_struct(double* const p,
                      const int n, const int ncol, const size_t sizeof_struct);
} // namespace np_array

#endif
