//
// Create numpy array from C++ vector
//
#include "Python.h"
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "numpy/arrayobject.h"

#include <iostream> // DEBUG!!!

#include <vector>
#include <type_traits>

#include "error.h"
#include "np_array.h"

using std::vector;
using std::string;

namespace {
//
// Return NPY data type
// numpy_dtype<int> => NPY_INT
//
template<typename T>
int dtype()
{
  if (typeid(T) == typeid(int))
    return NPY_INT;
  else if (typeid(T) == typeid(long))
    return NPY_LONG;
  else if (typeid(T) == typeid(float))
    return NPY_FLOAT;
  else if (typeid(T) == typeid(double))
    return NPY_DOUBLE;
  
  throw TypeError();

  return 0;
}


template<typename T> PyObject* copy_from_vector_template(const vector<T>& v)
{
  // Args:
  //   v <T>
  // Returns:
  //   numpy array of dtype T
  npy_intp len = static_cast<npy_intp>(v.size());
  
  PyObject* const py_arr = PyArray_SimpleNew(1, &len, dtype<T>());
  T* const arr =
    reinterpret_cast<T*>(PyArray_DATA((PyArrayObject*) py_arr));

  for(size_t i=0; i<v.size(); ++i)
    arr[i] = v[i];

  return py_arr;
}

  
//
// Create an np.array pointing to exising memeory
// 
template<typename T>
PyObject* view_from_vector_template(vector<T>& v)
{
  // Note: np.array would not work after vector v deallocated
  const int ndim = 1;
  npy_intp len = static_cast<npy_intp>(v.size());

  npy_intp dims[]= {len};

  return PyArray_SimpleNewFromData(ndim, dims, dtype<T>(), v.data());
}


//
// Array view from Vector<struct>
//  T: int/float/double
template<typename T>
PyObject* view_from_vector_struct_template(T* const p,
               const int n, const int ncol, const size_t sizeof_struct)
{
  // Return vector<S> as an np.array
  // Args:
  //   p: pointer to the first data
  //   n: length of the vector
  //   ncol: number of elements <T> to put in array
  //   sizeof_struct: sizeof<S>
  //
  // Example:
  //   struct S {
  //     double x[3];
  //     ...
  //   };
  //   vector<S> v; // length n
  // Provide a view as an n x 3 array
  // array_view_from_vector_struct<double>(v.front().x, 2, {v.size(), 3},
  // sizeof(S)}
  //
  // Note: np.array would not work after vector v deallocated

  int ndim = 2;
  if(ncol == 1)
    ndim = 1;
  
  npy_intp dims[] = {n, ncol};
    
  npy_intp strides[]= {static_cast<npy_intp>(sizeof_struct),
                       sizeof(T)};

  return PyArray_New(&PyArray_Type, ndim, dims, dtype<T>(), strides,
                     p, 0, NPY_ARRAY_WRITEABLE, 0);
}


} // unnamed namespace


namespace np_array {
//
// numpy initialiser
//
PyMODINIT_FUNC module_init()
{
  import_array();

  return NULL;
}



int dtype(const string name)
{
  if (name == "int")
    return NPY_INT;
  else if (name == "long")
    return NPY_LONG;
  else if (name == "float")
    return NPY_FLOAT;
  else if (name == "double")
    return NPY_DOUBLE;  
  
  throw TypeError();

  return 0;
}



//
// low-tech explicit instantation
//

PyObject* copy_from_vector(const vector<int>& v)
{
  return copy_from_vector_template<int>(v);
}

PyObject* copy_from_vector(const vector<long>& v)
{
  return copy_from_vector_template<long>(v);
}

PyObject* copy_from_vector(const vector<float>& v)
{
  return copy_from_vector_template<float>(v);
}

PyObject* copy_from_vector(const vector<double>& v)
{
  return copy_from_vector_template<double>(v);
}

// Don't know how to do this correctly; not working.
//template PyObject* array_from_vector<int>(const vector<int>& v);


  
// array_vew_from_vector naive instantation
PyObject* view_from_vector(vector<int>& v)
{
  return view_from_vector_template(v);
}


PyObject* view_from_vector(vector<long>& v)
{
  return view_from_vector_template(v);
}


PyObject* view_from_vector(vector<double>& v)
{
  return view_from_vector_template(v);
}

  

PyObject* view_from_vector_struct(int* const p,
                                        const int n, const int ncol,
                                        const size_t sizeof_struct)
{
  return view_from_vector_struct_template(p, n, ncol, sizeof_struct);
}


PyObject* view_from_vector_struct(long* const p,
                                        const int n, const int ncol,
                                        const size_t sizeof_struct)
{
  return view_from_vector_struct_template(p, n, ncol, sizeof_struct);
}

PyObject* view_from_vector_struct(float* const p,
                                        const int n, const int ncol,
                                        const size_t sizeof_struct)
{
  return view_from_vector_struct_template(p, n, ncol, sizeof_struct);
}

PyObject* view_from_vector_struct(double* const p,
                                        const int n, const int ncol,
                                        const size_t sizeof_struct)
{
  return view_from_vector_struct_template(p, n, ncol, sizeof_struct);
}

} // namespace array


