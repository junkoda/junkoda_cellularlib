#ifndef PY_UTIL_H
#define PY_UTIL_H 1

#include <queue>
#include <complex>
#include <cassert>

#include "Python.h"

//void util_sequence_as_queue(PyObject* seq, std::queue<double>& q);

//
// Return template variable type as a string
//
template<typename T>
const char* util_type_name()
{
  if (typeid(T) == typeid(int))
    return "int";
  else if (typeid(T) == typeid(long))
    return "long";
  else if (typeid(T) == typeid(long long))
    return "long long";
  else if (typeid(T) == typeid(float))
    return "float";
  else if (typeid(T) == typeid(double))
    return "double";
  else if (typeid(T) == typeid(std::complex<float>))
    return "complex<float>";
  else if (typeid(T) == typeid(std::complex<double>))
    return "complex<double>";
  else if (typeid(T) == typeid(char))
    return "char";
  else if (typeid(T) == typeid(unsigned char))
    return "unsigned char";

  //else {
  //  PyErr_SetString(PyExc_TypeError, "Unknown type for util_type_name");
  // }

  throw typeid(T).name();

  return nullptr;
}


//
// Convert sequence to vector
//
template<typename T>
void util_sequence_as_vector(PyObject* seq, std::vector<T>& v)
{
  //
  // Add the contens of the Python squence to C++ vector v
  // T must be either int or double
  //
  if(PySequence_Check(seq) == false) {
    PyErr_SetString(PyExc_TypeError, "sequence is not provide");
  }

  size_t n = (size_t) PySequence_Length(seq);
  v.reserve(v.size() + n);

  for(size_t i=0; i<n; ++i) {
    PyObject* elem = PySequence_GetItem(seq, i); assert(elem);
    // new reference
    if(std::is_same<T, int>::value)
      v.push_back(PyLong_AsLong(elem));
    else if(std::is_same<T, double>::value)
      v.push_back(PyFloat_AsDouble(elem));
    else {
      PyErr_SetString(PyExc_TypeError,
		      "The value in the sequence is not the type expected");
    }
    Py_DECREF(elem);
  }
}

#endif
