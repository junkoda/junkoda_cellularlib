#ifndef BUFFER_H
#define BUFFER_H 1

// Resource management class for Py_Buffer
// Python Buffer protocol is handles np.array-like data structures

#include <string>
#include <complex>
#include <vector>
#include <typeinfo>


#include "Python.h"

#include "error.h"
#include "py_util.h"

template <typename T>
class Buffer {
 public:
  // Buffer buf(py_array);
  // Args:
  //   py_array: Buffer object, e.g., numpy array
  Buffer();
  Buffer(PyObject* const py_array, char const * const name_=nullptr);
  ~Buffer();

  // delete assignment and copy
  Buffer(const Buffer&) = delete;
  Buffer& operator=(Buffer const&) = delete;

  T operator[](const size_t i) const {
    return buf[stride.back()*i];
  }
  T& operator[](const size_t i) {
    return buf[stride.back()*i];
  }
  
  // 1D
  T operator()(const size_t ix) const {
    size_t i = stride[0]*ix;
    return buf[i];
  }
  T& operator()(const size_t ix) {
    size_t i = stride[0]*ix;
    return buf[i];
  }

  // 2D
  T operator()(const size_t ix, const size_t iy) const {
    size_t i = stride[0]*ix + stride[1]*iy;
    return buf[i];
  }
  T& operator()(const size_t ix, const size_t iy) {
    size_t i = stride[0]*ix + stride[1]*iy;
    return buf[i];
  }
  
  void assign(PyObject* const py_array);
  void release();
  T* data(){ return buf; }
  T const * data() const { return buf; }

  int ndim;
  std::vector<size_t> shape;
  std::vector<size_t> stride; // in units of elements converted from bytes
  std::string str_format;
  T* buf;
 private:
  Py_buffer pybuf;
  char const * const name;
  void type_error(const char error_msg[]) const;
};

//
// Default constructre
//
template<typename T>
Buffer<T>::Buffer() : ndim(0), buf(nullptr)
{
}

//
// Contruct with an array
//
template<typename T>
Buffer<T>::Buffer(PyObject* const py_array, char const * const name_) :
ndim(0), buf(nullptr), name(name_)
{
  // May throw TypeError
  assign(py_array);
}

//
// Deconstrucor: automatically realease buffer object
//
template<typename T>
Buffer<T>::~Buffer()
{
  release();
}

template<typename T>
void Buffer<T>::assign(PyObject* const py_array)
{
  // Args:
  //   py_array: Python buffer object
  //
  // Throws: TypeError
  //
  release();
  
  char msg[128];

  if(PyObject_GetBuffer(py_array, &pybuf,
			PyBUF_FORMAT | PyBUF_STRIDED) == -1) {

    type_error("Expected a contiguous or strided buffer protocol");
  }

  //
  // Check data type
  //
  const char c = pybuf.format[0];

  // first character may indicate the byteorder
  if(c == '@' || c == '=' || c == '<') {
    str_format = pybuf.format + 1;
  }
  else {
    str_format = pybuf.format;
  }

  if(str_format == "i") {
    if(typeid(T) != typeid(int)) {
      sprintf(msg, "Expected an array of int but it is %.4s",
	      pybuf.format);
      type_error(msg);
    }
  }  
  else if(str_format == "l") {
    if(typeid(T) != typeid(long)) {
      sprintf(msg, "Expected an array of %.16s but type is %.4s",
	      util_type_name<T>(), pybuf.format);
      type_error(msg);
    }
  }  
  else if(str_format == "q") {
    if(typeid(T) != typeid(long long)) {
      sprintf(msg, "Expected an array of %.16s but type is %.4s",
	      util_type_name<T>(), pybuf.format);
      type_error(msg);
    }
  }
  else if(str_format == "f") {
    if(typeid(T) != typeid(float)) {
      sprintf(msg, "Expected an array of %.16s but type is %.4s",
	      util_type_name<T>(), pybuf.format);
      type_error(msg);
    }
  }
  else if(str_format == "d") {
    if(typeid(T) != typeid(double)) {
      sprintf(msg, "Expected an array of %.16s but type is %.4s",
	      util_type_name<T>(), pybuf.format);
      type_error(msg);
    }
  }
  else if(str_format == "Zd") {
    if(typeid(T) != typeid(std::complex<double>)) {
      sprintf(msg, "Expected an array of %.16s but type is %.4s",
	      util_type_name<T>(), pybuf.format);
      type_error(msg);
    }
  }
  else if(str_format == "Zf") {
    if(typeid(T) != typeid(std::complex<float>)) {
      sprintf(msg, "Expected an array of %.16s but type is %.4s",
	      util_type_name<T>(), pybuf.format);
      type_error(msg);
    }
  }
  else if(str_format == "?") { // bool
    if(typeid(T) != typeid(bool)) {
      sprintf(msg, "Expected an array of char but it is %.4s",
	      pybuf.format);
      type_error(msg);
    }
  }
  else if(str_format == "b") { // np.int8
    if(typeid(T) != typeid(char)) {
      sprintf(msg, "Expected an array of char but it is %.4s",
	      pybuf.format);
      type_error(msg);
    }
  }  
  else if(str_format == "B") {
    if(typeid(T) != typeid(unsigned char)) {
      sprintf(msg, "Expected an array of unsigned char but it is %.4s",
	      pybuf.format);
      type_error(msg);
    }
  }
  else {
    sprintf(msg, "Array has an unknown dtype %.4s.", str_format.c_str());
    type_error(msg);
  }

  buf = reinterpret_cast<T*>(pybuf.buf);

  
  //
  // ndim, shape, and stride
  //
  ndim = pybuf.ndim;

  for(int i=0; i<ndim; ++i) {
    shape.push_back(pybuf.shape[i]);
  }


  //PyBuffer_IsContiguous(Py_buffer *view, char order)¶
  //PyBuffer_FillContiguousStrides(ndim, pybuf

  //
  // Stride
  //
  if(pybuf.strides) {
    for(int i=0; i<ndim; ++i) {
      assert(pybuf.strides[i] % sizeof(T) == 0);
      stride.push_back(pybuf.strides[i] / sizeof(T));
    }
  }
  else {
    assert(false);
  }
}


//
// Release buffer
//
template<typename T>
void Buffer<T>::release()
{
  if(buf) {
    PyBuffer_Release(&pybuf);
  }
  buf = nullptr;
  ndim = 0;
}

template<typename T>
void Buffer<T>::type_error(const char error_msg[]) const
{
  // Set Python exception message and then thow TypeError
  std::string msg(error_msg);
  
  if(name)
    msg += " (" + std::string(name) + ").";

  PyErr_SetString(PyExc_TypeError, msg.c_str());

  throw TypeError();
}


#endif
