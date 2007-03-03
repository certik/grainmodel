%module(directors="1") libmeshpy

%{
#define SWIG_FILE_WITH_INIT
#include "libmeshpy.h"
%}

%include "std_string.i"

/* Get the Numeric typemaps */
%include "numpy.i"

%init %{
  import_array();
%}

/* Apply the Numeric typemaps for 2D input arrays */
%apply (int*    IN_ARRAY2, int DIM1, int DIM2) {(int*    matrix, int rows, int cols)};
%apply (double* IN_ARRAY2, int DIM1, int DIM2) {(double* matrix, int rows, int cols)};

/* Apply the Numeric typemaps for 1D input/output arrays */
%apply (int*    INPLACE_ARRAY1, int DIM1) {(int*    array,   int size)};
%apply (double* INPLACE_ARRAY1, int DIM1) {(double* array,   int size)};

%apply (double* IN_ARRAY1, int DIM1) {(double* x,   int xsize)};
%apply (double* IN_ARRAY1, int DIM1) {(double* y,   int ysize)};
%apply (double* IN_ARRAY1, int DIM1) {(double* z,   int zsize)};
%apply (double* INPLACE_ARRAY1, int DIM1) {(double* gx,   int gxsize)};
%apply (double* INPLACE_ARRAY1, int DIM1) {(double* gy,   int gysize)};
%apply (double* INPLACE_ARRAY1, int DIM1) {(double* gz,   int gzsize)};

%apply (double* IN_ARRAY1, int DIM1) {(double* bvalues,   int vsize)};
%apply (int* IN_ARRAY1, int DIM1) {(int* bidx,   int isize)};

%feature("director") Updater;
/* Include the header file to be wrapped */
%include "libmeshpy.h"
