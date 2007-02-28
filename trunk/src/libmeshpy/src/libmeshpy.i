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

/* Apply the Numeric typemaps for 1D input arrays */
%apply (short*  IN_ARRAY1, int DIM1) {(short*  series, int size)};
%apply (int*    IN_ARRAY1, int DIM1) {(int*    series, int size)};
%apply (long*   IN_ARRAY1, int DIM1) {(long*   series, int size)};
%apply (float*  IN_ARRAY1, int DIM1) {(float*  series, int size)};
%apply (double* IN_ARRAY1, int DIM1) {(double* series, int size)};

/* Apply the Numeric typemaps for 2D input arrays */
%apply (int*    IN_ARRAY2, int DIM1, int DIM2) {(int*    matrix, int rows, int cols)};
%apply (double* IN_ARRAY2, int DIM1, int DIM2) {(double* matrix, int rows, int cols)};

/* Apply the Numeric typemaps for 1D input/output arrays */
%apply (int*    INPLACE_ARRAY1, int DIM1) {(int*    array,   int size)};
%apply (double* INPLACE_ARRAY1, int DIM1) {(double* array,   int size)};

%apply (double* IN_ARRAY1, int DIM1) {(double* x,   int xsize)};
%apply (double* INPLACE_ARRAY1, int DIM1) {(double* g,   int gsize)};

%feature("director") Updater;
/* Include the header file to be wrapped */
%include "libmeshpy.h"
