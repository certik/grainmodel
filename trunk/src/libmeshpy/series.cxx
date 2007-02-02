#include <stdlib.h>
#include <math.h>
#include <iostream>
#include "series.h"

double doubleSum(double* series, int size) {
  double result = 0.0;
  for (int i=0; i<size; ++i) result += series[i];
  return result;
}

void doubleOnes(double* array, int size) {
  for (int i=0; i<size; ++i) array[i] = 1.0;
}

#include "libmesh.h"
#include "mesh.h"

void mesh()
{
    std::cout << "starting..." << std::endl;

    int argc=1; char *p="./lmesh\n"; char **argv=&p;
    libMesh::init (argc, argv);
    {    
        Mesh mesh(3);
        mesh.print_info();
    }
    libMesh::close();
}


