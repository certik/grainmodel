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
    std::cout << "hello" << std::endl;

    int argc=0;
    char **argv=NULL;
    libMesh::init (argc, argv);
    {    
        const unsigned int dim = 3;
        Mesh mesh(dim);
        //mesh.read (argv[3]);
        mesh.print_info();
        //mesh.write (argv[4]);
    }
    libMesh::close();
}


