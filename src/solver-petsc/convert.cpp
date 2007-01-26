
static char help[] = "Reads a PETSc matrix and vector from a file and solves the normal equations.\n\n";
/*T
   Concepts: KSP^solving a linear system
   Concepts: Normal equations
   Processors: n
T*/

/* 
  Include "petscksp.h" so that we can use KSP solvers.  Note that this file
  automatically includes:
     petsc.h       - base PETSc routines   petscvec.h - vectors
     petscsys.h    - system routines       petscmat.h - matrices
     petscis.h     - index sets            petscksp.h - Krylov subspace methods
     petscviewer.h - viewers               petscpc.h  - preconditioners
*/
#include "petscksp.h"
#include <iostream>


#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **args)
{
  Vec            x;          /* exact solution, RHS, temp */
  PetscViewer    fd;               /* viewer */
  char           file[PETSC_MAX_PATH_LEN];     /* input file name */
  char           file2[PETSC_MAX_PATH_LEN];     /* input file name */
  PetscErrorCode ierr,ierrp;

  PetscInitialize(&argc,&args,(char *)0,help);

  ierr = PetscPrintf(PETSC_COMM_WORLD,"initializing...\n");CHKERRQ(ierr);


  ierr = PetscOptionsGetString(PETSC_NULL,"-f1",file,PETSC_MAX_PATH_LEN-1,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetString(PETSC_NULL,"-f2",file2,PETSC_MAX_PATH_LEN-1,PETSC_NULL);CHKERRQ(ierr);

    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,file,FILE_MODE_READ,&fd);CHKERRQ(ierr);
    ierrp = VecLoad(fd,PETSC_NULL,&x);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(fd);CHKERRQ(ierr);

    PetscScalar       *a;
    VecGetArray(x,&a);
    PetscInt n;
    VecGetSize(x,&n);
    for (int i=0;i<n;i++) 
    {
        std::cout << a[i] << " ";
    }

    ierr = VecDestroy(x);CHKERRQ(ierr);

  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}

