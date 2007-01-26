
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


#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **args)
{
  KSP            ksp;             /* linear solver context */
  Mat            A;                /* matrix */
  Vec            x,b,u;          /* exact solution, RHS, temp */
  PetscViewer    fd;               /* viewer */
  char           file[PETSC_MAX_PATH_LEN];     /* input file name */
  char           file2[PETSC_MAX_PATH_LEN];     /* input file name */
  char           file3[PETSC_MAX_PATH_LEN];     /* input file name */
  PetscErrorCode ierr,ierrp;
  PetscInt       its;
  PetscReal      norm;

  PetscInitialize(&argc,&args,(char *)0,help);

  ierr = PetscPrintf(PETSC_COMM_WORLD,"initializing...\n");CHKERRQ(ierr);


  ierr = PetscOptionsGetString(PETSC_NULL,"-f1",file,PETSC_MAX_PATH_LEN-1,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetString(PETSC_NULL,"-f2",file2,PETSC_MAX_PATH_LEN-1,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetString(PETSC_NULL,"-f3",file3,PETSC_MAX_PATH_LEN-1,PETSC_NULL);CHKERRQ(ierr);

    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,file,FILE_MODE_READ,&fd);CHKERRQ(ierr);
    ierr  = MatLoad(fd,MATMPIAIJ,&A);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(fd);CHKERRQ(ierr);

    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,file2,FILE_MODE_READ,&fd);CHKERRQ(ierr);
    ierrp = VecLoad(fd,PETSC_NULL,&b);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(fd);CHKERRQ(ierr);

    ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);

    ierr = KSPSetOperators(ksp,A,A,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
    ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);

    ierr = VecDuplicate(b,&x);CHKERRQ(ierr);

    ierr = PetscPrintf(PETSC_COMM_WORLD,"solver start...\n");CHKERRQ(ierr);
    ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"solver end...\n");CHKERRQ(ierr);

    ierr = VecDuplicate(b,&u);CHKERRQ(ierr);
    ierr = MatMult(A,x,u);CHKERRQ(ierr);
    ierr = VecAXPY(u,-1.0,b);CHKERRQ(ierr);
    ierr = VecNorm(u,NORM_2,&norm);CHKERRQ(ierr);
    ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Number of iterations = %3D\n",its);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Residual norm %A\n",norm);CHKERRQ(ierr);


    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, file3, FILE_MODE_WRITE, 
            &fd);CHKERRQ(ierr);
    ierr = PetscViewerSetFormat (fd, PETSC_VIEWER_BINARY_DEFAULT);CHKERRQ(ierr);
    ierr = VecView (x, fd);CHKERRQ(ierr);
    ierr = PetscViewerDestroy (fd);CHKERRQ(ierr);

    /* 
       Free work space.  All PETSc objects should be destroyed when they
       are no longer needed.
    */
    ierr = MatDestroy(A);CHKERRQ(ierr); ierr = VecDestroy(b);CHKERRQ(ierr);
    ierr = VecDestroy(u);CHKERRQ(ierr); ierr = VecDestroy(x);CHKERRQ(ierr);
    ierr = KSPDestroy(ksp);CHKERRQ(ierr); 

  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}

