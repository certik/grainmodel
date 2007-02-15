#include <stdlib.h>
#include <math.h>
#include <iostream>
#include "libmesh_mod.h"

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
        mesh.read("../../tmp/in.xda");
        mesh.print_info();
    }
    libMesh::close();
}

#include <fstream>

#include "mesh.h"
#include "mesh_data.h"
#include "mesh_tools.h"
#include "boundary_info.h"
//#include "eigen_system.h"
#include "linear_implicit_system.h"
#include "equation_systems.h"

#include "fe.h"
#include "quadrature_gauss.h"

#include "sparse_matrix.h"
#include "petsc_vector.h"
#include "petsc_matrix.h"
#include "dense_matrix.h"
#include "dense_vector.h"
#include "elem.h"

#include "dof_map.h"

//#include "configfile.h"

//extern Config config;

class BC
{
public:
    BC(const char *fname)
    {
        std::ifstream f(fname);
        f >> NN;
        elements=new std::vector<unsigned int>[NN];
        sides=new std::vector<unsigned int>[NN];
        for (int k=0;k<NN;k++)
        {
            int n;
            f >> n;
            int count;
            f >> count;
            for (int j=0;j<count;j++)
            {
                int i;
                f >> i;
                elements[k].push_back(i);
                f >> i;
                sides[k].push_back(i);
            }
//            std::cout << "HE:" << count << " " << elements[k].size();
        }
    }
    ~BC()
    {
        //this doesn't work - but it's not necessary after all
        //delete nodes;
    }
    bool find(int id,int *b, int *s)
    //beware: libmesh counts elements, nodes and sides from 0, but in this
    //routine we count from 1.
    {
/*        if (id==60) {
            *b=1;
            *s=2;
            return true;
        }
        if (id==93) {
            *b=1;
            *s=4;
            return true;
        }
        if (id==197) {
            *b=1;
            *s=3;
            return true;
        }
        if (id==203) {
            *b=2;
            *s=1;
            return true;
        }
        if (id==214) {
            *b=2;
            *s=4;
            return true;
        }
        if (id==216) {
            *b=2;
            *s=4;
            return true;
        }
        if (id==234) {
            *b=2;
            *s=1;
            return true;
        }
        return false;
        */
        for (int k=0;k<NN;k++)
            for (int i=0;i<elements[k].size();i++)
                if (elements[k][i]==id)
                {
                    *b=k+1;
                    *s=sides[k][i];
//                    std::cout << id << " " << *b << " " << *s << " SS " << i <<
//                        " " << k << " " << elements[k].size() << std::endl;
                    return true;
                }
        return false;
    }

int NN;
std::vector<unsigned int> *elements;
std::vector<unsigned int> *sides;
};

inline int get_local_id(const Elem *elem, int i)
{
    for (int s=0;s<elem->n_nodes();s++)
        if (elem->node(s)==i) return s;
    std::cout << "node " << i << " not found in elem!" << std::endl;
    error();
}

void save_matrix(std::ostream& f, DenseMatrix<Number>& M)
{
    M.print_scientific(f);
}

void save_vector(NumericVector<Number>& M, const char *fname)
{
//    M.print_matlab(fname);
/*    
    int ierr=0; 
    PetscViewer petsc_viewer;

    M.close();
    ierr = PetscViewerCreate (libMesh::COMM_WORLD, &petsc_viewer);
    CHKERRABORT(libMesh::COMM_WORLD,ierr);

    ierr = PetscViewerBinaryOpen( libMesh::COMM_WORLD, fname, FILE_MODE_WRITE,
            &petsc_viewer);
    CHKERRABORT(libMesh::COMM_WORLD,ierr);

    ierr = PetscViewerSetFormat (petsc_viewer, PETSC_VIEWER_BINARY_DEFAULT);
    CHKERRABORT(libMesh::COMM_WORLD,ierr);

    Vec mat = ((PetscVector<Number>&)(M)).vec();

    ierr = VecView (mat, petsc_viewer);
    CHKERRABORT(libMesh::COMM_WORLD,ierr);

    ierr = PetscViewerDestroy (petsc_viewer);
    CHKERRABORT(libMesh::COMM_WORLD,ierr);
    */
}

void save_sparse_matrix(SparseMatrix<Number>& M, const char *fname)
{
//    M.print_matlab(fname);

/*    int ierr=0; 
    PetscViewer petsc_viewer;

    M.close();
    ierr = PetscViewerCreate (libMesh::COMM_WORLD, &petsc_viewer);
    CHKERRABORT(libMesh::COMM_WORLD,ierr);

    ierr = PetscViewerBinaryOpen( libMesh::COMM_WORLD, fname, FILE_MODE_WRITE,
            &petsc_viewer);
    CHKERRABORT(libMesh::COMM_WORLD,ierr);

    ierr = PetscViewerSetFormat (petsc_viewer, PETSC_VIEWER_BINARY_DEFAULT);
    CHKERRABORT(libMesh::COMM_WORLD,ierr);

    Mat mat = ((PetscMatrix<Number>&)(M)).mat();

    ierr = MatView (mat, petsc_viewer);
    CHKERRABORT(libMesh::COMM_WORLD,ierr);

    ierr = PetscViewerDestroy (petsc_viewer);
    CHKERRABORT(libMesh::COMM_WORLD,ierr);
    */
}

void save_elem(std::ostream& f, const Elem* el,
        std::vector<unsigned int>& dof_indices)
{
    f << el->id() << ": ";
    //indices in the global matrix
    for (int i=0;i<dof_indices.size();i++)
        //the same as: el->get_node(i)->dof_number(0,0,0)
        f << dof_indices[i] << " ";
    f << "|| ";
    //original node numbers
    for (int i=0;i<el->n_nodes();i++)
        f << el->node(i) << " ";
    f << std::endl;
}

//mapping from orig ids to libmesh's ids.
void save_node_map(const char* fname, const Mesh& mesh)
{
    std::ofstream f(fname);
    for (int i=0;i<mesh.n_nodes();i++)
        f<<mesh.node(i).dof_number(0,0,0) << " ";
    f << std::endl;
}

void load_zeronodes(const char* fname, std::vector<unsigned int>& nodes)
{
    std::ifstream f(fname);
    while (!f.eof())
    {
        int i;
        f >> i;
        nodes.push_back(i);
    }
}

void assemble_poisson(EquationSystems& es,
                      const std::string& system_name)
{
	assert (system_name == "Poisson");

    std::cout << "assembling..." << std::endl;
	PerfLog perf("Matrix Assembly");
	const Mesh& mesh = es.get_mesh();
	const MeshData& mesh_data = es.get_mesh_data();
	const unsigned int dim = mesh.mesh_dimension();
	LinearImplicitSystem& system=es.get_system<LinearImplicitSystem>("Poisson");
	const DofMap& dof_map = system.get_dof_map();
	FEType fe_type = dof_map.variable_type(0);
	AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));
	QGauss qrule (dim, FIFTH);
	fe->attach_quadrature_rule (&qrule);
	AutoPtr<FEBase> fe_face (FEBase::build(dim, fe_type));
	QGauss qface(dim-1, FIFTH);
	fe_face->attach_quadrature_rule (&qface);
	const std::vector<Real>& JxW = fe->get_JxW();
	const std::vector<Point>& q_point = fe->get_xyz();
	const std::vector<std::vector<Real> >& phi = fe->get_phi();
	const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();

	DenseMatrix<Number> Ke;
	DenseVector<Number> Fee;
	std::vector<unsigned int> dof_indices;

    SparseMatrix<Number>&  matrix_A = *system.matrix;
    NumericVector<Number>&  vector_F = *system.rhs;
//    PetscVector<Number> vector_F;

    QGauss qquad(2,FIFTH);
    //QGauss qquad(2,FORTYTHIRD);
    qquad.init(QUAD4);
    std::vector<std::vector<const Elem *> > nodes_to_elem_map;
    MeshTools::build_nodes_to_elem_map(mesh,nodes_to_elem_map);

    std::cout << "nodes: " << mesh.n_nodes() << "; elements: " 
        << mesh.n_elem() << std::endl;

    std::ofstream f("tmp/matrices.dat");
    std::ofstream ftopo("tmp/topo.dat");
    save_node_map("tmp/nodemap.libmesh", mesh);
//	std::vector<unsigned int> zeronodes;
//    load_zeronodes("tmp/zeronodes.gmsh", zeronodes);
    BC bc("tmp/t12.boundaries");

	MeshBase::const_element_iterator       el     = mesh.elements_begin();
	const MeshBase::const_element_iterator end_el = mesh.elements_end();
    std::cout << "Start" << std::endl;
	for ( ; el != end_el ; ++el)
	{
		perf.start_event("elem init");
		const Elem* elem = *el;
        if (elem->id() % 10000 == 0)
            std::cout << 100.0*elem->id()/mesh.n_elem() << "%" << std::endl;
		dof_map.dof_indices (elem, dof_indices);
		fe->reinit (elem);
		Ke.resize (dof_indices.size(), dof_indices.size());
		Fee.resize (dof_indices.size());
		perf.stop_event("elem init");

		perf.start_event("Ke");
		for (unsigned int qp=0; qp<qrule.n_points(); qp++)
		{
//            Real hbar=config.hbar;
//            Real m=config.m;
            Real lambda=1.0;
			for (unsigned int i=0; i<phi.size(); i++)
			for (unsigned int j=0; j<phi.size(); j++)
				Ke(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp])*lambda;
		} 
		perf.stop_event("Ke");

		{
		perf.start_event("Fe");
        int b,s;
        if (bc.find(elem->id()+1,&b,&s))
		for (unsigned int side=0; side<elem->n_sides(); side++)
			if (side+1==s)
		{
 /*           std::cout << b << " " << s << " "<<  elem->neighbor(side) << std::endl;
            std::cout << elem->id() << ": ";
            for (int i=0;i<dof_indices.size();i++)
                std::cout << dof_indices[i] << " ";
            std::cout << "|| ";
            for (int i=0;i<elem->n_nodes();i++)
                std::cout << elem->node(i) << " ";
            std::cout << std::endl;
            */
            if (elem->neighbor(side) != NULL) error();
			const std::vector<std::vector<Real> >&  phi_face=fe_face->get_phi();
			const std::vector<Real>& JxW_face = fe_face->get_JxW();
			const std::vector<Point >& qface_point = fe_face->get_xyz();
			fe_face->reinit(elem, side);

			Real value;
            if (b==1) value=1.0;
            else if (b==2) value=0.0;
            else error()

/*            const Real penalty = 1.e10;
				for (unsigned int i=0; i<phi_face.size(); i++)
				for (unsigned int j=0; j<phi_face.size(); j++)
					Ke(i,j) += penalty;
//                    Ke.print_scientific(std::cout);
//                    std::cout << "||";

				for (unsigned int i=0; i<phi_face.size(); i++)
					Fee(i) += penalty*value;
                    */

			for (unsigned int qp=0; qp<qface.n_points(); qp++)
			{
				const Real xf = qface_point[qp](0);
				const Real yf = qface_point[qp](1);
				const Real penalty = 1.e10;

				for (unsigned int i=0; i<phi_face.size(); i++)
				for (unsigned int j=0; j<phi_face.size(); j++)
					Ke(i,j) += JxW_face[qp]*
						penalty*phi_face[i][qp]*phi_face[j][qp];

				for (unsigned int i=0; i<phi_face.size(); i++)
					Fee(i) += JxW_face[qp]*penalty*value*phi_face[i][qp];
			} 

		}
		perf.stop_event("Fe");
        }

		perf.start_event("matrix insertion");

//        matrix_A.add_matrix (Ke, dof_indices);
//        vector_F.add_vector (Fee, dof_indices);
		perf.stop_event("matrix insertion");
	} //for element

/*    //assing the value 0 to every node from "zeronodes" using a penalty method
    for (int i=0;i<zeronodes.size();i++)
    {
        unsigned int nd =mesh.node(zeronodes[i]-1).dof_number(0,0,0);
        //std::cout << nd << " ";
        const Real penalty = 1.e10;
        matrix_A.add(nd,nd,penalty*1001.5);
        matrix_B.add(nd,nd,penalty);
    }
    */

    //print matrices A and M
    std::cout << "saving matrix to tmp/matA.petsc" << std::endl;
    save_sparse_matrix(matrix_A,"tmp/matA.petsc");
    std::cout << "saving vector to tmp/vecF.petsc" << std::endl;
    save_vector(vector_F,"tmp/vecF.petsc");
//    save_sparse_matrix(matrix_B,"tmp/matM.petsc");
    //matrix_A.print_matlab("tmp/matA.matlab");
    //matrix_B.print_matlab("tmp/matM.matlab");
//    if (!config.printlog) perf.clear();
    std::cout << "done." << std::endl;
}
