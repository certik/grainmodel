#include <stdlib.h>
#include <math.h>
#include <iostream>
#include "libmeshpy.h"

double doubleSum(double* series, int size) {
  double result = 0.0;
  for (int i=0; i<size; ++i) result += series[i];
  return result;
}

void doubleOnes(double* array, int size) {
  for (int i=0; i<size; ++i) array[i] = 1.0;
}

#include "libmesh.h"


#include <fstream>

#include "mesh.h"
#include "mesh_data.h"
#include "mesh_tools.h"
#include "boundary_info.h"
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
        }
    }
    ~BC()
    {
        //this doesn't work - but it's not necessary after all
        //delete nodes;
    }
    bool find(unsigned int id,int *b, int *s)
    //beware: libmesh counts elements, nodes and sides from 0, but in this
    //routine we count from 1.
    {
        for (int k=0;k<NN;k++)
            for (unsigned int i=0;i<elements[k].size();i++)
                if (elements[k][i]==id)
                {
                    *b=k+1;
                    *s=sides[k][i];
                    return true;
                }
        return false;
    }

int NN;
std::vector<unsigned int> *elements;
std::vector<unsigned int> *sides;
};

void write_int(std::ostream& f, unsigned int i)
{
    f.write((char *) &i,sizeof(unsigned int));
}

void write_float(std::ostream& f, double d)
{
    f.write((char *) &d,sizeof(double));
}

double read_float(std::istream& f)
{
    double d;
    f.read((char *)&d, sizeof(double));
    return d;
}

unsigned int read_int(std::istream& f)
{
    unsigned int i;
    f.read((char *)&i, sizeof(unsigned int));
    return i;
}

void save_matrix(std::ostream& f, DenseMatrix<Number>& M)
{
    for (unsigned int i=0;i<M.m();i++)
        for (unsigned int j=0;j<M.n();j++)
            write_float(f,M.el(i,j));
}

void save_vector(std::ostream& f, DenseVector<Number>& V)
{
    for (unsigned int i=0;i<V.size();i++)
        write_float(f,V.el(i));
}

unsigned int finda(unsigned int a, const Mesh& mesh)
{
    for (unsigned int i=0;i<mesh.n_nodes();i++)
        if (mesh.node(i).dof_number(0,0,0)==a) return i;
    error();
}

void save_vector_int(std::ostream& f, std::vector<unsigned int>& V, 
        const Mesh& mesh)
{
    for (unsigned int i=0;i<V.size();i++)
    {
        unsigned int a=V[i];
        a=finda(a,mesh);
        write_int(f,a);
    }
}

loadmatrices::loadmatrices(const std::string& fname)
{
    f=new std::ifstream(fname.c_str());
}

loadmatrices::~loadmatrices()
{
    delete f;
}

double loadmatrices::readfloat()
{
    return read_float(*f);
}

unsigned int loadmatrices::readint()
{
    return read_int(*f);
}

class matrices
//currently Ax=F
//supports construction of A and F
{
    public:
        matrices(const char *fname)
        {
            f=new std::ofstream(fname);
        }
        ~matrices()
        {
            delete f;
        }
        void setsize(int nn, int ne)
        {
            write_int(*f,nn);
            write_int(*f,ne);
        }
        void addtoA(DenseMatrix<Number> Ae,
                std::vector<unsigned int>& dof_indices,
                const Mesh& mesh)
        {
            save_matrix(*f,Ae);
            save_vector_int(*f,dof_indices,mesh);
        }
        void addtoF(DenseVector<Number> Fe,
                std::vector<unsigned int>& dof_indices,
                const Mesh& mesh)
        {
            save_vector(*f,Fe);
            save_vector_int(*f,dof_indices,mesh);
        }
    private:
        std::ofstream *f;
};

inline int get_local_id(const Elem *elem, unsigned int i)
{
    for (unsigned int s=0;s<elem->n_nodes();s++)
        if (elem->node(s)==i) return s;
    std::cout << "node " << i << " not found in elem!" << std::endl;
    error();
}

void save_elem(std::ostream& f, const Elem* el,
        std::vector<unsigned int>& dof_indices)
{
    f << el->id() << ": ";
    //indices in the global matrix
    for (unsigned int i=0;i<dof_indices.size();i++)
        //the same as: el->get_node(i)->dof_number(0,0,0)
        f << dof_indices[i] << " ";
    f << "|| ";
    //original node numbers
    for (unsigned int i=0;i<el->n_nodes();i++)
        f << el->node(i) << " ";
    f << std::endl;
}

//mapping from orig ids to libmesh's ids.
void save_node_map(const char* fname, const Mesh& mesh)
{
    std::ofstream f(fname);
    for (unsigned int i=0;i<mesh.n_nodes();i++)
        f<<mesh.node(i).dof_number(0,0,0) << " ";
    f << std::endl;
}

void assemble_poisson(EquationSystems& es)
{
    std::cout << "assembling..." << std::endl;
	PerfLog perf("Matrix Assembly");
	const Mesh& mesh = es.get_mesh();
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
	const std::vector<std::vector<Real> >& phi = fe->get_phi();
	const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();

	DenseMatrix<Number> Ke;
	DenseVector<Number> Fee;
	std::vector<unsigned int> dof_indices;

    QGauss qquad(2,FIFTH);
    qquad.init(QUAD4);

    std::cout << "nodes: " << mesh.n_nodes() << "; elements: " 
        << mesh.n_elem() << std::endl;

    BC bc("../../tmp/t12.boundaries");
    matrices mymatrices("../../tmp/matrices");
    mymatrices.setsize(mesh.n_nodes(),mesh.n_elem());

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
			if (side+1==(unsigned int)s)
		{
            if (elem->neighbor(side) != NULL) error();
			const std::vector<std::vector<Real> >&  phi_face=fe_face->get_phi();
			const std::vector<Real>& JxW_face = fe_face->get_JxW();
			fe_face->reinit(elem, side);

			Real value;
            if (b==1) value=1.0;
            else if (b==2) value=0.0;
            else error()

			for (unsigned int qp=0; qp<qface.n_points(); qp++)
			{
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

        mymatrices.addtoA(Ke,dof_indices,mesh);
        mymatrices.addtoF(Fee,dof_indices,mesh);
		perf.stop_event("matrix insertion");
	} //for element

    perf.clear();
    std::cout << "done." << std::endl;
}

void mesh(const std::string& meshfile)
{
    std::cout << "starting..." << std::endl;

    int argc=1; char *p="./lmesh\n"; char **argv=&p;
    libMesh::init (argc, argv);
    {    
        Mesh mesh(3);
        mesh.read(meshfile);
        mesh.find_neighbors();
        mesh.print_info();
        EquationSystems equation_systems (mesh);
        equation_systems.add_system<LinearImplicitSystem> ("Poisson");
        equation_systems.get_system("Poisson").add_variable("u", FIRST);
        equation_systems.init();
        assemble_poisson(equation_systems);
    }
    libMesh::close();
}
