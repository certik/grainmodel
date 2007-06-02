#include <stdlib.h>
#include <math.h>
#include <iostream>

#include "libmeshpy.h"

/*
 * some notes about node numbering:
 * mesh.node(i)->id()   .... the i (and the returned value) is the node number
 * as read from the in.xda file
 * mesh.node(i)->dof_number(0,0,0) returns a number, which is used in the global
 * matrix construction, which is returned in dof_indices, returned by
 * dof_map.dof_indices (elem, dof_indices);
 *
 * The policy is to export only the original numbers. That libmesh also uses
 * some other numberings should be hidden from the user.
*/

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
    BC(const char *fname, double *bvalues, int* bidx, int bsize)
    {
        this->bidx=bidx;
        this->bvalues=bvalues;
        this->bsize=bsize;
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
    bool isin(int i, int* bidx, int isize, int *k)
    {
        for (int j=0;j<isize;j++)
            if (bidx[j]==i) 
            {
                *k=j;
                return true;
            }
        return false;
    }
    bool find2(unsigned int id,int *b, int *s, double *bval)
    //beware: libmesh counts elements, nodes and sides from 0, but in this
    //routine we count from 1.
    {
        if (find(id,b,s)) 
        {
            int j;
            if (isin(*b,bidx,bsize,&j))
            {
                *bval=bvalues[j];
                return true;
            }
        }
        return false;
    }

int NN;
std::vector<unsigned int> *elements;
std::vector<unsigned int> *sides;
double *bvalues;
int *bidx;
int bsize;
};

double myabs(double d)
{
    if (d<0) return -d;
    else return d;
}

void write_int(std::ostream& f, unsigned int i)
{
    f.write((char *) &i,sizeof(unsigned int));
}

void write_float(std::ostream& f, double d)
{
    f.write((char *) &d,sizeof(double));
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

void save_vector_int(std::ostream& f, std::vector<unsigned int>& V, 
        unsigned int *map)
{
    for (unsigned int i=0;i<V.size();i++)
    {
        write_int(f,map[V[i]]);
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
    double d;
    (*f).read((char *)&d, sizeof(double));
    return d;
}

unsigned int loadmatrices::readint()
{
    unsigned int i;
    (*f).read((char *)&i, sizeof(unsigned int));
    return i;
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
        void setsize(int nn, int ne, int linear)
        {
            write_int(*f,nn);
            write_int(*f,ne);
            write_int(*f,linear);
        }
        void addtoA(DenseMatrix<Number> Ae,
                std::vector<unsigned int>& dof_indices,
                unsigned int *map)
        {
            save_vector_int(*f,dof_indices,map);
            save_matrix(*f,Ae);
        }
        void addtoF(DenseVector<Number> Fe)
        {
            save_vector(*f,Fe);
        }
    private:
        std::ofstream *f;
};

void mesh(const std::string& fmesh, const std::string& fmatrices,
    const std::string& fboundaries,
    double* bvalues, int vsize,
    int* bidx, int isize,
    double* lambda, int lsize,
    Updater *up)
{
    //int argc=1; char *p="./lmesh\n"; char **argv=&p;
    //char *p[3]={"./lmesh","--disable-mpi","--disable-petsc"}; 
    //int argc=2; char **argv=p;
    char *p[1]={"./lmesh"}; 
    int argc=1; char **argv=p;
    libMesh::init(argc, argv);
    {    
        Mesh mesh(3);
        mesh.read(fmesh);
        mesh.find_neighbors();
        int linear=mesh.elem(0)->type()==TET4;
        EquationSystems equation_systems (mesh);
        equation_systems.add_system<LinearImplicitSystem> ("Poisson");
        if (linear)
            equation_systems.get_system("Poisson").add_variable("u", FIRST);
        else
            equation_systems.get_system("Poisson").add_variable("u", SECOND);
        //this needs the MPI initialized:
        equation_systems.init();

        const unsigned int dim = mesh.mesh_dimension();
        LinearImplicitSystem& system=
            equation_systems.get_system<LinearImplicitSystem>("Poisson");
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

        BC bc(fboundaries.c_str(),bvalues,bidx,isize);
        matrices mymatrices(fmatrices.c_str());
        mymatrices.setsize(mesh.n_nodes(),mesh.n_elem(), linear);

        unsigned int nodemap[mesh.n_nodes()];
        for (unsigned int i=0;i<mesh.n_nodes();i++)
            nodemap[mesh.node(i).dof_number(0,0,0)]=i;

        MeshBase::const_element_iterator       el     = mesh.elements_begin();
        const MeshBase::const_element_iterator end_el = mesh.elements_end();
        if (up!=NULL) up->init(mesh.n_elem()-1);
        for ( ; el != end_el ; ++el)
        {
            const Elem* elem = *el;
            if (up!=NULL) up->update(elem->id());
            dof_map.dof_indices (elem, dof_indices);
            //std::cout << dof_indices.size() << " " <<
            //    elem->type() << std::endl;
            fe->reinit (elem);
            Ke.resize (dof_indices.size(), dof_indices.size());
            Fee.resize (dof_indices.size());

            for (unsigned int qp=0; qp<qrule.n_points(); qp++)
            {
                Real lam=lambda[elem->id()];
                for (unsigned int i=0; i<phi.size(); i++)
                for (unsigned int j=0; j<phi.size(); j++)
                    Ke(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp])*lam;
            } 

            {
            int b,s;
            double bval;
            if (bc.find2(elem->id()+1,&b,&s,&bval))
            for (unsigned int side=0; side<elem->n_sides(); side++)
                if ((side+1==(unsigned int)s) )
            {
                if (elem->neighbor(side) != NULL) error();
                const std::vector<std::vector<Real> >&  phi_face=fe_face->get_phi();
                const std::vector<Real>& JxW_face = fe_face->get_JxW();
                fe_face->reinit(elem, side);

                Real value=bval;

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
            }


            mymatrices.addtoA(Ke,dof_indices,nodemap);
            mymatrices.addtoF(Fee);
        } //for element
    }
    libMesh::close();
}

void grad(const std::string& meshfile, double* x, int xsize, 
        double* gx, int gxsize,
        double* gy, int gysize,
        double* gz, int gzsize,
        Updater *up)
{
    int argc=1; char *p="./lmesh\n"; char **argv=&p;
    libMesh::init (argc, argv);
    {    
        Mesh mesh(3);
        mesh.read(meshfile);
        mesh.find_neighbors();
        EquationSystems equation_systems (mesh);
        equation_systems.add_system<LinearImplicitSystem> ("Poisson");
        equation_systems.get_system("Poisson").add_variable("u", FIRST);
        equation_systems.init();
        
	const unsigned int dim = mesh.mesh_dimension();
	LinearImplicitSystem& system=equation_systems.get_system<LinearImplicitSystem>("Poisson");
	const DofMap& dof_map = system.get_dof_map();
	FEType fe_type = dof_map.variable_type(0);
	AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));
	QGauss qrule (dim, FIFTH);
	fe->attach_quadrature_rule (&qrule);
	AutoPtr<FEBase> fe_face (FEBase::build(dim, fe_type));
	QGauss qface(dim-1, FIFTH);
	fe_face->attach_quadrature_rule (&qface);
	const std::vector<std::vector<Real> >& phi = fe->get_phi();
	const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();

	DenseMatrix<Number> Ke;
	DenseVector<Number> Fee;
	std::vector<unsigned int> dof_indices;

	MeshBase::const_element_iterator       el     = mesh.elements_begin();
	const MeshBase::const_element_iterator end_el = mesh.elements_end();
    if (up!=NULL) up->init(mesh.n_elem()-1);
	for ( ; el != end_el ; ++el)
	{
		const Elem* elem = *el;
        if (up!=NULL) up->update(elem->id());
		dof_map.dof_indices (elem, dof_indices);
		fe->reinit (elem);
		Ke.resize (dof_indices.size(), dof_indices.size());
		Fee.resize (dof_indices.size());

        RealGradient gra=0;
		for (unsigned int qp=0; qp<qrule.n_points(); qp++)
		{
            RealGradient gr=0;
			for (unsigned int i=0; i<phi.size(); i++)
            {
                unsigned int nod=elem->get_node(i)->id();
				gr+=x[nod]*dphi[i][qp];
            }
            //g=average(gr,gr,gr,gr....) for gr in all gauss points
            gra=(gra*((Real)qp)+gr)/(qp+1);
		} 
        gx[elem->id()] = gra(0);
        gy[elem->id()] = gra(1);
        gz[elem->id()] = gra(2);

	} 

    }
    libMesh::close();
}

double integ(const std::string& meshfile, const std::string& fboundaries, 
        double* x, int xsize,
        double* y, int ysize,
        double* z, int zsize,
        int b, Updater *up) 
{
    double S=0.0;

    int argc=1; char *p="./lmesh\n"; char **argv=&p;
    libMesh::init (argc, argv);
    {    
        Mesh mesh(3);
        mesh.read(meshfile);
        mesh.find_neighbors();
        EquationSystems equation_systems (mesh);
        equation_systems.add_system<LinearImplicitSystem> ("Poisson");
        equation_systems.get_system("Poisson").add_variable("u", FIRST);
        equation_systems.init();
        

	const unsigned int dim = mesh.mesh_dimension();
	LinearImplicitSystem& system=equation_systems.get_system<LinearImplicitSystem>("Poisson");
	const DofMap& dof_map = system.get_dof_map();
	FEType fe_type = dof_map.variable_type(0);
	AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));
	QGauss qrule (dim, FIFTH);
	fe->attach_quadrature_rule (&qrule);
	AutoPtr<FEBase> fe_face (FEBase::build(dim, fe_type));
	QGauss qface(dim-1, FIFTH);
	fe_face->attach_quadrature_rule (&qface);

	DenseMatrix<Number> Ke;
	DenseVector<Number> Fee;
	std::vector<unsigned int> dof_indices;

    BC bc(fboundaries.c_str());

	MeshBase::const_element_iterator       el     = mesh.elements_begin();
	const MeshBase::const_element_iterator end_el = mesh.elements_end();
    if (up!=NULL) up->init(mesh.n_elem()-1);
	for ( ; el != end_el ; ++el)
	{
		const Elem* elem = *el;
        if (up!=NULL) up->update(elem->id());
		dof_map.dof_indices (elem, dof_indices);
		fe->reinit (elem);
		Ke.resize (dof_indices.size(), dof_indices.size());
		Fee.resize (dof_indices.size());

        {
        int bb,s;
        if (bc.find(elem->id()+1,&bb,&s))
		for (unsigned int side=0; side<elem->n_sides(); side++)
			if ((side+1==(unsigned int)s) and (bb==b))
		{
            if (elem->neighbor(side) != NULL) error();
			const std::vector<std::vector<Real> >&  phi_face=fe_face->get_phi();
			const std::vector<Real>& JxW_face = fe_face->get_JxW();
			const std::vector<Point>& normals = fe_face->get_normals();
			fe_face->reinit(elem, side);

//            std::cout << normals[0] << std::endl;
			for (unsigned int qp=0; qp<qface.n_points(); qp++)
			{
                double scalar;
                scalar=x[elem->id()]*normals[qp](0)+
                    y[elem->id()]*normals[qp](1)+
                    z[elem->id()]*normals[qp](2);
                //scalar=y[elem->id()];
                //scalar=myabs(scalar);
                //scalar=1;
				for (unsigned int i=0; i<phi_face.size(); i++)
					S += scalar * JxW_face[qp]*phi_face[i][qp];
			} 

		}
        }

	} 

    }
    libMesh::close();
    return S;
}
