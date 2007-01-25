/* Quantum */
/* Copyright (C) 2006  Ondrej Certik */

/* This library is free software; you can redistribute it and/or */
/* modify it under the terms of the GNU Lesser General Public */
/* License as published by the Free Software Foundation; either */
/* version 2.1 of the License, or (at your option) any later version. */

/* This library is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU */
/* Lesser General Public License for more details. */

/* You should have received a copy of the GNU Lesser General Public */
/* License along with this library; if not, write to the Free Software */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */

#include <fstream>

#include "mesh.h"
#include "mesh_data.h"
#include "mesh_tools.h"
#include "boundary_info.h"
#include "eigen_system.h"
#include "equation_systems.h"

#include "fe.h"
#include "quadrature_gauss.h"

#include "sparse_matrix.h"
#include "petsc_matrix.h"
#include "dense_matrix.h"
#include "elem.h"

#include "dof_map.h"

#include "configfile.h"

extern Config config;
Real pi=2*asin(1);

class pseudopotentials
{
public:
    pseudopotentials(const char* fname)
    {
        std::ifstream *f;
        f=new std::ifstream(fname);
        int n=500;
#define load(vec) vec=new double[n]; for (int i=0;i<n;i++) *f >> vec[i];
        load(r);
        load(U3);
        load(Vloc);
        load(V0);
        load(V1);
        load(V2);
        delete f;
    };
    ~pseudopotentials()
    {
        delete r;
        delete U3;
        delete Vloc;
        delete V0;
        delete V1;
        delete V2;
    };
    double getVloc(double r)
    {
        return Vloc[getr(r)];
    };
    double getVhVxc(double r)
    {
        return U3[getr(r)];
    };
    double getVl(int l,double r)
    {
        switch (l)
        {
            case 0: return V0[getr(r)];
            case 1: return V1[getr(r)];
            case 2: return V2[getr(r)];
            default: error();
        }
    };
    int getr(double r)
    {
        int imax=499;
        int imin=0;
        int i;
        assert(r>this->r[imin]);
        assert(r<this->r[imax]);
        while (imax-imin>1)
        {
            i=(imin+imax)/2;
            if (r>this->r[i]) 
                imin=i;
            else
                imax=i;
        }
        return i;
    };
    int getr2(double r)
    {
        int i=0;
        while ((i<500-1) && (r> this->r[i])) i++;
        return i;
    };
private:
    double *r,*U3,*Vloc,*V0,*V1,*V2;
};

pseudopotentials Vpseudo("data/potentials.dat");

//Legendre polynomials
inline Real P(int l, Real x)
{
    switch (l)
    {
        case 0: return 1.0;
        case 1: return x;
        case 2: return 0.5*(3*x*x-1);
        case 3: return 0.5*x*(5*x*x-3);
        case 4: return 1.0/8*(35*pow(x,4)-30*x*x+3);
        case 5: return 1.0/8*x*(63*pow(x,4)-70*x*x+15);
        case 6: return 1.0/16*(231*pow(x,6)-315*pow(x,4)+105*x*x-5);
        default: 
        {
            std::cout << "Legendre polynomial for l="<<l<<" not implemented!"
                << std::endl;
            error();
        }
    }
}

inline int get_local_id(const Elem *elem, int i)
{
    for (int s=0;s<elem->n_nodes();s++)
        if (elem->node(s)==i) return s;
    std::cout << "node " << i << " not found in elem!" << std::endl;
    error();
}

//convert a point "p" from (x,y) to (ksi,eta), ie. to a reference element
inline Real phi_value(int i, Point p, AutoPtr<FEBase> &fe,
    const Elem *elem_orig, 
    std::vector<std::vector<const Elem *> > &nodes_to_elem_map)
{
    //FEInterface::inverse_map
    
    //is "p" on the element "elem_orig"?
    const Elem *elem=elem_orig;
    Point p_ref=FE<3,LAGRANGE>::inverse_map(elem,p,0.01,true);
    if (FEBase::on_reference_element(p_ref,elem->type()))
        //yes -- return phi_i(p)
        return (FE<3,LAGRANGE>::shape(elem->type(),fe->get_order(),i,p_ref));

    //comment out to try other elements as well.
    //return 0.0;

    //no -- try all other elements containing the global shape function "i".
    for (int s=0;s<nodes_to_elem_map[elem_orig->node(i)].size();s++)
    {
        elem=nodes_to_elem_map[elem_orig->node(i)][s];
        //we already tried elem_orig, skip it
        if (elem!=elem_orig)
        {
            p_ref=FE<3,LAGRANGE>::inverse_map(elem,p,0.01,true);
            if (FEBase::on_reference_element(p_ref,elem->type()))
                return (FE<3,LAGRANGE>::shape(elem->type(),fe->get_order(),
                            get_local_id(elem,elem_orig->node(i)),p_ref));
        }
    }
    return 0.0;
}

//computes phi and theta of the point "p"
inline void compute_polar(Point p, Real *phi, Real *theta)
{
    Real x,y,z;
    x=p(0);y=p(1);z=p(2);
    /*
    Real r2=sqrt(x*x+y*y);
    if (y==0.0)
        *phi=0;
    else
        *phi=2*atan((r2-x)/y);
    if (*phi<0) *phi+=2*pi;
    //Real r=sqrt(x*x+y*y+z*z);
    // *theta=acos(z/r);
    *theta=pi/2-atan(z/r2);
    */
    *phi=atan(y/x);
    *theta=atan(sqrt(x*x+y*y)/z); //   =acos(z/r); 
}

inline void compute_box(std::vector<const Elem *> &elems,
        Real *phi1, Real *phi2, Real *theta1, Real *theta2)
{
    *phi1=2*pi;
    *phi2=0;
    *theta1=pi;
    *theta2=0;
    for (unsigned int s=0; s<elems.size(); s++)
        for (int n=0;n<elems[s]->n_nodes();n++)
        {
            Real phi,theta;
            compute_polar(elems[s]->point(n),&phi,&theta);
            if (phi<*phi1) *phi1=phi;
            if (phi>*phi2) *phi2=phi;
            if (theta<*theta1) *theta1=theta;
            if (theta>*theta2) *theta2=theta;
        }
}

inline void compute_box2(const Elem *elem,
        Real *phi1, Real *phi2, Real *theta1, Real *theta2)
{
    *phi1=2*pi;
    *phi2=0;
    *theta1=pi;
    *theta2=0;
    for (int n=0;n<elem->n_nodes();n++)
    {
        Real phi,theta;
        compute_polar(elem->point(n),&phi,&theta);
        if (phi<*phi1) *phi1=phi;
        if (phi>*phi2) *phi2=phi;
        if (theta<*theta1) *theta1=theta;
        if (theta>*theta2) *theta2=theta;
    }
}

//Computes the value of V*phi_i at the quadrature point qp 
//  for various potentials according to config.potential
inline Real Vphi(int i, int qp, AutoPtr<FEBase> &fe, const Elem *elem, 
        QGauss &gauss,
        std::vector<std::vector<const Elem *> > &nodes_to_elem_map,
        Point localp)
{
    //global coordinates of the quadrature point qp
    Real x=fe->get_xyz()[qp](0); 
    Real y=fe->get_xyz()[qp](1); 
    Real z=fe->get_xyz()[qp](2); 
    Real r=sqrt(x*x+y*y+z*z);
    Real Vxpsi=0.;
    switch (config.potential) {
        case POT_WELL: break;
        case POT_LHO: {
            Real w=config.lhofrequency;
            Real m=config.m;
            Vxpsi=0.5*m*w*w*r*r; 
            break;
        }
        case POT_HYDROGEN: {
            if (r!=0.0)
            {
                r=r*config.scale;
                Vxpsi=-1.0/r;
            }
            break;
        }
        case POT_DFT: {
            Vxpsi= Vpseudo.getVloc(r)+Vpseudo.getVhVxc(r);
            break;
        }
        default: 
        {
            std::cout << "Potential type not implemented!" << std::endl;
            error();
        }
    }
    Vxpsi+=config.Eshift;

    Vxpsi*= FE<3,LAGRANGE>::shape(elem->type(),fe->get_order(),i,localp);
    //for scalar potentials, we are already done
    if (config.potential!=POT_DFT) return Vxpsi; 

    return Vxpsi;

    Real phi1,phi2,theta1,theta2;
    bool boxcomputed=false;

#define phi_xyz(i,x,y,z) phi_value(i,Point(x,y,z),fe,elem,nodes_to_elem_map)
#define integrand(phi,theta) \
(P(l,(x*sin(theta)*cos(phi)+y*sin(theta)*sin(phi)+z*cos(theta))/r)\
*phi_xyz(i,r*sin(theta)*cos(phi),r*sin(theta)*sin(phi),r*cos(theta))\
*pow(sin(theta),2))

    for (int l=0;l<3;l+=1)
    {
        Real V=Vpseudo.getVl(l,r);
        if (std::abs(V)<1e-3) continue;
        if (!boxcomputed)
        {
            compute_box(nodes_to_elem_map[elem->node(i)],
                    &phi1,&phi2,&theta1,&theta2);
            //compute_box2(elem, &phi1,&phi2,&theta1,&theta2);
            boxcomputed=true;
        }
        Real I=0.0;
        for (unsigned int q=0; q<gauss.n_points(); q++)
        {
            Real phi=(gauss.qp(q)(0)+1)*(phi2-phi1)/2.0+phi1;
            Real theta=(gauss.qp(q)(1)+1)*(theta2-theta1)/2.0+theta1;
            I+=gauss.w(q)*integrand(phi,theta);
            //I+=gauss.w(q)*1.0;
        } 
        I*=(phi2-phi1)*(theta2-theta1)/4.0;
        Vxpsi+=V*I*(2*l+1)/(4*pi);
    }

    return Vxpsi;
}

void save_matrix(std::ostream& f, DenseMatrix<Number>& M)
{
    M.print_scientific(f);
}

void save_sparse_matrix(SparseMatrix<Number>& M, const char *fname)
{
    //M.print_matlab(fname);

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

    Mat mat = ((PetscMatrix<Number>&)(M)).mat();

    ierr = MatView (mat, petsc_viewer);
    CHKERRABORT(libMesh::COMM_WORLD,ierr);

    ierr = PetscViewerDestroy (petsc_viewer);
    CHKERRABORT(libMesh::COMM_WORLD,ierr);
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

void save_topo(std::ostream& f, const Elem* el,
        std::vector<unsigned int>& dof_indices)
{
    f << el->id() << " ";
    for (int i=0;i<el->n_nodes();i++)
        //save the original node numbers
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
	EigenSystem& system=es.get_system<EigenSystem>("Poisson");
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
	DenseMatrix<Number> Ve;
	DenseMatrix<Number> Me;
	std::vector<unsigned int> dof_indices;

    SparseMatrix<Number>&  matrix_A = *system.matrix_A;
    SparseMatrix<Number>&  matrix_B = *system.matrix_B;

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
	std::vector<unsigned int> zeronodes;
    load_zeronodes("tmp/zeronodes.gmsh", zeronodes);

	MeshBase::const_element_iterator       el     = mesh.elements_begin();
	const MeshBase::const_element_iterator end_el = mesh.elements_end();
	for ( ; el != end_el ; ++el)
	{
		perf.start_event("elem init");
		const Elem* elem = *el;
        if (elem->id() % 100 == 0)
            std::cout << 100.0*elem->id()/mesh.n_elem() << "%" << std::endl;
		dof_map.dof_indices (elem, dof_indices);
		fe->reinit (elem);
		Ke.resize (dof_indices.size(), dof_indices.size());
		Ve.resize (dof_indices.size(), dof_indices.size());
		Me.resize (dof_indices.size(), dof_indices.size());
		perf.stop_event("elem init");

		perf.start_event("Ke");
		for (unsigned int qp=0; qp<qrule.n_points(); qp++)
		{
            Real hbar=config.hbar;
            Real m=config.m;
            Real C=hbar*hbar/(2.0*m);
			for (unsigned int i=0; i<phi.size(); i++)
			for (unsigned int j=0; j<phi.size(); j++)
				Ke(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp])*C;
		} 
		perf.stop_event("Ke");

		perf.start_event("Ve");
		for (unsigned int qp=0; qp<qrule.n_points(); qp++)
		{
			for (unsigned int i=0; i<phi.size(); i++)
			for (unsigned int j=0; j<phi.size(); j++)
				//Ve(i,j) += JxW[qp]*(phi[i][qp]*phi[j][qp])*V;
				Ve(i,j) += JxW[qp]*phi[i][qp]*Vphi(j,qp,fe,elem,qquad,
                    nodes_to_elem_map,qrule.qp(qp));
		} 
		perf.stop_event("Ve");

		perf.start_event("Me");
		for (unsigned int qp=0; qp<qrule.n_points(); qp++)
		{
			for (unsigned int i=0; i<phi.size(); i++)
			for (unsigned int j=0; j<phi.size(); j++)
				Me(i,j) += JxW[qp]*(phi[i][qp]*phi[j][qp]);
		} 
		perf.stop_event("Me");

		perf.start_event("matrix insertion");
        matrix_A.add_matrix (Ke, dof_indices);
        matrix_A.add_matrix (Ve, dof_indices);
        matrix_B.add_matrix (Me, dof_indices);
		perf.stop_event("matrix insertion");

		perf.start_event("matrix saving");
        //save matrices A and M to "matrices.dat"
        /*
        save_matrix(f,Ke);
        save_matrix(f,Ve);
        save_matrix(f,Me);
        */
        //save the topology to "topo.dat"
        save_topo(ftopo,elem,dof_indices);
		perf.stop_event("matrix saving");
	} //for element

    //assing the value 0 to every node from "zeronodes" using a penalty method
    for (int i=0;i<zeronodes.size();i++)
    {
        unsigned int nd =mesh.node(zeronodes[i]-1).dof_number(0,0,0);
        //std::cout << nd << " ";
        const Real penalty = 1.e10;
        matrix_A.add(nd,nd,penalty*1001.5);
        matrix_B.add(nd,nd,penalty);
    }

    //print matrices A and M
    std::cout << "saving matrices..." << std::endl;
    save_sparse_matrix(matrix_A,"tmp/matA.petsc");
    save_sparse_matrix(matrix_B,"tmp/matM.petsc");
    //matrix_A.print_matlab("tmp/matA.matlab");
    //matrix_B.print_matlab("tmp/matM.matlab");
    if (!config.printlog) perf.clear();
    std::cout << "done." << std::endl;
}
