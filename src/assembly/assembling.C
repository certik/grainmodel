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
				Ve(i,j) += JxW[qp]*(phi[i][qp]*phi[j][qp]);
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
