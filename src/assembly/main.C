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

#include "mesh.h"
#include "mesh_data.h"
#include "gmv_io.h"
//#include "eigen_system.h"
#include "linear_implicit_system.h"
#include "equation_systems.h"

#include "configfile.h"
#include "assembling.h"

Config config("tmp/ex3.xml");

int main (int argc, char** argv)
{
	libMesh::init (argc, argv);
	{
		std::cout << "Running " << argv[0];
		for (int i=1; i<argc; i++)
			std::cout << " " << argv[i];
		std::cout << std::endl << std::endl;

        config.load();
        config.print();

		PerfLog perf("Main Program");
		perf.start_event("program init");
		Mesh mesh (config.dim);
		MeshData mesh_data(mesh);
		//mesh_data.activate();
		mesh_data.enable_compatibility_mode();

		mesh.read("tmp/in.xda",&mesh_data);
        //mesh.find_neighbors();
		//mesh_data.read("data.xta");
		//mesh.print_info();

		EquationSystems equation_systems (mesh,&mesh_data);
        LinearImplicitSystem & eigen_system =
	    	equation_systems.add_system<LinearImplicitSystem> ("Poisson");
	    	//equation_systems.add_system<EigenSystem> ("Poisson");
		equation_systems.get_system("Poisson").add_variable("u", FIRST);
		equation_systems.get_system("Poisson").attach_assemble_function (assemble_poisson);
        unsigned int nev = config.eigs;
        equation_systems.parameters.set<unsigned int>("eigenpairs")    = nev;
        equation_systems.parameters.set<unsigned int>("basis vectors") = nev*3;
//        eigen_system.eigen_solver-> set_eigensolver_type(ARNOLDI);
        //eigen_system.eigen_solver-> set_eigensolver_type(SUBSPACE);
        //eigen_system.eigen_solver-> set_eigensolver_type(POWER);
//        eigen_system.eigen_solver-> set_eigensolver_type(LANCZOS);
//        eigen_system.set_eigenproblem_type(GHEP);
//        eigen_system.eigen_solver->set_position_of_spectrum(SMALLEST_MAGNITUDE);
        //eigen_system.eigen_solver->set_position_of_spectrum(LARGEST_MAGNITUDE);
        equation_systems.parameters.set<Real>("linear solver tolerance") = 
            pow(TOLERANCE, 5./3.);
        equation_systems.parameters.set<unsigned int>
	        ("linear solver maximum iterations") = 1000;
		equation_systems.init();
		//equation_systems.print_info();
		perf.stop_event("program init");

        if (config.assemble) {
            assemble_poisson(equation_systems, "Poisson");
        }
        if (!config.printlog) perf.clear();
	}
	return libMesh::close();
}
