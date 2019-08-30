/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Markus Breit
 * Creation date: 2013-03
 *
 * This file is part of NeuroBox, which is based on UG4.
 *
 * NeuroBox and UG4 are free software: You can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 3
 * (as published by the Free Software Foundation) with the following additional
 * attribution requirements (according to LGPL/GPL v3 §7):
 *
 * (1) The following notice must be displayed in the appropriate legal notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 *
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 *
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating PDE based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * "Stepniewski, M., Breit, M., Hoffer, M. and Queisser, G.
 *   NeuroBox: computational mathematics in multiscale neuroscience.
 *   Computing and visualization in science (2019).
 * "Breit, M. et al. Anatomically detailed and large-scale simulations studying
 *   synapse loss and synchrony using NeuroBox. Front. Neuroanat. 10 (2016), 8"
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */


#ifndef UG__PLUGINS__NERNST_PLANCK__UTIL__EDL_1D_H
#define UG__PLUGINS__NERNST_PLANCK__UTIL__EDL_1D_H

#include <cstddef>                                                    // for size_t
#include <string>                                                     // for string
#include <vector>                                                     // for vector

#include "common/types.h"                                             // for number
#include "lib_algebra/small_algebra/small_matrix/densematrix.h"       // for DenseMatrix
#include "lib_algebra/small_algebra/small_matrix/densevector.h"       // for DenseVector
#include "lib_algebra/small_algebra/storage/variable_array.h"         // for VariableArray2, Var...


namespace ug {
namespace nernst_planck {


/**
 * This class implements a simple FD discretization of a rotationally symmetrical
 * PNP problem (cylindrical cable).
 * Domain: - cylinder with radius a: cytosol
 *         - radius d around it: membrane
 *         - radius A around that: outer medium
 * Problem in 1D:
 * c_i' = -z_i*F/(RT)*c_i*phi'
 * phi'' + 1/r*phi' = -sum_i(z_i*F*c_i / (eps_0*eps_r))
 * where r is the cylindrical radius coordinate.
**/

class EDLSimulation
{
	public:
		typedef DenseVector< VariableArray1<number> > Vector;
		typedef DenseMatrix< VariableArray2<number> > Matrix;

		const number R;		// universal gas constant
		const number T;		// temperature
		const number F;		// Faraday constant
		const number eps_0;	// vacuum permettivity
		const number rho_f;	// charge density membrane

		EDLSimulation(number rhof = -0.03) : R(8.314), T(310.0), F(96485.0), eps_0(8.854e-12), rho_f(rhof) {};

		void set_geom_specs(number a, number d, number A, number h);
		void set_bnd_cond(std::vector<number> inner, std::vector<number> outer);
		void set_constants(std::vector<number> val, std::vector<number> perm);

		// set verbosity level
		void set_verbosity_level(size_t lvl);

		// set minimal reduction to be reached for defect in Newton iteration
		void set_min_reduction(number mr);

		// set defect norm minimum (if defect is lower,
		// then Newton iteration is considered converged regardless of the defect reduction)
		void set_min_defect(number md);

		// set maximum number of iterations for Newton method
		void set_max_iter(size_t mi);

		// set output file name
		void set_output_file(std::string file);


		void compute_solution();

		void calc_EDL_ions();

	protected:
		// assembles the defect vector
		void defect(Vector& d, Vector& u);

		// assembles the Jacobian matrix
		void jacobian(Matrix& J, Vector& u);

	protected:
		size_t m_verbLvl;			// verbosity level
		number m_minRed;				// minimal defect reduction Newton
		number m_minDefect;				// minimal defect norm Newton
		size_t m_maxIt;			// maximum iterations Newton

		std::string m_outFile;			// output file name

	private:
		number a;						// radius dendrite
		number d;						// thickness membrane
		number A;						// radius outer medium
		size_t N1, N2, N3;
		Vector u;						// unknowns
		std::vector<number> grid;
		std::vector<number> valency;
		number eps_r1;					// permettivity dendrite
		number eps_r2;					// permettivity outer medium
		number eps_r3;					// permettivity membrane
};


} // namespace nernst_planck
} // namespace ug


#endif // UG__PLUGINS__NERNST_PLANCK__UTIL__EDL_1D_H
