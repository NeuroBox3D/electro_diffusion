/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Markus Breit
 * Creation date: 2016-07-11
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

#ifndef UG__PLUGINS__NERNST_PLANCK__PNP_UPWIND_H
#define UG__PLUGINS__NERNST_PLANCK__PNP_UPWIND_H

#include "common/types.h"                                          // for number
#include "common/math/math_vector_matrix/math_matrix.h"            // for MathMatrix
#include "common/math/math_vector_matrix/math_vector.h"            // for MathVector
#include "lib_disc/spatial_disc/disc_util/conv_shape_interface.h"  // for IConvectionShapes


namespace ug {

// forward declaration
template <int N> struct Int2Type;

namespace nernst_planck {


/** @brief Upwind scheme for the PNP problem
 *
 * For the SCVF integration point at coordinates s,
 * the convection shape for the shape function i,
 * which is 1 at position p_i, is calculated as:
 *
 *     \frac{\exp{\left( \beta \frac{v^T (s - p_i)}{(\|v\| + \alpha)}\right)}}
 *          {\sum\limits_j \exp{\left( \beta \frac{v^T (s - p_j)}{(\|v\| + \alpha)}\right)}}
 *     v^T n
 *
 * where v is the velocity at the SCVF integration point
 * and n the normal on the SCVF in that point.
 *
 * The parameter \alpha is used to prevent problems for very small v,
 * it is a problem-specific constant that should be chosen much smaller
 * than important velocity norms (at the appropriate scaling), the
 * default is 1e-3.
 * The larger the value is chosen, the more the upwinding will tend to use
 * concentrations from the center of elements.
 * The parameter \beta is a constant chosen by the developer at \beta = 20.
 *
 * The rationale is to use concentration values that will
 * cause only minimal stationary fluxes at charged surfaces.
 *
 * @note: The exponential might be expensive (number of calls!).
 *        Maybe use (x+1)^n with n = 2^k + 1 instead of \exp{(\beta x)}
 *        to create a strongly nonlinear mapping from (-1, 1) -> (0, \infty).
 */
template <int worldDim>
class PNPUpwind
	: public IConvectionShapes<worldDim>
{
	public:
		typedef IConvectionShapes<worldDim> base_type;
		typedef PNPUpwind<worldDim> this_type;

	protected:
		// explicitly forward some functions
		using base_type::set_non_zero_deriv_diffusion_flag;
		using base_type::conv_shape;
		using base_type::D_vel;
		using base_type::conv_shape_diffusion;
		using base_type::non_zero_deriv_diffusion;
		using base_type::register_update_func;

	public:
		/// constructor
		PNPUpwind();

		/// destructor
		virtual ~PNPUpwind();

		/// set alpha parameter
		void set_alpha(number a);

		/// update convection shape values for FV1Geometry
		template <typename TFVGeom>
		bool update
		(
			const TFVGeom* geo,
			const MathVector<worldDim>* Velocity,
			const MathMatrix<worldDim, worldDim>* DiffDisp,
			bool computeDeriv
		);

	private:
		void register_func(Int2Type<1>);
		void register_func(Int2Type<2>);
		void register_func(Int2Type<3>);

		template <typename TElem>
		void register_func_for_elem();

		template <int refDim>
		void register_func_for_refDim();

	protected:
		number m_alpha;

	private:
		std::vector<MathVector<worldDim> > m_vDir;
		std::vector<number> m_vDirVelProdExp;
};


} // namespace nernst_planck
} // namespace ug

#endif // UG__PLUGINS__NERNST_PLANCK__PNP_UPWIND_H

