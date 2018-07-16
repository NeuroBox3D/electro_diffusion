/*
 * pnp_upwind.h
 *
 *  Created on: 11.07.2016
 *      Author: mbreit
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

