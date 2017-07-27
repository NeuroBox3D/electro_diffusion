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


template <int worldDim>
class PNPUpwind
	: public IConvectionShapes<worldDim>
{
	public:
	///	Base class
		typedef IConvectionShapes<worldDim> base_type;

	///	This class
		typedef PNPUpwind<worldDim> this_type;

	protected:
	//	explicitly forward some function
		using base_type::set_non_zero_deriv_diffusion_flag;
		using base_type::conv_shape;
		using base_type::D_vel;
		using base_type::conv_shape_diffusion;
		using base_type::non_zero_deriv_diffusion;
		using base_type::register_update_func;

	public:
	///	constructor
		PNPUpwind();

	/// destructor
		virtual ~PNPUpwind();

	/// set exponential weighing factor
		void set_exp_factor(number a);

	///	update of values for FV1Geometry
		template <typename TFVGeom>
		bool update(const TFVGeom* geo,
					const MathVector<worldDim>* Velocity,
					const MathMatrix<worldDim, worldDim>* DiffDisp,
		            bool computeDeriv);

	private:
		void register_func(Int2Type<1>);
		void register_func(Int2Type<2>);
		void register_func(Int2Type<3>);

		template <typename TElem>
		void register_func_for_elem();

		template <int refDim>
		void register_func_for_refDim();

	protected:
		number m_expFact;
};


} // namespace nernst_planck
} // namespace ug

#endif // UG__PLUGINS__NERNST_PLANCK__PNP_UPWIND_H

