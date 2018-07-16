/*
 * refinement_error_estimator.h
 *
 *  Created on: 2018-06-13
 *      Author: mbreit
 */

#ifndef UG__PLUGINS__NERNST_PLANCK__REFINEMENT_ERROR_ESTIMATOR_H_
#define UG__PLUGINS__NERNST_PLANCK__REFINEMENT_ERROR_ESTIMATOR_H_


#include <cstddef>  // for size_t

#include "common/util/smart_pointer.h"  // for SmartPtr
#include "lib_grid/grid_objects/grid_dim_traits.h"  // for grid_dim_traits
#include "lib_grid/grid/grid.h"  // for AttachmentAccessor
#include "lib_grid/refinement/refiner_interface.h"  // for IRefiner
#include "lib_disc/dof_manager/dof_distribution.h"  // for DoFDistribution
#include "lib_disc/function_spaces/grid_function.h"  // for GridFunction
#include "lib_disc/function_spaces/error_elem_marking_strategy.h"  // for IElementMarkingStrategy

namespace ug {
namespace nernst_planck {

template <typename TDomain, typename TAlgebra>
class RefinementErrorEstimator
{
	public:
		static const int dim = TDomain::dim;
		typedef GridFunction<TDomain, TAlgebra> TGridFunction;

	public:
		RefinementErrorEstimator();
		virtual ~RefinementErrorEstimator();

		number compute_elementwise_errors
		(
			SmartPtr<TGridFunction> uFine,
			SmartPtr<TGridFunction> uCoarse,
			size_t cmp
		);

		SmartPtr<TGridFunction> error_grid_function(SmartPtr<TDomain> dom);

		void mark_with_strategy
		(
			IRefiner& refiner,
			SmartPtr<IElementMarkingStrategy<TDomain > > markingStrategy
		);


	private:
		typedef typename grid_dim_traits<dim>::element_type elem_type;
		typedef MultiGrid::AttachmentAccessor<elem_type, Attachment<number> > aa_type;

		ConstSmartPtr<DoFDistribution> m_spDD;

		Attachment<number> m_aError;
		aa_type m_aaError;

		bool m_bErrorsCalculated;
};

} // namespace nernst_planck
} // namespace ug

#endif // UG__PLUGINS__NERNST_PLANCK__REFINEMENT_ERROR_ESTIMATOR_H_
