/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Markus Breit
 * Creation date: 2018-06-13
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

#ifndef UG__PLUGINS__NERNST_PLANCK__UTIL__REFINEMENT_ERROR_ESTIMATOR_H
#define UG__PLUGINS__NERNST_PLANCK__UTIL__REFINEMENT_ERROR_ESTIMATOR_H


#include <cstddef>  // for size_t

#include "common/util/smart_pointer.h"  // for SmartPtr
#include "lib_algebra/cpu_algebra_types.h"  // for CPUAlgebra
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
		typedef GridFunction<TDomain, CPUAlgebra> TErrorGridFunction;

	public:
		RefinementErrorEstimator();
		virtual ~RefinementErrorEstimator();

		number compute_elementwise_errors
		(
			SmartPtr<TGridFunction> uFine,
			SmartPtr<TGridFunction> uCoarse,
			size_t cmp
		);

		SmartPtr<TErrorGridFunction> error_grid_function(SmartPtr<TDomain> dom);

		void mark_with_strategy
		(
			IRefiner& refiner,
			SmartPtr<IElementMarkingStrategy<TDomain > > markingStrategy
		);


	private:
		typedef typename grid_dim_traits<dim>::element_type elem_type;
		typedef MultiGrid::AttachmentAccessor<elem_type, Attachment<number> > aa_type;

		ConstSmartPtr<DoFDistribution> m_spDD;

		IMultigridElementIndicators<TDomain> m_mgElemErrors;  // holding the indicators

		bool m_bErrorsCalculated;
};

} // namespace nernst_planck
} // namespace ug

#endif // UG__PLUGINS__NERNST_PLANCK__UTIL__REFINEMENT_ERROR_ESTIMATOR_H
