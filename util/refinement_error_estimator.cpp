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

#include "refinement_error_estimator.h"

#include "lib_algebra/cpu_algebra_types.h"  // for CPUAlgebra etc.
#include "lib_disc/domain.h"  // for Domain1d etc.
#include "lib_disc/function_spaces/integrate.h"  // for Integrate
//#include "lib_disc/domain_util.h"  // for CollectCornerCoordinates
//#include "lib_disc/reference_element/reference_mapping_provider.h"  // for ReferenceMappingProvider
//#include "lib_disc/quadrature/quadrature_provider.h"  // for QuadratureProvider

#include <vector>


namespace ug {
namespace nernst_planck {


template <typename TDomain, typename TAlgebra>
RefinementErrorEstimator<TDomain, TAlgebra>::RefinementErrorEstimator()
: m_spDD(SPNULL), m_bErrorsCalculated(false)
{}


template <typename TDomain, typename TAlgebra>
RefinementErrorEstimator<TDomain, TAlgebra>::~RefinementErrorEstimator()
{}


template <typename TDomain, typename TAlgebra>
number RefinementErrorEstimator<TDomain, TAlgebra>::
compute_elementwise_errors(SmartPtr<TGridFunction> uFine, SmartPtr<TGridFunction> uCoarse, size_t cmp)
{
	SmartPtr<MultiGrid> spMG = uFine->domain()->grid();
	if (!spMG->has_attachment<elem_type>(m_aError))
		spMG->template attach_to_dv<elem_type>(m_aError, -1.0);
	m_aaError = aa_type(*spMG, m_aError);

	// prepare error computation
	//L2DistIntegrand<TGridFunction> integrand(*uFine, cmp, *uCoarse, cmp);
	H1DistIntegrand<TGridFunction> integrand(*uFine, cmp, *uCoarse, cmp);
	int quadOrder = 3;

	typename TDomain::position_accessor_type& aaPos = uFine->domain()->position_accessor();

	m_spDD = uFine->dd();
	typedef typename DoFDistribution::dim_traits<dim>::const_iterator elem_it;
	elem_it it = m_spDD->begin<elem_type>();
	elem_it itEnd = m_spDD->end<elem_type>();

	// compute errors
	number diffNormSq = 0.0;
	try {diffNormSq = Integrate<dim, dim>(it, itEnd, aaPos, integrand, quadOrder, "best", &m_aaError);}
	UG_CATCH_THROW("Failed trying to calculate element errors via Integrate().");

	m_bErrorsCalculated = true;

	return sqrt(diffNormSq);


	// Integrate assigns (instead of adds) values to m_aaError.
	// So in order to compute a norm for all components, we might need the code below again.
#if 0
	for (; it != itEnd; ++it)
	{
		elem_type* elem = *it;

	// compute error on this element //
		// get reference object id
		ReferenceObjectID roid = elem->reference_object_id();

		// get quadrature rule for ROID and order
		const QuadratureRule<dim>* quadRule;
		try {quadRule = &QuadratureRuleProvider<dim>::get(roid, 3);}
		UG_CATCH_THROW("Could not get a quadrature rule for roid " << roid << " and order 3.");

		const size_t numIP = quadRule->size();

		// get reference element mapping by reference object id
		DimReferenceMapping<dim, dim>* mapping;
		try {mapping = &ReferenceMappingProvider::get<dim, dim>(roid);}
		UG_CATCH_THROW("Could not get a reference mapping for roid " << roid << ".");

		// update the reference mapping for the corners
		CollectCornerCoordinates(vCorner, *elem, *uFine->domain(), true);
		mapping->update(&vCorner[0]);

		// compute global integration points and transformation matrices
		vGlobIP.resize(numIP);
		mapping->local_to_global(&(vGlobIP[0]), quadRule->points(), numIP);
		vJT.resize(numIP);
		mapping->jacobian_transposed(&(vJT[0]), quadRule->points(), numIP);

		// compute integrand values at integration points
		vValue.resize(numIP);
		try
		{
			integrand.values(&(vValue[0]), &(vGlobIP[0]), elem, &vCorner[0],
				quadRule->points(), &(vJT[0]), numIP);
		}
		UG_CATCH_THROW("Unable to compute values of integrand at integration points.");

		// now perform the quadrature
		number& error = m_aaError[elem];
		error = 0;
		for (size_t ip = 0; ip < numIP; ++ip)
		{
			const number weightIP = quadRule->weight(ip);
			const number det = SqrtGramDeterminant(vJT[ip]);
			error += vValue[ip] * weightIP * det;
		}

		// add to global sum
		h1DiffNormSq += intValElem;
		if (aaElemContribs.valid())
				aaElemContribs[pElem] = intValElem;
	}
#endif
}


template <typename TDomain, typename TAlgebra>
SmartPtr<typename RefinementErrorEstimator<TDomain, TAlgebra>::TGridFunction>
RefinementErrorEstimator<TDomain, TAlgebra>::error_grid_function(SmartPtr<TDomain> dom)
{
	// check that error indicators have been calculated
	UG_COND_THROW(!m_bErrorsCalculated, "Error indicators have to be "
		"calculated first by a call to 'compute_elementwise_errors()'.");

	typedef typename DoFDistribution::dim_traits<dim>::const_iterator elem_it;

	SmartPtr<ApproximationSpace<TDomain> > spApprox =
		make_sp(new ApproximationSpace<TDomain>(dom, AlgebraType(AlgebraType::CPU, 1)));
	spApprox->add("eta2", "piecewise-constant");
	SmartPtr<TGridFunction> u_vtk = make_sp(new TGridFunction(spApprox));
	u_vtk->set(0.0);

	std::vector<DoFIndex> vDI;

	// map attachments to grid function
	ConstSmartPtr<DoFDistribution> dd = u_vtk->dd();
	elem_it it = dd->begin<elem_type>();
	elem_it itEnd = dd->end<elem_type>();
	for (; it != itEnd; ++it)
	{
		// get element DoF index
		dd->inner_dof_indices(*it, 0, vDI, true);
		UG_COND_THROW(vDI.size() != 1, "Not exactly 1 DoF index on full-dim element "
			"for piecewise constant approximation space.")

		// assign error value
		DoFRef(*u_vtk, vDI[0]) = m_aaError[*it];
	}

	return u_vtk;
}


template <typename TDomain, typename TAlgebra>
void RefinementErrorEstimator<TDomain, TAlgebra>::mark_with_strategy
(
	IRefiner& refiner,
	SmartPtr<IElementMarkingStrategy<TDomain> > markingStrategy
)
{
	// check that error indicators have been calculated
	UG_COND_THROW(!m_bErrorsCalculated, "Error indicators have to be "
		"calculated first by a call to 'compute_elementwise_errors()'.");

	// mark elements for coarsening
	UG_COND_THROW(!markingStrategy.valid(), "No valid marking strategy given.")
	markingStrategy->mark(m_aaError, refiner, m_spDD);

	m_bErrorsCalculated = false;
}



// explicit template specializations
#ifdef UG_CPU_1
	#ifdef UG_DIM_1
		template class RefinementErrorEstimator<Domain1d, CPUAlgebra>;
	#endif
	#ifdef UG_DIM_2
		template class RefinementErrorEstimator<Domain2d, CPUAlgebra>;
	#endif
	#ifdef UG_DIM_3
		template class RefinementErrorEstimator<Domain3d, CPUAlgebra>;
	#endif
#endif
#ifdef UG_CPU_5
	#ifdef UG_DIM_1
		template class RefinementErrorEstimator<Domain1d, CPUBlockAlgebra<5> >;
	#endif
	#ifdef UG_DIM_2
		template class RefinementErrorEstimator<Domain2d, CPUBlockAlgebra<5> >;
	#endif
	#ifdef UG_DIM_3
		template class RefinementErrorEstimator<Domain3d, CPUBlockAlgebra<5> >;
	#endif
#endif
#ifdef UG_CPU_6
	#ifdef UG_DIM_1
		template class RefinementErrorEstimator<Domain1d, CPUBlockAlgebra<6> >;
	#endif
	#ifdef UG_DIM_2
		template class RefinementErrorEstimator<Domain2d, CPUBlockAlgebra<6> >;
	#endif
	#ifdef UG_DIM_3
		template class RefinementErrorEstimator<Domain3d, CPUBlockAlgebra<6> >;
	#endif
#endif


} // namespace nernst_planck
} // namespace ug

