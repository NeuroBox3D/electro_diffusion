/*
 * order.h
 *
 *  Created on: 12.08.2015
 *      Author: mbreit
 */

#ifndef __UG__PLUGINS__NERNST_PLANCK__ORDER_H__
#define __UG__PLUGINS__NERNST_PLANCK__ORDER_H__

#include <vector>

#include "lib_disc/function_spaces/approximation_space.h"


namespace ug {
namespace nernst_planck {


template <typename TBaseElem>
void collect_obj_indices
(
	std::vector<int>& vObjInd,
	SmartPtr<DoFDistribution> dd,
	ConstSmartPtr<MGSubsetHandler> sh,
	const SubsetGroup& csg
);

template <typename TDomain>
void reorder_dofs(SmartPtr<ApproximationSpace<TDomain> > approxSpace, const char* constrained);


/**
 * This function implements algorithm LEX M from
 * Rose et al.: "Algorithmic aspects of vertex elimination on graphs"
**/
void lex_order
(
	std::vector<size_t>& newIndOut,
	const std::vector<std::vector<size_t> >& vAdj,
	bool preserveConsec = true
);

template <typename TDomain>
void reorder_dof_distros_lex(SmartPtr<ApproximationSpace<TDomain> > approx);



} // namespace nernst_planck
} // namespace ug

#endif // __UG__PLUGINS__NERNST_PLANCK__ORDER_H__
