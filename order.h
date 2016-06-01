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


} // namespace nernst_planck
} // namespace ug


#endif // __UG__PLUGINS__NERNST_PLANCK__ORDER_H__
