/*
 * nernst_planck_util.h
 *
 *  Created on: 17.07.2014
 *      Author: mbreit
 */

#ifndef UG__PLUGINS__EXPERIMENTAL__NERNST_PLANCK__NERNST_PLANCK_UTIL_H
#define UG__PLUGINS__EXPERIMENTAL__NERNST_PLANCK__NERNST_PLANCK_UTIL_H


#include "common/common.h"
#include "lib_disc/function_spaces/approximation_space.h"
#include "lib_disc/function_spaces/grid_function.h"

#include <iostream>
#include <fstream>


namespace ug{
namespace nernst_planck{


template <typename TGridFunction>
number writeResidualsToFile
(
	SmartPtr<TGridFunction> sol1,
	SmartPtr<TGridFunction> sol2,
	const char* cmp,
	const char* fileName
);


/// adjusts interface after a refinement of the geometry
/**
 * When a geometry with an interface is refined, the interface node on the full-dimensional
 * side will no longer be correctly located on the top level (one layer too far from the interface).
 * We therefore need to re-locate this interface node every time we refine the interface (globally).
 */
template <typename TDomain>
void adjust_geom_after_refinement
(
	SmartPtr<ApproximationSpace<TDomain> > approx,
	const char* innerSubset,
	const char* fullDimIntfNodeSubset,
	const char* lowDimIntfNodeSubset
);



} // namspace calciumDynamics
} // namespace ug

#include "nernst_planck_util_impl.h"

#endif // UG__PLUGINS__EXPERIMENTAL__NERNST_PLANCK__NERNST_PLANCK_UTIL_H
