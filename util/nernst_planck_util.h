/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Markus Breit
 * Creation date: 2014-07-17
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

#ifndef UG__PLUGINS__NERNST_PLANCK__UTIL__NERNST_PLANCK_UTIL_H
#define UG__PLUGINS__NERNST_PLANCK__UTIL__NERNST_PLANCK_UTIL_H

#include <vector>                       // for vector

#include "common/types.h"               // for number
#include "common/util/smart_pointer.h"  // for SmartPtr, ConstSmartPtr


namespace ug {

// forward declarations
class IRefiner;
template <typename TDomain> class ApproximationSpace;

namespace nernst_planck {


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
	const char* fullDimIntfNodeSubset,
	const char* lowDimIntfNodeSubset
);




} // namspace nernst_planck
} // namespace ug


#endif // UG__PLUGINS__NERNST_PLANCK__UTIL__NERNST_PLANCK_UTIL_H
