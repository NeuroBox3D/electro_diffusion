/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Markus Breit
 * Creation date: 2016-12-18
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


#ifndef UG__PLUGINS__NERNST_PLANCK__UTIL__CHOOSE_FVGEOM_H
#define UG__PLUGINS__NERNST_PLANCK__UTIL__CHOOSE_FVGEOM_H

#include "lib_disc/spatial_disc/disc_util/fv1_geom.h"                // for FV1Geometry
#include "lib_disc/spatial_disc/disc_util/fvho_geom.h"               // for FVGeometry, DimFVGeometry
#include "lib_disc/spatial_disc/disc_util/hfv1_geom.h"               // for HFV1Geometry

namespace ug {
namespace nernst_planck {


// Choose FVGeom as in ConvectionDiffusionFV.

// standard FVGeom
template <int dim, typename TElem, bool hanging, int order, int quadOrder = order+1>
struct ChooseProperFVGeom
{
	typedef DimFVGeometry<dim> FVGeom;
};

// for order 1 and !hanging: use FV1Geom
template <int dim, typename TElem, int quadOrder>
struct ChooseProperFVGeom<dim, TElem, false, 1, quadOrder>
{
	typedef FV1Geometry<TElem, dim> FVGeom;
};

// for order 1 and hanging: use HFV1Geom
template <int dim, typename TElem, int quadOrder>
struct ChooseProperFVGeom<dim, TElem, true, 1, quadOrder>
{
	typedef HFV1Geometry<TElem, dim> FVGeom;
};

// take FVGeom instead if dim > 1 and 1 < order < 4 and quadOrder == order+1
template <typename TElem, bool hanging>
struct ChooseProperFVGeom<2, TElem, hanging, 2>
{
	typedef FVGeometry<2, TElem, 2> FVGeom;
};
template <typename TElem, bool hanging>
struct ChooseProperFVGeom<2, TElem, hanging, 3>
{
	typedef FVGeometry<3, TElem, 2> FVGeom;
};
template <typename TElem, bool hanging>
struct ChooseProperFVGeom<3, TElem, hanging, 2>
{
	typedef FVGeometry<2, TElem, 3> FVGeom;
};
template <typename TElem, bool hanging>
struct ChooseProperFVGeom<3, TElem, hanging, 3>
{
	typedef FVGeometry<3, TElem, 3> FVGeom;
};

// take DimFVGeometry again instead if dim == 3 and TElem == Pyramid, Octahedron
template <bool hanging>
struct ChooseProperFVGeom<3, Pyramid, hanging, 2>
{
	typedef DimFVGeometry<3> FVGeom;
};
template <bool hanging>
struct ChooseProperFVGeom<3, Pyramid, hanging, 3>
{
	typedef DimFVGeometry<3> FVGeom;
};
template <bool hanging>
struct ChooseProperFVGeom<3, Octahedron, hanging, 2>
{
	typedef DimFVGeometry<3> FVGeom;
};
template <bool hanging>
struct ChooseProperFVGeom<3, Octahedron, hanging, 3>
{
	typedef DimFVGeometry<3> FVGeom;
};




} // namespace nernst_planck
} // namespace ug

#endif // UG__PLUGINS__NERNST_PLANCK__UTIL__CHOOSE_FVGEOM_H
