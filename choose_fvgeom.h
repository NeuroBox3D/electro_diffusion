/*
 * choose_fvgeom.h
 *
 *  Created on: 2016-12-18
 *      Author: mbreit
 */

#ifndef UG__PLUGINS__NERNST_PLANCK__CHOOSE_FVGEOM_H
#define UG__PLUGINS__NERNST_PLANCK__CHOOSE_FVGEOM_H

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

#endif // UG__PLUGINS__NERNST_PLANCK__FLUX_EXPORTER_H
