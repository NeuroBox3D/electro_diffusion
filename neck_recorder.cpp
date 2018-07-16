/*
 * NeckRecorder.cpp
 *
 *  Created on: 2017-11-29
 *      Author: mbreit
 */

#include "neck_recorder.h"

#include <boost/mpl/for_each.hpp>  // for_each

#include "common/util/provider.h"  // Provider
#include "common/util/string_util.h"  // FilenameAndPathWithoutExtension, GetFilenameExtension
#include "lib_algebra/cpu_algebra_types.h" // CPUAlgebra
#include "lib_disc/domain.h"  // Domain3d
#include "lib_disc/domain_util.h"  // CollectCornerCoordinates
#include "lib_disc/local_finite_element/lagrange/lagrangep1.h"  // LagrangeP1
#include "lib_disc/quadrature/quadrature.h"  // QuadratureRule
#include "lib_disc/quadrature/quadrature_provider.h"  // QuadratureRuleProvider
#include "lib_disc/reference_element/reference_mapping_provider.h"  // DimReferenceMapping
#include "lib_disc/spatial_disc/disc_util/conv_shape.h"  // ConvectionShapesNoUpwind
#include "lib_disc/spatial_disc/disc_util/geom_provider.h"  // GeomProvider
#include "lib_disc/spatial_disc/disc_util/fv1_geom.h"  // FV1Geometry
#include "lib_grid/algorithms/debug_util.h"  // ElementDebugInfo
#include "lib_grid/algorithms/element_side_util.h"  // vertexGroupsMatch
#include "lib_grid/algorithms/subset_dim_util.h"  // SubsetIsRegularGrid
#include "lib_grid/grid/grid_base_objects.h"  // Volume
#include "lib_grid/grid_objects/grid_objects_3d.h"  // geometry_traits<Tetrahedron>
#include "lib_grid/tools/subset_group.h"  // SubsetGroup

#include "lib_grid/file_io/file_io.h"  // SaveGridToFile  (only for testing, remove later)
#include "choose_fvgeom.h"  // ChooseProperFVGeom


//#define NECK_RECORDER_DEBUG 1

namespace ug {
namespace nernst_planck {


#if NECK_RECORDER_DEBUG
static Grid testGrid;

template <typename TDomain>
struct TestAttachment
{
	static typename TDomain::position_attachment_type aPosTest;
};

template <> Domain1d::position_attachment_type TestAttachment<Domain1d>::aPosTest = Domain1d::position_attachment_type();
template <> Domain2d::position_attachment_type TestAttachment<Domain2d>::aPosTest = Domain2d::position_attachment_type();
template <> Domain3d::position_attachment_type TestAttachment<Domain3d>::aPosTest = Domain3d::position_attachment_type();
#endif

enum coState
{
	CS_BELOW = 1,
	CS_ABOVE = 1 << 1  // equality is to be understood as above
};

static int cornerState(number pos, number thresh)
{
	if (pos < thresh)
		return CS_BELOW;
	return CS_ABOVE;  // equality is regarded as being above
}


template <typename TDomain>
static bool check_for_intersection
(
	IVertexGroup* elem,
	number thresh,
	typename TDomain::position_accessor_type& aaPos,
	int cornerStates[domain_traits<TDomain::dim>::MaxNumVerticesOfElem] = NULL
)
{
	// intersection exists if two of the volumes' vertices are located on opposite sides
	int check = 0;

	size_t nVrt = elem->num_vertices();
	for (size_t v = 0; v < nVrt; ++v)
	{
		const int coState = cornerState(aaPos[elem->vertex(v)][TDomain::dim-1], thresh);
		check |= coState;
		if (cornerStates)
			cornerStates[v] = coState;
	}

	// no intersection -> continue
	return (check & (CS_BELOW | CS_ABOVE)) == (CS_BELOW | CS_ABOVE);
}



template <typename TDomain>
static int hangingCornerState
(
	Vertex* vrt,
	SmartPtr<TDomain> dom,
	number thresh
)
{
	typename TDomain::position_accessor_type& aaPos = dom->position_accessor();

	// normal vertex
	if (!vrt->is_constrained())
		return cornerState(aaPos[vrt][TDomain::dim-1], thresh);

	// hanging vertex: bitwise-or of constrainer states
	int state = 0;
	typedef typename MultiGrid::traits<Edge>::secure_container edge_list;
	edge_list el;
	dom->grid()->associated_elements(el, vrt);
	const size_t ne = el.size();
	for (size_t e = 0; e < ne; ++e)
	{
		Edge* edge = el[e];
		if (edge->is_constrained())
		{
			Vertex* other = GetOpposingSide(*dom->grid(), edge, vrt);
			if (!other->is_constrained())
				state |= cornerState(aaPos[other][TDomain::dim-1], thresh);
		}
	}

	return state;
}


template <typename TDomain>
static int hangingCornerStateCoarse
(
	const typename TDomain::position_type& coords,
	number thresh,
	typename grid_dim_traits<TDomain::dim>::element_type* elem,
	ConstSmartPtr<TDomain> dom,
	number refLenSq
)
{
	// find constraingin edge
	const typename TDomain::position_accessor_type aaPos = dom->position_accessor();
	size_t ne = elem->num_edges();
	EdgeDescriptor ed;
	size_t e = 0;
	typename TDomain::position_type eCenter;
	for (; e < ne; ++e)
	{
		ed = elem->edge_desc(e);
		VecScaleAdd(eCenter, 0.5, aaPos[ed.vertex(0)], 0.5, aaPos[ed.vertex(1)]);
		if (VecDistanceSq(eCenter, coords) < 1e-6 * refLenSq)
			break;
	}
	UG_COND_THROW(e == ne, "No edge center is located at the hanging node's position.");

	// now investigate its vertices
	int state = cornerState(aaPos[ed.vertex(0)][TDomain::dim-1], thresh);
	state |= cornerState(aaPos[ed.vertex(1)][TDomain::dim-1], thresh);

	return state;
}

#ifdef UG_DIM_1
template <>
int hangingCornerStateCoarse<Domain1d>
(
	const Domain1d::position_type& coords,
	number thresh,
	grid_dim_traits<Domain1d::dim>::element_type* elem,
	ConstSmartPtr<Domain1d> dom,
	number refLenSq
)
{
	UG_THROW("There can be no 1d hanging nodes.");
}
#endif

/// Checks whether a volume is in the extended intersection, i.e.,
/// if it is intersected itself or of it contains a hanging node
/// which is constrained by an intersected edge.
template <typename TDomain>
static bool check_for_hanging_intersection
(
	typename domain_traits<TDomain::dim>::element_type* vol,
	SmartPtr<TDomain> dom,
	number thresh,
	int cytSI,
	int cornerStates[domain_traits<TDomain::dim>::MaxNumVerticesOfElem]
)
{
	int check = 0;

	const size_t nVrt = vol->num_vertices();
	for (size_t v = 0; v < nVrt; ++v)
	{
		Vertex* vrt = vol->vertex(v);

		// check for normal intersection
		check |= cornerStates[v] = hangingCornerState(vrt, dom, thresh);
	}

	return check == (CS_BELOW | CS_ABOVE);
}



template <typename TElem>
static void elem_edge_connections(std::vector<std::vector<size_t> >& vConn)
{
	UG_THROW("This unspecialized implementation is never to be used.");
}

template <>
void elem_edge_connections<Tetrahedron>(std::vector<std::vector<size_t> >& vConn)
{
	vConn[0].resize(3);
	vConn[0][0] = 1;
	vConn[0][1] = 2;
	vConn[0][2] = 3;

	vConn[1].resize(2);
	vConn[1][0] = 2;
	vConn[1][1] = 3;

	vConn[2].resize(1);
	vConn[2][0] = 3;
}

template <>
void elem_edge_connections<Prism>(std::vector<std::vector<size_t> >& vConn)
{
	vConn[0].resize(3);
	vConn[0][0] = 1;
	vConn[0][1] = 2;
	vConn[0][2] = 3;

	vConn[1].resize(2);
	vConn[1][0] = 2;
	vConn[1][1] = 4;

	vConn[2].resize(1);
	vConn[2][0] = 5;

	vConn[3].resize(2);
	vConn[3][0] = 4;
	vConn[3][1] = 5;

	vConn[4].resize(1);
	vConn[4][0] = 5;
}

template <>
void elem_edge_connections<Triangle>(std::vector<std::vector<size_t> >& vConn)
{
	vConn[0].resize(1);
	vConn[0][0] = 1;

	vConn[1].resize(1);
	vConn[1][0] = 2;

	vConn[2].resize(1);
	vConn[2][0] = 0;
}

template <>
void elem_edge_connections<Quadrilateral>(std::vector<std::vector<size_t> >& vConn)
{
	vConn[0].resize(1);
	vConn[0][0] = 1;

	vConn[1].resize(1);
	vConn[1][0] = 2;

	vConn[2].resize(1);
	vConn[2][0] = 3;

	vConn[3].resize(1);
	vConn[3][0] = 0;
}


template <typename TDomain>
static void calculate_intersection_corners
(
	std::vector<typename TDomain::position_type>& vIsecCorners,
	typename domain_traits<TDomain::dim>::element_type* vol,
	number thresh,
	int cornerState[domain_traits<TDomain::dim>::MaxNumVerticesOfElem],
	typename TDomain::position_accessor_type& aaPos,
	MultiGrid& mg
)
{
	const int dim = TDomain::dim;

	size_t nVrt = vol->num_vertices();
	std::vector<std::vector<size_t> > vConn(nVrt);

	// determine exact element type and get edge connections
	if (dynamic_cast<Tetrahedron*>(vol))
		elem_edge_connections<Tetrahedron>(vConn);
	else if (dynamic_cast<Prism*>(vol))
		elem_edge_connections<Prism>(vConn);
	else if (dynamic_cast<Triangle*>(vol))
		elem_edge_connections<Triangle>(vConn);
	else if (dynamic_cast<Quadrilateral*>(vol))
		elem_edge_connections<Quadrilateral>(vConn);
	else
	{
		UG_THROW("Volume is neither tetrahedron or prism nor triangle or quadrilateral, "
			"but other volume types are not supported by this method.");
	}

	CustomVertexGroup vrtGrp;
	for (size_t v1 = 0; v1 < nVrt; ++v1)
	{
		typename TDomain::position_type& p1 = aaPos[vol->vertex(v1)];

		// a volume corner is also intersection corner
		if (cornerState[v1] == 3)
		{
			vIsecCorners.push_back(p1);
			vrtGrp.push_back(vol->vertex(v1));
			continue;
		}

		// intersection corner on volume edge
		std::vector<size_t>& vConnVrt = vConn[v1];
		size_t nConn = vConnVrt.size();
		for (size_t i = 0; i < nConn; ++i)
		{
			size_t v2 = vConnVrt[i];
			if (cornerState[v1] + cornerState[v2] == 3)
			{
				typename TDomain::position_type& p2 = aaPos[vol->vertex(v2)];
				number lambda = (p2[dim-1] - thresh) / (p2[dim-1] - p1[dim-1]);
				vIsecCorners.resize(vIsecCorners.size()+1);
				typename TDomain::position_type& isecPos = vIsecCorners.back();
				VecScaleAdd(isecPos, lambda, p1, 1.0-lambda, p2);
			}
		}
	}

#ifdef UG_PARALLEL
	if (pcl::NumProcs() <= 1)
		return;

	// check whether intersection is a side
	// and if it is, check whether it is a vmaster (in that case, drop it)
	if (!vrtGrp.size()) // most of the cases
		return;

	const DistributedGridManager& dgm = *mg.distributed_grid_manager();
	typedef typename domain_traits<TDomain::dim>::side_type side_type;
	typedef typename MultiGrid::traits<side_type>::secure_container side_list;
	side_list sl;
	mg.associated_elements(sl, vol);
	size_t slsz = sl.size();
	for (size_t s = 0; s < slsz; ++s)
	{
		side_type* side = sl[s];
		if (VertexGroupsMatch(side, vrtGrp))
		{
			// check for VMASTER
			if (dgm.contains_status(side, INT_V_MASTER))
				vIsecCorners.clear();

			return;
		}
	}
#endif
}


template <typename TDomain>
static bool postprocess_intersection_polygon
(
	ReferenceObjectID& roid,
	std::vector<typename TDomain::position_type>& vIsecCorners
)
{
	size_t nCo = vIsecCorners.size();
	if (TDomain::dim == 3)
	{
		// lower-dim manifold elements can be ignored (zero measure)
		if (nCo <= 2)
			return false;

		if (nCo == 3)
			return true;

		// four- and five-corner polygons are divided into two and three triangles, resp.
		// to that end, find the correct order to make polygon convex (we suppose that this is possible)
		if (nCo == 4 || nCo == 5)
		{
			// sort using selection sort on x
			for (size_t i = 0; i < nCo; ++i)
			{
				number min = vIsecCorners[i][0];
				number minInd = i;
				for (size_t j = i+1; j < nCo; ++j)
				{
					if (vIsecCorners[j][0] < min)
					{
						min = vIsecCorners[j][0];
						minInd = j;
					}
				}
				std::swap(vIsecCorners[i], vIsecCorners[minInd]);
			}

			// now use x-order to sort ccw
			std::vector<typename TDomain::position_type> tmp(vIsecCorners.size());
			tmp[0] = vIsecCorners[0];

			typename TDomain::position_type e0;
			VecSubtract(e0, vIsecCorners[nCo-1], vIsecCorners[0]);

			size_t a = 1; // current next from start
			size_t e = nCo-1; // current next from end
			for (size_t i = 1; i < nCo-1; ++i)
			{
				typename TDomain::position_type e1;
				VecSubtract(e1, vIsecCorners[i], vIsecCorners[0]);
				bool upperSide = e0[0]*e1[1] - e0[1]*e1[0] >= 0;
				if (upperSide)
				{
					tmp[e] = vIsecCorners[i];
					--e;
				}
				else
				{
					tmp[a] = vIsecCorners[i];
					++a;
				}
			}
			tmp[a] = vIsecCorners[nCo-1];
			UG_COND_THROW(a != e, "Something went wrong here.");

			vIsecCorners.swap(tmp);

			return true;
		}

		// this location should not be reached
		UG_THROW("Number of corners of intersection polygon is " << nCo
			<< " (>5), this should be impossible with convex elements and indicates "
				"the presence of non-convex prisms in the geometry.");
	}

	if (TDomain::dim == 2)
	{
		if (nCo <= 1)
			return false;

		if (nCo == 2)
			return true;

		// only edges should occur in 2d
		UG_COND_THROW(nCo > 2, "Number of corners of intersection polygon is " << nCo
			<< ", this should be impossible with convex elements and indicates "
				"the presence of non-convex quads in the geometry.");
	}

	UG_THROW("Only dimensions 2 and 3 are supported.");
}


template <typename TDomain, typename TElem>
static void quadrature_weights_and_points
(
	NeckRecorderBase::IntegrationVolume& iv,
	std::vector<MathVector<TDomain::dim> >& vQuadPts,
	const std::vector<typename TDomain::position_type>& vIsecCorners
)
{
	const int dim = TDomain::dim;
	ReferenceObjectID roid = geometry_traits<TElem>::REFERENCE_OBJECT_ID;

	try
	{
		const size_t order = 3;
		const QuadratureRule<dim-1>& rQuadRule
			= QuadratureRuleProvider<dim-1>::get(roid, order, GetQuadratureType("best")); // "gauss"?

		// get reference element mapping by reference object id
		DimReferenceMapping<dim-1, dim>& mapping
			= ReferenceMappingProvider::get<dim-1, dim>(roid);

		// number of integration points
		const size_t numIP = rQuadRule.size();

		// update the reference mapping for the corners
		mapping.update(&vIsecCorners[0]);

		// compute global integration points
		size_t qpsz = vQuadPts.size();
		vQuadPts.resize(qpsz + numIP);
		mapping.local_to_global(&(vQuadPts[qpsz]), rQuadRule.points(), numIP);

		// compute transformation matrices
		std::vector<MathMatrix<dim-1, dim> > vJT(numIP);
		mapping.jacobian_transposed(&(vJT[0]), rQuadRule.points(), numIP);

		// save all computed data in intermediate structure
		std::vector<NeckRecorderBase::IPData>& vIPData = iv.vIPData;
		size_t ipdsz = vIPData.size();
		vIPData.resize(ipdsz + numIP);
		for (size_t ip = 0; ip < numIP; ++ip)
		{
			NeckRecorderBase::IPData& ipd = vIPData[ip+ipdsz];
			ipd.weight = rQuadRule.weight(ip);
			ipd.detJ = SqrtGramDeterminant(vJT[ip]);

			// it can happen that the triangle is almost an edge, then the det is practically 0,
			// but due to minuscule computation errors, the determinant of J^TJ is negative
			// (which is mathematically impossible) and thus we assign a NaN value here
			if (ipd.detJ != ipd.detJ)
				ipd.detJ = 0.0;
		}
	}
	UG_CATCH_THROW("Error during quadrature or reference mapping.");
}


template <typename TDomain>
static void quadrature_weights_and_points
(
	NeckRecorderBase::IntegrationVolume& iv,
	std::vector<MathVector<TDomain::dim> >& vQuadPts,
	const std::vector<typename TDomain::position_type>& vIsecCorners
)
{
	UG_THROW("This unspecialized implementation is never to be used.");
}


#ifdef UG_DIM_3
template <>
void quadrature_weights_and_points<Domain3d>
(
	NeckRecorderBase::IntegrationVolume& iv,
	std::vector<MathVector<Domain3d::dim> >& vQuadPts,
	const std::vector<Domain3d::position_type>& vIsecCorners
)
{
	size_t nCo = vIsecCorners.size();

	if (nCo == 3)
		quadrature_weights_and_points<Domain3d, Triangle>(iv, vQuadPts, vIsecCorners);
	else
	{
		// split polygon into triangles
		std::vector<Domain3d::position_type> threeCorners(3);
		threeCorners[0] = vIsecCorners[0];
		for (size_t i = 2; i < nCo; ++i)
		{
			threeCorners[1] = vIsecCorners[i-1];
			threeCorners[2] = vIsecCorners[i];
			quadrature_weights_and_points<Domain3d, Triangle>(iv, vQuadPts, threeCorners);
		}
	}
}
#endif

#ifdef UG_DIM_2
template <>
void quadrature_weights_and_points<Domain2d>
(
	NeckRecorderBase::IntegrationVolume& iv,
	std::vector<MathVector<Domain2d::dim> >& vQuadPts,
	const std::vector<Domain2d::position_type>& vIsecCorners
)
{
	quadrature_weights_and_points<Domain2d, Edge>(iv, vQuadPts, vIsecCorners);
}
#endif


template <typename TDomain, typename TElem, bool hanging>
static void integration_points_of_approximating_scvfs_elem
(
	NeckRecorderBase::IntegrationVolume& iv,
	TElem* elem,
	number thresh,
	int cornerState[domain_traits<TDomain::dim>::MaxNumVerticesOfElem],
	ConstSmartPtr<TDomain> dom
)
{
#if NECK_RECORDER_DEBUG
	typename TDomain::position_accessor_type aaPosTest(testGrid, TestAttachment<TDomain>::aPosTest);
#endif
	typedef typename ChooseProperFVGeom<TDomain::dim, TElem, hanging, 1>::FVGeom fvgeom_type;

	// get FV1 geometry for element
	std::vector<typename TDomain::position_type> vCornerCoords;
	CollectCornerCoordinates(vCornerCoords, elem, dom->position_accessor(), false);

	static fvgeom_type& geo = GeomProvider<fvgeom_type>::get();
	try {geo.update(elem, &vCornerCoords[0], dom->subset_handler().get());}
	UG_CATCH_THROW("Failed creating FV1Geometry for element.");

	// find the edges that are cut by the integration manifold
	// if the manifold happens to go exactly through a corner, take all adjacent edges,
	// but only with half the integration weight (0.5 for FV1) at the corresponding SCVF IP

	// loop SCVFs
	size_t nscvf = geo.num_scvf();
	for (size_t i = 0; i < nscvf; ++i)
	{
		const typename fvgeom_type::SCVF& scvf = geo.scvf(i);
		size_t from = scvf.from();
		size_t to = scvf.to();

		const bool fromIsHangingInCoarse = from >= elem->num_vertices();
		const bool toIsHangingInCoarse = to >= elem->num_vertices();

		const typename TDomain::position_type& coCoFrom = geo.global_node_position(from);
		const typename TDomain::position_type& coCoTo = geo.global_node_position(to);
		const number refLenSq = VecDistanceSq(coCoFrom, coCoTo);

		int coStFrom = 0;
		int coStTo = 0;

		if (!fromIsHangingInCoarse)
			// corner state is already known
			coStFrom = cornerState[from];
		else
			// corner state must be composed from constrainers
			coStFrom = hangingCornerStateCoarse(coCoFrom, thresh, elem, dom, refLenSq);

		if (!toIsHangingInCoarse)
			// corner state is already known
			coStTo = cornerState[to];
		else
			// corner state must be composed from constrainers
			coStTo = hangingCornerStateCoarse(coCoTo, thresh, elem, dom, refLenSq);


		number weight = 0.0;
		if (coStFrom == CS_BELOW && coStTo == CS_ABOVE)
			weight = 1.0;
		else if (coStFrom == CS_ABOVE && coStTo == CS_BELOW)
			weight = -1.0;
		else if ((coStFrom == CS_BELOW && coStTo == (CS_BELOW | CS_ABOVE))
			|| (coStFrom == (CS_BELOW | CS_ABOVE) && coStTo == CS_ABOVE))
			weight = 0.5;
		else if ((coStFrom == (CS_BELOW | CS_ABOVE) && coStTo == CS_BELOW)
			|| (coStFrom == CS_ABOVE && coStTo == (CS_BELOW | CS_ABOVE)))
			weight = -0.5;

		if (weight == 0.0)
			continue;


#if 0
		// both normal
		if (coStFrom == CS_BELOW && coStTo == CS_ABOVE)
		{
			number weight = 1.0;
		}
		else if (coStFrom == CS_ABOVE && coStTo == CS_BELOW)
		{
			number weight = -1.0;
		}


		// from normal, to hanging
		if (coStFrom == CS_BELOW && coStTo == (CS_BELOW | CS_ABOVE | CS_HANGING))
		{
			number weight = 0.5;
		}
		else if (coStFrom == CS_ABOVE && coStTo == (CS_BELOW | CS_ABOVE | CS_HANGING))
		{
			number weight = -0.5;
		}
		else if (coStFrom == CS_BELOW && coStTo == CS_ABOVE)
		{
			number weight = 1.0;
		}
		else if (coStFrom == CS_ABOVE && coStTo == CS_BELOW)
		{
			number weight = -1.0;
		}


		// from hanging, to normal
		if (coStFrom == (CS_BELOW | CS_ABOVE | CS_HANGING) && coStTo == CS_ABOVE)
		{
			number weight = 0.5;
		}
		else if (coStFrom == (CS_BELOW | CS_ABOVE | CS_HANGING) && coStTo == CS_BELOW)
		{
			number weight = -0.5;
		}
		else if (coStFrom == CS_BELOW && coStTo == CS_ABOVE)
		{
			number weight = 1.0;
		}
		else if (coStFrom == CS_ABOVE && coStTo == CS_BELOW)
		{
			number weight = -1.0;
		}


		// both hanging
		if (coStFrom == CS_BELOW && coStTo == (CS_BELOW | CS_ABOVE | CS_HANGING))
		{
			number weight = 0.5;
		}
		else if (coStFrom == CS_BELOW && coStTo == CS_ABOVE)
		{
			number weight = 1.0;
		}
		else if (coStFrom == (CS_BELOW | CS_ABOVE | CS_HANGING) && coStTo == CS_BELOW)
		{
			number weight = -0.5;
		}
		else if (coStFrom == (CS_BELOW | CS_ABOVE | CS_HANGING) && coStTo == CS_ABOVE)
		{
			number weight = 0.5;
		}
		else if (coStFrom == CS_ABOVE && coStTo == CS_BELOW)
		{
			number weight = -1.0;
		}
		else if (coStFrom == CS_ABOVE && coStTo == (CS_BELOW | CS_ABOVE | CS_HANGING))
		{
			number weight = -0.5;
		}
#endif


#if NECK_RECORDER_DEBUG
if (scvf.num_corners() == 4)
{
	Vertex* cornerVrts[4];
	for (size_t j = 0; j < 4; ++j)
	{
		cornerVrts[j] = *testGrid.create<RegularVertex>();
		aaPosTest[cornerVrts[j]] = scvf.global_corner(j);
	}
	testGrid.create<Quadrilateral>(QuadrilateralDescriptor(cornerVrts[0], cornerVrts[1], cornerVrts[2], cornerVrts[3]));
}
else if (scvf.num_corners() == 3)
{
	Vertex* cornerVrts[3];
	for (size_t j = 0; j < 3; ++j)
	{
		cornerVrts[j] = *testGrid.create<RegularVertex>();
		aaPosTest[cornerVrts[j]] = scvf.global_corner(j);
	}
	testGrid.create<Triangle>(TriangleDescriptor(cornerVrts[0], cornerVrts[1], cornerVrts[2]));
}
else if (scvf.num_corners() == 2)
{
	Vertex* cornerVrts[2];
	cornerVrts[0] = *testGrid.create<RegularVertex>();
	cornerVrts[1] = *testGrid.create<RegularVertex>();
	aaPosTest[cornerVrts[0]] = scvf.global_corner(0);
	aaPosTest[cornerVrts[1]] = scvf.global_corner(1);

	testGrid.create<RegularEdge>(EdgeDescriptor(cornerVrts[0], cornerVrts[1]));
}
else UG_THROW("Unexpected number of corners in SCVF: " << scvf.num_corners() << " instead of 2 or 3 or 4"
		" in " << ElementDebugInfo(*dom->grid(), elem) << ".");
#endif

		// get integration point on scvf (for generalization to higher-order schemes, this would have to be adapted)
		iv.vIPData.resize(iv.vIPData.size() + 1);
		NeckRecorderBase::IPData& ipd = iv.vIPData.back();

		ipd.weight = weight;
		ipd.detJ = 1.0; // more correctly, this would be the area of the scvf, but it's already contained in scvf.normal()

		// orientation of normal must be from below manifold to above it; invert it if this is not the case
		typename TDomain::position_type normal = scvf.normal();

		size_t nSh = scvf.num_sh();
		ipd.vShapes.resize(nSh);
		ipd.vGradZ.resize(nSh);
		for (size_t sh = 0; sh < nSh; ++sh)
		{
			ipd.vShapes[sh] = scvf.shape(sh);
			ipd.vGradZ[sh] = VecDot(scvf.global_grad(sh), normal);
		}
		MathVector<TDomain::dim> pos(0.0);
		VecScaleAdd(pos, 0.5, scvf.global_corner(0), 0.5, scvf.global_corner(1));
		ipd.pos[0] = pos[0];
		ipd.pos[1] = pos[1];
	}
}


template <typename TDomain, bool hanging>
struct wrap_ipoas
{
	wrap_ipoas
	(
		NeckRecorderBase::IntegrationVolume& _iv,
		typename domain_traits<TDomain::dim>::element_type* _elem,
		number _thresh,
		int _cornerState[domain_traits<TDomain::dim>::MaxNumVerticesOfElem],
		ConstSmartPtr<TDomain> _dom
	)
	: iv(_iv), elem(_elem), thresh(_thresh), cornerState(_cornerState), dom(_dom) {}

	template <typename TElem>
	void operator() (TElem&)
	{
		TElem* castElem = dynamic_cast<TElem*>(elem);
		if (castElem)
			integration_points_of_approximating_scvfs_elem<TDomain, TElem, hanging>
				(iv, castElem, thresh, cornerState, dom);
	}

	NeckRecorderBase::IntegrationVolume& iv;
	typename domain_traits<TDomain::dim>::element_type* elem;
	number thresh;
	int* cornerState;
	ConstSmartPtr<TDomain> dom;
};


template <typename TDomain, bool hanging>
static void integration_points_of_approximating_scvfs
(
	NeckRecorderBase::IntegrationVolume& iv,
	typename domain_traits<TDomain::dim>::element_type* elem,
	number thresh,
	int cornerState[domain_traits<TDomain::dim>::MaxNumVerticesOfElem],
	ConstSmartPtr<TDomain> dom
)
{
	typedef typename domain_traits<TDomain::dim>::DimElemList ElemList;
	boost::mpl::for_each<ElemList>(wrap_ipoas<TDomain, hanging>(iv, elem, thresh, cornerState, dom));
}



template <typename TDomain, typename TElem>
void eval_shapes_and_grads
(
	NeckRecorderBase::IntegrationVolume& iv,
	const std::vector<MathVector<TDomain::dim> >& vPts,
	typename domain_traits<TDomain::dim>::element_type* vol,
	typename TDomain::position_accessor_type& aaPos
)
{
	const int dim = TDomain::dim;

	// collect volume corner coords
	size_t nCo = vol->num_vertices();
	std::vector<typename TDomain::position_type> coco(nCo);
	for (size_t co = 0; co < nCo; ++co)
		coco[co] = aaPos[vol->vertex(co)];

	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
	typedef LagrangeP1<ref_elem_type> local_shape_fct_set_type;
	local_shape_fct_set_type& lsfs = Provider<LagrangeP1<ref_elem_type> >::get();

	// prepare global-to-local mapping
	DimReferenceMapping<dim, dim>& mapping = ReferenceMappingProvider::get<dim, dim>
		(geometry_traits<TElem>::REFERENCE_OBJECT_ID, coco);

	MathMatrix<dim, dim> jacInvT;

	size_t nSh = lsfs.num_sh();
	size_t numIP = vPts.size();
	for (size_t ip = 0; ip < numIP; ++ip)
	{
		NeckRecorderBase::IPData& ipd = iv.vIPData[ip];

		// get local coords of IP and J^-T
		MathVector<dim> locPt;
		try
		{
			mapping.global_to_local(locPt, vPts[ip]);
			mapping.jacobian_transposed_inverse(jacInvT, locPt);
		}
		UG_CATCH_THROW("Global-local mapping not successful.");

		// get shapes and gradients in IP
		ipd.vShapes.resize(nSh);
		ipd.vGradZ.resize(nSh);
		for (size_t sh = 0; sh < nSh; ++sh)
		{
			ipd.vShapes[sh] = lsfs.shape(sh, locPt);
			MathVector<dim> locGrad, globGrad;
			lsfs.grad(locGrad, sh, locPt);
			MatVecMult(globGrad, jacInvT, locGrad);
			ipd.vGradZ[sh] = globGrad[dim-1];
		}
	}
}

template <typename TDomain>
static void eval_shapes_and_grads
(
	NeckRecorderBase::IntegrationVolume& iv,
	const std::vector<MathVector<TDomain::dim> >& vPts,
	typename domain_traits<TDomain::dim>::element_type* vol,
	typename TDomain::position_accessor_type& aaPos
)
{
	UG_THROW("This unspecialized implementation is never to be used.");
}

#ifdef UG_DIM_3
template <>
void eval_shapes_and_grads<Domain3d>
(
	NeckRecorderBase::IntegrationVolume& iv,
	const std::vector<MathVector<3> >& vPts,
	domain_traits<Domain3d::dim>::element_type* vol,
	Domain3d::position_accessor_type& aaPos
)
{
	if (dynamic_cast<Tetrahedron*>(vol))
		eval_shapes_and_grads<Domain3d, Tetrahedron>(iv, vPts, vol, aaPos);
	else if (dynamic_cast<Prism*>(vol))
		eval_shapes_and_grads<Domain3d, Prism>(iv, vPts, vol, aaPos);
	else
	{
		UG_THROW("Volume is neither tetrahedron nor prism, "
			"but other volume types are not supported by this method.");
	}
}
#endif

#ifdef UG_DIM_2
template <>
void eval_shapes_and_grads<Domain2d>
(
	NeckRecorderBase::IntegrationVolume& iv,
	const std::vector<MathVector<2> >& vPts,
	domain_traits<Domain2d::dim>::element_type* vol,
	Domain2d::position_accessor_type& aaPos
)
{
	if (dynamic_cast<Triangle*>(vol))
		eval_shapes_and_grads<Domain2d, Triangle>(iv, vPts, vol, aaPos);
	else if (dynamic_cast<Quadrilateral*>(vol))
		eval_shapes_and_grads<Domain2d, Quadrilateral>(iv, vPts, vol, aaPos);
	else
	{
		UG_THROW("Face is neither triangle nor quadrilateral, "
			"but other face types are not supported by this method.");
	}
}
#endif


template <typename TDomain>
static void dof_indices_on_element
(
	NeckRecorderBase::IntegrationVolume& iv,
	typename domain_traits<TDomain::dim>::element_type* elem,
	ConstSmartPtr<DoFDistribution> dd
)
{
	// order of unknowns in approx space
	size_t fctMap[5];
	fctMap[NeckRecorderBase::_PHI_] = dd->fct_id_by_name("Phi");
	fctMap[NeckRecorderBase::_K_] = dd->fct_id_by_name("K");
	fctMap[NeckRecorderBase::_NA_] = dd->fct_id_by_name("Na");
	fctMap[NeckRecorderBase::_CL_] = dd->fct_id_by_name("Cl");
	fctMap[NeckRecorderBase::_A_] = dd->fct_id_by_name("A");

	// get dof indices for functions on element
	std::vector<std::vector<DoFIndex> >& vvDofIndices = iv.vDofIndex;
	size_t nVrt = elem->num_vertices();
	vvDofIndices.resize(nVrt);
	std::vector<DoFIndex> vDI;
	for (size_t co = 0; co < nVrt; ++co)
	{
		vvDofIndices[co].resize(5);
		Vertex* vrt = elem->vertex(co);
		for (size_t i = 0; i < 5; ++i)
		{
			size_t fct = fctMap[i];
			dd->inner_dof_indices(vrt, fct, vDI, true);
			UG_COND_THROW(vDI.size() != 1, "Not exactly one DoF of fct " << fct
				<< " located on " << ElementDebugInfo(*dd->multi_grid(), vrt));
			vvDofIndices[co][i] = vDI[0];
		}
	}
}



NeckRecorderBase::NeckRecorderBase()
: F(96485.0), R(8.31451)
{};


template <typename TDomain, typename TAlgebra>
NeckRecorder<TDomain, TAlgebra>::NeckRecorder(SmartPtr<ApproximationSpace<TDomain> > approx)
: m_spApprox(approx),
  m_spConvShape(new ConvectionShapesNoUpwind<dim>),
  m_siCyt(-1),
  m_vDiffConst(4),
  m_vConvConst(4),
  m_temp(298.15),
  m_bIndividualCurrents(false),
  bPreparedCurrent(false),
  bPreparedPot(false),
  m_bSsIsRegular(false)
{
	m_vDiffConst[_K_] = 1.96e-9;
	m_vDiffConst[_NA_] = 1.33e-9;
	m_vDiffConst[_CL_] = 2.03e-9;
	m_vDiffConst[_A_] = 2.0e-9;

	m_vConvConst[_K_] = -m_vDiffConst[_K_] * F/(R*m_temp);
	m_vConvConst[_NA_] = -m_vDiffConst[_NA_] * F/(R*m_temp);
	m_vConvConst[_CL_] = m_vDiffConst[_CL_] * F/(R*m_temp);
	m_vConvConst[_A_] = m_vDiffConst[_A_] * F/(R*m_temp);
}


template <typename TDomain, typename TAlgebra>
void NeckRecorder<TDomain, TAlgebra>::add_measurement_zone(number zCoord, const std::string& name)
{
	m_vMeasCoords.push_back(zCoord);
	m_vMeasZoneNames.push_back(name);
	m_vIntegrationDataCurrent.resize(m_vIntegrationDataCurrent.size()+1);
	m_vIntegrationDataPot.resize(m_vIntegrationDataPot.size()+1);
	bPreparedCurrent = false;
	bPreparedPot = false;
}


template <typename TDomain, typename TAlgebra>
void NeckRecorder<TDomain, TAlgebra>::set_cytosolic_subset(const std::string& ss)
{
	SubsetGroup ssg(m_spApprox->domain()->subset_handler());
	try {ssg.add(ss);}
	UG_CATCH_THROW("Cannot use subset named '" << ss << "' as cytosolic subset.");

	UG_COND_THROW(ssg.size() < 1, "No cytosolic subset name contained in '" << ss << "'.");
	UG_COND_THROW(ssg.size() > 1, "Please only provide one subset as cytosolic subset.");

	m_siCyt = ssg[0];
	bPreparedCurrent = false;
	bPreparedPot = false;
}


template <typename TDomain, typename TAlgebra>
void NeckRecorder<TDomain, TAlgebra>::
set_upwind(SmartPtr<IConvectionShapes<TDomain::dim> > spUpwind)
{
	m_spConvShape = spUpwind;
}


template <typename TDomain, typename TAlgebra>
void NeckRecorder<TDomain, TAlgebra>::set_diffusion_constants(const std::vector<number>& diffConsts)
{
	m_vDiffConst = diffConsts;
}


template <typename TDomain, typename TAlgebra>
void NeckRecorder<TDomain, TAlgebra>::set_convection_constants(const std::vector<number>& convConsts)
{
	m_vConvConst = convConsts;
}


template <typename TDomain, typename TAlgebra>
void NeckRecorder<TDomain, TAlgebra>::set_temperature(number t)
{
	UG_COND_THROW(t <= 0, "Temperature must be positive.");
	m_temp = t;
}


template <typename TDomain, typename TAlgebra>
void NeckRecorder<TDomain, TAlgebra>::set_record_individual_currents(bool indivCurr)
{
	m_bIndividualCurrents = indivCurr;
}


template <typename TDomain, typename TAlgebra>
void NeckRecorder<TDomain, TAlgebra>::prepare
(
	std::vector<std::vector<IntegrationVolume> >& integData,
	bool useSCVFMode
)
{
	typedef typename domain_traits<dim>::element_type elem_type;

	// surface dof distro
	SmartPtr<DoFDistribution> dd = m_spApprox->dof_distribution(GridLevel());

	// check that cytosolic subset is set
	UG_COND_THROW(m_siCyt == -1, "Cytosolic subset name not specified. "
		"Do so using 'set_cytosolic_subset()'.");

	// access to position attachment
	typename TDomain::position_accessor_type& aaPos = m_spApprox->domain()->position_accessor();

	// loop all volumes in cytosolic subset
	typedef typename DoFDistribution::traits<elem_type>::const_iterator vol_iter_type;
	vol_iter_type it = dd->template begin<elem_type>(m_siCyt);
	vol_iter_type it_end = dd->template end<elem_type>(m_siCyt);
	for (; it != it_end; ++it)
	{
		elem_type* vol = *it;

		// loop measurement zones
		size_t nMZ = m_vMeasCoords.size();
		for (size_t mz = 0; mz < nMZ; ++mz)
		{
			number thresh = m_vMeasCoords[mz];

			// check for intersection with measurement zone manifold;
			int cornerState[domain_traits<dim>::MaxNumVerticesOfElem]; // can be re-used later on
			bool elemIntersected = false;
			if (useSCVFMode)
			{
				elemIntersected = check_for_hanging_intersection<TDomain>(
					vol, m_spApprox->domain(), thresh, m_siCyt, cornerState);
			}
			else
				elemIntersected = check_for_intersection<TDomain>(vol, thresh, aaPos, cornerState);

			if (!elemIntersected)
				continue;


			// prepare intermediate structure
			// that allows fast computation of flux once grid function is given
			integData[mz].resize(integData[mz].size() + 1);
			IntegrationVolume& iv = integData[mz].back();
			iv.elem = vol;

			if (!useSCVFMode)
			{
				// determine corners of intersection polygon
				std::vector<typename TDomain::position_type> vIsecCorners;
				try {calculate_intersection_corners<TDomain>(vIsecCorners, vol, thresh, cornerState,
					aaPos, *m_spApprox->domain()->grid());}
				UG_CATCH_THROW("Exception during calculation of intersection corners.");

				// post-processing (notably, quads are checked for correct corner order)
				ReferenceObjectID roid;
				bool goOn;
				try {goOn = postprocess_intersection_polygon<TDomain>(roid, vIsecCorners);}
				UG_CATCH_THROW("Exception during post-processing of intersection polygon.");

				// only process intersection integral if dimensions are correct (otherwise zero measure)
				if (!goOn)
				{
					integData[mz].resize(integData[mz].size() - 1); // needs to be removed
					continue;
				}

				// determine quadrature weights and points in element
				std::vector<MathVector<dim> > vQuadPts;
				try {quadrature_weights_and_points<TDomain>(iv, vQuadPts, vIsecCorners);}
				UG_CATCH_THROW("Exception during determination of quadrature weights and points.");

				// evaluate shape functions in global IPs
				try {eval_shapes_and_grads<TDomain>(iv, vQuadPts, vol, aaPos);}
				UG_CATCH_THROW("Exception during evaluation of shape functions.");
			}
			else
			{
				// determine quadrature weights and points in element
				m_bSsIsRegular = SubsetIsRegularGrid(*m_spApprox->domain()->subset_handler(), m_siCyt);
				try
				{
					if (m_bSsIsRegular)
						integration_points_of_approximating_scvfs<TDomain, false>(iv, vol, thresh, cornerState, m_spApprox->domain());
					else
						integration_points_of_approximating_scvfs<TDomain, true>(iv, vol, thresh, cornerState, m_spApprox->domain());
				}
				UG_CATCH_THROW("Exception during determination of integration points on approximating SCVFs.");
			}

			// get dof indices for functions on element
			try {dof_indices_on_element<TDomain>(iv, vol, dd);}
			UG_CATCH_THROW("Exception during acquisition of DoF indices.");
		}
	}
}



template <typename TDomain, typename TAlgebra>
template <bool hanging>
NeckRecorder<TDomain, TAlgebra>::wrap_ccs<hanging>::wrap_ccs
(
	std::vector<std::vector<number> >* _convShapes,
	const NeckRecorderBase::IntegrationVolume& _iv,
	number _thresh,
	ConstSmartPtr<GridFunction<TDomain, TAlgebra> > _u,
	NeckRecorder<TDomain, TAlgebra>* _nr
)
: convShapes(_convShapes), iv(_iv), thresh(_thresh), u(_u), nr(_nr)
{}

template <typename TDomain, typename TAlgebra>
template <bool hanging>
template <typename TElem>
void NeckRecorder<TDomain, TAlgebra>::wrap_ccs<hanging>::operator() (TElem&)
{
	TElem* castElem = dynamic_cast<TElem*>(iv.elem);
	if (castElem)
		nr->template compute_convection_shapes_elem<TElem, hanging>(convShapes, iv, thresh, castElem, u);
}


template <typename TDomain, typename TAlgebra>
void NeckRecorder<TDomain, TAlgebra>::compute_convection_shapes
(
	std::vector<std::vector<number> >* convShapes,
	const IntegrationVolume& iv,
	number thresh,
	ConstSmartPtr<GridFunction<TDomain, TAlgebra> > u
)
{
	typedef typename domain_traits<TDomain::dim>::DimElemList ElemList;
	if (!m_bSsIsRegular)
		boost::mpl::for_each<ElemList>(wrap_ccs<true>(convShapes, iv, thresh, u, this));
	else
		boost::mpl::for_each<ElemList>(wrap_ccs<false>(convShapes, iv, thresh, u, this));
}


template <typename TDomain, typename TAlgebra>
template <typename TElem, bool hanging>
void NeckRecorder<TDomain, TAlgebra>::compute_convection_shapes_elem
(
	std::vector<std::vector<number> >* convShapes,
	const IntegrationVolume& iv,
	number thresh,
	TElem* elem,
	ConstSmartPtr<GridFunction<TDomain, TAlgebra> > u
)
{
	ConstSmartPtr<TDomain> dom = m_spApprox->domain();

	// get FV1 geometry for element
	std::vector<typename TDomain::position_type> vCornerCoords;
	CollectCornerCoordinates(vCornerCoords, elem, dom->position_accessor(), false);

	typedef typename ChooseProperFVGeom<TDomain::dim, TElem, hanging, 1>::FVGeom fvgeom_type;
	static fvgeom_type& geo = GeomProvider<fvgeom_type>::get();
	try {geo.update(elem, &vCornerCoords[0], dom->subset_handler().get());}
	UG_CATCH_THROW("Failed creating FV1Geometry for element.");

	// init convection shapes on volume
	UG_COND_THROW(!m_spConvShape->template set_geometry_type<fvgeom_type>(geo),
		"Cannot init upwind scheme for element type.");

	// compute velocity at IPs
	const size_t nSCVF = geo.num_scvf_ips();
	std::vector<MathVector<dim> > velocityAtIPs(nSCVF);
	for (size_t s = 0; s < nSCVF; ++s)
	{
		const typename fvgeom_type::SCVF& scvf = geo.scvf(s);
		MathVector<dim>& vel = velocityAtIPs[s] = 0.0;

		const size_t nSh = scvf.num_sh();
		UG_COND_THROW(nSh != iv.vDofIndex.size(),
			"Number of shape functions for SCVF (" << nSh << ") does not match number of DoFs "
			"for potential function in IntegrationVolume (" << iv.vDofIndex.size() << ").")
		for (size_t sh = 0; sh < nSh; ++sh)
			VecScaleAppend(vel, DoFRef(*u, iv.vDofIndex[sh][_PHI_]), scvf.global_grad(sh));
	}

	// compute convection shapes
	std::vector<MathVector<dim> > scaledVelocityAtIPs(nSCVF);
	for (size_t f = 0; f < 4; ++f)
	{
		// scale velocity according to ion species
		for (size_t s = 0; s < nSCVF; ++s)
			VecScale(scaledVelocityAtIPs[s], velocityAtIPs[s], m_vConvConst[f]);

		// compute the actual convection shapes
		UG_COND_THROW(!m_spConvShape->update(&geo, &scaledVelocityAtIPs[0], NULL, false),
				"Cannot compute convection shapes.");

		// copy into given structure
		std::vector<std::vector<number> >& functionConvShapes = convShapes[f];
		if (functionConvShapes.size() < nSCVF)
			functionConvShapes.resize(nSCVF);
		size_t nCutEdges = 0;
		for (size_t s = 0; s < nSCVF; ++s)
		{
			const typename fvgeom_type::SCVF& scvf = geo.scvf(s);
			size_t from = scvf.from();
			size_t to = scvf.to();

			// find out whether current scvf is part of approximating manifold;
			// this is the case if the corners of the corresponding edge are on different sides

			const bool fromIsHangingInCoarse = from >= elem->num_vertices();
			const bool toIsHangingInCoarse = to >= elem->num_vertices();

			const typename TDomain::position_type& coCoFrom = geo.global_node_position(from);
			const typename TDomain::position_type& coCoTo = geo.global_node_position(to);
			const number refLenSq = VecDistanceSq(coCoFrom, coCoTo);

			int coStFrom = 0;
			int coStTo = 0;

			if (!fromIsHangingInCoarse)
				// natural corner of the elem
				coStFrom = hangingCornerState(elem->vertex(from), m_spApprox->domain(), thresh);
			else
				// hanging corner in coarse elem
				coStFrom = hangingCornerStateCoarse(coCoFrom, thresh, elem, dom, refLenSq);

			if (!toIsHangingInCoarse)
				// natural corner of the elem
				coStTo = hangingCornerState(elem->vertex(to), m_spApprox->domain(), thresh);
			else
				// hanging corner in coarse elem
				coStTo = hangingCornerStateCoarse(coCoTo, thresh, elem, dom, refLenSq);


			number weight = 0.0;
			if (coStFrom == CS_BELOW && coStTo == CS_ABOVE)
				weight = 1.0;
			else if (coStFrom == CS_ABOVE && coStTo == CS_BELOW)
				weight = -1.0;
			else if ((coStFrom == CS_BELOW && coStTo == (CS_BELOW | CS_ABOVE))
				|| (coStFrom == (CS_BELOW | CS_ABOVE) && coStTo == CS_ABOVE))
				weight = 0.5;
			else if ((coStFrom == (CS_BELOW | CS_ABOVE) && coStTo == CS_BELOW)
				|| (coStFrom == CS_ABOVE && coStTo == (CS_BELOW | CS_ABOVE)))
				weight = -0.5;

			if (weight == 0.0)
				continue;


			std::vector<number>& scvfConvShapes = functionConvShapes[nCutEdges];
			++nCutEdges;

			size_t nSh = scvf.num_sh();
			if (scvfConvShapes.size() < nSh)
				scvfConvShapes.resize(nSh);
			for (size_t sh = 0; sh < nSh; ++sh)
				scvfConvShapes[sh] = m_spConvShape->operator()(s, sh);
		}
	}
}


template <typename TDomain, typename TAlgebra>
void NeckRecorder<TDomain, TAlgebra>::record_current
(
	const std::string& fileName,
	number time,
	ConstSmartPtr<GridFunction<TDomain, TAlgebra> > u,
	number scale
)
{
#if NECK_RECORDER_DEBUG
	if (!testGrid.has_vertex_attachment(TestAttachment<TDomain>::aPosTest))
		testGrid.attach_to_vertices(TestAttachment<TDomain>::aPosTest);
#endif

	const std::string fnWoExt = FilenameAndPathWithoutExtension(fileName);
	std::string ext = GetFilenameExtension(fileName);
	if (ext == "")
		ext = "dat";

	// prepare if not yet done
	if (!bPreparedCurrent)
	{
		prepare(m_vIntegrationDataCurrent, true);
		bPreparedCurrent = true;
	}

#if NECK_RECORDER_DEBUG
	// construct outFile name for fluxes
	std::ostringstream offluxes(fnWoExt, std::ios::app);
	offluxes << "_fluxes.csv";

	// open file
	std::ofstream outFileFluxes(offluxes.str().c_str(), std::ios_base::out | std::ios_base::app);
	UG_COND_THROW(!outFileFluxes.is_open(), "File '" << offluxes.str() << "' could not be opened for writing.");

	// write header
	std::string header = std::string("coordX");
	if (dim >= 2) header.append(std::string(",coordY"));
	if (dim >= 3) header.append(std::string(",coordZ"));
	header.append(std::string(",fluxX"));
	if (dim >= 2) header.append(std::string(",fluxY"));
	if (dim >= 3) header.append(std::string(",fluxZ"));
	outFileFluxes << header << std::endl;
#endif

	// loop measurement zones
	std::vector<std::vector<std::vector<number> > > convShapes(4);
	const size_t nMZ = m_vIntegrationDataCurrent.size();
	for (size_t mz = 0; mz < nMZ; ++mz)
	{
		number currentK = 0.0;
		number currentNa = 0.0;
		number currentCl = 0.0;
		number currentA = 0.0;

		// loop intersecting volumes
		const size_t nIV = m_vIntegrationDataCurrent[mz].size();
		for (size_t v = 0; v < nIV; ++v)
		{
			const IntegrationVolume& iv = m_vIntegrationDataCurrent[mz][v];
			const size_t nCo = iv.vDofIndex.size();

			// compute upwinding convection shapes
			compute_convection_shapes(&convShapes[0], iv, m_vMeasCoords[mz], u);

			// loop IPs
			const size_t nIP = iv.vIPData.size();
			UG_COND_THROW(nIP > convShapes[_K_].size(),
				"Number of IPs in IntegrationVolume (" << nIP << ") does not match "
				"number of SCVF IPs of convection shapes (" << convShapes[_K_].size()
				<< ") on " << ElementDebugInfo(*m_spApprox->domain()->grid(), iv.elem) << "\b.");

			for (size_t ip = 0; ip < nIP; ++ip)
			{
#if NECK_RECORDER_DEBUG
				// get coordinates and normal
				typename TDomain::position_accessor_type::ValueType center;
				center[0] = iv.vIPData[ip].pos[0];
				center[1] = iv.vIPData[ip].pos[1];
#endif

				const IPData& ipd = iv.vIPData[ip];

				// construct function values and gradients at IP (by looping corners)
				number gradK = 0.0;
				number gradNa = 0.0;
				number gradCl = 0.0;
				number gradA = 0.0;

				number elFluxK = 0.0;
				number elFluxNa = 0.0;
				number elFluxCl = 0.0;
				number elFluxA = 0.0;

				for (size_t co = 0; co < nCo; ++co)
				{
					const std::vector<DoFIndex>& di = iv.vDofIndex[co];
					const number grad = ipd.vGradZ[co];

					gradK   += grad * DoFRef(*u, di[_K_]);
					gradNa  += grad * DoFRef(*u, di[_NA_]);
					gradCl  += grad * DoFRef(*u, di[_CL_]);
					gradA   += grad * DoFRef(*u, di[_A_]);

					elFluxK  += convShapes[_K_][ip][co] * DoFRef(*u, di[_K_]);
					elFluxNa += convShapes[_NA_][ip][co] * DoFRef(*u, di[_NA_]);
					elFluxCl += convShapes[_CL_][ip][co] * DoFRef(*u, di[_CL_]);
					elFluxA  += convShapes[_A_][ip][co] * DoFRef(*u, di[_A_]);
				}

				// construct current density
				const number currK = - m_vDiffConst[_K_] * gradK + elFluxK;
				const number currNa = - m_vDiffConst[_NA_] * gradNa + elFluxNa;
				const number currCl = - m_vDiffConst[_CL_] * gradCl + elFluxCl;
				const number currA = - m_vDiffConst[_A_] * gradA + elFluxA;

				currentK += ipd.weight * ipd.detJ * currK;
				currentNa += ipd.weight * ipd.detJ * currNa;
				currentCl -= ipd.weight * ipd.detJ * currCl;
				currentA -= ipd.weight * ipd.detJ * currA;

#ifdef NECK_RECORDER_DEBUG
				for (size_t d = 0; d < dim; ++d)
					outFileFluxes << center[d] << ",\t";

				outFileFluxes << "0.0,\t" << scale * ipd.weight * ipd.detJ * (currK) << "\n";
				//outFileFluxes << scale * normal[dim-1] * ipd.weight * ipd.detJ * (currNa) << ",\t";
				//outFileFluxes << scale * normal[dim-1] * ipd.weight * ipd.detJ * (-currCl) << ",\t";
				//outFileFluxes << scale * normal[dim-1] * ipd.weight * ipd.detJ * (-currA) << "\n";
#endif
			}
		}

		currentK *= scale;
		currentNa *= scale;
		currentCl *= scale;
		currentA *= scale;

#ifdef UG_PARALLEL
		// sum over processes
		if (pcl::NumProcs() > 1)
		{
			pcl::ProcessCommunicator com;
			number local[4] = {currentK, currentNa, currentCl, currentA};
			number global[4] = {0.0, 0.0, 0.0, 0.0};
			com.allreduce(&local, &global, 4, PCL_DT_DOUBLE, PCL_RO_SUM);

			currentK = global[0];
			currentNa = global[1];
			currentCl = global[2];
			currentA = global[3];
		}

		// check if this proc is output proc
		if (GetLogAssistant().is_output_process())
#endif
		{
			number current = currentK + currentNa + currentCl + currentA;

			// construct outFile name
			std::ostringstream ofnss(fnWoExt, std::ios::app);
			ofnss << "_" << m_vMeasZoneNames[mz] << "." << ext;

			// open file
			std::ofstream outFile(ofnss.str().c_str(), std::ios_base::out | std::ios_base::app);
			UG_COND_THROW(!outFile.is_open(), "File '" << ofnss.str() << "' could not be opened for writing.");

			// write record
			try
			{
				outFile << time << "\t" << current;
				if (m_bIndividualCurrents)
				{
					outFile << "\t" << currentK;
					outFile << "\t" << currentNa;
					outFile << "\t" << currentCl;
					outFile << "\t" << currentA;
				}
				outFile << "\n";
			}
			UG_CATCH_THROW("Output file " << ofnss.str() << " could not be written to.");

			// close file
			outFile.close();
		}
	}

#if NECK_RECORDER_DEBUG
	outFileFluxes.close();

	SubsetHandler sh(testGrid);
	std::ostringstream ofnss(fnWoExt, std::ios::app);
	ofnss << "_testGrid";
#ifdef UG_PARALLEL
	if (pcl::NumProcs() > 1)
		ofnss << "_p" << pcl::ProcRank();
#endif
	ofnss << ".ugx";
	UG_COND_THROW(!SaveGridToFile(testGrid, sh, ofnss.str().c_str(), TestAttachment<TDomain>::aPosTest),
		"Saving grid to file failed.");
#endif
}



template <typename TDomain, typename TAlgebra>
void NeckRecorder<TDomain, TAlgebra>::record_potential
(
	const std::string& fileName,
	number time,
	ConstSmartPtr<GridFunction<TDomain, TAlgebra> > u,
	number scale
)
{
	const std::string fnWoExt = FilenameAndPathWithoutExtension(fileName);
	std::string ext = GetFilenameExtension(fileName);
	if (ext == "")
		ext = "dat";

	// prepare if not yet done
	if (!bPreparedPot)
	{
		prepare(m_vIntegrationDataPot, false);
		bPreparedPot = true;
	}

	// loop measurement zones
	const size_t nMZ = m_vIntegrationDataPot.size();
	for (size_t mz = 0; mz < nMZ; ++mz)
	{
		number potInt = 0.0;
		number area = 0.0;

		// loop intersecting volumes
		const size_t nIV = m_vIntegrationDataPot[mz].size();
		for (size_t v = 0; v < nIV; ++v)
		{
			const IntegrationVolume& iv = m_vIntegrationDataPot[mz][v];
			const size_t nCo = iv.vDofIndex.size();

			// loop IPs
			const size_t nIP = iv.vIPData.size();
			for (size_t ip = 0; ip < nIP; ++ip)
			{
				const IPData& ipd = iv.vIPData[ip];

				// construct function values and gradients at IP (by looping corners)
				number phi = 0.0;

				for (size_t co = 0; co < nCo; ++co)
					phi += ipd.vShapes[co] * DoFRef(*u, iv.vDofIndex[co][_PHI_]);

				potInt += ipd.weight * ipd.detJ * phi;
				area += ipd.weight * ipd.detJ;
			}
		}

		potInt *= scale;

#ifdef UG_PARALLEL
		// sum over processes
		if (pcl::NumProcs() > 1)
		{
			pcl::ProcessCommunicator com;
			number local[2] = {potInt, area};
			number global[2];
			com.allreduce(local, global, 2, PCL_DT_DOUBLE, PCL_RO_SUM);
			potInt = global[0];
			area = global[1];
		}
		//UG_LOGN("zone " << mz << ": " << potInt << ", " << area);

		// check if this proc is output proc
		if (GetLogAssistant().is_output_process())
#endif
		{
			// construct outFile name
			std::ostringstream ofnss(fnWoExt, std::ios::app);
			ofnss << "_" << m_vMeasZoneNames[mz] << "." << ext;

			// open file
			std::ofstream outFile(ofnss.str().c_str(), std::ios_base::out | std::ios_base::app);
			UG_COND_THROW(!outFile.is_open(), "File '" << ofnss.str() << "' could not be opened for writing.");

			// write record
			try {outFile << time << "\t" << (potInt / area) << "\n";}
			UG_CATCH_THROW("Output file " << ofnss.str() << " could not be written to.");

			// close file
			outFile.close();
		}
	}
}





// explicit template specializations
#ifdef UG_DIM_1
	template class NeckRecorder<Domain1d, CPUAlgebra>;
	template class NeckRecorder<Domain1d, CPUBlockAlgebra<5> >;
#endif
#ifdef UG_DIM_2
	template class NeckRecorder<Domain2d, CPUAlgebra>;
	template class NeckRecorder<Domain2d, CPUBlockAlgebra<5> >;
#endif
#ifdef UG_DIM_3
	template class NeckRecorder<Domain3d, CPUAlgebra>;
	template class NeckRecorder<Domain3d, CPUBlockAlgebra<5> >;
#endif

} // namespace nernst_planck
} // namespace ug
