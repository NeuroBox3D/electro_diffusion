/*
 * morpho_gen.cpp
 *
 *  Created on: 07.09.2016
 *      Author: mbreit
 */

#include "morpho_gen.h"
#include "lib_grid/grid/geometry.h" // MakeGeometry3d
#include "lib_grid/file_io/file_io_ugx.h"  // GridWriterUGX
#include "lib_grid/algorithms/subset_util.h"  // AssignSubsetColors
#include "common/util/file_util.h"     // FindDirInStandardPaths
#include "common/util/string_util.h"   // GetFilenameExtension etc.
#include "lib_grid/algorithms/remeshing/delaunay_triangulation.h" // QualityGridGeneration
#include "lib_grid/algorithms/extrusion/extrusion.h" // Extrude
#include "lib_grid/algorithms/smoothing/manifold_smoothing.h"	// Laplacian smooth
#include "lib_grid/algorithms/element_side_util.h" // vertexGroupsMatch
#include "lib_grid/algorithms/debug_util.h"	// ElemDebugInfo
#include "lib_algebra/small_algebra/small_algebra.h" // Invert
#include "lib_grid/refinement/projectors/cylinder_projector.h"	// CylinderProjector
#include "lib_grid/refinement/projectors/sphere_projector.h"	// SphereProjector
#include "lib_grid/refinement/regular_refinement.h"	// Refine
#include "lib_grid/algorithms/geom_obj_util/volume_util.h"	// FixOrientation
#include "lib_grid/algorithms/geom_obj_util/face_util.h"	// Triangulate, CalculateFaceNormals, FixFaceOrientation
#include "lib_grid/algorithms/remove_duplicates_util.h"	// RemoveDuplicates
#include "lib_grid/algorithms/grid_generation/triangle_fill_sweep_line.h"	// TriangleFillSweepLine
#include "lib_grid/algorithms/selection_util.h"  // SelectLinkedFlatFaces
#include "lib_grid/refinement/global_multi_grid_refiner.h"	// GlobalMultiGridRefiner
#include "lib_disc/domain_util.h"	// LoadDomain
#include "lib_grid/file_io/file_io.h"  // SaveGridHierarchyTransformed
#include "lib_grid/algorithms/grid_generation/icosahedron.h"  // icosahedron generation

#include <cstdlib>  // rand
#include <limits>	// numeric_limits

#include "tetgen_config.h"

#ifdef TETGEN_15_ENABLED
	#include "tetgen.h"
#endif


namespace ug {
namespace nernst_planck {

MorphoGen::MorphoGen()
: m_grid(), m_sh(m_grid), m_sel(m_grid),
  m_shProj(m_grid), m_projHandler(&m_shProj),
  m_tmpHeadHeight(0.0),
  m_bFilAnisotropic(false),
  DENDRITE_LENGTH(1200.0),
  DENDRITE_RADIUS(200.0),
  MEMBRANE_RADIUS(10.0),
  MEMBRANE_ENVELOPE_RADIUS(20.0),
  DENDRITE_RIM_VERTICES(16),
  SPINE_NECK_RADIUS(100.0),
  SPINE_NECK_LENGTH(450.0),
  SPINE_HEAD_RADIUS(175.0),
  FILAMENT_WIDTH(10.0),
  FILAMENT_RIM_VERTICES(6),
  FILAMENT_ENVELOPE_RADIUS(10.0),
  EXTENSION_LENGTH(1.0e5),
  EXTENSION_COMPARTMENT_LENGTH(1.0e3),
  BOX_MARGIN(150.0),
  SHORT_EDGES_FACTOR(0.5),
  NECK_FILAMENTS(5),
  NUM_FILAMENTS(20),
  TETRAHEDRALIZATION_QUALITY(20.0),
  INNER_SI(0),
  OUTER_SI(1),
  MEM_SI(2),
  FIL_NECK_SI(3),
  MEM_INNER_BND_SI(4),
  MEM_OUTER_BND_SI(5),
  MEM_NOFLUX_BND_SI(6),
  SURF_CH_BND(7),
  NOFLUX_BND(8),
  DIRI_BND(9),
  PSD_INNER_SI(10),
  PSD_OUTER_SI(11),
  EXT_LEFT_SI(12),
  EXT_RIGHT_SI(13),
  INTF_LEFT_CONSTRD_SI(14),
  INTF_RIGHT_CONSTRD_SI(15),
  INTF_LEFT_NODE1D_SI(16),
  INTF_RIGHT_NODE1D_SI(17),
  INTF_LEFT_NODEHD_SI(18),
  INTF_RIGHT_NODEHD_SI(19),
  EXT_LEFT_BND_SI(20),
  EXT_RIGHT_BND_SI(21),
  USELESS_SI(22),
  INNER_FRONT_TMP_SI(14),
  OUTER_FRONT_TMP_SI(15),
  ENVELOPE_END_TMP_SI(16)
{
	// handle attachments and accessors
	m_grid.attach_to_vertices(aPosition);
	m_grid.attach_to_faces(aNormal);

	m_aaPos = AAPosition(m_grid, aPosition);
	m_aaNorm = AANormal(m_grid, aNormal);

	m_sh.set_default_subset_index(0);

	// set geometry in refinement projector
	m_projHandler.set_geometry(MakeGeometry3d(m_grid, aPosition));
}

void MorphoGen::set_num_neck_filaments(size_t nFil)
{
	NECK_FILAMENTS = nFil;
}

void MorphoGen::set_num_filaments(size_t nFil)
{
	NUM_FILAMENTS = nFil;
}


void MorphoGen::set_fil_anisotropic(bool filAniso)
{
	m_bFilAnisotropic = filAniso;
}



void MorphoGen::create_circle(const vector3& center, const vector3& axis, number radius, size_t numRimVertices)
{
	m_sel.clear();
	bool autoselEnabled = m_sel.autoselection_enabled();
	m_sel.enable_autoselection(true);

	// make sure axis is normed
	UG_COND_THROW(fabs(1.0-VecLength(axis)) > 1e-4, "Calling create_circle, axis must be normed.");

	// calculate trafo matrix
	vector3 axis2;
	if (axis[2] != 1.0)
	{
		axis2[2] = 0.0;
		axis2[1] = -axis[0];
		axis2[0] = axis[1];
	}
	else
	{
		axis2[0] = 1.0;
		axis2[1] = axis2[2] = 0.0;
	}
	axis2 /= VecLength(axis2);

	vector3 axis3;
	VecCross(axis3, axis, axis2);

	matrix33 trafo;
	trafo(0,0) = axis[0];
	trafo(1,0) = axis[1];
	trafo(2,0) = axis[2];

	trafo(0,1) = axis2[0];
	trafo(1,1) = axis2[1];
	trafo(2,1) = axis2[2];

	trafo(0,2) = axis3[0];
	trafo(1,2) = axis3[1];
	trafo(2,2) = axis3[2];


	// create edges by going round in the circle
	Vertex* firstVrt = *m_grid.create<RegularVertex>();
	MatVecMult(m_aaPos[firstVrt], trafo, vector3(0.0, radius, 0.0));
	VecAdd(m_aaPos[firstVrt], m_aaPos[firstVrt], center);

	number arc_length = 2.0 * PI / (number) numRimVertices;

	Vertex* lastVrt = firstVrt;
	for (size_t i = 1; i < numRimVertices; ++i)
	{
		// create a new vertex
		Vertex* vNew = *m_grid.create<RegularVertex>();
		MatVecMult(m_aaPos[vNew], trafo, vector3(0.0, cos(i*arc_length), -sin(i*arc_length)));
		VecScale(m_aaPos[vNew], m_aaPos[vNew], radius);
		VecAdd(m_aaPos[vNew], m_aaPos[vNew], center);

		// create a new edge
		m_grid.create<RegularEdge>(EdgeDescriptor(lastVrt, vNew));

		lastVrt = vNew;
	}

	// one edge is still missing
	m_grid.create<RegularEdge>(EdgeDescriptor(lastVrt, firstVrt));

//	restore selector
	m_sel.enable_autoselection(autoselEnabled);
}


void MorphoGen::create_shaft()
{
	m_sel.clear();
	bool autoselEnabled = m_sel.autoselection_enabled();
	m_sel.enable_autoselection(true);

	// create left end circle (new elems are auto-selected)
	vector3 left_end_center(0.0);
	left_end_center.coord(0) = -DENDRITE_LENGTH / 2.0;
	number radius = DENDRITE_RADIUS;

	create_circle(left_end_center, vector3(1,0,0), radius, DENDRITE_RIM_VERTICES);

	// assign subset
	m_sh.assign_subset(m_sel.begin<Vertex>(), m_sel.end<Vertex>(), MEM_OUTER_BND_SI);
	m_sh.assign_subset(m_sel.begin<Edge>(), m_sel.end<Edge>(), MEM_OUTER_BND_SI);

	// extrude to build shaft
	std::vector<Vertex*> vrts;
	vrts.assign(m_sel.vertices_begin(), m_sel.vertices_end());
	std::vector<Edge*> edges;
	edges.assign(m_sel.edges_begin(), m_sel.edges_end());

	// make extruded faces as regular as possible:
	// length of extrusion should be as near as possible to current edge length in circle
	number numExtrudes = floor(DENDRITE_LENGTH/(2*radius*sin(PI/DENDRITE_RIM_VERTICES)));
	size_t nExtrudes = std::max((size_t) numExtrudes, (size_t) 1);
	number extrude_length = DENDRITE_LENGTH / nExtrudes;

	vector3 extrudeDir(0.0);
	extrudeDir.coord(0) = extrude_length;

	for (size_t i = 0; i < nExtrudes; ++i)
		Extrude(m_grid, &vrts, &edges, NULL, extrudeDir, m_aaPos, EO_CREATE_FACES, NULL);

	// retriangulate
	//QualityGridGeneration(m_grid, m_sel.begin<Triangle>(), m_sel.end<Triangle>(),
	//					  m_aaPos, 30);

	// restore selector
	m_sel.enable_autoselection(autoselEnabled);
}


void MorphoGen::smooth_branching_point(number alpha, int numIterations, const vector3& spine_anchor)
{
	Grid::edge_traits::secure_container edges;
	Grid::face_traits::secure_container faces;
	Grid::volume_traits::secure_container vols;

	vector2 toOpp;
	vector2 toAnchor;

	for (int i = 0; i < numIterations; ++i)
	{
		std::queue<vector3> newPos;

		// iterate through all vertices
		geometry_traits<Vertex>::iterator it = m_sel.begin<Vertex>();
		geometry_traits<Vertex>::iterator it_end = m_sel.end<Vertex>();

		for (; it != it_end; ++it)
		{
			Vertex* vrt = *it;
			const vector3& vPos = m_aaPos[vrt];

			vector3 v;
			VecSet(v, 0.0);
			number weight = 0.0;

			// calculate smoothing vector relative to neighbors
			m_grid.associated_elements(edges, vrt);
			for (size_t ei = 0; ei < edges.size(); ++ei)
			{
				const vector3& opp = m_aaPos[GetConnectedVertex(edges[ei], vrt)];

				toOpp.coord(0) = opp.coord(0) - vPos.coord(0);
				toOpp.coord(1) = opp.coord(1) - vPos.coord(1);
				toAnchor.coord(0) = spine_anchor.coord(0) - 0.5*(vPos.coord(0) + opp.coord(0));
				toAnchor.coord(1) = spine_anchor.coord(1) - 0.5*(vPos.coord(1) + opp.coord(1));

				// smooth radially
				number prodLength = VecLength(toAnchor) * VecLength(toOpp);
				number w = prodLength ? fabs(VecProd(toAnchor, toOpp)) / prodLength : 0.0;

				//UG_LOGN("Vertex " << vrt << "at " << m_aaPos[vrt] << ": Add " << opp << " with weight " << w << ".");

				VecScaleAdd(v, 1.0, v, w, opp);
				weight += w;				
			}
/*
			m_grid.associated_elements(faces, vrt);
			for (size_t fi = 0; fi < faces.size(); ++fi)
			{
				Face* f = faces[fi];
				const vector3& opp = CalculateGridObjectCenter(
									 m_grid.get_opposing_object(vrt, f), m_aaPos);

				toOpp.coord(0) = opp.coord(0) - vPos.coord(0);
				toOpp.coord(1) = opp.coord(1) - vPos.coord(1);
				toAnchor.coord(0) = spine_anchor.coord(0) - 0.5*(vPos.coord(0) + opp.coord(0));
				toAnchor.coord(1) = spine_anchor.coord(1) - 0.5*(vPos.coord(1) + opp.coord(1));

				// smooth radially
				number prodLength = VecLength(toAnchor) * VecLength(toOpp);
				number w = prodLength ? fabs(VecProd(toAnchor, toOpp)) / prodLength : 0.0;

				VecScaleAdd(v, 1.0, v, w, opp);
				weight += w;
			}
*/

			newPos.push(m_aaPos[vrt]);
			if (weight > 0)
			{
				VecScale(v, v, 1.0 / weight);
				VecSubtract(v, v, m_aaPos[vrt]);
				VecScale(v, v, alpha);
				VecAdd(newPos.back(), m_aaPos[vrt], v);
			}
		}

		// assign new positions
		for (it = m_sel.begin<Vertex>(); it != it_end; ++it)
		{
			m_aaPos[*it] = newPos.front();
			newPos.pop();
		}
	}
}



bool MorphoGen::is_outside_spine_shaft(Vertex* v, const vector3& center, number radius)
{
	vector3* coord = &m_aaPos[v];
	number subtractZDistSq = center.coord(2) - coord->coord(2);
	subtractZDistSq = subtractZDistSq*subtractZDistSq;
	number distSqXY1 = VecDistanceSq(center, *coord) - subtractZDistSq;
	return coord->coord(2) < 0 || distSqXY1 > radius*radius;
}


bool MorphoGen::is_cut_by_spine_shaft
(
	Face* f,
	const vector3& center,
	number radius,
	std::vector<Edge*>& outCutEdges,
	std::vector<size_t>& outCutEdgeIndices,
	std::vector<vector3>& outCutPositions
)
{
	outCutEdges.clear();
	outCutEdgeIndices.clear();
	outCutPositions.clear();

	size_t nEdges = f->num_edges();

	for (size_t i = 0; i < nEdges; ++i)
	{
		const EdgeDescriptor& ed = f->edge_desc(i);
		bool outside = is_outside_spine_shaft(ed.vertex(0), center, radius);
		if (outside != is_outside_spine_shaft(ed.vertex(1), center, radius))
		{
			// this edge is cut, find concrete edge in grid
			typedef Grid::traits<Edge>::secure_container edge_list_type;
			edge_list_type el;
			m_grid.associated_elements(el, f);
			size_t el_sz = el.size();
			for (size_t e = 0; e < el_sz; ++e)
			{
				if (vertexGroupsMatch(el[e], ed))
				{
					outCutEdges.push_back(el[e]);
					outCutEdgeIndices.push_back(i);
					break;
				}
			}

			// calculate cut position
			// a quadratic equation needs to be solved here:
			// || v1 + a*(v2-v1) - x ||^2 = r^2
			vector2 v2minusv1;
			v2minusv1.coord(0) = m_aaPos[ed.vertex(1)].coord(0) - m_aaPos[ed.vertex(0)].coord(0);
			v2minusv1.coord(1) = m_aaPos[ed.vertex(1)].coord(1) - m_aaPos[ed.vertex(0)].coord(1);

			vector2 v1minusCenter;
			v1minusCenter.coord(0) = m_aaPos[ed.vertex(0)].coord(0) - center.coord(0);
			v1minusCenter.coord(1) = m_aaPos[ed.vertex(0)].coord(1) - center.coord(1);

			number scProd = VecProd(v2minusv1, v1minusCenter);
			number dist12Sq = VecProd(v2minusv1, v2minusv1);
			number dist1CenterSq = VecProd(v1minusCenter, v1minusCenter);

			number localCoord = - scProd;
			if (localCoord < 0)
				localCoord += sqrt(scProd*scProd - (dist1CenterSq-radius*radius)*dist12Sq);
			else
				localCoord -= sqrt(scProd*scProd - (dist1CenterSq-radius*radius)*dist12Sq);
			localCoord /= dist12Sq;

			outCutPositions.push_back(vector3());
			VecScaleAdd(outCutPositions.back(), (1.0-localCoord), m_aaPos[ed.vertex(0)],
						localCoord, m_aaPos[ed.vertex(1)]);
		}
	}

	return outCutEdges.size();
}


void MorphoGen::graft_spine()
{
	typedef Grid::traits<Edge>::secure_container edge_list;
	typedef Grid::traits<Face>::secure_container face_list;

	m_sel.clear();

	// find the face nearest do the anchor point of the spine
	vector3 anchorPoint(0.0);
	anchorPoint.coord(2) = DENDRITE_RADIUS;
	Face* anchorFace = FindClosestByCoordinate<Face>(anchorPoint, m_grid.begin<Face>(), m_grid.end<Face>(), m_aaPos);
	UG_COND_THROW(!anchorFace, "Could not find face near anchor point for spine.");

	// now search neighbors until one is found that is cut by the spine shaft's circumcircle
	m_grid.begin_marking();
	std::queue<Face*> searchQueue;
	m_grid.mark(anchorFace);
	searchQueue.push(anchorFace);
	anchorFace = NULL;

	std::vector<Edge*> cutEdges;
	std::vector<size_t> cutEdgeIndices;
	std::vector<vector3> cutPositions;
	while (!searchQueue.empty())
	{
		Face* f = searchQueue.front();
		searchQueue.pop();
		if (is_cut_by_spine_shaft(f, anchorPoint, SPINE_NECK_RADIUS, cutEdges, cutEdgeIndices, cutPositions))
		{
			anchorFace = f;
			while (!searchQueue.empty()) searchQueue.pop();
			break;
		}

		// push neighboring faces to queue
		edge_list el;
		m_grid.associated_elements(el, f);
		size_t el_sz = el.size();
		for (size_t e = 0; e < el_sz; ++e)
		{
			face_list fl;
			m_grid.associated_elements(fl, el[e]);
			size_t fl_sz = fl.size();
			for (size_t f1 = 0; f1 < fl_sz; ++f1)
			{
				if (!m_grid.is_marked(fl[f1]))
				{
					m_grid.mark(fl[f1]);
					searchQueue.push(fl[f1]);
				}
			}
		}
	}
	m_grid.end_marking();

	UG_COND_THROW(!anchorFace, "No face intersecting with spine cylinder could be found.");

	// iterating along the cutting border:
	// replace cut face by proper elements realizing the cutting border;
	int oldDefSI = m_sh.get_default_subset_index();
	m_sh.set_default_subset_index(MEM_OUTER_BND_SI);
	UG_COND_THROW(cutEdges.size() != 2, "Number of cut edges is not exactly 2.");
	Edge* firstEdge = cutEdges[0];
	Face* nextFace = NULL;
	Vertex* firstVertex;
	Vertex* v5; // new vertex backwards
	Vertex* v6; // new vertex forwards
	while (true)
	{
		bool firstStep = cutEdges[0] == firstEdge;
		bool lastStep = cutEdges[1] == firstEdge;

		// save next face
		face_list fl;
		m_grid.associated_elements(fl, cutEdges[1]);
		size_t fl_sz = fl.size();
		for (size_t f = 0; f < fl_sz; ++f)
		{
			if (fl[f] != anchorFace)
			{
				nextFace = fl[f];
				break;
			}
		}

		// get vertices (to ensure they are in the rifght order: use descriptors)
		Vertex* v1 = anchorFace->edge_desc(cutEdgeIndices[0]).vertex(0);
		Vertex* v2 = anchorFace->edge_desc(cutEdgeIndices[0]).vertex(1);
		Vertex* v3 = anchorFace->edge_desc(cutEdgeIndices[1]).vertex(0);
		Vertex* v4 = anchorFace->edge_desc(cutEdgeIndices[1]).vertex(1);

		// position new vertices
		if (firstStep)
		{
			firstVertex = v5 = *m_grid.create<RegularVertex>();
			m_aaPos[v5] = cutPositions[0];
		}
		else v5 = v6;

		if (!lastStep)
		{
			v6 = *m_grid.create<RegularVertex>();
			m_aaPos[v6] = cutPositions[1];
		}
		else v6 = firstVertex;

		// create new faces
		Face* circleFace;	// is to point to a face at the circle rim
		switch ((cutEdgeIndices[1] - cutEdgeIndices[0] + 4) % 4)
		{
			case 2: // opposing sides
			{
				// cut into two quadrilaterals
				m_grid.create<Quadrilateral>(QuadrilateralDescriptor(v1, v5, v6, v4));
				circleFace = *m_grid.create<Quadrilateral>(QuadrilateralDescriptor(v5, v2, v3, v6));
				break;
			}
			case 3: // edges must be 0 and 3
			{
				// cut into four triangles
				Vertex* v7 = anchorFace->vertex((cutEdgeIndices[0]+2)%4); // missing vertex of quad
				circleFace = *m_grid.create<Triangle>(TriangleDescriptor(v1, v5, v6));
				m_grid.create<Triangle>(TriangleDescriptor(v2, v7, v5));
				m_grid.create<Triangle>(TriangleDescriptor(v7, v3, v6));
				m_grid.create<Triangle>(TriangleDescriptor(v5, v7, v6));
				break;
			}
			case 1:
			{
				// cut into four triangles
				Vertex* v7 = anchorFace->vertex((cutEdgeIndices[1]+2)%4); // missing vertex of quad
				m_grid.create<Triangle>(TriangleDescriptor(v1, v5, v7));
				circleFace = *m_grid.create<Triangle>(TriangleDescriptor(v3, v6, v5));
				m_grid.create<Triangle>(TriangleDescriptor(v4, v7, v6));
				m_grid.create<Triangle>(TriangleDescriptor(v5, v6, v7));
				break;
			}
			default: UG_THROW("Invalid cut edge indices: " << cutEdgeIndices[0]
							  << ", " << cutEdgeIndices[1] << ".");
		}

		// select circle edge
		EdgeDescriptor ed(v5, v6);
		typedef Grid::traits<Edge>::secure_container edge_list_type;
		edge_list_type el;
		m_grid.associated_elements(el, circleFace);
		size_t el_sz = el.size();
		for (size_t e = 0; e < el_sz; ++e)
		{
			if (vertexGroupsMatch(el[e], ed))
			{
				m_sel.select(el[e]);
				break;
			}
		}

		// remove old face and edge
		m_grid.erase(anchorFace);
		if (!firstStep)
			m_grid.erase(cutEdges[0]);

		// remove first edge in last step; then break
		if (lastStep)
		{
			m_grid.erase(cutEdges[1]);
			break;
		}

		// prepare next step
		Edge* connectingEdge = cutEdges[1];

		UG_COND_THROW(!nextFace, "Next face could not be found.");
		anchorFace = nextFace;
		is_cut_by_spine_shaft(anchorFace, anchorPoint, SPINE_NECK_RADIUS, cutEdges, cutEdgeIndices, cutPositions);

		// it should not be possible to have other than 2 cut edges
		UG_COND_THROW(cutEdges.size() != 2, "Number of cut edges is not exactly 2.");

		// reorder cut edges if need be
		if (cutEdges[0] != connectingEdge)
		{
			std::swap(cutEdges[0], cutEdges[1]);
			std::swap(cutEdgeIndices[0], cutEdgeIndices[1]);
			std::swap(cutPositions[0], cutPositions[1]);
		}
	}
	m_sh.set_default_subset_index(oldDefSI);

	// retain rim edges
	std::vector<Edge*> rimEdges(m_sel.edges_begin(), m_sel.edges_end());

	// remove faces from inside spine shaft
	anchorFace = FindClosestByCoordinate<Face>(anchorPoint, m_grid.begin<Face>(), m_grid.end<Face>(), m_aaPos);
	UG_COND_THROW(!anchorFace, "Could not find face near anchor point for spine.");
	m_sel.select(anchorFace);

	SelectionFill<Face>(m_sel, IsSelected(m_sel));
	m_sel.deselect(m_sel.edges_begin(), m_sel.edges_end());
	SelectInnerSelectionEdges(m_sel);
	SelectInnerSelectionVertices(m_sel);

	m_grid.erase(m_sel.begin<Face>(), m_sel.end<Face>());
	m_grid.erase(m_sel.begin<Edge>(), m_sel.end<Edge>());
	m_grid.erase(m_sel.begin<Vertex>(), m_sel.end<Vertex>());


	// regularize rim (kind of Laplacian smoothing along rim)
	m_sel.select(rimEdges.begin(), rimEdges.end());
	SelectAssociatedVertices(m_sel, m_sel.begin<Edge>(), m_sel.end<Edge>());
	std::vector<Vertex*> rimVerts(m_sel.begin<Vertex>(), m_sel.end<Vertex>());
	std::vector<vector3> rimCoordsNew(rimVerts.size());
	size_t nVrt = rimVerts.size();

	edge_list el;

	for (size_t i = 0; i < 8; ++i)
	{
		for (size_t rvi = 0; rvi < nVrt; ++rvi)
		{
			Vertex* v = rimVerts[rvi];

			// find connecting rim edges
			Edge* connEdge[2];
			size_t foundConnEdges = 0;
			m_grid.associated_elements(el, v);
			size_t el_sz = el.size();
			for (size_t e = 0; e < el_sz; ++e)
			{
				if (m_sel.is_selected(el[e]))
				{
					connEdge[foundConnEdges] = el[e];
					++foundConnEdges;
				}
			}
			UG_COND_THROW(foundConnEdges != 2, "Not exactly two connecting rim edges for rim vertex.");

			// compare edge lengths
			number el0 = EdgeLength(connEdge[0], m_aaPos);
			number el1 = EdgeLength(connEdge[1], m_aaPos);

			// move vertex into direction of longer edge
			Edge* moveTo = el0 >= el1 ? connEdge[0] : connEdge[1];
			number fac = el0 >= el1 ? (el0-el1)/(4.0*el0) : (el1-el0)/(4.0*el1);
			Vertex* opp = GetOpposingSide(m_grid, moveTo, v);
			VecScaleAdd(rimCoordsNew[rvi], 1.0-fac, m_aaPos[v], fac, m_aaPos[opp]);

			// shift radially to achieve correct (x-y) distance
			el0 = rimCoordsNew[rvi].coord(0) - anchorPoint.coord(0);
			el1 = rimCoordsNew[rvi].coord(1) - anchorPoint.coord(1);
			el0 = sqrt(el0*el0 + el1*el1);
			fac = SPINE_NECK_RADIUS / el0;
			VecScaleAdd(rimCoordsNew[rvi], 1.0-fac, anchorPoint, fac, rimCoordsNew[rvi]);

			// finally, shift z coordinate to ensure vertex is inside original face
			// to do that solve: v + k*e1 + l*e2 = w + m*(0, 0, -1)
			// with e1, e2being  the two vectors defining the face associated to long edge
			// and w being the current new position of v
			face_list fl;
			m_grid.associated_elements(fl, moveTo);
			UG_COND_THROW(fl.size() != 1, "Not exactly one face associated to rim edge.")
			Face* rimFace = fl[0];

			const EdgeDescriptor& ed1 = rimFace->edge_desc(0);
			const EdgeDescriptor& ed2 = rimFace->edge_desc(1);
			vector3 e1, e2;
			VecScaleAdd(e1, 1.0, m_aaPos[ed1.vertex(1)], -1.0, m_aaPos[ed1.vertex(0)]);
			VecScaleAdd(e2, 1.0, m_aaPos[ed2.vertex(1)], -1.0, m_aaPos[ed2.vertex(0)]);
			if (fabs(e2.coord(1)) < 1e-4*fabs(e2.coord(0))) // ensure diagonal is non-zero
				std::swap(e1, e2);
			vector3 diff;
			VecScaleAdd(diff, 1.0, rimCoordsNew[rvi], -1.0, m_aaPos[v]);

			number k = (e2.coord(1) * diff.coord(0) - e2.coord(0) * diff.coord(1))
					   / (e2.coord(1) * e1.coord(0) - e2.coord(0) * e1.coord(1));
			number l = (diff.coord(1) - k*e1.coord(1)) / e2.coord(1);
			number m = diff.coord(2) - k*e1.coord(2) - l*e2.coord(2);

			rimCoordsNew[rvi].coord(2) -= m;
		}

		// copy newly calculated positions to attachments
		for (size_t rvi = 0; rvi < nVrt; ++rvi)
			m_aaPos[rimVerts[rvi]] = rimCoordsNew[rvi];
	}

	// remove short edges
	// TODO: replace this calculation by member: m_elemLength
	number arc_length = 2.0 * PI / (number) DENDRITE_RIM_VERTICES;
	number numExtrudes = floor(DENDRITE_LENGTH/(2*DENDRITE_RADIUS*sin(arc_length/2.0)));
	size_t nExtrudes = std::max((size_t) numExtrudes, (size_t) 1);
	number extrude_length = DENDRITE_LENGTH / nExtrudes;
	number min_length = SHORT_EDGES_FACTOR * extrude_length;
	for (size_t rvi = 0; rvi < nVrt; ++rvi)
	{
		Vertex* v = rimVerts[rvi];
		m_grid.associated_elements(el, v);
		size_t el_sz = el.size();
		for (size_t e = 0; e < el_sz; ++e)
		{
			if (!m_sel.is_selected(el[e]))
			{
				if (EdgeLength(el[e], m_aaPos) < min_length)
				{
					Vertex* opp = GetOpposingSide(m_grid, el[e], v);
					UG_COND_THROW(!opp, "Opposing side not found.");
					MergeVertices(m_grid, v, opp);

					// edge list might be invalid after merge
					m_grid.associated_elements(el, v);
					el_sz = el.size();
					e = 0;
				}
			}
		}
	}

	// extrude twice (with reduced length)
	std::vector<Vertex*> vrts;
	vrts.assign(m_sel.vertices_begin(), m_sel.vertices_end());
	std::vector<Edge*> edges;
	edges.assign(m_sel.edges_begin(), m_sel.edges_end());

	UG_COND_THROW(!edges.size(), "No rim edges present.");
	extrude_length = EdgeLength(edges[0], m_aaPos);
	numExtrudes = floor(SPINE_NECK_LENGTH / extrude_length);
	nExtrudes = std::max((size_t) numExtrudes, (size_t) 1);
	extrude_length = SPINE_NECK_LENGTH / nExtrudes;
	vector3 extrudeDir(0.0);
	extrudeDir.coord(2) = 0.8*extrude_length;

	Extrude(m_grid, &vrts, &edges, NULL, extrudeDir, m_aaPos, EO_CREATE_FACES, NULL);
	Extrude(m_grid, &vrts, &edges, NULL, extrudeDir, m_aaPos, EO_CREATE_FACES, NULL);

	// save current spine vertices for later usage in other methods
	m_tmpSpineEdges = edges;

	// project to highest z coordinate in selection
	size_t sz = vrts.size();
	number maxZ = 0.0;
	for (size_t i = 0; i < sz; ++i)
		maxZ = std::max(maxZ, m_aaPos[vrts[i]].coord(2));
	for (size_t i = 0; i < sz; ++i)
		m_aaPos[vrts[i]].coord(2) = maxZ;


	// smoothing of branching point
	m_sel.deselect(vrts.begin(), vrts.end());
	smooth_branching_point(0.25, 8, anchorPoint);


	// remove short edges again, but this time merge at other vertex
	// TODO: replace this calculation by member: m_elemLength
	for (size_t rvi = 0; rvi < nVrt; ++rvi)
	{
		Vertex* v = rimVerts[rvi];
		m_grid.associated_elements(el, v);
		size_t el_sz = el.size();
		for (size_t e = 0; e < el_sz; ++e)
		{
			if (!m_sel.is_selected(el[e]))
			{
				if (EdgeLength(el[e], m_aaPos) < min_length)
				{
					Vertex* opp = GetOpposingSide(m_grid, el[e], v);
					UG_COND_THROW(!opp, "Opposing side not found.");
					MergeVertices(m_grid, opp, v);
					v = opp;

					// edge list might be invalid after merge
					m_grid.associated_elements(el, v);
					el_sz = el.size();
					e = 0;
				}
			}
		}
	}


	// extrude to full spine neck length
	extrude_length = EdgeLength(edges[0], m_aaPos);
	number rest_length = SPINE_NECK_LENGTH - (maxZ - DENDRITE_RADIUS);
	numExtrudes = floor(rest_length / extrude_length);
	nExtrudes = std::max((size_t) numExtrudes, (size_t) 1);
	extrude_length = rest_length / nExtrudes;
	extrudeDir.coord(2) = extrude_length;

	for (size_t i = 0; i < nExtrudes; ++i)
		Extrude(m_grid, &vrts, &edges, NULL, extrudeDir, m_aaPos, EO_CREATE_FACES, NULL);


	// make spine head
	// calculate center, current angle
	number curr_z = DENDRITE_RADIUS + SPINE_NECK_LENGTH;
	number z_center = curr_z + sqrt(SPINE_HEAD_RADIUS*SPINE_HEAD_RADIUS - SPINE_NECK_RADIUS*SPINE_NECK_RADIUS);
	number curr_angle = asin(SPINE_NECK_RADIUS/SPINE_HEAD_RADIUS);
	number step_angle = (PI - curr_angle) / 5.0;
	number curr_radius = SPINE_NECK_RADIUS;

	for (size_t i =  0; i < 5; ++i)
	{
		curr_angle += step_angle;
		number dz = -curr_z;
		curr_z = z_center - SPINE_HEAD_RADIUS*cos(curr_angle);
		dz += curr_z;
		number fac_rad = 1.0 / curr_radius;
		curr_radius = SPINE_HEAD_RADIUS*sin(curr_angle);
		fac_rad *= curr_radius;

		// extrude in z-direction
		extrudeDir.coord(2) = dz;
		Extrude(m_grid, &vrts, &edges, NULL, extrudeDir, m_aaPos, EO_CREATE_FACES, NULL);

		// scale around center of selection
		vector3 center = CalculateCenter(vrts.begin(), vrts.end(), m_aaPos);

		size_t nVrt = vrts.size();
		for (size_t i = 0; i < nVrt; ++i)
		{
			vector3& v = m_aaPos[vrts[i]];
			VecSubtract(v, v, center);
			VecScale(v, v, fac_rad);
			VecAdd(v, v, center);
		}
	}

	// close spine head by merging top vertices
	Vertex* v = MergeMultipleVertices(m_grid, vrts.begin(), vrts.end());
	m_tmpHeadHeight = m_aaPos[v].z() - (SPINE_NECK_LENGTH + DENDRITE_RADIUS);

	// assign PSD subset
	m_sel.clear<Vertex>();
	m_sel.clear<Edge>();
	m_sel.clear<Face>();
	m_sel.select(v);
	ExtendSelection(m_sel, 1);
	SelectAssociatedEdges(m_sel, m_sel.begin<Face>(), m_sel.end<Face>());
	SelectAssociatedVertices(m_sel, m_sel.begin<Edge>(), m_sel.end<Edge>());

	m_sh.assign_subset(m_sel.begin<Vertex>(), m_sel.end<Vertex>(), PSD_OUTER_SI);
	m_sh.assign_subset(m_sel.begin<Edge>(), m_sel.end<Edge>(), PSD_OUTER_SI);
	m_sh.assign_subset(m_sel.begin<Face>(), m_sel.end<Face>(), PSD_OUTER_SI);
}

void MorphoGen::defect_for_filament_distribution
(
	const DenseVector<VariableArray1<number> >& bndPos,
	const DenseVector<VariableArray1<number> >& solNew,
	const DenseVector<VariableArray1<number> >& solOld,
	DenseVector<VariableArray1<number> >& defOut,
	number dt
) const
{
	// reset defect
	defOut = 0.0;

	size_t nFil = solNew.size() / 4;
	size_t nBnd = bndPos.size() / 2;
	for (size_t i = 0; i < nFil; ++i)
	{
		const number& xi1 = solNew[4*i];
		const number& xi2 = solNew[4*i+1];
		const number& vi1 = solNew[4*i+2];
		const number& vi2 = solNew[4*i+3];

		number& dxi1 = defOut[4*i];
		number& dxi2 = defOut[4*i+1];
		number& dvi1 = defOut[4*i+2];
		number& dvi2 = defOut[4*i+3];

		// force, coupling with other filaments
		for (size_t j = 0; j < nFil; ++j)
		{
			if (i == j) continue;

			const number& xj1 = solNew[4*j];
			const number& xj2 = solNew[4*j+1];

			// distance
			number dist3 = 0.0;
			number diffX = xi1-xj1;
			dist3 += diffX*diffX;
			number diffY = xi2-xj2;
			dist3 += diffY*diffY;
			dist3 = sqrt(dist3);
			dist3 = 1.0 / (dist3*dist3*dist3);

			dvi1 += diffX * dist3;
			dvi2 += diffY * dist3;
		}

		// force, coupling with boundary
		for (size_t k = 0; k < nBnd; ++k)
		{
			// distance
			number dist3 = 0.0;
			number diffX = xi1 - bndPos[2*k];
			dist3 += diffX*diffX;
			number diffY = xi2 - bndPos[2*k+1];
			dist3 += diffY*diffY;
			dist3 = sqrt(dist3);
			dist3 = 1.0 / (dist3*dist3*dist3);

			dvi1 += diffX * dist3;
			dvi2 += diffY * dist3;
		}

		// circular pseudo-force (to force filaments to the inner)
		number r = sqrt(xi1*xi1 + xi2*xi2);
		dvi1 -= 2*nFil * xi1 * r;
		dvi2 -= 2*nFil * xi2 * r;

		// dampening (for energy dissipation); the stronger the more charges
		dvi1 -= 10*nFil*vi1;
		dvi2 -= 10*nFil*vi2;

		// add stiffness terms for position defects
		dxi1 += vi1;
		dxi2 += vi2;
	}

	// multiply all stiffness by (-dt)
	defOut *= -dt;

	// mass terms
	for (size_t i = 0; i < nFil; ++i)
	{
		defOut[4*i] += solNew[4*i] - solOld[4*i];
		defOut[4*i+1] += solNew[4*i+1] - solOld[4*i+1];
		defOut[4*i+2] += solNew[4*i+2] - solOld[4*i+2];
		defOut[4*i+3] += solNew[4*i+3] - solOld[4*i+3];
	}
}


void MorphoGen::jacobian_for_filament_distribution
(
	const DenseVector<VariableArray1<number> >& bndPos,
	const DenseVector<VariableArray1<number> >& solNew,
	DenseMatrix<VariableArray2<number> >& jacOut,
	number dt
) const
{
	// reset jacobian
	jacOut = 0.0;

	size_t nFil = solNew.size() / 4;
	size_t nBnd = bndPos.size() / 2;
	for (size_t i = 0; i < nFil; ++i)
	{
		// solution
		const number& xi1 = solNew[4*i];
		const number& xi2 = solNew[4*i+1];

		// diagonal entries v/x
		number& jvi1xi1 = jacOut(4*i+2, 4*i);
		number& jvi1xi2 = jacOut(4*i+2, 4*i+1);
		number& jvi2xi1 = jacOut(4*i+3, 4*i);
		number& jvi2xi2 = jacOut(4*i+3, 4*i+1);

		// diagonal entries x/v
		number& jxi1vi1 = jacOut(4*i, 4*i+2);
		number& jxi2vi2 = jacOut(4*i+1, 4*i+3);

		// diagonal entries v/v
		number& jvi1vi1 = jacOut(4*i+2, 4*i+2);
		number& jvi2vi2 = jacOut(4*i+3, 4*i+3);

		for (size_t j = 0; j < nFil; ++j)
		{
			const number& xj1 = solNew[4*j];
			const number& xj2 = solNew[4*j+1];

			// off-diagonal entries
			number& jvi1xj1 = jacOut(4*i+2, 4*j);
			number& jvi1xj2 = jacOut(4*i+2, 4*j+1);
			number& jvi2xj1 = jacOut(4*i+3, 4*j);
			number& jvi2xj2 = jacOut(4*i+3, 4*j+1);

			// diagonal coupling assembled at off-diagonal entries
			// but boundary terms only for diagonal
			if (i == j)
			{

				// coupling with boundary
				for (size_t k = 0; k < nBnd; ++k)
				{
					// distance
					number dist = 0.0;
					number diffX = xi1-bndPos[2*k];
					dist += diffX*diffX;
					number diffY = xi2-bndPos[2*k+1];
					dist += diffY*diffY;
					dist = sqrt(dist);
					number dist3 = 1.0 / (dist*dist*dist);
					number dist5 = dist3 / (dist*dist);

					jvi1xi1 += dist3;
					jvi2xi2 += dist3;

					number add;
					add = 3.0 * diffX * diffX * dist5;
					jvi1xi1 -= add;
					add = 3.0 * diffX * diffY * dist5;
					jvi1xi2 -= add;
					jvi2xi1 -= add;
					add = 3.0 * diffY * diffY * dist5;
					jvi2xi2 -= add;
				}

				// circular pseudo-force (to force filaments to the inner)
				// At the circular bnd, this creates a force like that between two
				// charges at 0.1*SPINE_RADIUS apart. Exponential decay towards center.
				number r = sqrt(xi1*xi1 + xi2*xi2);
				if (r)
				{
					jvi1xi1 -= 2*nFil * (r + xi1*xi1 / r);
					jvi1xi2 -= 2*nFil * xi1*xi2 / r;
					jvi2xi1 -= 2*nFil * xi2*xi1 / r;
					jvi2xi2 -= 2*nFil * (r + xi2*xi2 / r);
				}

				// dampening (for energy dissipation)
				jvi1vi1 -= 10*nFil;
				jvi2vi2 -= 10*nFil;
			}
			else
			{
				// distance
				number dist = 0.0;
				number diffX = xi1-xj1;
				dist += diffX*diffX;
				number diffY = xi2-xj2;
				dist += diffY*diffY;
				dist = sqrt(dist);
				number dist3 = 1.0 / (dist*dist*dist);
				number dist5 = dist3 / (dist*dist);

				jvi1xj1 -= dist3;	jvi1xi1 += dist3;
				jvi2xj2 -= dist3;	jvi2xi2 += dist3;

				number add;
				add = 3.0 * diffX * diffX * dist5;
				jvi1xj1 += add; 		jvi1xi1 -= add;
				add = 3.0 * diffX * diffY * dist5;
				jvi1xj2 += add; 		jvi1xi2 -= add;
				jvi2xj1 += add; 		jvi2xi1 -= add;
				add = 3.0 * diffY * diffY * dist5;
				jvi2xj2 += add; 		jvi2xi2 -= add;
			}
		}

		// add stiffness terms for position defects
		jxi1vi1 += 1.0;
		jxi2vi2 += 1.0;
	}

	// multiply all stiffness by (-dt)
	jacOut *= -dt;

	// mass terms
	for (size_t i = 0; i < nFil; ++i)
	{
		jacOut(4*i, 4*i) += 1.0;
		jacOut(4*i+1, 4*i+1) += 1.0;
		jacOut(4*i+2, 4*i+2) += 1.0;
		jacOut(4*i+3, 4*i+3) += 1.0;
	}
}


void MorphoGen::distribute_filaments(const DenseVector<VariableArray1<number> >& bndPos, DenseVector<VariableArray1<number> >& posOut)
{
	// We get the bounding vertices (all in one plane)
	// and want to position nFil filaments with maximal
	// distance from each other as well as from the bounding vertices.
	// To achieve this, filament positions and bounding vertices
	// are treated as point charges. The final positions are given as
	// solution to the problem of finding a configuration where all
	// effective forces on the (filament) charges are zero.
	// This solution is found by solving the instationary problem given
	// by Newton's law of motion until no motion occurs any more.
	// To prevent solutions with filaments outside of the spine,
	// a corrective force effective only near the spine radius is added.
	// To prevent oscillations, a friction term is introduced.

	// The problem:
	// \f[
	//
	//     \dot{\vec{v}}_{i} & = \sum\limits _{j\neq i}\frac{\vec{x}_{i}-\vec{x}_{j}}{\|\vec{x}_{i}-\vec{x}_{j}\|^{3}}\;
	//                           + \;\sum\limits _{k}\frac{\vec{x}_{i}-\vec{b}_{k}}{\|\vec{x}_{i}-\vec{b}_{k}\|^{3}}\;
	//                           - \;2N\left\Vert \vec{x}_{i}\right\Vert \vec{x}_{i}\;
	//                           - \;10N\left\Vert \vec{v}_{i}\right\Vert \vec{v}_{i}\\
	//     \dot{\vec{x}}_{i} & = \vec{v}_{i}
	// \F]

	// This is far from being the best possible way to do this, but it seems to work.

	// start solution
	number minX = 0.0;
	number maxX = 0.0;
	number minY = 0.0;
	number maxY = 0.0;
	size_t nBnd = bndPos.size() / 2;
	for (size_t i = 0; i < nBnd; ++i)
	{
		minX = std::min(minX, bndPos[2*i]);
		maxX = std::max(maxX, bndPos[2*i]);
		minY = std::min(minY, bndPos[2*i+1]);
		maxY = std::max(maxY, bndPos[2*i+1]);
	}
	number diff = maxX - minX;
	minX += 0.25*diff;
	maxX -= 0.25*diff;
	diff = maxY - minY;
	minY += 0.25*diff;
	maxY-= 0.25*diff;

	size_t nFil = posOut.size() / 2;
	size_t perDim = (size_t) ceil(sqrt((number)nFil));
	size_t nSteps = perDim == 1 ? 1 : perDim - 1;

	// start solution with start positions on lattice, start velocity 0
	DenseVector<VariableArray1<number> > x_old, x_new;
	x_old.resize(4*nFil);
	x_new.resize(4*nFil);
	for (size_t i = 0; i < nFil; ++i)
	{
		x_new[4*i] = x_old[4*i] = minX + (i % perDim)*(maxX-minX) / nSteps;
		x_new[4*i+1] = x_old[4*i+1] = minY + (i / perDim)*(maxY-minY) / nSteps;
		x_new[4*i+2] = x_old[4*i+2] = 0.0;
		x_new[4*i+3] = x_old[4*i+3] = 0.0;
	}

	number dt = 0.05 * std::min((maxX-minX) / nSteps, (maxY-minY) / nSteps);


	// create algebra elements
	DenseVector<VariableArray1<number> > def;
	def.resize(4*nFil);
	DenseMatrix<VariableArray2<number> > jac;
	jac.resize(4*nFil, 4*nFil);
	DenseVector<VariableArray1<number> > cor;
	cor.resize(4*nFil);

	size_t maxIt = 200;

	// time stepping
	size_t tstep = 1;
	size_t maxStep = 200;
	while (tstep <= maxStep)
	{
//UG_LOGN("\nStep " << tstep);

		// start defect
		defect_for_filament_distribution(bndPos, x_new, x_old, def, dt);
		number normSqInit = VecNormSquared(def);
		number normSq = normSqInit;
		number newNormSq;
//UG_LOGN("dt = " << dt);
//UG_LOGN("iter 0" << "   " << normSq);

		size_t it = 1;
		while (normSq > 1e-16*normSqInit && it <= maxIt)
		{
			// assemble jacobian
			jacobian_for_filament_distribution(bndPos, x_new, jac, dt);

//UG_LOGN(jac);

			// invert
			UG_COND_THROW(!Invert(jac), "Failed to invert Jacobian for filament position calculation.")

			// correction
			cor = jac*def;

			// perform line search to prevent over-stepping boundary
			number lambda = 1.0;
			do
			{
				if (lambda == 1.0)
					x_new -= cor;
				else
					x_new += lambda*cor;

				defect_for_filament_distribution(bndPos, x_new, x_old, def, dt);
				newNormSq = VecNormSquared(def);

				lambda *= 0.5;
			}
			while (newNormSq >= 10*normSq && lambda > 1e-3);


//UG_LOGN(posOut);

			// new defect
			defect_for_filament_distribution(bndPos, x_new, x_old, def, dt);
			normSq = VecNormSquared(def);

//			number newNormSq = VecNormSquared(def);
//			number red = newNormSq / normSq;
//			normSq = newNormSq;
//UG_LOGN("iter " << it << "   " << normSq << "   " << red << "   " << normSq / normSqInit);
//UG_LOGN(def);

			++it;
		}

		UG_COND_THROW(normSq != normSq, "Newton iteration ended with NaN defect in time step " << tstep << ".");
		UG_COND_THROW(normSq > 1e-16*normSqInit, "Newton iteration did not converge in time step " << tstep << ".");

//UG_LOGN(x_new);
		// calculate current maximal force on and velocity of charges
		VecSubtract(def, def, x_new);
		VecAdd(def, def, x_old);
		VecScale(def, def, 1.0/dt);

		number maxFSq = 0.0;
		number maxVSq = 0.0;
		for (size_t i = 0; i < nFil; ++i)
		{
			maxVSq = std::max(maxVSq, def[4*i]*def[4*i] + def[4*i+1]*def[4*i+1]);
			maxFSq = std::max(maxFSq, def[4*i+2]*def[4*i+2] + def[4*i+3]*def[4*i+3]);
		}

		// if velocity small then break
		if (maxFSq < 1e-8) break;

		// set new dt according to current velocity
		dt = std::min(1.0, std::min(10.0/maxFSq, 0.01/maxVSq));

		x_old = x_new;
		++tstep;
	}

	UG_COND_THROW(tstep > maxStep, "Time stepping did not finish in " << maxStep << " steps.");

	// write positions to output
	for (size_t i = 0; i < nFil; ++i)
	{
		posOut[2*i] = x_new[4*i];
		posOut[2*i+1] = x_new[4*i+1];
	}
}



void MorphoGen::make_neck_filaments()
{
	// calculate filament positions
	size_t nBnd = m_tmpSpineEdges.size();
	size_t nFil = NECK_FILAMENTS;

	DenseVector<VariableArray1<number> > bndPos;
	bndPos.resize(2*nBnd);
	DenseVector<VariableArray1<number> > filPos;
	filPos.resize(2*nFil);

	for (size_t i = 0; i < nBnd; ++i)
	{
		vector3 edgeCenter = CalculateGridObjectCenter(m_tmpSpineEdges[i], m_aaPos);
		bndPos[2*i] = edgeCenter[0] / SPINE_NECK_RADIUS;
		bndPos[2*i+1] = edgeCenter[1] / SPINE_NECK_RADIUS;
	}

	try {distribute_filaments(bndPos, filPos);}
	UG_CATCH_THROW("Could not calculate filament positions.");

	for (size_t i = 0; i < nFil; ++i)
	{
		vector3 basePos;
		basePos.coord(0) = filPos[2*i] * (SPINE_NECK_RADIUS-MEMBRANE_RADIUS);
		basePos.coord(1) = filPos[2*i+1] * (SPINE_NECK_RADIUS-MEMBRANE_RADIUS);
		basePos.coord(2) = DENDRITE_RADIUS;
		m_tmpFilPos.push_back(vector2(basePos.coord(0), basePos.coord(1)));

//UG_LOGN("filament base position: " << basePos);

		// create circular cross-section for filaments at calculated positions
		// newly created elems are auto-selected
		create_circle(basePos, vector3(0,0,1), FILAMENT_WIDTH/2.0, FILAMENT_RIM_VERTICES);

		// assign subset
		m_sh.assign_subset(m_sel.begin<Vertex>(), m_sel.end<Vertex>(), SURF_CH_BND);
		m_sh.assign_subset(m_sel.begin<Edge>(), m_sel.end<Edge>(), SURF_CH_BND);

		// assign projection handler subset
		int proj_ssi_shaft = m_shProj.num_subsets();
		m_shProj.assign_subset(m_sel.begin<Vertex>(), m_sel.end<Vertex>(), proj_ssi_shaft);
		m_shProj.assign_subset(m_sel.begin<Edge>(), m_sel.end<Edge>(), proj_ssi_shaft);
		SmartPtr<CylinderProjector> shaftProj(new CylinderProjector(basePos, vector3(0.0, 0.0, 1.0), -1, -1));
		m_projHandler.set_projector(proj_ssi_shaft, shaftProj);

		// extrude head-wards
		std::vector<Vertex*> vrts;
		vrts.assign(m_sel.vertices_begin(), m_sel.vertices_end());
		std::vector<Edge*> edges;
		edges.assign(m_sel.edges_begin(), m_sel.edges_end());

		// retain base edges and vertices for later closing of filament
		std::vector<Vertex*> baseEndVrts = vrts;
		std::vector<Edge*> baseEndEdges = edges;

		// make extruded faces as regular as possible:
		// length of extrusion should be as near as possible to current edge length in circle
		number numExtrudes = floor(SPINE_NECK_LENGTH/(4*FILAMENT_WIDTH*sin(PI/FILAMENT_RIM_VERTICES)));
		size_t nExtrudes = std::max((size_t) numExtrudes, (size_t) 1);
		number extrude_length = SPINE_NECK_LENGTH / nExtrudes;

		vector3 extrudeDir(0.0);
		extrudeDir.coord(2) = extrude_length;

		for (size_t i = 0; i < nExtrudes; ++i)
			Extrude(m_grid, &vrts, &edges, NULL, extrudeDir, m_aaPos, EO_CREATE_FACES, NULL);

		// close both ends, upper one first
		int proj_ssi_upperEnd = m_shProj.num_subsets();
		m_shProj.assign_subset(vrts.begin(), vrts.end(), proj_ssi_upperEnd);
		m_shProj.assign_subset(edges.begin(), edges.end(), proj_ssi_upperEnd);
		SmartPtr<SphereProjector> upperEndProj(new SphereProjector(CalculateCenter(vrts.begin(), vrts.end(), m_aaPos), -1, -1));
		m_projHandler.set_projector(proj_ssi_upperEnd, upperEndProj);

		extrudeDir.coord(2) = 0.25*FILAMENT_WIDTH * sqrt(2.0);
		Extrude(m_grid, &vrts, &edges, NULL, extrudeDir, m_aaPos, EO_CREATE_FACES, NULL);

		vector3 center = CalculateCenter(vrts.begin(), vrts.end(), m_aaPos);
		size_t nVrt = vrts.size();
		for (size_t i = 0; i < nVrt; ++i)
		{
			vector3& v = m_aaPos[vrts[i]];
			VecSubtract(v, v, center);
			VecScale(v, v, 0.5*sqrt(2.0));
			VecAdd(v, v, center);
		}

		extrudeDir.coord(2) = 0.5*FILAMENT_WIDTH * (1.0 - 0.5*sqrt(2.0));
		Extrude(m_grid, &vrts, &edges, NULL, extrudeDir, m_aaPos, EO_CREATE_FACES, NULL);
		center = CalculateCenter(vrts.begin(), vrts.end(), m_aaPos);
		for (size_t i = 0; i < nVrt; ++i)
			m_aaPos[vrts[i]] = center;
		MergeMultipleVertices(m_grid, vrts.begin(), vrts.end());

		// close both ends, now lower one
		m_sel.clear();
		bool autoselEnabled = m_sel.autoselection_enabled();
		m_sel.enable_autoselection(true);

		int proj_ssi_lowerEnd = m_shProj.num_subsets();
		m_shProj.assign_subset(baseEndVrts.begin(), baseEndVrts.end(), proj_ssi_lowerEnd);
		m_shProj.assign_subset(baseEndEdges.begin(), baseEndEdges.end(), proj_ssi_lowerEnd);
		SmartPtr<SphereProjector> lowerEndProj(new SphereProjector(basePos, -1, -1));
		m_projHandler.set_projector(proj_ssi_lowerEnd, lowerEndProj);

		extrudeDir.coord(2) = -0.25*FILAMENT_WIDTH * sqrt(2.0);
		Extrude(m_grid, &baseEndVrts, &baseEndEdges, NULL, extrudeDir, m_aaPos, EO_CREATE_FACES, NULL);

		center = CalculateCenter(baseEndVrts.begin(), baseEndVrts.end(), m_aaPos);
		nVrt = baseEndVrts.size();
		for (size_t i = 0; i < nVrt; ++i)
		{
			vector3& v = m_aaPos[baseEndVrts[i]];
			VecSubtract(v, v, center);
			VecScale(v, v, 0.5*sqrt(2.0));
			VecAdd(v, v, center);
		}

		extrudeDir.coord(2) = -0.5*FILAMENT_WIDTH * (1.0 - 0.5*sqrt(2.0));
		Extrude(m_grid, &baseEndVrts, &baseEndEdges, NULL, extrudeDir, m_aaPos, EO_CREATE_FACES, NULL);
		center = CalculateCenter(baseEndVrts.begin(), baseEndVrts.end(), m_aaPos);
		for (size_t i = 0; i < nVrt; ++i)
			m_aaPos[baseEndVrts[i]] = center;
		MergeMultipleVertices(m_grid, baseEndVrts.begin(), baseEndVrts.end());

		// face orientation of lower closing faces needs to be inverted
		InvertOrientation(m_grid, m_sel.begin<Face>(), m_sel.end<Face>());

		m_sel.enable_autoselection(autoselEnabled);
	}
}


void MorphoGen::MinDistCalculator::calculate_minDist
(
	number& minDistOut,
	size_t& minIndOut,
	size_t ind,
	const std::vector<vector3>& vPoints,
	const std::vector<size_t>& vInd,
	size_t nPts
)
{
	const vector3& test_i = vPoints[ind];
	minDistOut = std::numeric_limits<number>::infinity();

	// distance to other mid points
	for (size_t j = 0; j < nPts; ++j)
	{
		if (vInd[j] == ind) continue;

		number dist = VecDistanceSq(test_i, vPoints[vInd[j]]);
		if (dist < minDistOut)
		{
			minDistOut = dist;
			minIndOut = vInd[j];
		}
	}

	// check boundary (boundary has no thickness, in contrast to the filaments)
	number dist;
	if (test_i.z() > neck_end_z)
	{
		// the situation is a slight bit complicated for the head:
		// if we are in the "opening cone", calculate distance to rim,
		// otherwise just take radius
		number theta = acos((head_center_z - test_i.z()) / VecDistance(test_i, vector3(0,0,head_center_z)));
		if (theta > head_opening_theta)
		{
			dist = head_radius - VecDistance(test_i, vector3(0,0,head_center_z)) + filament_width;
			dist = dist*dist;
		}
		else
		{
			number phi;
			if (fabs(test_i.x()) < fabs(1e-8*test_i.y()))
				phi = test_i.y() > 0 ? 0.5*PI : -0.5*PI;
			else
			{
				phi = atan(test_i.y() / test_i.x());
				if (test_i.x() < 0)
					phi += PI;
			}
			vector3 nearestRimPt(head_radius*sin(head_opening_theta)*cos(phi),
								 head_radius*sin(head_opening_theta)*sin(phi),
								 head_center_z - head_radius*cos(head_opening_theta));
			dist = VecDistanceSq(test_i, nearestRimPt);
		}
	}
	else
	{
		dist = neck_radius - sqrt(test_i.x()*test_i.x() + test_i.y()*test_i.y()) + filament_width;
		dist = dist*dist;
	}

	if (dist < minDistOut)
	{
		minDistOut = dist;
		minIndOut = nPts;  // index code for bnd
	}
}


void MorphoGen::make_spherical_filaments()
{
	// We will distribute a number of spheres representing anisotropic filament.
	// They are to be equally distributed in the head and the neck and as far apart
	// from both each other and the enclosing membrane.
	// This is achieved by generating a lot more spheres than necessary and then
	// iteratively removing one realizing the smallest distance to another or the membrane.

	// some geometrical parameters
	number neck_radius = SPINE_NECK_RADIUS - MEMBRANE_RADIUS;
	number head_radius = SPINE_HEAD_RADIUS - MEMBRANE_RADIUS;
	number neck_start_z = DENDRITE_RADIUS;
	number neck_end_z = DENDRITE_RADIUS + SPINE_NECK_LENGTH;
	number head_opening_theta = asin(neck_radius / head_radius);
	number dist_head_center_neck_end = head_radius * cos(head_opening_theta);
	number head_center_z = neck_end_z + dist_head_center_neck_end;
	number head_vol = 2.0/3.0*PI*(1+cos(head_opening_theta))*head_radius*head_radius*head_radius;
	number neck_vol = PI*neck_radius*neck_radius*SPINE_NECK_LENGTH;
	number head_vol_rel = head_vol / (head_vol + neck_vol);

	// distribute mid points
	size_t nTest = 10*NUM_FILAMENTS;
	std::vector<vector3> vTest;
	vTest.reserve(nTest);
	for (size_t i = 0; i < nTest; ++i)
	{
		number random = (number) std::rand() / RAND_MAX;

		// head
		if (random < head_vol_rel)
		{
			number r, theta, phi;
			bool allowed = false;
			while (!allowed)
			{
				r = (number) std::rand() / RAND_MAX;
				r = r*r*r; // distribution density must be quadratic in r
				r = head_radius * (1-r);

				theta = PI * (number) std::rand() / RAND_MAX;
				if (theta > head_opening_theta || r*cos(theta) < dist_head_center_neck_end)
					allowed = true;
			}

			phi = 2.0*PI * (number) std::rand() / RAND_MAX;

			vTest.push_back(vector3(r*sin(theta)*cos(phi),
				                    r*sin(theta)*sin(phi),
									head_center_z - r*cos(theta)));
		}

		// neck
		else
		{
			number r, phi, z;
			r = (number) std::rand() / RAND_MAX;
			r = r*r; // distribution density must be linear in r
			r = neck_radius * (1-r);
			phi = 2.0*PI * (number) std::rand() / RAND_MAX;
			z = SPINE_NECK_LENGTH * (number) std::rand() / RAND_MAX;

			vTest.push_back(vector3(r*cos(phi),
				                    r*sin(phi),
									neck_start_z + z));
		}
	}


	// find minimal distances
	std::vector<number> vMinDist(nTest, std::numeric_limits<number>::infinity());
	std::vector<size_t> vNN(nTest, 0);
	std::vector<size_t> vInd(nTest);
	for (size_t i = 0; i < nTest; ++i) vInd[i] = i;

	MinDistCalculator mdc;
	mdc.neck_radius = neck_radius;
	mdc.neck_end_z = neck_end_z;
	mdc.head_radius = head_radius;
	mdc.head_center_z = head_center_z;
	mdc.head_opening_theta = head_opening_theta;
	mdc.filament_width = FILAMENT_WIDTH;
	for (size_t i = 0; i < nTest; ++i)
		mdc.calculate_minDist(vMinDist[i], vNN[i], i, vTest, vInd, nTest);


	// iteratively take away points verifying current minDist
	IndexCompare cmp(vMinDist);
	std::vector<bool> vStillExists(nTest, true);

	for (size_t i = 0; i < nTest - NUM_FILAMENTS; ++i)
	{
		// sort according to minDist
		std::sort(vInd.begin(), vInd.end()-i, cmp);

		size_t minInd = vInd[0];
		UG_COND_THROW(!vStillExists[minInd], "Something went wrong during distribution of filaments.");

		// remove point
		vStillExists[minInd] = false;
		vMinDist[minInd] = std::numeric_limits<number>::infinity();
		vNN[minInd] = -1;
		std::swap(vInd[0], vInd[nTest-1-i]);

		// update minDists
		for (size_t j = 0; j < nTest-i-1; ++j)
			if (vNN[vInd[j]] == minInd)
				mdc.calculate_minDist(vMinDist[vInd[j]], vNN[vInd[j]], vInd[j], vTest, vInd, nTest-i-1);
	}

	// check sizes
	if (NUM_FILAMENTS)
	{
		vInd.resize(NUM_FILAMENTS);
		std::sort(vInd.begin(), vInd.end(), cmp);
		UG_COND_THROW(sqrt(vMinDist[vInd[0]]) < 2*FILAMENT_WIDTH,
			"Filaments could not be positioned far enough apart to exclude intersections.");
		UG_COND_THROW(sqrt(vMinDist[vInd[0]]) < 2*FILAMENT_WIDTH+FILAMENT_ENVELOPE_RADIUS+MEMBRANE_ENVELOPE_RADIUS,
			"Filaments could not be positioned far enough apart to exclude intersections of their envelopes.");
	}


	// create spheres
	bool autoselectionEnabled = m_sel.autoselection_enabled();
	m_sel.enable_autoselection(true);

	for (size_t i = 0; i < NUM_FILAMENTS; ++i)
	{
		//UG_LOGN("Filament sphere center " << i << " at " << vTest[vInd[i]] << ".");

		// deselect currently selected elements
		m_sel.deselect(m_sel.vertices_begin(), m_sel.vertices_end());
		m_sel.deselect(m_sel.edges_begin(), m_sel.edges_end());
		m_sel.deselect(m_sel.faces_begin(), m_sel.faces_end());

		// create the coarse ico
		GenerateIcosahedron(m_grid, vTest[vInd[i]], FILAMENT_WIDTH, aPosition);

		// assign subset
		m_sh.assign_subset(m_sel.begin<Vertex>(), m_sel.end<Vertex>(), SURF_CH_BND);
		m_sh.assign_subset(m_sel.begin<Edge>(), m_sel.end<Edge>(), SURF_CH_BND);
		m_sh.assign_subset(m_sel.begin<Face>(), m_sel.end<Face>(), SURF_CH_BND);

		// set sphere projector
		int proj_ssi = m_shProj.num_subsets();
		m_shProj.assign_subset(m_sel.vertices_begin(), m_sel.vertices_end(), proj_ssi);
		m_shProj.assign_subset(m_sel.edges_begin(), m_sel.edges_end(), proj_ssi);
		m_shProj.assign_subset(m_sel.faces_begin(), m_sel.faces_end(), proj_ssi);
		SmartPtr<SphereProjector> lowerEndProj(new SphereProjector(vTest[vInd[i]], -1, -1));
		m_projHandler.set_projector(proj_ssi, lowerEndProj);
	}

	m_sel.enable_autoselection(autoselectionEnabled);
}



template <typename TIterator>
void MorphoGen::replace_subsets
(
	const TIterator& begin,
	const TIterator& end,
	const std::vector<int>& oldSI,
	const std::vector<int>& newSI,
	bool deselect
)
{
	size_t sz = oldSI.size();
	UG_COND_THROW(sz != newSI.size(), "Old and new subset index vectors do not have the same sizes.");

	for (TIterator it = begin; it != end; ++it)
	{
		for (size_t k = 0; k < sz; ++k)
		{
			if (m_sh.get_subset_index(*it) == oldSI[k])
			{
				m_sh.assign_subset(*it, newSI[k]);
				if (deselect) m_sel.deselect(*it);
				break;
			}
		}
	}
}



void MorphoGen::create_envelope(const std::vector<int>& vExtrudeSI, number offset, int newVolSI,  const std::vector<int>& vNewFrontSI)
{
	// have all newly created elements be selected
	m_sel.clear();
	bool autoselEnabled = m_sel.autoselection_enabled();
	m_sel.enable_autoselection(true);

	// extrude by zero
	std::vector<Vertex*> vrts;
	std::vector<Edge*> edges;
	std::vector<Face*> faces;

	for (size_t i = 0; i < vExtrudeSI.size(); ++i)
	{
		vrts.insert(vrts.end(), m_sh.begin<Vertex>(vExtrudeSI[i]), m_sh.end<Vertex>(vExtrudeSI[i]));
		edges.insert(edges.end(), m_sh.begin<Edge>(vExtrudeSI[i]), m_sh.end<Edge>(vExtrudeSI[i]));
		faces.insert(faces.end(), m_sh.begin<Face>(vExtrudeSI[i]), m_sh.end<Face>(vExtrudeSI[i]));
	}
	Extrude(m_grid, &vrts, &edges, &faces, vector3(0,0,0), m_aaPos, EO_CREATE_FACES | EO_CREATE_VOLUMES, NULL);

	// fix face orientation (might be wrong)
	FixFaceOrientation(m_grid, faces.begin(), faces.end());

	// normal move
	CalculateFaceNormals(m_grid, m_sel.begin<Face>(), m_sel.end<Face>(), aPosition, aNormal);
	for (std::vector<Vertex*>::iterator iter = vrts.begin(); iter != vrts.end(); ++iter)
	{
		vector3& pos = m_aaPos[*iter];

		// calculate vertex normal by averaging face normals
		vector3 normal(0.0, 0.0, 0.0);
		Grid::face_traits::secure_container assFaces;
		m_grid.associated_elements(assFaces, *iter);
		size_t afsz = assFaces.size();
		for (size_t i = 0; i < afsz; ++i)
		{
			VecAdd(normal, normal, m_aaNorm[assFaces[i]]);
		}
		VecNormalize(normal, normal);

		VecScaleAdd(pos, 1.0, pos, offset, normal);
	}

	m_sel.enable_autoselection(autoselEnabled);

	// adjust subset index for new front
	replace_subsets(vrts.begin(), vrts.end(), vExtrudeSI, vNewFrontSI, true);
	replace_subsets(edges.begin(), edges.end(), vExtrudeSI, vNewFrontSI, true);
	replace_subsets(faces.begin(), faces.end(), vExtrudeSI, vNewFrontSI, true);

	// adjust subset index for created volumes (and associated)
	std::vector<int> vNewVolSI(vExtrudeSI.size(), newVolSI);
	replace_subsets(m_sel.begin<Vertex>(), m_sel.end<Vertex>(), vExtrudeSI, vNewVolSI);
	replace_subsets(m_sel.begin<Edge>(), m_sel.end<Edge>(), vExtrudeSI, vNewVolSI);
	replace_subsets(m_sel.begin<Face>(), m_sel.end<Face>(), vExtrudeSI, vNewVolSI);
	replace_subsets(m_sel.begin<Volume>(), m_sel.end<Volume>(), vExtrudeSI, vNewVolSI);
}


void MorphoGen::create_membrane_and_envelopes()
{
	// first of all, check that distances between filaments and spine membrane are sufficient
	// also check that distances between filaments are sufficient
	// this is only necessary for the anisotropic case
	if (m_bFilAnisotropic)
	{
		size_t nFil = m_tmpFilPos.size();

		number safety = FILAMENT_WIDTH;
		for (size_t i = 0; i < nFil; ++i)
		{
			number r = VecLength(m_tmpFilPos[i]);

			if (r + FILAMENT_WIDTH/2.0 + FILAMENT_ENVELOPE_RADIUS + safety
				> SPINE_NECK_RADIUS - MEMBRANE_RADIUS - MEMBRANE_ENVELOPE_RADIUS)
			{
				UG_THROW("Envelopped filaments will (almost) intersect with envelopped membrane.\n"
					"Choose smaller envelope widths, fewer filaments or a bigger spine radius.");
			}

			number minDist = 2*SPINE_NECK_RADIUS;
			for (size_t j = i+1; j < nFil; ++j)
				minDist = std::min(minDist, VecDistance(m_tmpFilPos[i], m_tmpFilPos[j]));

			if (minDist - FILAMENT_WIDTH - 2*FILAMENT_ENVELOPE_RADIUS - safety < 0)
			{
				UG_THROW("Envelopped filaments will (almost) intersect with each other.\n"
					"Choose smaller envelope widths, fewer filaments or a bigger spine radius.");
			}
		}
	}

	std::vector<int> vExtrudeSI;
	std::vector<int> vNewFrontSI;

	// create envelope around filaments
	vExtrudeSI.push_back(SURF_CH_BND);
	vNewFrontSI.push_back(INNER_FRONT_TMP_SI);
	create_envelope(vExtrudeSI, FILAMENT_ENVELOPE_RADIUS, INNER_SI, vNewFrontSI);

	// create envelope around outer membrane surface
	vExtrudeSI[0] = MEM_OUTER_BND_SI;
	vExtrudeSI.push_back(PSD_OUTER_SI);
	vNewFrontSI[0] = OUTER_FRONT_TMP_SI;
	vNewFrontSI.push_back(OUTER_FRONT_TMP_SI);
	create_envelope(vExtrudeSI, MEMBRANE_ENVELOPE_RADIUS, OUTER_SI, vNewFrontSI);

	// create membrane
	vNewFrontSI[0] = MEM_INNER_BND_SI;
	vNewFrontSI[1]= PSD_INNER_SI;
	create_envelope(vExtrudeSI, -MEMBRANE_RADIUS, MEM_SI, vNewFrontSI);

	// create envelope around inner membrane surface
	vExtrudeSI = vNewFrontSI;
	vNewFrontSI[0] = INNER_FRONT_TMP_SI;
	vNewFrontSI[1] = INNER_FRONT_TMP_SI;
	create_envelope(vExtrudeSI, -MEMBRANE_ENVELOPE_RADIUS, INNER_SI, vNewFrontSI);
}


void MorphoGen::find_ring_at_dendritic_ends
(
	int si,
	int sideMultiplier,
	std::vector<Vertex*>& vrts,
	std::vector<Edge*>& edges
)
{
	typedef Grid::traits<Edge>::secure_container edge_list;

	// find left/right-most vertex in inner subset
	number optX = 0;
	Vertex* optVrt = NULL;
	VertexIterator it = m_sh.begin<Vertex>(si);
	VertexIterator it_end = m_sh.end<Vertex>(si);
	for (; it != it_end; ++it)
	{
		if (m_aaPos[*it].x()*sideMultiplier < optX*sideMultiplier)
		{
			optX = m_aaPos[*it].x();
			optVrt = *it;
		}
	}
	UG_COND_THROW(!optVrt, "No start vertex found for closing of dendritic end.");

	vector2 minVec(m_aaPos[optVrt].y(), m_aaPos[optVrt].z());
	number rSq = VecLengthSq(minVec);

	// fill vectors with all vertices and edges in left/right rim (counter-clockwise)
	vrts.push_back(optVrt);
	Vertex* currVrt = optVrt;
	Vertex* nextVrt = NULL;
	edge_list el;
	while (true)
	{
		m_grid.associated_elements(el, currVrt);
		size_t szel = el.size();
		size_t e = 0;
		for (; e < szel; ++e)
		{
			nextVrt = GetOpposingSide(m_grid, el[e], currVrt);
			vector2 nextVec(m_aaPos[nextVrt].y(), m_aaPos[nextVrt].z());
			number orient = optX*(m_aaPos[currVrt].y()*nextVec[1] - m_aaPos[currVrt].z()*nextVec[0]);
			if (fabs(VecLengthSq(nextVec) - rSq) < 1e-2 * rSq				// correct radius
				&& fabs(m_aaPos[nextVrt].x() - optX) < 1e-2 * fabs(optX)	// correct x coord
				&&  orient > 0) 											// correct orientation
			{
				break;
			}
		}
		UG_COND_THROW(e == szel, "No closed circle of delimiting edges found for closing of dendritic end.");

		edges.push_back(el[e]);
		currVrt = nextVrt;
		if (currVrt == optVrt) break;
		vrts.push_back(currVrt);
	}
}


void MorphoGen::close_ends()
{
	typedef Grid::traits<Edge>::secure_container edge_list;

	for (size_t i = 0; i < 2; ++i)
	{
		int mult = 1;
		if (i == 1) mult = -1;

		std::vector<Vertex*> vrts;
		std::vector<Edge*> edges;
		find_ring_at_dendritic_ends(INNER_FRONT_TMP_SI, mult, vrts, edges);

		EdgeSelector edgeSel(m_grid);
		edgeSel.select(edges.begin(), edges.end());

		// create faces using center vertex
		m_sel.clear();
		m_sh.set_default_subset_index(INNER_FRONT_TMP_SI);
		vector3 center = CalculateBarycenter(vrts.begin(), vrts.end(), m_aaPos);
		Vertex* cv = *m_grid.create<RegularVertex>();
		m_aaPos[cv] = center;

		size_t nVrt = vrts.size();
		for (size_t v = 0; v < nVrt; ++v)
		{
			Face* f = *m_grid.create<Triangle>(TriangleDescriptor(vrts[v], vrts[(v+1)%nVrt], cv));
			m_sel.select(f);
		}

		// enhance quality
		QualityGridGeneration(m_grid, m_sel.begin<Face>(), m_sel.end<Face>(), m_aaPos, 30.0, IsSelected(edgeSel));
	}
}


void MorphoGen::create_mark_straight_edge_connection(Vertex* v0, Vertex* v1)
{
	bool autoselEnabled = m_sel.autoselection_enabled();
	m_sel.enable_autoselection(true);

	const vector3& vec0 = m_aaPos[v0];
	const vector3& vec1 = m_aaPos[v1];

	number length = VecDistance(vec0, vec1);
	number numExtrudes = floor(length / (4*DENDRITE_RADIUS*sin(PI/DENDRITE_RIM_VERTICES)));
	size_t nExtrudes = std::max((size_t) numExtrudes, (size_t) 1);

	std::vector<Vertex*> vrts(1, v0);
	vector3 dir;
	VecSubtract(dir, vec1, vec0);
	VecScale(dir, dir, 1.0/nExtrudes);
	for (size_t j = 0; j < nExtrudes; ++j)
		Extrude(m_grid, &vrts, NULL, NULL, dir, m_aaPos, 0, NULL);

	MergeVertices(m_grid, v1, vrts[0]);

	m_sel.select(v0);
	m_sel.select(v1);

	m_sel.enable_autoselection(autoselEnabled);
}


void MorphoGen::mark_straight_edge_connection(Vertex* v0, Vertex* v1)
{
	typedef Grid::traits<Edge>::secure_container edge_list;

	vector3 dir;
	const vector3& vec0 = m_aaPos[v0];
	const vector3& vec1 = m_aaPos[v1];
	VecSubtract(dir, vec1, vec0);
	VecNormalize(dir, dir);
	m_sel.select(v0);
	Vertex* curr = v0;
	Vertex* next = NULL;
	vector3 edge_center;

	while (curr != v1)
	{
		edge_list el;
		m_grid.associated_elements(el, curr);
		size_t elsz = el.size();
		size_t k = 0;
		for (; k < elsz; ++k)
		{
			Edge* e = el[k];
			next = GetOpposingSide(m_grid, e, curr);
			vector3 edge_dir;
			VecSubtract(edge_dir, m_aaPos[next], m_aaPos[curr]);
			VecNormalize(edge_dir, edge_dir);
			if (VecProd(edge_dir, dir) < 0.999)
				continue;

			curr = next;
			m_sel.select(e);
			m_sel.select(curr);
			break;
		}
		UG_COND_THROW(k == elsz, "Failed to find straight edge connection.");
	}
}


void MorphoGen::create_bounding_box()
{
	// check margin is big enough
	UG_COND_THROW(BOX_MARGIN < 2 * MEMBRANE_ENVELOPE_RADIUS,
		"Box margin is not big enough, should be at least " <<
		MEMBRANE_ENVELOPE_RADIUS << ".");

	// calculate bounding box
	number z_low = -(DENDRITE_RADIUS + MEMBRANE_ENVELOPE_RADIUS + BOX_MARGIN);
	number z_high = (DENDRITE_RADIUS + SPINE_NECK_LENGTH + m_tmpHeadHeight
					 + MEMBRANE_ENVELOPE_RADIUS + BOX_MARGIN);
	number y_low = z_low;
	number y_high = -y_low;
	number x_low = -0.5*DENDRITE_LENGTH;

	// enable auto_selection
	bool autoselEnabled = m_sel.autoselection_enabled();
	m_sel.enable_autoselection(true);

	std::vector<Vertex*> vCorner[2];

	number numExtrudes = floor(y_high / (2*DENDRITE_RADIUS*sin(PI/DENDRITE_RIM_VERTICES)));
	size_t nExtrudes = std::max((size_t) numExtrudes, (size_t) 1);
	number extrude_length = 2*y_high / nExtrudes;
	numExtrudes = floor((z_high-y_low) / (4*DENDRITE_RADIUS*sin(PI/DENDRITE_RIM_VERTICES)));
	size_t nExtrudesUp = std::max((size_t) numExtrudes, (size_t) 1);
	number extrude_length_up = (z_high - y_low) / nExtrudesUp;

	// first, triangulate partial areas around dendritic ends
	// 0 encodes left; 1 encodes right
	AInt aInt;
	m_grid.attach_to_vertices(aInt);
	m_sh.set_default_subset_index(NOFLUX_BND);
	for (size_t i = 0; i < 2; ++i)
	{
		m_sel.clear();

		int mult = 1;
		if (i == 1) mult = -1;

		// create start vertices left and right
		Vertex* btm_lft = *m_grid.create<RegularVertex>();
		Vertex* top_lft = *m_grid.create<RegularVertex>();

		m_aaPos[btm_lft] = vector3(mult*x_low, mult*y_low, z_low);
		m_aaPos[top_lft] = vector3(mult*x_low, mult*y_low, z_high);

		// extrude verts in y-direction
		std::vector<Vertex*> vrts(2);
		vrts[0] = btm_lft;
		vrts[1] = top_lft;

		vector3 dir(0.0, mult*extrude_length, 0.0);

		for (size_t j = 0; j < nExtrudes; ++j)
			Extrude(m_grid, &vrts, NULL, NULL, dir, m_aaPos, 0, NULL);

		Vertex* btm_rgt = vrts[0];
		Vertex* top_rgt = vrts[1];

		// extrude verts in z-direction
		vrts[1] = btm_lft;
		dir.y() = 0.0; dir.z() = extrude_length_up;

		for (size_t j = 0; j < nExtrudesUp; ++j)
			Extrude(m_grid, &vrts, NULL, NULL, dir, m_aaPos, 0, NULL);

		// remove doubles
		MergeVertices(m_grid, top_lft, vrts[1]);
		MergeVertices(m_grid, top_rgt, vrts[0]);

		// save corners for later use
		if (i == 0)
		{
			vCorner[i].push_back(btm_rgt);
			vCorner[i].push_back(top_rgt);
			vCorner[i].push_back(top_lft);
			vCorner[i].push_back(btm_lft);
		}
		else
		{
			vCorner[i].push_back(btm_lft);
			vCorner[i].push_back(top_lft);
			vCorner[i].push_back(top_rgt);
			vCorner[i].push_back(btm_rgt);
		}

		// select outer envelope ring
		vrts.clear();
		std::vector<Edge*> edges;
		try {find_ring_at_dendritic_ends(OUTER_FRONT_TMP_SI, mult, vrts, edges);}
		UG_CATCH_THROW("Could not find bounding envelope ring at dendritic end.");

		m_sel.select(vrts.begin(), vrts.end());
		m_sel.select(edges.begin(), edges.end());

		// triangulate
		m_sel.enable_autoselection(false);

		FaceSelector faceSel(m_grid);
		faceSel.enable_autoselection(true);

		UG_COND_THROW(!TriangleFill_SweepLine(m_grid, m_sel.edges_begin(), m_sel.edges_end(),
											  aPosition, aInt, &m_sh, NOFLUX_BND),
			"TriangleFill_SweepLine failed in bounding box creation of dendritic sides.");

		// remove faces created inside dendrite
		Selector eraseSel(m_grid);
		for (FaceIterator fit = faceSel.begin<Face>(); fit != faceSel.end<Face>(); ++fit)
		{
			vector3 c = CalculateCenter(*fit, m_aaPos);
			if (sqrt(c.y()*c.y() + c.z()*c.z()) < DENDRITE_RADIUS + MEMBRANE_ENVELOPE_RADIUS)
				eraseSel.select(*fit);
		}
		SelectInnerSelectionEdges(eraseSel);
		m_grid.erase(eraseSel.begin<Face>(), eraseSel.end<Face>());
		m_grid.erase(eraseSel.begin<Edge>(), eraseSel.end<Edge>());

		// better quality
		QualityGridGeneration(m_grid, faceSel.begin<Face>(), faceSel.end<Face>(), m_aaPos,
			30.0, IsSelected(m_sel));

		// invert face orientation for left side
		if (i == 0)
			InvertOrientation(m_grid, faceSel.begin<Face>(), faceSel.end<Face>());

		// correction for envelope subset
		// This subset has been wrongly assigned by envelope extrusion.
		// We cannot simply assign NOFLUX_BND, as this will interfere with tetrahedralization
		// of extracellular space. Instead, we assign a temporary subset
		// which we will change in the end.
		VertexIterator vit = m_sel.begin<Vertex>();
		VertexIterator vit_end = m_sel.end<Vertex>();
		for (; vit != vit_end; ++vit)
		{
			if (m_sh.get_subset_index(*vit) != OUTER_FRONT_TMP_SI)
				continue;

			m_sh.assign_subset(*vit, NOFLUX_BND);

			Grid::AssociatedEdgeIterator edgeIt = m_grid.associated_edges_begin(*vit);
			Grid::AssociatedEdgeIterator edge_end = m_grid.associated_edges_end(*vit);
			for (; edgeIt != edge_end; ++edgeIt)
			{
				vector3 ec = CalculateCenter(*edgeIt, m_aaPos);
				if (fabs(ec.x() - mult*x_low) < 1e-4*fabs(x_low))
				{
					int oldSi = m_sh.get_subset_index(*edgeIt);
					if (oldSi == OUTER_FRONT_TMP_SI || oldSi == NOFLUX_BND)
						m_sh.assign_subset(*edgeIt, NOFLUX_BND);
					else
						m_sh.assign_subset(*edgeIt, ENVELOPE_END_TMP_SI);

				}
			}

			Grid::AssociatedFaceIterator faceIt = m_grid.associated_faces_begin(*vit);
			Grid::AssociatedFaceIterator face_end = m_grid.associated_faces_end(*vit);
			for (; faceIt != face_end; ++faceIt)
			{
				vector3 fc = CalculateCenter(*faceIt, m_aaPos);
				if (fabs(fc.x() - mult*x_low) < 1e-4*fabs(x_low))
				{
					if (m_sh.get_subset_index(*faceIt) != NOFLUX_BND)
						m_sh.assign_subset(*faceIt, ENVELOPE_END_TMP_SI);
				}
			}
		}

		m_sel.enable_autoselection(true);
	}

	// now, go around the x-axis and triangulate the rest
	std::vector<Vertex*> vrts(1);
	vector3 dir;
	m_sh.set_default_subset_index(DIRI_BND);
	for (size_t i = 0; i < 4; ++i)
	{
		m_sel.clear();

		// create/mark edges in x-dir for first face side
		if (i == 0)
		{
			try {create_mark_straight_edge_connection(vCorner[0][i], vCorner[1][i]);}
			UG_CATCH_THROW("Failed to create-mark straight edge connection.");
		}
		else
		{
			try {mark_straight_edge_connection(vCorner[0][i], vCorner[1][i]);}
			UG_CATCH_THROW("Failed to mark straight edge connection.");
		}

		// create edges in x-dir for second face side
		if (i != 3)
		{
			try {create_mark_straight_edge_connection(vCorner[1][i+1], vCorner[0][i+1]);}
			UG_CATCH_THROW("Failed to create-mark straight edge connection.");

		}
		else
		{
			try {mark_straight_edge_connection(vCorner[1][0], vCorner[0][0]);}
			UG_CATCH_THROW("Failed to mark straight edge connection.");
		}

		// mark the two other edges
		try
		{
			mark_straight_edge_connection(vCorner[1][i], vCorner[1][(i+1)%4]);
			mark_straight_edge_connection(vCorner[0][(i+1)%4], vCorner[0][i]);
		}
		UG_CATCH_THROW("Failed to mark straight edge connection.");

		// first, assign subset
		m_sh.assign_subset(m_sel.begin<Vertex>(), m_sel.end<Vertex>(), DIRI_BND);
		m_sh.assign_subset(m_sel.begin<Edge>(), m_sel.end<Edge>(), DIRI_BND);

		// then triangulate
		m_sel.enable_autoselection(false);

		FaceSelector faceSel(m_grid);
		faceSel.enable_autoselection(true);

		UG_COND_THROW(!TriangleFill_SweepLine(m_grid, m_sel.edges_begin(), m_sel.edges_end(),
											  aPosition, aInt, &m_sh, DIRI_BND),
				"TriangleFill_SweepLine failed in bounding box creation of side " << i << ".");

		// better quality
		QualityGridGeneration(m_grid, faceSel.begin<Face>(), faceSel.end<Face>(), m_aaPos,
			30.0, IsSelected(m_sel));

		// invert bottom and front side
		if (i >= 2)
			InvertOrientation(m_grid, faceSel.begin<Face>(), faceSel.end<Face>());

		m_sel.enable_autoselection(true);
	}

	m_sel.enable_autoselection(autoselEnabled);
}


void MorphoGen::tetrahedralize_selection(int tetSI, const std::vector<vector3>* pvHoles)
{
#ifdef TETGEN_15_ENABLED
	// attach an index to the vertices
	AInt aInd;
	m_grid.attach_to_vertices(aInd);
	Grid::VertexAttachmentAccessor<AInt> aaInd(m_grid, aInd);

	// data structures used to communicate with tetgen
	tetgenio in, out;

	// set up points
	in.numberofpoints = m_sel.num<Vertex>();
	in.pointlist = new REAL[in.numberofpoints*3];

	int counter = 0;
	VertexIterator vit = m_sel.begin<Vertex>();
	VertexIterator vit_end = m_sel.end<Vertex>();
	for (; vit != vit_end; ++vit, ++counter)
	{
		aaInd[*vit] = counter;
		vector3& v = m_aaPos[*vit];
		in.pointlist[counter*3] = (REAL) v.x();
		in.pointlist[counter*3+1] = (REAL) v.y();
		in.pointlist[counter*3+2] = (REAL) v.z();
	}

	// set up facets
	counter = 0;
	in.numberoffacets = m_sel.num<Face>();
	in.facetlist = new tetgenio::facet[in.numberoffacets];
	in.facetmarkerlist = new int[in.numberoffacets];

	FaceIterator fit = m_sel.begin<Face>();
	FaceIterator fit_end = m_sel.end<Face>();
	for (; fit != fit_end; ++fit, ++counter)
	{
		Face* f = *fit;

		// copy the face to the facetlist
		tetgenio::facet* tf = &in.facetlist[counter];
		tf->numberofpolygons = 1;
		tf->polygonlist = new tetgenio::polygon[tf->numberofpolygons];
		tf->numberofholes = 0;
		tf->holelist = NULL;
		tetgenio::polygon* p = &tf->polygonlist[0];
		p->numberofvertices = f->num_vertices();
		p->vertexlist = new int[p->numberofvertices];
		for (int i = 0; i < p->numberofvertices; ++i)
			p->vertexlist[i] = aaInd[f->vertex(i)];

		// set the face mark
		in.facetmarkerlist[counter] = m_sh.get_subset_index(f);
	}

	// the aInd attachment is no longer required
	m_grid.detach_from_vertices(aInd);
	aaInd.invalidate();

	// set up holes
	if (pvHoles && pvHoles->size())
	{
		int nh = in.numberofholes = pvHoles->size();
		in.holelist = new REAL[in.numberofholes*3];
		for (int i = 0; i < nh; ++i)
		{
			const vector3& hv = (*pvHoles)[i];
			in.holelist[3*i] = (REAL) hv.x();
			in.holelist[3*i+1] = (REAL) hv.y();
			in.holelist[3*i+2] = (REAL) hv.z();
		}
	}

	// call tetgen
	std::stringstream ss;
	ss << "p"; // "piecewise linear complex" (polygonal mesh)
	ss << "Y"; // preservation of bounds
	ss << "q," << TETRAHEDRALIZATION_QUALITY; // quality setting
	ss << "Q"; // quiet
	try {::tetrahedralize(const_cast<char*>(ss.str().c_str()), &in, &out);}
	catch (int errCode)
	{UG_THROW("TetGen tetrahedralization failed with error code " << errCode);}
	UG_CATCH_THROW("TetGen tetrahedralization failed.");

	// add new vertices to the grid, store all vertices in a vector
	std::vector<Vertex*> vVrts(out.numberofpoints);
	counter = 0;

	// add the old ones to the vector
	vit = m_sel.begin<Vertex>();
	for (; vit != vit_end && counter < out.numberofpoints; ++vit, ++counter)
	{
		m_aaPos[*vit].x() = out.pointlist[counter*3];
		m_aaPos[*vit].y() = out.pointlist[counter*3+1];
		m_aaPos[*vit].z() = out.pointlist[counter*3+2];
		vVrts[counter] = *vit;
	}
	UG_COND_THROW(vit != vit_end, "Not all vertices of the input are present in the output.");

	// create new ones and add them to the vector
	for (; counter < out.numberofpoints; ++counter)
	{
		RegularVertex* v = *m_grid.create<RegularVertex>();
		m_aaPos[v].x() = out.pointlist[counter*3];
		m_aaPos[v].y() = out.pointlist[counter*3+1];
		m_aaPos[v].z() = out.pointlist[counter*3+2];
		vVrts[counter] = v;
	}

	// add new volumes
	UG_COND_THROW(out.numberoftetrahedra < 1, "Not a single tetrahedron was created.");
	for (int i = 0; i < out.numberoftetrahedra; ++i)
	{
		Tetrahedron* tet = *m_grid.create<Tetrahedron>(
			TetrahedronDescriptor(vVrts[out.tetrahedronlist[i*4]],
								  vVrts[out.tetrahedronlist[i*4+1]],
								  vVrts[out.tetrahedronlist[i*4+2]],
								  vVrts[out.tetrahedronlist[i*4+3]]));
		m_sh.assign_subset(tet, tetSI);
	}

	// freeing dynamically allocated memory is done within tetgen!
#else
	TETRAHEDRALIZATION_QUALITY += 0;	// suppress -Wunused-private-field
	UG_THROW("TetGen is not available in the current build.\n"
			 "Recompile with TetGen 1.5 support to use tetrahedralization.");

#endif
}


void MorphoGen::tetrahedralize()
{
	// tetrahedralize the filaments
	m_sel.clear();
	m_sel.select(m_sh.begin<Face>(SURF_CH_BND), m_sh.end<Face>(SURF_CH_BND));
	SelectAssociatedEdges(m_sel, m_sel.begin<Face>(), m_sel.end<Face>());
	SelectAssociatedVertices(m_sel, m_sel.begin<Edge>(), m_sel.end<Edge>());
	m_sh.set_default_subset_index(FIL_NECK_SI);
	try {tetrahedralize_selection(FIL_NECK_SI);}
	UG_CATCH_THROW("Tetrahedralization of filaments failed.");

	// tetrahedralize the cytosol
	m_sel.clear();
	m_sel.select(m_sh.begin<Face>(INNER_FRONT_TMP_SI), m_sh.end<Face>(INNER_FRONT_TMP_SI));
	SelectAssociatedEdges(m_sel, m_sel.begin<Face>(), m_sel.end<Face>());
	SelectAssociatedVertices(m_sel, m_sel.begin<Edge>(), m_sel.end<Edge>());

	size_t nFil = m_tmpFilPos.size();
	std::vector<vector3> holes(nFil);
	for (size_t i = 0; i < nFil; ++i)
	{
		holes[i].x() = m_tmpFilPos[i].x();
		holes[i].y() = m_tmpFilPos[i].y();
		holes[i].z() = DENDRITE_RADIUS + SPINE_NECK_LENGTH / 2.0;
	}

	m_sh.set_default_subset_index(INNER_SI);
	try {tetrahedralize_selection(INNER_SI, &holes);}
	UG_CATCH_THROW("Tetrahedralization of intracellular space failed.");

	// change temporal inner front subset to inner
	m_sel.clear();
	m_sel.select(m_sh.begin<Vertex>(INNER_FRONT_TMP_SI), m_sh.end<Vertex>(INNER_FRONT_TMP_SI));
	m_sel.select(m_sh.begin<Edge>(INNER_FRONT_TMP_SI), m_sh.end<Edge>(INNER_FRONT_TMP_SI));
	m_sel.select(m_sh.begin<Face>(INNER_FRONT_TMP_SI), m_sh.end<Face>(INNER_FRONT_TMP_SI));
	m_sh.assign_subset(m_sel.begin<Vertex>(), m_sel.end<Vertex>(), INNER_SI);
	m_sh.assign_subset(m_sel.begin<Edge>(), m_sel.end<Edge>(), INNER_SI);
	m_sh.assign_subset(m_sel.begin<Face>(), m_sel.end<Face>(), INNER_SI);

	// tetrahedralize extracellular space
	m_sel.clear();
	m_sel.select(m_sh.begin<Face>(OUTER_FRONT_TMP_SI), m_sh.end<Face>(OUTER_FRONT_TMP_SI));
	m_sel.select(m_sh.begin<Face>(DIRI_BND), m_sh.end<Face>(DIRI_BND));
	m_sel.select(m_sh.begin<Face>(NOFLUX_BND), m_sh.end<Face>(NOFLUX_BND));
	SelectAssociatedEdges(m_sel, m_sel.begin<Face>(), m_sel.end<Face>());
	SelectAssociatedVertices(m_sel, m_sel.begin<Edge>(), m_sel.end<Edge>());

	m_sh.set_default_subset_index(OUTER_SI);
	holes.resize(1);
	holes[0] = vector3(0.0, 0.0, 0.0);
	try {tetrahedralize_selection(OUTER_SI, &holes);}
	UG_CATCH_THROW("Tetrahedralization of extracellular space failed.");

	// change temporal outer front subset to outer
	m_sel.clear();
	m_sel.select(m_sh.begin<Vertex>(OUTER_FRONT_TMP_SI), m_sh.end<Vertex>(OUTER_FRONT_TMP_SI));
	m_sel.select(m_sh.begin<Edge>(OUTER_FRONT_TMP_SI), m_sh.end<Edge>(OUTER_FRONT_TMP_SI));
	m_sel.select(m_sh.begin<Face>(OUTER_FRONT_TMP_SI), m_sh.end<Face>(OUTER_FRONT_TMP_SI));
	m_sh.assign_subset(m_sel.begin<Vertex>(), m_sel.end<Vertex>(), OUTER_SI);
	m_sh.assign_subset(m_sel.begin<Edge>(), m_sel.end<Edge>(), OUTER_SI);
	m_sh.assign_subset(m_sel.begin<Face>(), m_sel.end<Face>(), OUTER_SI);
}



void MorphoGen::fix_dendritic_ends()
{
	// correct envelope end subset index previously set to temporary subset
	// for tetrahedralization of extracellular space
	m_sel.clear();
	m_sel.select(m_sh.begin<Edge>(ENVELOPE_END_TMP_SI), m_sh.end<Edge>(ENVELOPE_END_TMP_SI));
	m_sel.select(m_sh.begin<Face>(ENVELOPE_END_TMP_SI), m_sh.end<Face>(ENVELOPE_END_TMP_SI));
	m_sh.assign_subset(m_sel.begin<Edge>(), m_sel.end<Edge>(), NOFLUX_BND);
	m_sh.assign_subset(m_sel.begin<Face>(), m_sel.end<Face>(), NOFLUX_BND);

	// correct wrong subset at rim of membrane
	// select outer envelope ring
	for (size_t i = 0; i < 2; ++i)
	{
		int mult = i == 0 ? 1 : -1;

		m_sel.clear();

		std::vector<Vertex*> vrts;
		std::vector<Edge*> edges;
		try {find_ring_at_dendritic_ends(MEM_OUTER_BND_SI, mult, vrts, edges);}
		UG_CATCH_THROW("Could not find ring of outer membrane at dendritic end.");

		m_sel.select(vrts.begin(), vrts.end());
		m_sel.select(edges.begin(), edges.end());

		VertexIterator vit = m_sel.begin<Vertex>();
		VertexIterator vit_end = m_sel.end<Vertex>();
		for (; vit != vit_end; ++vit)
		{
			number xCoord = m_aaPos[*vit].x();

			// every associated edge in membrane subset needs to be changed
			Grid::AssociatedEdgeIterator edgeIt = m_grid.associated_edges_begin(*vit);
			Grid::AssociatedEdgeIterator edge_end = m_grid.associated_edges_end(*vit);
			for (; edgeIt != edge_end; ++edgeIt)
			{
				if (m_sh.get_subset_index(*edgeIt) == MEM_SI)
					m_sh.assign_subset(*edgeIt, MEM_NOFLUX_BND_SI);
			}

			// change every associated face in membrane subset and with boundary coords
			Grid::AssociatedFaceIterator faceIt = m_grid.associated_faces_begin(*vit);
			Grid::AssociatedFaceIterator face_end = m_grid.associated_faces_end(*vit);
			for (; faceIt != face_end; ++faceIt)
			{
				vector3 fc = CalculateCenter(*faceIt, m_aaPos);
				if (fabs(fc.x() - xCoord) < 1e-4*fabs(xCoord))
				{
					if (m_sh.get_subset_index(*faceIt) == MEM_SI)
						m_sh.assign_subset(*faceIt, MEM_NOFLUX_BND_SI);
				}
			}
		}
	}

// fix face orientation
		// find one membrane bnd face on each side and make sure they are correctly oriented
	Face* left = NULL;
	Face* right = NULL;
	FaceIterator fit = m_sh.begin<Face>(MEM_NOFLUX_BND_SI);
	FaceIterator fit_end = m_sh.end<Face>(MEM_NOFLUX_BND_SI);
	for (; fit != fit_end; ++fit)
	{
		Face* f = *fit;

		// find out side
		if (m_aaPos[f->vertex(0)].x() < 0)
		{
			if (left) continue;
			left = f;

			// ensure correct orientation
			vector3 v1, v2;
			VecSubtract(v1, m_aaPos[left->vertex(1)], m_aaPos[left->vertex(0)]);
			VecSubtract(v2, m_aaPos[left->vertex(2)], m_aaPos[left->vertex(0)]);
			VecCross(v1, v1, v2);
			if (v1.z() > 0)
				m_grid.flip_orientation(left);

			if (right) break;
		}
		else
		{
			if (right) continue;
			right = f;

			// ensure correct orientation
			vector3 v1, v2;
			VecSubtract(v1, m_aaPos[right->vertex(1)], m_aaPos[right->vertex(0)]);
			VecSubtract(v2, m_aaPos[right->vertex(2)], m_aaPos[right->vertex(0)]);
			VecCross(v1, v1, v2);
			if (v1.z() < 0)
				m_grid.flip_orientation(right);

			if (left) break;
		}
	}
	UG_COND_THROW(!left || !right, "Left and/or right membrane boundary face not found.");

	// select all remaining dendritic side faces via SelectLinkedFlatFaces
	m_sel.clear();
	m_sel.select(left);
	m_sel.select(right);
	SelectLinkedFlatFaces(m_sel, 1, true, false,  aPosition);

	// fix orientation; i.e.: invert faces that are not oriented in the same direction
	// as the corresponding left or right face first selected here
	FixFaceOrientation(m_grid, m_sel.begin<Face>(), m_sel.end<Face>());
}



void MorphoGen::create_interfaces()
{
	typedef Grid::traits<Face>::secure_container face_list;

	for (size_t i = 0; i < 2; ++i)
	{
		int mult = i == 0 ? 1 : -1;

		int constrd_si = i == 0 ? INTF_LEFT_CONSTRD_SI : INTF_RIGHT_CONSTRD_SI;

		// select inner membrane bnd side
		std::vector<Vertex*> vrts;
		std::vector<Edge*> edges;
		try {find_ring_at_dendritic_ends(MEM_INNER_BND_SI, mult, vrts, edges);}
		UG_CATCH_THROW("Could not find ring of inner membrane boundary, "
			<< (i==0 ? "left" : "right") << " side.");

		m_sel.clear();
		m_sel.select(vrts.begin(), vrts.end());
		m_sel.select(edges.begin(), edges.end());

		// select one face from within
		typedef Grid::traits<Face>::secure_container face_list;
		face_list fl;
		Edge* e = *edges.begin();
		number ex = m_aaPos[e->vertex(0)].x();
		m_grid.associated_elements(fl, e);
		size_t flsz = fl.size();
		size_t f = 0;
		for (; f < flsz; ++f)
		{
			vector3 fc = CalculateCenter(fl[f], m_aaPos);
			if (sqrt(fc.y()*fc.y() + fc.z()*fc.z())
				< 0.99*(DENDRITE_RADIUS - MEMBRANE_RADIUS)
				&& fabs(fc.x() - ex) < 1e-4 * fabs(ex))
			{
				m_sel.select(fl[f]);
				break;
			}
		}
		UG_COND_THROW(f == flsz, "No inner dendrite boundary face found.");

		// select the rest of the inner boundary faces and close
		SelectLinkedFlatFaces(m_sel, 1, false, true, aPosition);
		SelectAssociatedEdges(m_sel, m_sel.begin<Face>(), m_sel.end<Face>());
		SelectAssociatedVertices(m_sel, m_sel.begin<Edge>(), m_sel.end<Edge>());

		// extrude
		std::vector<Face*> faces;
		vrts.assign(m_sel.begin<Vertex>(), m_sel.end<Vertex>());
		edges.assign(m_sel.begin<Edge>(), m_sel.end<Edge>());
		faces.assign(m_sel.begin<Face>(), m_sel.end<Face>());
		vector3 dir(-mult*MEMBRANE_ENVELOPE_RADIUS, 0.0, 0.0);
		m_sh.set_default_subset_index(INNER_SI);
		Extrude(m_grid, &vrts, &edges, &faces, dir, m_aaPos,
				EO_CREATE_FACES | EO_CREATE_VOLUMES, NULL);

		// fix orientation
		FixFaceOrientation(m_grid, faces.begin(), faces.end());

		// assign constrained subset
		m_sh.assign_subset(vrts.begin(), vrts.end(), constrd_si);
		m_sh.assign_subset(edges.begin(), edges.end(), constrd_si);
		m_sh.assign_subset(faces.begin(), faces.end(), constrd_si);

		// assign interface node (selected vrt with minimal radius)
		std::vector<Vertex*>::const_iterator vit = vrts.begin();
		std::vector<Vertex*>::const_iterator vit_end = vrts.end();
		Vertex* minVrt = NULL;
		number minR2 = DENDRITE_RADIUS*DENDRITE_RADIUS + DENDRITE_LENGTH*DENDRITE_LENGTH;
		for (; vit != vit_end; ++vit)
		{
			number r2 = VecLength(m_aaPos[*vit]);
			if (r2 < minR2)
			{
				minR2 = r2;
				minVrt = *vit;
			}
		}
		m_sh.assign_subset(minVrt, i == 0 ? INTF_LEFT_NODE1D_SI : INTF_RIGHT_NODE1D_SI);

		// assign inner interface node
		Grid::AssociatedEdgeIterator edgeIt = m_grid.associated_edges_begin(minVrt);
		Grid::AssociatedEdgeIterator edge_end = m_grid.associated_edges_end(minVrt);
		for (; edgeIt != edge_end; ++edgeIt)
		{
			if (m_sh.get_subset_index(*edgeIt) == INNER_SI)
			{
				Vertex* intf_hd = GetOpposingSide(m_grid, *edgeIt, minVrt);
				m_sh.assign_subset(intf_hd, i == 0 ? INTF_LEFT_NODEHD_SI : INTF_RIGHT_NODEHD_SI);
				break;
			}
		}
		UG_COND_THROW(edgeIt == edge_end, "Unable to find high-dimensional interface node.");

		// fix orientation of lateral extrusion faces
		faces.clear();
		edges.assign(m_sh.begin<Edge>(constrd_si), m_sh.end<Edge>(constrd_si));
		size_t nEdge = edges.size();
		for (size_t j = 0; j < nEdge; ++j)
		{
			face_list fl;
			m_grid.associated_elements(fl, edges[j]);
			size_t flsz = fl.size();
			if (flsz != 2) continue; // only rim edges
			size_t fi = 0;
			for (; fi < flsz; ++fi)
			{
				if (m_sh.get_subset_index(fl[fi]) == MEM_INNER_BND_SI)
					break;
			}
			UG_COND_THROW(fi == flsz, "No lateral extrusion face for edge in extrusion rim.");
			faces.push_back(fl[fi]);
		}

		Face* face = faces[0];
		vector3 v1, v2;
		VecSubtract(v1, m_aaPos[face->vertex(1)], m_aaPos[face->vertex(0)]);
		VecSubtract(v2, m_aaPos[face->vertex(2)], m_aaPos[face->vertex(0)]);
		VecCross(v1, v1, v2);
		v2 = m_aaPos[face->vertex(0)];
		v2.x() = 0.0;
		if (VecProd(v1, v2) < 0)
			m_grid.flip_orientation(face);

		FixFaceOrientation(m_grid, faces.begin(), faces.end());
	}
}


void MorphoGen::create_extensions()
{
	m_sh.set_default_subset_index(USELESS_SI);

	for (size_t i = 0; i < 2; ++i)
	{
		int mult = i == 0 ? 1 : -1;
		int ext_si = i == 0 ? EXT_LEFT_SI : EXT_RIGHT_SI;

		// get 1d intf node
		std::vector<Vertex*> v(8);
		v[0] = *m_sh.begin<Vertex>(i == 0 ? INTF_LEFT_NODE1D_SI : INTF_RIGHT_NODE1D_SI);

		// construct 1 extension element
		v[1] = *m_grid.create<RegularVertex>();	// next in extension
		v[2] = *m_grid.create<RegularVertex>();
		v[3] = *m_grid.create<RegularVertex>();
		v[4] = *m_grid.create<RegularVertex>();
		v[5] = *m_grid.create<RegularVertex>();
		v[6] = *m_grid.create<RegularVertex>();
		v[7] = *m_grid.create<RegularVertex>();

		number fac = DENDRITE_RADIUS-MEMBRANE_RADIUS;
		vector3 dir(-mult*MEMBRANE_ENVELOPE_RADIUS, 0.0, 0.0);
		VecAdd(m_aaPos[v[1]], m_aaPos[v[0]], dir);
		VecAdd(m_aaPos[v[2]], m_aaPos[v[1]], vector3(-0.01*mult*fac, -mult*0.5*fac, 0.5*fac));
		VecAdd(m_aaPos[v[3]], m_aaPos[v[0]], vector3(-0.01*mult*fac, -mult*0.5*fac, 0.5*fac));
		for (size_t j = 0; j < 4; ++j)
			VecAdd(m_aaPos[v[j+4]], m_aaPos[v[j]], vector3(-0.01*mult*fac, mult*0.5*fac, 0.5*fac));

		Edge* e = *m_grid.create<RegularEdge>(EdgeDescriptor(v[0], v[1]));
		m_grid.create<Hexahedron>(HexahedronDescriptor(v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7]));

		// assign subsets
		m_sh.assign_subset(e, ext_si);
		m_sh.assign_subset(v[1], ext_si);

		// find outer face and mark
		Grid::AssociatedFaceIterator fit = m_grid.associated_faces_begin(v[1]);
		Grid::AssociatedFaceIterator fit_end = m_grid.associated_faces_end(v[1]);
		for (; fit != fit_end; ++fit)
		{
			if ((*fit)->vertex((*fit)->get_opposing_object(v[1]).second) == v[6])
				break;
		}
		UG_COND_THROW(fit == fit_end, "Front face for extension extrusion not found.");

		m_sel.clear();
		m_sel.select(*fit);
		SelectAssociatedEdges(m_sel, m_sel.begin<Face>(), m_sel.end<Face>());
		SelectAssociatedVertices(m_sel, m_sel.begin<Edge>(), m_sel.end<Edge>());

		// double extrusion length until compartment length reached
		number nExtrusion = floor(log2(EXTENSION_COMPARTMENT_LENGTH / MEMBRANE_ENVELOPE_RADIUS));
		size_t nExtr = (size_t) std::max(0.0, nExtrusion);

		std::vector<Vertex*> vrts;
		std::vector<Edge*> edges;
		std::vector<Face*> faces;
		vrts.assign(m_sel.begin<Vertex>(), m_sel.end<Vertex>());
		edges.assign(m_sel.begin<Edge>(), m_sel.end<Edge>());
		faces.assign(m_sel.begin<Face>(), m_sel.end<Face>());
		for (size_t j = 0; j < nExtr; ++j)
		{
			dir.x() *= 2.0;
			Extrude(m_grid, &vrts, &edges, &faces, dir, m_aaPos,
					EO_CREATE_FACES | EO_CREATE_VOLUMES, NULL);
		}

		// make rest of compartments
		dir.x() = -mult*EXTENSION_COMPARTMENT_LENGTH;
		nExtrusion = floor(EXTENSION_LENGTH / EXTENSION_COMPARTMENT_LENGTH) - 1;
		nExtr = (size_t) std::max(0.0, nExtrusion);
		for (size_t j = 0; j < nExtr; ++j)
		{
			Extrude(m_grid, &vrts, &edges, &faces, dir, m_aaPos,
					EO_CREATE_FACES | EO_CREATE_VOLUMES, NULL);
		}

		// mark last vertex as boundary
		size_t sz = vrts.size();
		size_t k = 0;
		for (; k < sz; ++k)
		{
			if (m_sh.get_subset_index(vrts[k]) == ext_si)
				break;
		}
		UG_COND_THROW(k == sz, "Boundary vertex of " << (i == 0 ? "left" : "right")
					  << " extension not found.");

		m_sh.assign_subset(vrts[k], i == 0 ? EXT_LEFT_BND_SI : EXT_RIGHT_BND_SI);
	}
}


void MorphoGen::create_dendrite(const std::string& filename)
{
	// check that filename ends in ".ugx"
	std::string useFileName = filename;
	if (GetFilenameExtension(filename) != std::string("ugx"))
	{
		UG_LOGN("File name extension needs to be '.ugx' - appending extension.")
		useFileName.append(".ugx");
	}
	std::string filePath = FindDirInStandardPaths(PathFromFilename(useFileName).c_str());
	if (filePath.empty())
		UG_THROW("Directory '" << PathFromFilename(useFileName) << "' could not be located. "
				"The file cannot be written there.");

	// start off by creating the left end of the dendrite
	try {create_shaft();}
	UG_CATCH_THROW("Could not create dendritic shaft morphology.");

	// graft spine onto shaft
	try {graft_spine();}
	UG_CATCH_THROW("Could not graft spine onto dendritic shaft.");

	// insert filaments
	if (m_bFilAnisotropic)
	{
		try {make_neck_filaments();}
		UG_CATCH_THROW("Could not create anisotropic neck filaments.");
	}
	else
	{
		try {make_spherical_filaments();}
		UG_CATCH_THROW("Could not create isotropic filaments.");
	}

	// triangulate the whole surface
	try {Triangulate(m_grid, m_grid.begin<Quadrilateral>(), m_grid.end<Quadrilateral>(), &m_aaPos);}
	UG_CATCH_THROW("Triangulation of surface failed.")

	// create membrane and prisms envelope around charged surfaces
	try {create_membrane_and_envelopes();}
	UG_CATCH_THROW("Could not create membrane and charge envelopes.");

	// close dendritic ends
	try {close_ends();}
	UG_CATCH_THROW("Could not close dendritic ends.");

	// create box bounding extracellular space
	try {create_bounding_box();}
	UG_CATCH_THROW("Could not create bounding box.");
/*
	// tetrahedralize filaments, inner and outer
	try {tetrahedralize();}
	UG_CATCH_THROW("Tetrahedralization failed.");

	// subset correction
	try {fix_dendritic_ends();}
	UG_CATCH_THROW("Subset correction after tetrahedralization failed.");

	// create interfaces
	try {create_interfaces();}
	UG_CATCH_THROW("Interface creation failed.");

	// create extensions
	try {create_extensions();}
	UG_CATCH_THROW("Extension creation failed.");

	// fix volume orientation
	FixOrientation(m_grid, m_grid.begin<Volume>(), m_grid.end<Volume>(), m_aaPos);
*/
	// name and colorize subsets
	AssignSubsetColors(m_sh);
	m_sh.set_subset_name("Intracellular_Domain", 0);
	m_sh.set_subset_name("Extracellular_Domain", 1);
	m_sh.set_subset_name("Membrane_Domain", 2);
	//m_sh.set_subset_name("HeadFil_Domain", 3);
	m_sh.set_subset_name("Fil_Domain", 3);
	m_sh.set_subset_name("MembraneInner_Bnd", 4);
	m_sh.set_subset_name("MembraneOuter_Bnd", 5);
	m_sh.set_subset_name("MemNoFlx_Bnd", 6);
	m_sh.set_subset_name("SurfCh_Bnd", 7);
	m_sh.set_subset_name("NoFlx_Bnd", 8);
	m_sh.set_subset_name("Diri_Bnd", 9);
	m_sh.set_subset_name("PSDInner_Bnd", 10);
	m_sh.set_subset_name("PSDOuter_Bnd", 11);
	m_sh.set_subset_name("ExtensionLeft_Domain", 12);
	m_sh.set_subset_name("ExtensionRight_Domain", 13);
	m_sh.set_subset_name("InterfaceLeft_Constrained", 14);
	m_sh.set_subset_name("InterfaceRight_Constrained", 15);
	m_sh.set_subset_name("InterfaceLeft_Node1D", 16);
	m_sh.set_subset_name("InterfaceRight_Node1D", 17);
	m_sh.set_subset_name("InterfaceLeft_NodehD", 18);
	m_sh.set_subset_name("InterfaceRight_NodehD", 19);
	m_sh.set_subset_name("ExtensionLeft_Bnd", 20);
	m_sh.set_subset_name("ExtensionRight_Bnd", 21);
	m_sh.set_subset_name("useless", 22);

	// save to .ugx file
	std::string fileName = FilenameWithoutPath(useFileName);
	GridWriterUGX ugxWriter;
	ugxWriter.add_grid(m_grid, "defGrid", aPosition);
	ugxWriter.add_subset_handler(m_sh, "defSH", 0);
	ugxWriter.add_subset_handler(m_shProj, "projSH", 0);
	ugxWriter.add_projection_handler(m_projHandler, "defPH", 0);
	if (!ugxWriter.write_to_file((filePath+fileName).c_str()))
		UG_THROW("Grid could not be written to file '" << filePath << fileName << "'.");


	// test loading/refining
	Domain3d dom;
	dom.create_additional_subset_handler("projSH");
	try {LoadDomain(dom, (filePath+fileName).c_str());}
	UG_CATCH_THROW("Failed loading domain from '" << filePath << fileName << "'.");

	GlobalMultiGridRefiner ref(*dom.grid(), dom.refinement_projector());
	ref.refine();

	fileName = fileName.substr(0, fileName.size()-4).append("_refined.ugx");
	number offset = 2*DENDRITE_RADIUS + 2*MEMBRANE_RADIUS + 2*MEMBRANE_ENVELOPE_RADIUS
		+ m_tmpHeadHeight + SPINE_NECK_LENGTH + 4*BOX_MARGIN;
	try {SaveGridHierarchyTransformed(*dom.grid(), *dom.subset_handler(), (filePath+fileName).c_str(), offset);}
	UG_CATCH_THROW("Grid could not be written to file '" << filePath << fileName << "'.");

}


void MorphoGen::create_dendrite_2d(const std::string& filename, number radius)
{
	number ELEM_LENGTH = 250.0;
	number MEM_RADIUS = 10.0;
	size_t NUM_ELEMS = 200;

	Grid g;
	SubsetHandler sh(g);
	sh.set_default_subset_index(0);
	g.attach_to_vertices(aPosition);
	AAPosition aaPos = AAPosition(g, aPosition);

	// create vertex at zero
	Vertex* v_start = *g.create<RegularVertex>();
	aaPos[v_start] = vector3(0,0,0);
	sh.assign_subset(v_start, 0);

	// extrude
	std::vector<Vertex*> vrts;
	vrts.push_back(v_start);
	vector3 extrudeDir(0.0);
	extrudeDir.coord(0) = ELEM_LENGTH;

	for (size_t i = 0; i < NUM_ELEMS; ++i)
		Extrude(g, &vrts, NULL, NULL, extrudeDir, aaPos, 0, NULL);

	// extrude inner domain
	Vertex* v_end = vrts[0];
	sh.assign_subset(v_start, 6);
	sh.assign_subset(v_end, 8);
	vrts.assign(g.begin<Vertex>(), g.end<Vertex>());
	std::vector<Edge*> edges;
	edges.assign(g.begin<Edge>(), g.end<Edge>());

	std::vector<Vertex*> oldVrts;
	oldVrts = vrts;
	std::vector<Edge*> oldEdges;
	oldEdges = edges;

	extrudeDir.coord(0) = 0.0;
	extrudeDir.coord(1) = radius;
	Extrude(g, &vrts, &edges, NULL, extrudeDir, aaPos, EO_CREATE_FACES, NULL);

	// assign proper subsets
	sh.assign_subset(oldEdges.begin(), oldEdges.end(), 9);
	sh.assign_subset(oldVrts.begin(), oldVrts.end(), 9);
	sh.assign_subset(v_start, 11);
	sh.assign_subset(v_end, 12);
	v_start = *vrts.begin();
	v_end = *(vrts.end()-1);
	oldVrts = vrts;
	oldEdges = edges;
	sh.assign_subset(edges.begin(), edges.end(), 2);
	sh.assign_subset(vrts.begin(), vrts.end(), 2);
	sh.assign_subset(v_start, 10);
	sh.assign_subset(v_end, 10);

	// extrude membrane
	extrudeDir.coord(1) = MEM_RADIUS;
	Extrude(g, &vrts, &edges, NULL, extrudeDir, aaPos, EO_CREATE_FACES, NULL);

	// assign proper subsets
	sh.assign_subset(oldEdges.begin(), oldEdges.end(), 3);
	sh.assign_subset(oldVrts.begin(), oldVrts.end(), 3);
	sh.assign_subset(v_start, 3);
	sh.assign_subset(v_end, 3);
	v_start = *vrts.begin();
	v_end = *(vrts.end()-1);
	oldVrts = vrts;
	oldEdges = edges;
	sh.assign_subset(edges.begin(), edges.end(), 1);
	sh.assign_subset(vrts.begin(), vrts.end(), 1);
	sh.assign_subset(v_start, 7);
	sh.assign_subset(v_end, 7);

	// extrude extracellular domain
	extrudeDir.coord(1) = radius;
	Extrude(g, &vrts, &edges, NULL, extrudeDir, aaPos, EO_CREATE_FACES, NULL);

	// assign proper subsets
	sh.assign_subset(oldEdges.begin(), oldEdges.end(), 4);
	sh.assign_subset(oldVrts.begin(), oldVrts.end(), 4);
	sh.assign_subset(v_start, 4);
	sh.assign_subset(v_end, 4);
	v_start = *vrts.begin();
	v_end = *(vrts.end()-1);
	oldVrts = vrts;
	oldEdges = edges;
	sh.assign_subset(edges.begin(), edges.end(), 5);
	sh.assign_subset(vrts.begin(), vrts.end(), 5);
	sh.assign_subset(v_start, 5);
	sh.assign_subset(v_end, 5);

	// invert wrong orientation of all faces
	InvertOrientation(g, g.begin<Face>(), g.end<Face>());

	// subset names
	AssignSubsetColors(sh);
	sh.set_subset_name("intracellular", 0);
	sh.set_subset_name("extracellular", 1);
	sh.set_subset_name("membrane", 2);
	sh.set_subset_name("mem_surf_inner", 3);
	sh.set_subset_name("mem_surf_outer", 4);
	sh.set_subset_name("bnd_dirichlet_outer", 5);
	sh.set_subset_name("bnd_influx_dendrite", 6);
	sh.set_subset_name("bnd_noflux", 7);
	sh.set_subset_name("bnd_dirichlet_dendrite", 8);
	sh.set_subset_name("bnd_inner", 9);
	sh.set_subset_name("bnd_mem_noflux", 10);
	sh.set_subset_name("bnd_influx_dendrite1d", 11);
	sh.set_subset_name("bnd_dirichlet_dendrite1d", 12);

	// save to .ugx file
	// check that filename ends in ".ugx"
	std::string useFileName = filename;
	if (GetFilenameExtension(filename) != std::string("ugx"))
	{
		UG_LOGN("File name extension needs to be '.ugx' - appending extension.")
		useFileName.append(".ugx");
	}
	std::string filePath = FindDirInStandardPaths(PathFromFilename(useFileName).c_str());
	if (filePath.empty())
		UG_THROW("Directory '" << PathFromFilename(useFileName) << "' could not be located. "
				"The file cannot be written there.");

	std::string fileName = FilenameWithoutPath(useFileName);
	GridWriterUGX ugxWriter;
	ugxWriter.add_grid(g, "defGrid", aPosition);
	ugxWriter.add_subset_handler(sh, "defSH", 0);
	if (!ugxWriter.write_to_file((filePath+fileName).c_str()))
		UG_THROW("Grid could not be written to file '" << filePath << fileName << "'.");
}


void MorphoGen::create_dendrite_1d(const std::string& filename)
{
	number ELEM_LENGTH = 250.0;
	size_t NUM_ELEMS = 200;

	Grid g;
	SubsetHandler sh(g);
	sh.set_default_subset_index(0);
	g.attach_to_vertices(aPosition);
	AAPosition aaPos(g, aPosition);

	// create vertex at zero
	Vertex* v = *g.create<RegularVertex>();
	aaPos[v] = vector3(0,0,0);
	sh.assign_subset(v, 0);

	// extrude
	std::vector<Vertex*> vrts;
	vrts.push_back(v);
	vector3 extrudeDir(0.0);
	extrudeDir.coord(0) = ELEM_LENGTH;

	for (size_t i = 0; i < NUM_ELEMS; ++i)
		Extrude(g, &vrts, NULL, NULL, extrudeDir, aaPos, 0, NULL);

	// assign proper subsets
	sh.assign_subset(v, 1);
	sh.assign_subset(vrts[0], 2);

	// subset names
	AssignSubsetColors(sh);
	sh.set_subset_name("bnd_inner", 0);
	sh.set_subset_name("bnd_influx_dendrite1d", 1);
	sh.set_subset_name("bnd_dirichlet_dendrite1d", 2);

	// save to .ugx file
	// check that filename ends in ".ugx"
	std::string useFileName = filename;
	if (GetFilenameExtension(filename) != std::string("ugx"))
	{
		UG_LOGN("File name extension needs to be '.ugx' - appending extension.")
		useFileName.append(".ugx");
	}
	std::string filePath = FindDirInStandardPaths(PathFromFilename(useFileName).c_str());
	if (filePath.empty())
		UG_THROW("Directory '" << PathFromFilename(useFileName) << "' could not be located. "
				"The file cannot be written there.");

	std::string fileName = FilenameWithoutPath(useFileName);
	GridWriterUGX ugxWriter;
	ugxWriter.add_grid(g, "defGrid", aPosition);
	ugxWriter.add_subset_handler(sh, "defSH", 0);
	if (!ugxWriter.write_to_file((filePath+fileName).c_str()))
		UG_THROW("Grid could not be written to file '" << filePath << fileName << "'.");
}


} // namespace nernst_planck
} // namespace ug
