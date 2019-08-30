/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Markus Breit
 * Creation date: 2016-04-29
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


#include "intf_refMarkAdjuster.h"

#include <cstddef>                                               // for size_t

#include "common/error.h"                                         // for UG_COND_THROW
#include "lib_grid/algorithms/element_side_util.h"                // for GetOpposingSide
#include "lib_grid/grid/grid.h"                                   // for Grid, Grid::traits, Grid::traits<...
#include "lib_grid/grid/grid_base_objects.h"                      // for Edge, Vertex, Face, Volume (ptr o...
#include "lib_grid/refinement/hanging_node_refiner_multi_grid.h"  // for HangingNodeRefiner_MultiGrid
#include "lib_grid/refinement/refiner_interface.h"                // for IRefiner, RefinementMark::RM_ANIS...
#include "lib_grid/multi_grid.h"                                  // for MultiGrid


namespace ug {
namespace nernst_planck {


void InterfaceRefMarkAdjuster::ref_marks_changed
(
	IRefiner& ref,
	const std::vector<Vertex*>& vrts,
	const std::vector<Edge*>& edges,
	const std::vector<Face*>& faces,
	const std::vector<Volume*>& vols
)
{
	UG_COND_THROW(!m_ssh.valid(), "No SubsetHandler set in InterfaceRefMarkAdjuster.");

	if (!ref.grid()) return;
	Grid& grid = *ref.grid();

	Grid::edge_traits::secure_container		assEdges;
	Grid::face_traits::secure_container		assFaces;
	Grid::volume_traits::secure_container 	assVols;

	/*
	// find out grid dimension
	size_t dim = 0;
	if (grid.num_volumes()) dim = 3;
	else if (grid.num_faces()) dim = 2;
	else if (grid.num_edges()) dim = 1;

	UG_COND_THROW(dim < 2, "This refinement mark adjuster only works for grids of dimension 2 or 3.");
	 */


	// EDGES //
	size_t sz = edges.size();
	for (size_t i = 0; i < sz; ++i)
	{
		// if one vertex is in an interface subset and the other is its constrainer
		// then this edge must be marked anisotropic
		Edge* e  = edges[i];

		for (size_t l = 0; l < m_vIntf.size(); ++l)
			if ((m_ssh->get_subset_index(e->vertex(0)) == m_vIntf[l]->constrained_subset_index()
					&& e->vertex(1) == m_vIntf[l]->get_constrainer_object(e->vertex(0)))
				||
				(m_ssh->get_subset_index(e->vertex(1)) == m_vIntf[l]->constrained_subset_index()
					&& e->vertex(0) == m_vIntf[l]->get_constrainer_object(e->vertex(1)))
				||
				(m_ssh->get_subset_index(e->vertex(0)) == m_vIntf[l]->intf_node_1d_subset_index()
					&& e->vertex(1) == m_vIntf[l]->get_constrainer_object(e->vertex(0)))
				||
				(m_ssh->get_subset_index(e->vertex(1)) == m_vIntf[l]->intf_node_1d_subset_index()
					&& e->vertex(0) == m_vIntf[l]->get_constrainer_object(e->vertex(1)))
				)
				ref.mark(e, RM_ANISOTROPIC);
	}

	// FACES //
	typedef MultiGrid::traits<Edge>::secure_container edge_list_type;
	sz = faces.size();
	for (size_t i = 0; i < sz; ++i)
	{
		Face* f = faces[i];

		edge_list_type el;
		grid.associated_elements(el, f);

		size_t el_sz = el.size();
		size_t k;
		for (k = 0; k < el_sz; ++k)
		{
			int edge_si = m_ssh->get_subset_index(el[k]);
			Edge* opp = GetOpposingSide(grid, f, el[k]);
			if (!opp) continue; // there might not be an opposing side

			// now check whether this face has an edge in one of the interface subsets
			// and if the opposing side is the corresponding constrainer
			for (size_t l = 0; l < m_vIntf.size(); ++l)
				if (edge_si == m_vIntf[l]->constrained_subset_index()
					&& opp == m_vIntf[l]->get_constrainer_object(el[k]))
					goto change_marks_face;
		}

		continue;

	change_marks_face:
		// mark the element for anisotropic refinement
		ref.mark(f, RM_ANISOTROPIC);

		// mark the interface side
		ref.mark(el[k], RM_REFINE);

		// find opposing side and mark that one as well
		Edge* opp = GetOpposingSide(grid, f, el[k]);
		ref.mark(opp, RM_REFINE);

		// mark other edges for anisotropic refinement
		for (size_t l = 0; l < el_sz; ++l)
		{
			if (l == k) continue;
			if (el[l] == opp) continue;

			ref.mark(el[l], RM_ANISOTROPIC);
		}
	}


	// VOLUMES //
	typedef MultiGrid::traits<Edge>::secure_container edge_list_type;
	typedef MultiGrid::traits<Face>::secure_container face_list_type;
	sz = vols.size();
	for (size_t i = 0; i < sz; ++i)
	{
		Volume* v = vols[i];

		// now check whether this volume has a face in one of the interface subsets
		face_list_type fl;
		grid.associated_elements(fl, v);

		size_t fl_sz = fl.size();
		size_t k;
		for (k = 0; k < fl_sz; ++k)
		{
			int face_si = m_ssh->get_subset_index(fl[k]);
			Face* opp = GetOpposingSide(grid, v, fl[k]);
			if (!opp) continue; // there might not be an opposing side

			for (size_t l = 0; l < m_vIntf.size(); ++l)
				if (face_si == m_vIntf[l]->constrained_subset_index()
					&& opp == m_vIntf[l]->get_constrainer_object(fl[k]))
					goto change_marks_volume;
		}

		continue;

	change_marks_volume:
		// mark the element for anisotropic refinement
		ref.mark(v, RM_ANISOTROPIC);

		// mark the interface side
		ref.mark(fl[k], RM_REFINE);

		// find opposing side and mark that one as well
		Face* opp = GetOpposingSide(grid, v, fl[k]);
		ref.mark(opp, RM_REFINE);

		// mark other faces for anisotropic refinement (do not forget their edges)
		for (size_t l = 0; l < fl_sz; ++l)
		{
			if (l == k) continue;
			if (fl[l] == opp) continue;

			ref.mark(fl[l], RM_ANISOTROPIC);

			// face edges (in case they are marked with RM_REFINE)
			edge_list_type el;
			grid.associated_elements(el, fl[l]);

			size_t el_sz = el.size();
			for (size_t m = 0; m < el_sz; ++m)
				ref.mark(el[m], RM_ANISOTROPIC);
		}

		// We have marked RM_ANISOTROPIC every single edge of the volume at this point;
		// now correct for the two opposing sides
		edge_list_type el;
		grid.associated_elements(el, fl[k]);
		size_t el_sz = el.size();
		for (size_t m = 0; m < el_sz; ++m)
			ref.mark(el[m], RM_REFINE);

		grid.associated_elements(el, opp);
		el_sz = el.size();
		for (size_t m = 0; m < el_sz; ++m)
			ref.mark(el[m], RM_REFINE);

	}
}



void add_interface_ref_mark_adjuster(IRefiner* ref, SmartPtr<InterfaceRefMarkAdjuster> irma)
{
	HangingNodeRefiner_MultiGrid* href = dynamic_cast<HangingNodeRefiner_MultiGrid*>(ref);
	UG_COND_THROW(!href, "An interface refinement mark adjuster can only be added to an instance of HangingNodeRefiner_MultiGrid.");

	href->add_ref_mark_adjuster(irma);
}




} // namespace nernst_planck
} // namespace ug
