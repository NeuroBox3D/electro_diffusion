/*
 * intf_refMarkAdjuster.cpp
 *
 *  Created on: 29.04.2016
 *      Author: mbreit
 */

#include "extension_refMarkAdjuster.h"

#include <cmath>                                                   // for fabs
#include <cstddef>                                                 // for size_t

#include "common/error.h"                                          // for UG_COND_...
#include "common/math/math_vector_matrix/math_vector_functions.h"  // for VecNorma...
#include "lib_disc/domain.h"                                       // for Domain1d, Domain2d, Dom...
#include "lib_grid/grid/grid.h"                                    // for Grid::associated_elements ...
#include "lib_grid/grid/grid_base_objects.h"                       // for Edge
#include "lib_grid/refinement/hanging_node_refiner_multi_grid.h"   // for HangingN...
#include "lib_grid/refinement/refiner_interface.h"                 // for IRefiner
#include "lib_grid/multi_grid.h"                                   // for MultiGrid


namespace ug {
namespace nernst_planck {

template <typename TDomain>
ExtensionRefMarkAdjuster<TDomain>::ExtensionRefMarkAdjuster
(
	SmartPtr<TDomain> dom,
	std::vector<number> dir,
	const std::string useless
)
: IRefMarkAdjuster(), m_si(-1)
{
	// set subset handler
	UG_COND_THROW(!dom.valid(), "Domain invalid.");
	m_ssh = dom->subset_handler();

	// find subset index for useless subset
	UG_COND_THROW(!m_ssh.valid(), "Subset handler not set. Please do so using set_subset_handler() method first.");
	m_si = m_ssh->get_subset_index(useless.c_str());
	UG_COND_THROW(m_si < 0, "Subset '" << useless << "' not found by subset handler.");

	// convert direction to MathVector
	UG_COND_THROW(dir.size() < (size_t) worldDim,
		"Expected direction with at least " << worldDim << " (worldDim) components.")

	for (size_t i = 0; i < (size_t) worldDim; ++i)
		m_direction[i] = dir[i];

	VecNormalize(m_direction, m_direction);

	// get position accessor from domain
	m_aaPos = dom->position_accessor();
}



template <typename TDomain>
void ExtensionRefMarkAdjuster<TDomain>::change_edge_marks(IRefiner& ref, Edge* e)
{
	// only edges in extension direction are refined, the rest is marked anisotropic
	MathVector<worldDim> edgeDir;
	VecSubtract(edgeDir, m_aaPos[e->vertex(1)], m_aaPos[e->vertex(0)]);
	VecNormalize(edgeDir, edgeDir);
	if (fabs(VecProd(edgeDir, m_direction)) < 0.95)	// that is about 18 degrees off
		ref.mark(e, RM_ANISOTROPIC);
	else
		ref.mark(e, RM_REFINE);
}


template <typename TDomain>
void ExtensionRefMarkAdjuster<TDomain>::change_face_marks(IRefiner& ref, Face* f)
{
	// all faces need to be marked anisotropic
	ref.mark(f, RM_ANISOTROPIC);

	// mark associated edges
	MultiGrid::traits<Edge>::secure_container el;
	ref.grid()->associated_elements(el, f);
	size_t el_sz = el.size();
	for (size_t k = 0; k < el_sz; ++k)
		change_edge_marks(ref, el[k]);
}


template <typename TDomain>
void ExtensionRefMarkAdjuster<TDomain>::change_vol_marks(IRefiner& ref, Volume* v)
{
	// all volumes need to be marked anisotropic
	ref.mark(v, RM_ANISOTROPIC);

	// mark associated faces
	MultiGrid::traits<Face>::secure_container fl;
	ref.grid()->associated_elements(fl, v);
	size_t fl_sz = fl.size();
	for (size_t k = 0; k < fl_sz; ++k)
		change_face_marks(ref, fl[k]);

}

template <typename TDomain>
void ExtensionRefMarkAdjuster<TDomain>::ref_marks_changed
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

	Grid::edge_traits::secure_container		assEdges;
	Grid::face_traits::secure_container		assFaces;
	Grid::volume_traits::secure_container 	assVols;


	// EDGES //
	size_t sz = edges.size();
	for (size_t i = 0; i < sz; ++i)
	{
		Edge* e = edges[i];

		// only change useless edges
		if (m_ssh->get_subset_index(e) != m_si)
			continue;

		change_edge_marks(ref, e);
	}

	// FACES //
	sz = faces.size();
	for (size_t i = 0; i < sz; ++i)
	{
		Face* f = faces[i];

		// only change useless faces
		if (m_ssh->get_subset_index(f) != m_si)
			continue;

		change_face_marks(ref, f);
	}

	// VOLUMES //
	typedef MultiGrid::traits<Edge>::secure_container edge_list_type;
	typedef MultiGrid::traits<Face>::secure_container face_list_type;
	sz = vols.size();
	for (size_t i = 0; i < sz; ++i)
	{
		Volume* v = vols[i];

		// only change useless volumes
		if (m_ssh->get_subset_index(v) != m_si)
			continue;

		change_vol_marks(ref, v);
	}
}


template <typename TDomain>
void add_extension_ref_mark_adjuster(IRefiner* ref, SmartPtr<ExtensionRefMarkAdjuster<TDomain> > erma)
{
	HangingNodeRefiner_MultiGrid* href = dynamic_cast<HangingNodeRefiner_MultiGrid*>(ref);
	UG_COND_THROW(!href, "An interface refinement mark adjuster can only be added to an instance of HangingNodeRefiner_MultiGrid.");

	href->add_ref_mark_adjuster(erma);
}



// explicit specializations
#ifdef UG_DIM_1
	template class ExtensionRefMarkAdjuster<Domain1d>;
	template void add_extension_ref_mark_adjuster<Domain1d>(IRefiner*, SmartPtr<ExtensionRefMarkAdjuster<Domain1d> >);
#endif
#ifdef UG_DIM_2
	template class ExtensionRefMarkAdjuster<Domain2d>;
	template void add_extension_ref_mark_adjuster<Domain2d>(IRefiner*, SmartPtr<ExtensionRefMarkAdjuster<Domain2d> >);
#endif
#ifdef UG_DIM_3
	template class ExtensionRefMarkAdjuster<Domain3d>;
	template void add_extension_ref_mark_adjuster<Domain3d>(IRefiner*, SmartPtr<ExtensionRefMarkAdjuster<Domain3d> >);
#endif


} // namespace nernst_planck
} // namespace ug
