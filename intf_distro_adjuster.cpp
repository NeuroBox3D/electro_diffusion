/*
 * intf_distro_adjuster.cpp
 *
 *  Created on: 19.05.2016
 *      Author: mbreit
 */

#include "intf_distro_adjuster.h"

#include <cstddef>                                          // for size_t
#include <set>                                              // for set, set<>::iterator

#include "common/error.h"                                   // for UG_COND_THROW
#include "lib_disc/domain.h"                                // for Domain1d, Domain2d, Domain3d
#include "lib_grid/grid/grid.h"                             // for Grid::traits, Grid::traits<>::secure_container
#include "lib_grid/grid/grid_base_object_traits.h"          // for geometry_traits, geometry_traits<>::const_i...
#include "lib_grid/grid/grid_base_objects.h"                // for Vertex, GridObject (ptr only)
#include "lib_grid/grid/neighborhood.h"                     // for CollectNeighbors
#include "lib_grid/multi_grid.h"                            // for MultiGrid
#ifdef UG_PARALLEL
	#include "lib_grid/parallelization/distributed_grid.h"  // for DistributedGridManager, ElementStatusTypes:...
#endif
#include "lib_grid/parallelization/distribution.h"          // for InterfaceStates::IS_VSLAVE, InterfaceStates...
#include "lib_grid/tools/grid_level.h"                      // for GridLevel
#include "lib_grid/tools/selector_multi_grid.h"             // for MGSelector
#include "lib_grid/tools/surface_view.h"                    // for SurfaceView, SurfaceView::SurfaceConstants:...


namespace ug {
namespace nernst_planck {

template <typename TDomain>
void PNPDistroManager<TDomain>::collect_neighbors
(
	std::vector<elem_type*>& neighborsOut,
	elem_type* elem
)
{
	// normal neighbors
	MultiGrid& mg = *m_dom->grid();
	CollectNeighbors(neighborsOut, elem, mg);

	// loop interfaces to find additional neighbors
	std::set<elem_type*> additionalNeighbors;
	size_t sz = m_vIntf.size();
	for (size_t i = 0; i < sz; ++i)
	{
		SmartPtr<IInterface1D> intf = m_vIntf[i];

		int constrdSI = intf->constrained_subset_index();
		int in1dSI = intf->intf_node_1d_subset_index();

		// check if elem is first elem of extension
		// that is the case if it has the 1d interface vertex, but no constrained side
		typedef typename MultiGrid::traits<Vertex>::secure_container vrt_list;
		vrt_list vl;
		mg.associated_elements(vl, elem);

		bool candidateElem = false;
		size_t sz = vl.size();
		for (size_t v = 0; v < sz; ++v)
		{
			if (m_sh->get_subset_index(vl[v]) == in1dSI)
			{
				candidateElem = true;
				break;
			}
		}
		if (!candidateElem) continue;

		typedef typename MultiGrid::traits<side_type>::secure_container side_list;
		side_list sl;
		mg.associated_elements(sl, elem);

		sz = sl.size();
		for (size_t s = 0; s < sz; ++s)
		{
			if (m_sh->get_subset_index(sl[s]) == constrdSI)
			{
				candidateElem = false;
				break;
			}
		}
		if (!candidateElem) continue;

		// add all (same-level) elems at interface to first extension element
		typedef typename geometry_traits<side_type>::const_iterator sh_it_type;
		typedef typename MultiGrid::traits<elem_type>::secure_container elem_list;

		int lvl = mg.get_level(elem);

		sh_it_type it = m_sh->begin<side_type>(constrdSI, lvl);
		sh_it_type it_end = m_sh->end<side_type>(constrdSI, lvl);
		for (; it != it_end; ++it)
		{
			elem_list el;
			mg.associated_elements(el, *it);

			UG_COND_THROW(el.size() != 1, "More (or less) than one element associated "
										  "to constrained side (" << el.size() << ").");

			UG_COND_THROW(!el[0], "NULL elem!");
			additionalNeighbors.insert(el[0]);
		}
	}

	// add additional neighbors to outgoing vector
	typename std::set<elem_type*>::iterator it = additionalNeighbors.begin();
	typename std::set<elem_type*>::iterator it_end = additionalNeighbors.end();
	for (; it != it_end; ++it)
		neighborsOut.push_back(*it);
}




template <typename TDomain>
PNPDistroManager<TDomain>::PNPDistroManager(SmartPtr<ApproximationSpace<TDomain> > approx)
:
#ifdef NPParmetis
parmetis::AnisotropyUnificator<TDomain, typename grid_dim_traits<TDomain::dim>::grid_base_object>(approx->domain()),
#endif
  m_approx(approx), m_dom(approx->domain()), m_sh(m_dom->subset_handler())
{
	// TODO: integrate this somehow (but not in AnisotropyUnificator)
	/*
	// set this class as its own neighbor collector
	SmartPtr<DualGraphNeighborCollector<elem_type> > spDGNC = make_sp(this);
	++*spDGNC.refcount_ptr(); // we do not want to accidentally destroy ourselves
	this->set_neighbor_collector(spDGNC);
	*/

	// no longer necessary!
	m_spGridCreationCallbackID = m_dom->grid()->message_hub()->register_class_callback(this,
		&ug::nernst_planck::PNPDistroManager<TDomain>::adjust_horizontal_interfaces);
}


template <typename TDomain>
void PNPDistroManager<TDomain>::adjust
(
	MGSelector& sel,
	bool partitionForLocalProc,
	bool createVerticalInterfaces
)
{
	typedef typename SurfaceView::traits<Vertex>::const_iterator sv_it_type;
	typedef typename geometry_traits<Vertex>::const_iterator sh_it_type;

	const MultiGrid& mg = *sel.multi_grid();
	const DistributedGridManager& dgm = *mg.distributed_grid_manager();
	const SurfaceView& sv = *m_approx->surface_view();

	// For any constrained node in the surface, a process needs the surface interface nodes.
	// For any constrained node in a level, a process needs the corresponding level interface nodes (if existing).

	size_t nLvl = mg.num_levels();
	size_t sz = m_vIntf.size();
	for (size_t i = 0; i < sz; ++i)
	{
		SmartPtr<IInterface1D> intf = m_vIntf[i];

		int constrdSI = intf->constrained_subset_index();
		int in1dSI = intf->intf_node_1d_subset_index();
		int inhdSI = intf->intf_node_hd_subset_index();

		std::vector<std::set<Vertex*> > levelVertices(nLvl);

	// 1. check whether process is to have any constrained nodes
		std::vector<size_t> nSelected(nLvl+1, 0);
		std::vector<size_t> nUnselected(nLvl+1, 0);

		// iterate surface constrained
		GridLevel gl;
		sv_it_type it = sv.begin<Vertex>(constrdSI, gl, SurfaceView::ALL_BUT_SHADOW_COPY);
		sv_it_type it_end = sv.end<Vertex>(constrdSI, gl, SurfaceView::ALL_BUT_SHADOW_COPY);

		size_t intfSrfLvl = 0;
		if (it != it_end)
		{
			sv_it_type it1d = sv.begin<Vertex>(in1dSI, gl, SurfaceView::ALL_BUT_SHADOW_COPY);
			sv_it_type it1d_end = sv.end<Vertex>(in1dSI, gl, SurfaceView::ALL_BUT_SHADOW_COPY);
			UG_COND_THROW(it1d == it1d_end, "Grid contains surface interface constrained, "
											"but not the necessary 1d interface node.");

			sv_it_type ithd = sv.begin<Vertex>(inhdSI, gl, SurfaceView::ALL_BUT_SHADOW_COPY);
			sv_it_type ithd_end = sv.end<Vertex>(inhdSI, gl, SurfaceView::ALL_BUT_SHADOW_COPY);
			UG_COND_THROW(ithd == ithd_end, "Grid contains surface interface constrained, "
											"but not the necessary hd interface node.");

			intfSrfLvl = mg.get_level(*it1d);
			UG_COND_THROW(intfSrfLvl != (size_t) mg.get_level(*ithd), "Surface level of 1d and hd interface nodes are not identical.\n"
																	 "This cannot be handled by the current implementation.");
		}

		for (; it != it_end; ++it)
		{
			if (sel.get_selection_status(*it) & (IS_NORMAL | IS_VSLAVE))
				++nSelected[intfSrfLvl];
			else
				++nUnselected[intfSrfLvl];
		}

		// now iterate level constrained
		for (size_t lv = 0; lv < nLvl; ++lv)
		{
			// iterate level constrained
			sh_it_type it = m_sh->begin<Vertex>(constrdSI, lv);
			sh_it_type it_end = m_sh->end<Vertex>(constrdSI, lv);

			for (; it != it_end; ++it)
			{
				if (sel.get_selection_status(*it) & (IS_NORMAL | IS_VSLAVE))
					++nSelected[lv];
				else
					++nUnselected[lv];
			}
		}

	// 2. add level interface nodes for level constrained
		for (size_t lv = 0; lv < nLvl; ++lv)
		{
			// if no surface constrained exist here: nothing to do
			if (!nSelected[lv]) continue;

			// find surface/level interface nodes
			sh_it_type it1d = m_sh->begin<Vertex>(in1dSI, lv);
			sh_it_type it1d_end = m_sh->end<Vertex>(in1dSI, lv);

			sh_it_type ithd = m_sh->begin<Vertex>(inhdSI, lv);
			sh_it_type ithd_end = m_sh->end<Vertex>(inhdSI, lv);

			// interface nodes might not exist on this level
			if (it1d != it1d_end)
			{
				sel.select(*it1d, sel.get_selection_status(*it1d) | IS_NORMAL);
				levelVertices[lv].insert(*it1d);
			}
			if (ithd != ithd_end)
			{
				sel.select(*ithd, sel.get_selection_status(*ithd) | IS_NORMAL);
				levelVertices[lv].insert(*ithd);
			}
		}

	// 3. treat root nodes
		if (partitionForLocalProc)
		{
			if (nUnselected[0])
			{
				sh_it_type it1d = m_sh->begin<Vertex>(in1dSI, 0);
				sh_it_type it1d_end = m_sh->end<Vertex>(in1dSI, 0);

				sh_it_type ithd = m_sh->begin<Vertex>(inhdSI, 0);
				sh_it_type ithd_end = m_sh->end<Vertex>(inhdSI, 0);

				if (it1d != it1d_end && !dgm.contains_status(*it1d, ES_V_SLAVE))
				{
//UG_LOGN("Selecting 1d node " << *it1d << " on level 0 as VMaster for local proc (unselected root).");
					sel.select(*it1d, sel.get_selection_status(*it1d) | IS_VMASTER);
				}
				if (ithd != ithd_end && !dgm.contains_status(*ithd, ES_V_SLAVE))
				{
//UG_LOGN("Selecting hd node " << *ithd << " on level 0 as VMaster for local proc (unselected root).");
					sel.select(*ithd, sel.get_selection_status(*ithd) | IS_VMASTER);
				}
			}
		}
		else
		{
			std::set<Vertex*>& vrtSet = levelVertices[0];
			std::set<Vertex*>::iterator it = vrtSet.begin();
			std::set<Vertex*>::iterator it_end = vrtSet.end();

			for (; it != it_end; ++it)
			{
				if (!dgm.contains_status(*it, ES_V_MASTER))
				{
//UG_LOGN("Selecting node " << *it << " on level 0 as VSlave for proc " << pcl::ProcRank() << " (selected root).");
					sel.select(*it, IS_VSLAVE);
				}
			}
		}

	// 4. assign vertical master and slave states on higher levels
		if (!createVerticalInterfaces) continue;

		// now for each level: mark appropriately
		for (size_t lv = nLvl-1; lv < nLvl; --lv)
		{
			std::set<Vertex*>& vrtSet = levelVertices[lv];
			std::set<Vertex*>::iterator it = vrtSet.begin();
			std::set<Vertex*>::iterator it_end = vrtSet.end();

			for (; it != it_end; ++it)
			{
				Vertex* vrt = *it;

				if ((sel.get_selection_status(vrt) & IS_VMASTER)
					|| (sel.get_selection_status(vrt) & IS_VSLAVE))
					continue;

				// assign vertical master states
				size_t numChildren = mg.num_children<Vertex>(vrt);
				GridObject* parent = mg.get_parent(vrt);
				bool parentIsSelected = false;
				if (parent)
					parentIsSelected = sel.is_selected(parent);

				if (numChildren)
				{
					Vertex* ch = mg.get_child<Vertex>(vrt, 0);
					if (!sel.is_selected(ch))
					{
//UG_LOGN("Selecting node " << ch << " on level " << lv+1 << " as VMaster (unselected child).");
						sel.select(ch, IS_VMASTER);
					}
				}
				else if (dgm.contains_status(vrt, ES_V_MASTER))
				{
					if (parentIsSelected || (partitionForLocalProc && (lv == 0)))
					{
//UG_LOGN("Selecting node " << vrt << " on level " << lv << " as VMaster (remaining vmaster).");
						sel.select(vrt, IS_VMASTER);
						continue;
					}
					else
					{
						//sel.deselect(vrt);	// TODO: think about whether this is correct!
						//continue;
					}
				}

				// assign slave states
				if (parent)
				{
					if (!sel.is_selected(parent))
					{
//UG_LOGN("Selecting node " << vrt << " on level " << lv << " as VSlave (unselected parent).");
						sel.select(vrt, IS_VSLAVE);
					}
				}
				else
				{
					if (dgm.contains_status(vrt, ES_V_SLAVE))
					{
//UG_LOGN("Selecting node " << vrt << " on level " << lv << " as VSlave (remaining vslave).");
						sel.select(vrt, IS_VSLAVE);
					}
				}
			}
		}
	}
}



template <typename TDomain>
void PNPDistroManager<TDomain>::adjust_horizontal_interfaces(const GridMessage_Creation& msg)
{
	if (msg.msg() != GMCT_CREATION_STOPS) return;

	// We only have to make sure that the interface node h-masters
	// are not located on such processes that do not have any adjacent elements
	// (this would be bad for smoothing).

	typedef typename geometry_traits<Vertex>::const_iterator sh_it_type;
	typedef typename MultiGrid::traits<elem_type>::secure_container elem_list;

	// determine subsets to treat (use std::set to exclude double entries)
	std::set<int> sis;
	const size_t nIntf = m_vIntf.size();
	for (size_t i = 0; i < nIntf; ++i)
	{
		SmartPtr<IInterface1D> intf = m_vIntf[i];
		sis.insert(intf->intf_node_1d_subset_index());
		sis.insert(intf->intf_node_hd_subset_index());
	}
	size_t nsi = sis.size();

	MultiGrid& mg = *m_dom->grid();

	// assign interface node states:
	//   state = 0 mod 3  does not exist on this proc
	//   state = 1 mod 3  exists but has no non-ghost elem neighbor
	//   state = 2 mod 3  exists and has non-ghost elem neighbor
	//   state / 3 encodes the element status (h-master or h-slave, etc.)
	const size_t nLvl = mg.num_levels();
	std::vector<int> localState(nLvl*nsi, 0);
	DistributedGridManager& dgm = *mg.distributed_grid_manager();

	std::set<int>::const_iterator itSI = sis.begin();
	std::set<int>::const_iterator itSIEnd = sis.end();
	for (size_t si = 0; itSI != itSIEnd; ++itSI, ++si)
	{
		for (size_t lv = 0; lv < nLvl; ++lv)
		{
			sh_it_type it = m_sh->begin<Vertex>(*itSI, lv);
			sh_it_type it_end = m_sh->end<Vertex>(*itSI, lv);
			if (it != it_end)
			{
				localState[nLvl*si + lv] = 1;

				// h-master needs to have at least one full-dim elem neighbor
				// AND that neighbor needs to be non-ghost
				elem_list el;
				mg.associated_elements(el, *it);
				const size_t elSz = el.size();
				for (size_t e = 0; e < elSz; ++e)
				{
					if (!dgm.is_ghost(el[e]))
					{
						localState[nLvl*si + lv] = 2;
						break;
					}
				}

				localState[nLvl*si + lv] += 3*dgm.get_status(*it);
			}
		}
	}

	// allgather interface node states
	size_t nProcs = pcl::NumProcs();
	std::vector<int> globalState(nProcs*nLvl*nsi);
	pcl::ProcessCommunicator pc;
	pc.allgather(GetDataPtr(localState), nLvl*nsi, PCL_DT_INT, GetDataPtr(globalState), nLvl*nsi, PCL_DT_INT);

	// loop interface nodes and levels and change HM where needed
	typedef typename GridLayoutMap::Types<Vertex>::Interface IntfType;
	GridLayoutMap& glm = dgm.grid_layout_map();
	CompareByAttachment<Vertex, AGeomObjID> gidCmp(*m_dom->grid(), aGeomObjID); // needed for interface entry search
	dgm.enable_interface_management(false);
	int locRank = pcl::ProcRank();

	size_t si = 0;
	for (itSI = sis.begin(); itSI != itSIEnd; ++itSI, ++si)
	{
		for (size_t lv = 0; lv < nLvl; ++lv)
		{
			// find out whether horizontal master has to be changed
			// and - if so - who shall be the new master
			int hMaster = -1;
			bool badHMaster = false;
			int minHMCandidate = -1;
			std::vector<int> hSlaveCandidates;
			int ownState = 0;

			for (size_t p = 0; p < nProcs; ++p)
			{
				const int& state = globalState[nLvl*(p*nsi + si) + lv];

				if (p == (size_t) locRank)
					ownState = state;

				if (state / 3 & ES_H_MASTER)
				{
					UG_COND_THROW(hMaster != -1, "Two h-masters!");
					hMaster = p;
					badHMaster = (state % 3 != 2);
				}

				// non-h-interface vertices are not considered for new interfaces either
				if (!(state / 3 & (ES_H_MASTER | ES_H_SLAVE)))
					continue;

				if (state % 3 == 1)
					hSlaveCandidates.push_back(p);
				else if (state % 3 == 2)
				{
					if (minHMCandidate == -1)
						minHMCandidate = p;
					else
						hSlaveCandidates.push_back(p);
				}
			}

			// if hmaster is already well chosen, do not change it
			if (!badHMaster)
				continue;

//UG_LOGN("  Bad choice of horizontal master for subset " << *itSI << " on level " << lv << ", adjusting.");

			UG_COND_THROW(minHMCandidate == -1, "No proc seems to have an interface node "
				"connected to a non-ghost element for subset " << *itSI << " on level " << lv << ".");

			// no changes to interfaces of procs that do not have the interface node
			// or of procs that did not participate in the horizontal interface before
			if ((ownState % 3 == 0) || !(ownState / 3 & (ES_H_MASTER | ES_H_SLAVE)))
				continue;

			// change interfaces
			// (1) remove old entries
			// (2) create new

			// (1) remove old entries
			Vertex* intfVrt = *m_sh->begin<Vertex>(*itSI, lv);
			size_t nSlaves = hSlaveCandidates.size();
			if (locRank == hMaster)
			{
				for (size_t i = 0; i < nSlaves; ++i)
				{
					int p = hSlaveCandidates[i];

					// instead of local proc (which surely is new HS)
					// treat the new HM which surely was HS before
					if (p == locRank)
						p = minHMCandidate;

					IntfType& ii = glm.get_layout<Vertex>(INT_H_MASTER).interface(p, lv);
					IntfType::iterator posIt = ii.find_insert_pos_sorted(intfVrt, gidCmp);
					UG_COND_THROW(posIt == ii.end() || gidCmp(intfVrt, ii.get_element(posIt)),
						"Interface entry to proc " << p << ", that is to be removed, does not exist.");
					ii.erase(posIt);
				}
			}
			else
			{
				IntfType& ii = glm.get_layout<Vertex>(INT_H_SLAVE).interface(hMaster, lv);
				IntfType::iterator posIt = ii.find_insert_pos_sorted(intfVrt, gidCmp);
				UG_COND_THROW(posIt == ii.end() || gidCmp(intfVrt, ii.get_element(posIt)),
					"Interface entry to proc " << hMaster << ", that is to be removed, does not exist.");
				ii.erase(posIt);
			}

			// (2) insert new entries
			if (locRank == minHMCandidate)
			{
				for (size_t i = 0; i < nSlaves; ++i)
				{
					int p = hSlaveCandidates[i];

					IntfType& ii = glm.get_layout<Vertex>(INT_H_MASTER).interface(p, lv);
					IntfType::iterator posIt = ii.find_insert_pos_sorted(intfVrt, gidCmp);
					UG_COND_THROW(posIt != ii.end() && !gidCmp(intfVrt, ii.get_element(posIt)),
						"Interface entry that is to be created already exists.");
					ii.insert(intfVrt, posIt);
				}
			}
			else
			{
				IntfType& ii = glm.get_layout<Vertex>(INT_H_SLAVE).interface(minHMCandidate, lv);
				IntfType::iterator posIt = ii.find_insert_pos_sorted(intfVrt, gidCmp);
				UG_COND_THROW(posIt != ii.end() && !gidCmp(intfVrt, ii.get_element(posIt)),
					"Interface entry that is to be created already exists.");
				ii.insert(intfVrt, posIt);
			}
		}
	}

	glm.remove_empty_interfaces();
	dgm.enable_interface_management(true);
	dgm.grid_layouts_changed(false);
}


#ifdef NPParmetis
template <typename TDomain>
void PNPDistroManager<TDomain>::unify
(
	MultiGrid* mg,
	int lvl,
	int localOffset,
	const Grid::AttachmentAccessor<elem_type, AElemIndex>& aaElemInd,
	const Grid::AttachmentAccessor<side_type, AElemIndices>& aaSideElemInd,
	std::vector<std::pair<int, int> >& unificationPairs
) const
{
	// call AnisotropyUnificator::unify() first
	parmetis::AnisotropyUnificator<TDomain, elem_type>::unify(mg, lvl, localOffset, aaElemInd, aaSideElemInd, unificationPairs);

	// add protection of the interfaces
	size_t nIntf = m_vIntf.size();
	for (size_t i = 0; i < nIntf; ++i)
	{
		SmartPtr<IInterface1D> intf = m_vIntf[i];
		const int nodeHdSI = intf->intf_node_1d_subset_index();
		const int constrainedSI = intf->constrained_subset_index();

		// get hd interface node
		typedef typename geometry_traits<Vertex>::const_iterator sh_it_type;
		Vertex* vrt = NULL;
		sh_it_type it = m_sh->begin<Vertex>(nodeHdSI, lvl);
		sh_it_type itEnd = m_sh->end<Vertex>(nodeHdSI, lvl);

		if (it == itEnd)
			continue;

		vrt = *it;
		UG_COND_THROW(++it != itEnd, "More than one hd interface node in interface on level " << lvl
			<< " for interface " << i << ".");

		// find all connected full-dim elements
		typedef typename MultiGrid::traits<elem_type>::secure_container elem_list;
		elem_list el;
		mg->associated_elements(el, vrt);

		// find the one that connects the 1d side
		// identifiable as the only element that has no constrained side
		elem_type* elem1d = NULL;
		size_t eSz = el.size();
		for (size_t e = 0; e < eSz; ++e)
		{
			typedef typename MultiGrid::traits<side_type>::secure_container side_list;
			side_list sl;
			mg->associated_elements(sl, el[e]);
			size_t sSz = sl.size();
			bool noConstr = true;
			for (size_t s = 0; s < sSz; ++s)
			{
				if (m_sh->get_subset_index(sl[s]) == constrainedSI)
				{
					noConstr = false;
					break;
				}
			}
			if (noConstr)
			{
				elem1d = el[e];
				break;
			}
		}

		if (!elem1d)
			continue;

		// no unification with ghosts
		if (aaElemInd[elem1d] == -1)
			continue;

		// get all other elems and create unification pairs
		UG_COND_THROW(eSz < 2, "Non-ghost 1d elem connected to 1d interface node, but no full-dim elem!");
		for (size_t otherInd = 0; otherInd < eSz; ++otherInd)
		{
			// no self-unification, no unification with ghosts
			if (el[otherInd] == elem1d || aaElemInd[el[otherInd]] == -1)
				continue;

			unificationPairs.push_back(std::make_pair(aaElemInd[elem1d] + localOffset,
													  aaElemInd[el[otherInd]] + localOffset));
		}
	}
}
#endif



template <typename TDomain>
void set_distro_adjuster(SmartPtr<TDomain> dom, SmartPtr<PNPDistroManager<TDomain> > adj)
{
#ifdef UG_PARALLEL
	DistributedGridManager& dgm = *dom->grid()->distributed_grid_manager();
	dgm.set_distro_adjuster(adj);
#endif
}



// explicit template specializations
#ifdef UG_DIM_1
	template class PNPDistroManager<Domain1d>;
	template void set_distro_adjuster<Domain1d>(SmartPtr<Domain1d> dom, SmartPtr<PNPDistroManager<Domain1d> > adj);
#endif
#ifdef UG_DIM_2
	template class PNPDistroManager<Domain2d>;
	template void set_distro_adjuster<Domain2d>(SmartPtr<Domain2d> dom, SmartPtr<PNPDistroManager<Domain2d> > adj);
#endif
#ifdef UG_DIM_3
	template class PNPDistroManager<Domain3d>;
	template void set_distro_adjuster<Domain3d>(SmartPtr<Domain3d> dom, SmartPtr<PNPDistroManager<Domain3d> > adj);
#endif




} // namespace nernst_planck
} // namespace ug4
