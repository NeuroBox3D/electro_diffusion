/*
 * intf_distro_adjuster.cpp
 *
 *  Created on: 19.05.2016
 *      Author: mbreit
 */

#include "lib_grid/parallelization/distribution.h" // for interface states
#include "intf_distro_adjuster.h"
#include "lib_grid/grid/grid_base_object_traits.h"
#include "lib_grid/grid/neighborhood.h"				// CollectNeighbors
#include "lib_grid/algorithms/attachment_util.h"

#include <set>

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
: parmetis::AnisotropyProtector<TDomain>(approx->domain()),
  m_approx(approx), m_dom(approx->domain()), m_sh(m_dom->subset_handler())
{
	// set this class as its own neighbor collector
	SmartPtr<DualGraphNeighborCollector<elem_type> > spDGNC = make_sp(this);
	++*spDGNC.refcount_ptr(); // we do not want to accidentally destroy ourselves
	this->set_neighbor_collector(spDGNC);

	// no longer necessary!
	//m_spGridCreationCallbackID = m_dom->grid()->message_hub()->register_class_callback(this,
	//	&ug::nernst_planck::InterfaceDistroAdjuster<TDomain>::adjust_horizontal_interfaces);
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

		size_t intfSrfLvl;
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

	// 2a. add level interface nodes for level constrained
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
/*
	// 2b. add surface interface nodes for surface constrained
		if (nSelected[nLvl])
		{
			sv_it_type it1d = sv.begin<Vertex>(in1dSI, gl, SurfaceView::ALL_BUT_SHADOW_COPY);
			sv_it_type it1d_end = sv.end<Vertex>(in1dSI, gl, SurfaceView::ALL_BUT_SHADOW_COPY);

			sv_it_type ithd = sv.begin<Vertex>(inhdSI, gl, SurfaceView::ALL_BUT_SHADOW_COPY);
			sv_it_type ithd_end = sv.end<Vertex>(inhdSI, gl, SurfaceView::ALL_BUT_SHADOW_COPY);

			// interface nodes MUST be present in surface
			if (it1d != it1d_end)
			{
				sel.select(*it1d, sel.get_selection_status(*it1d) | IS_NORMAL);
				levelVertices[mg.get_level(*it1d)].insert(*it1d);
			}
			if (ithd != ithd_end)
			{
				sel.select(*ithd, sel.get_selection_status(*ithd) | IS_NORMAL);
				levelVertices[mg.get_level(*ithd)].insert(*ithd);
			}
		}
*/
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

	// It seems to be enough to do this on root elements.

	// 4. assign vertical master and slave states:
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


#if 0
template <typename TDomain>
void PNPDistroManager<TDomain>::adjust_horizontal_interfaces(const GridMessage_Creation& msg)
{
	if (msg.msg() != GMCT_CREATION_STOPS) return;

	// We only have to make sure that the processes having constrained elements,
	// but where the interface nodes are pure V-MASTER, build h-interfaces for the interface nodes.
	// ALSO: It is possible that no constrained exists but an interface node is connected to a full-d
	// element which is NOT v-master. Another h-interface is needed in this situation.

	typedef typename geometry_traits<Vertex>::const_iterator sh_it_type;
	typedef typename SurfaceView::traits<Vertex>::const_iterator sv_it_type;
	typedef typename GridLayoutMap::Types<Vertex>::Interface IntfType;

	MultiGrid& mg = *m_dom->grid();
	DistributedGridManager& dgm = *mg.distributed_grid_manager();
	GridLayoutMap& glm = dgm.grid_layout_map();
	//const SurfaceView& sv = *m_approx->surface_view();

	size_t nLvl = mg.num_levels();

	dgm.enable_interface_management(false);

	CompareByAttachment<Vertex, AGeomObjID> gidCmp(*m_dom->grid(), aGeomObjID); // needed for interface entry search

	//size_t nLvl = m_dom->grid()->num_levels();
	int nProcs = pcl::NumProcs();
	int locRank = pcl::ProcRank();

	size_t sz = m_vIntf.size();
	for (size_t i = 0; i < sz; ++i)
	{
		SmartPtr<IInterface1D> intf = m_vIntf[i];

		int constrdSI = intf->constrained_subset_index();
		int i1dSI = intf->intf_node_1d_subset_index();
		int ihdSI = intf->intf_node_hd_subset_index();

		for (size_t lvl = 0; lvl < nLvl; ++lvl)
		{
			// find all non-vmaster procs that share a copy of the interface nodes
			int rankIN1d = nProcs;
			int rankINhd = nProcs;

			Vertex* vrtIN1d = NULL;
			Vertex* vrtINhd = NULL;

			sh_it_type it = m_sh->begin<Vertex>(i1dSI, lvl);
			sh_it_type it_end = m_sh->end<Vertex>(i1dSI, lvl);

			if (it != it_end)
			{
				if ((dgm.get_status(*it) & (ES_H_MASTER | ES_H_SLAVE | ES_V_SLAVE))
					|| (dgm.get_status(*it) == ES_NONE))
					rankIN1d = locRank;

				vrtIN1d = *it;
				//UG_LOGN("Intf " << i << ", lvl " << lvl << ": 1d inode state: " << dgm.get_status(*it));
			}
			bool isVMaster1d =  vrtIN1d && rankIN1d == nProcs;

			it = m_sh->begin<Vertex>(ihdSI, lvl);
			it_end = m_sh->end<Vertex>(ihdSI, lvl);

			if (it != it_end)
			{
				if ((dgm.get_status(*it) & (ES_H_MASTER | ES_H_SLAVE | ES_V_SLAVE))
					|| (dgm.get_status(*it) == ES_NONE))
					rankINhd = locRank;

				vrtINhd = *it;

				//GridObject* geomObj = mg.get_parent(*it);
				//UG_LOGN("Intf " << i << ", lvl " << lvl << ": hd inode state: "
				//		<< (unsigned int) dgm.get_status(*it) << ",   father: " << geomObj << ".");
			}
			bool isVMasterhd =  vrtINhd && rankINhd == nProcs;

			// now communicate who is to be h-master
			// (if h-master already exists, they coincide by construction)
			pcl::ProcessCommunicator comm;
			int rankIN1dGlob;
			comm.allreduce(&rankIN1d, &rankIN1dGlob, 1, PCL_RO_MIN);

			int rankINhdGlob;
			comm.allreduce(&rankINhd, &rankINhdGlob, 1, PCL_RO_MIN);

//UG_LOGN("rank_1d: " << rankIN1dGlob);

			// find out whether there are any constrained nodes
			bool constrdExist = false;
			it = m_sh->begin<Vertex>(constrdSI, lvl);
			it_end = m_sh->end<Vertex>(constrdSI, lvl);

			for (; it != it_end; ++it)
			{
				if ((dgm.get_status(*it) & ES_V_SLAVE) || (dgm.get_status(*it) == ES_NONE))
				{
					constrdExist = true;
					break;
				}
			}

			if (rankIN1dGlob < nProcs)
			{
				/*
				// find out whether 1d intf node is connected to non-v-slave element
				bool non_vmaster_connected = false;
				if (vrtIN1d)
				{
					typedef typename MultiGrid::traits<Edge>::secure_container edge_list;
					edge_list el;
					mg.associated_elements(el, vrtIN1d);
					size_t elsz = el.size();
					for (size_t ed = 0; ed < elsz; ++ed)
					{
						if (!dgm.contains_status(el[ed], ES_V_MASTER))
						{
							non_vmaster_connected = true;
							break;
						}
					}
				}
				*/

				// now send ranks of v-masters with constrained nodes to h-master
				int* recBuff = NULL;
				if (locRank == rankIN1dGlob) recBuff = new int[nProcs];

				if (isVMaster1d && (constrdExist /*|| non_vmaster_connected*/)) rankIN1d = locRank;
				else rankIN1d = nProcs;
				comm.gather(&rankIN1d, 1, PCL_DT_INT, recBuff, 1, PCL_DT_INT, rankIN1dGlob);

				// on h-master: add new slaves to h-master interface
				if (locRank == rankIN1dGlob)
				{
					for (size_t j = 0; j < (size_t) nProcs; ++j)
					{
						if (recBuff[j] != nProcs && j != (size_t) locRank)
						{
							IntfType& ii = glm.get_layout<Vertex>(INT_H_MASTER).interface(j, lvl);


							// push back elem if not already in interface
							IntfType::iterator it = ii.find_insert_pos_sorted(vrtIN1d, gidCmp);
							if (it == ii.end() || gidCmp(vrtIN1d, ii.get_element(it)))
							{
			//UG_LOG("Inserting 1d intf vertex " << vrtIN1d << " on level " << lvl << " as master for proc " << j << ".\n");
								ii.insert(vrtIN1d, it);

			//for (IntfType::iterator iter = ii.begin(); iter != ii.end(); ++iter)
			//UG_LOG(ii.get_element(iter) << "   " << ii.get_local_id(iter) << "\n");
							}
						}
					}
				}
				// on slaves: add master to slave interface
				else
				{
					if (rankIN1d != nProcs)
					{
						IntfType& ii = glm.get_layout<Vertex>(INT_H_SLAVE).interface(rankIN1dGlob, lvl);

						// push back elem if not already in interface
						IntfType::iterator it = ii.find_insert_pos_sorted(vrtIN1d, gidCmp);
						if (it == ii.end() || gidCmp(vrtIN1d, ii.get_element(it)))
						{
			//UG_LOG("Inserting 1d intf vertex " << vrtIN1d << " on level " << lvl << " as slave for proc " << rankIN1dGlob << ".\n");
							ii.insert(vrtIN1d, it);

			//for (IntfType::iterator iter = ii.begin(); iter != ii.end(); ++iter)
			//UG_LOG(ii.get_element(iter) << "   " << ii.get_local_id(iter) << "\n");
						}
					}
				}

				// delete receive buffer on h-master
				if (locRank == rankIN1dGlob) delete[] recBuff;
			}

			// same again for high-dim intf node
			if (rankIN1dGlob < nProcs)
			{
				/*
				// find out whether 1d intf node is connected to non-v-slave element
				bool non_vmaster_connected = false;
				if (vrtINhd)
				{
					typedef typename MultiGrid::traits<Edge>::secure_container edge_list;
					edge_list el;
					mg.associated_elements(el, vrtINhd);
					size_t elsz = el.size();
					for (size_t ed = 0; ed < elsz; ++ed)
					{
						if (!dgm.contains_status(el[ed], ES_V_MASTER))
						{
							non_vmaster_connected = true;
							break;
						}
					}
				}
				*/

				// now send ranks of v-masters with constrained nodes to h-master
				int* recBuff = NULL;
				if (locRank == rankINhdGlob) recBuff = new int[nProcs];
				if (isVMasterhd && (constrdExist /*|| non_vmaster_connected*/)) rankINhd = locRank;
				else rankINhd = nProcs;

				comm.gather(&rankINhd, 1, PCL_DT_INT, recBuff, 1, PCL_DT_INT, rankINhdGlob);

				// on h-master: add slaves to h-master interface
				if (locRank == rankINhdGlob)
				{
					for (size_t j = 0; j < (size_t) nProcs; ++j)
					{
						if (recBuff[j] != nProcs && j != (size_t) locRank)
						{
							IntfType& ii = glm.get_layout<Vertex>(INT_H_MASTER).interface(j, lvl);

							// push back elem if not already in interface
							IntfType::iterator it = ii.find_insert_pos_sorted(vrtINhd, gidCmp);
							if (it == ii.end() || gidCmp(vrtINhd, ii.get_element(it)))
							{
			//UG_LOG("Inserting hd intf vertex " << vrtINhd << " on level " << lvl << " as master for proc " << j << ".\n");
								ii.insert(vrtINhd, it);

			//for (IntfType::iterator iter = ii.begin(); iter != ii.end(); ++iter)
			//UG_LOG(ii.get_element(iter) << "   " << ii.get_local_id(iter) << "\n");
							}
						}
					}
				}
				// on slaves: add master to slave interface
				else
				{
					if (rankINhd != nProcs)
					{
						IntfType& ii = glm.get_layout<Vertex>(INT_H_SLAVE).interface(rankINhdGlob, lvl);

						// push back elem if not already in interface
						IntfType::iterator it = ii.find_insert_pos_sorted(vrtINhd, gidCmp);
						if (it == ii.end() || gidCmp(vrtINhd, ii.get_element(it)))
						{
			//UG_LOG("Inserting hd intf vertex " << vrtINhd << " on level " << lvl << " as slave for proc " << rankINhdGlob << ".\n");
							ii.insert(vrtINhd, it);

			//for (IntfType::iterator iter = ii.begin(); iter != ii.end(); ++iter)
			//UG_LOG(ii.get_element(iter) << "   " << ii.get_local_id(iter) << "\n");
						}
					}
				}

				// delete receive buffer on hmaster
				if (locRank == rankINhdGlob) delete[] recBuff;
			}
		}
	}

	glm.remove_empty_interfaces();
	dgm.enable_interface_management(true);
	dgm.grid_layouts_changed(false);
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
