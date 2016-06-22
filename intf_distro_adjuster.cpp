/*
 * intf_distro_adjuster.cpp
 *
 *  Created on: 19.05.2016
 *      Author: mbreit
 */

#include "lib_grid/parallelization/distribution.h" // for interface states
#include "intf_distro_adjuster.h"
#include "lib_grid/grid/grid_base_object_traits.h"
#include "lib_grid/algorithms/attachment_util.h"

namespace ug {
namespace nernst_planck {


template <typename TDomain>
InterfaceDistroAdjuster<TDomain>::InterfaceDistroAdjuster(SmartPtr<ApproximationSpace<TDomain> > approx)
: m_approx(approx), m_dom(approx->domain()), m_sh(m_dom->subset_handler())
{
	m_spGridCreationCallbackID = m_dom->grid()->message_hub()->register_class_callback(this,
		&ug::nernst_planck::InterfaceDistroAdjuster<TDomain>::adjust_horizontal_interfaces);
}


template <typename TDomain>
void InterfaceDistroAdjuster<TDomain>::adjust
(
	MGSelector& sel,
	bool partitionForLocalProc,
	bool createVerticalInterfaces
)
{
	const MultiGrid& mg = *sel.multi_grid();
	const DistributedGridManager& dgm = *mg.distributed_grid_manager();
	const SurfaceView& sv = *m_approx->surface_view();

	// check whether process is to have any constrained nodes
	//typedef typename MGSelector::traits<Vertex>::const_iterator sel_it_type;
	typedef typename SurfaceView::traits<Vertex>::const_iterator sv_it_type;
	//typedef typename geometry_traits<Vertex>::const_iterator sh_it_type;

	size_t sz = m_vIntf.size();
	for (size_t i = 0; i < sz; ++i)
	{
		SmartPtr<IInterface1D> intf = m_vIntf[i];

		int constrdSI = intf->constrained_subset_index();

		size_t nSelected = 0;
		size_t nUnselected = 0;

		// iterate surface constrained
		GridLevel gl;
		sv_it_type it = sv.begin<Vertex>(constrdSI, gl, SurfaceView::ALL_BUT_SHADOW_COPY);
		sv_it_type it_end = sv.end<Vertex>(constrdSI, gl, SurfaceView::ALL_BUT_SHADOW_COPY);

		for (; it != it_end; ++it)
		{
			if (sel.get_selection_status(*it) & (IS_NORMAL | IS_VSLAVE))
				++nSelected;
			else
				++nUnselected;
		}

		// if no surface constrained exist here: nothing to do
		if (! (nSelected + nUnselected)) continue;

//if (partitionForLocalProc)
//{
//	UG_LOGN("Interface " << i << ":");
//	UG_LOGN("    selected: " << nSelected);
//	UG_LOGN("  unselected: " << nUnselected);
//}

		// find surface interface nodes
		int in1dSI = intf->intf_node_1d_subset_index();
		int inhdSI = intf->intf_node_hd_subset_index();

		// select top-most interface nodes (surface)
		Vertex* i1d = NULL;
		Vertex* ihd = NULL;

		sv_it_type it1d = sv.begin<Vertex>(in1dSI, gl, SurfaceView::ALL_BUT_SHADOW_COPY);
		sv_it_type it1d_end = sv.end<Vertex>(in1dSI, gl, SurfaceView::ALL_BUT_SHADOW_COPY);

		sv_it_type ithd = sv.begin<Vertex>(inhdSI, gl, SurfaceView::ALL_BUT_SHADOW_COPY);
		sv_it_type ithd_end = sv.end<Vertex>(inhdSI, gl, SurfaceView::ALL_BUT_SHADOW_COPY);

		if (it1d != it1d_end)
		{
			i1d = *it1d;

			UG_COND_THROW(ithd == ithd_end,
				"Found 1d interface node, but not high-d interface node on the same level and proc.");

			ihd = *ithd;
		}

		UG_COND_THROW(!i1d || !ihd, "No interface nodes on proc distributing constrained nodes.");


//if (partitionForLocalProc && nUnselected)
//{
//	UG_LOGN("intf node 1d state = " << (size_t) dgm.get_status(i1d));
//	UG_LOGN("intf node hd state = " << (size_t) dgm.get_status(ihd));
//}

		// select and decide on v-interface status
		if (!createVerticalInterfaces)
		{
			if (!sel.is_selected(i1d)) sel.select(i1d, IS_NORMAL);
			if (!sel.is_selected(ihd)) sel.select(ihd, IS_NORMAL);

			continue;
		}

		// If we are processing selection for the local proc
		// and if any of the constrained vertices is unselected
		// (which means it is going somewhere else):
		// select Interface nodes as vMaster if not marked v-slave.
		// If they already are v-slaves, leave them as such.
		if (partitionForLocalProc)
		{
			if (nUnselected)
			{
				if (!(sel.get_selection_status(i1d) & IS_VSLAVE))
				{
					//if (mg.get_level(i1d) == 0)
					if (!dgm.contains_status(i1d, ES_V_SLAVE))
					{
//UG_LOGN("Setting intf 1d node on proc " << pcl::ProcRank() << " to v-master.");
						sel.select(i1d, IS_VMASTER);
					}
					//else if (selectedExist[i])
						//sel.select(i1d, IS_VSLAVE);
				}

				if (!(sel.get_selection_status(ihd) & IS_VSLAVE))
				{
					//if (mg.get_level(ihd) == 0)
					if (!dgm.contains_status(ihd, ES_V_SLAVE))
						sel.select(ihd, IS_VMASTER);
					//else if (selectedExist[i])
						//sel.select(ihd, IS_VSLAVE);
				}
			}

			// also ensure that vmaster status is kept if node is already vmaster
			else
			{
				if (dgm.contains_status(i1d, ES_V_MASTER)
					&& !(sel.get_selection_status(i1d) & IS_VSLAVE))
				{
//UG_LOGN("Setting intf 1d node on proc " << pcl::ProcRank() << " to v-master.");
					sel.select(i1d, sel.get_selection_status(i1d) | IS_VMASTER);
				}

				if (dgm.contains_status(ihd, ES_V_MASTER)
					&& !(sel.get_selection_status(ihd) & IS_VSLAVE))
					sel.select(ihd, sel.get_selection_status(ihd) | IS_VMASTER);
			}
		}

		// If we are processing selection for another proc
		// and if any of the constrained vertices is selected:
		// Select Interface nodes as vSlave.
		//else
		{
			if (nSelected)
			{
				if (!(sel.get_selection_status(i1d) & (IS_VMASTER | IS_NORMAL)))
					//if (mg.get_level(i1d) == 0 && !dgm.contains_status(i1d, ES_V_MASTER))
				{
//if (partitionForLocalProc)
//	{UG_LOGN("Setting intf 1d node on proc " << pcl::ProcRank() << " to v-slave.");}
//else
//	{UG_LOGN("Setting intf 1d node on some other proc to v-slave.");}
					sel.select(i1d, IS_VSLAVE);
				}
				if (!(sel.get_selection_status(ihd) & (IS_VMASTER | IS_NORMAL)))
					//if (mg.get_level(ihd) == 0 && !dgm.contains_status(ihd, ES_V_MASTER))
					sel.select(ihd, IS_VSLAVE);
			}
		}
	}
}


template <typename TDomain>
void InterfaceDistroAdjuster<TDomain>::adjust_horizontal_interfaces(const GridMessage_Creation& msg)
{
	if (msg.msg() != GMCT_CREATION_STOPS) return;

	// We only have to make sure that the processes having constrained elements,
	// but where the interface nodes are pure V-MASTER, build h-interfaces for the interface nodes.
	// (This will only be proc 0 if the interface nodes are not refined.)

	typedef typename geometry_traits<Vertex>::const_iterator sh_it_type;
	typedef typename SurfaceView::traits<Vertex>::const_iterator sv_it_type;
	typedef typename GridLayoutMap::Types<Vertex>::Interface IntfType;

	MultiGrid& mg = *m_dom->grid();
	DistributedGridManager& dgm = *mg.distributed_grid_manager();
	GridLayoutMap& glm = dgm.grid_layout_map();
	const SurfaceView& sv = *m_approx->surface_view();

	dgm.enable_interface_management(false);

	CompareByAttachment<Vertex, AGeomObjID> gidCmp(*m_dom->grid(), aGeomObjID); // needed for interface entry search

	//size_t nLvl = m_dom->grid()->num_levels();
	int nProcs = pcl::NumProcs();
	int locRank = pcl::ProcRank();

	size_t sz = m_vIntf.size();
	for (size_t i = 0; i < sz; ++i)
	{
		// find all non-vmaster procs that share a copy of the interface nodes
		SmartPtr<IInterface1D> intf = m_vIntf[i];

		int i1dSI = intf->intf_node_1d_subset_index();
		int ihdSI = intf->intf_node_hd_subset_index();

		int rankIN1d = nProcs;
		int rankINhd = nProcs;

		Vertex* vrtIN1d = NULL;
		Vertex* vrtINhd = NULL;
		int lvl1d = -1;
		int lvlhd = -1;

		//GridLevel surfWithGhosts = GridLevel(GridLevel::TOP, GridLevel::SURFACE, true); // also search for ghosts!

		// find surface level interface nodes
		uint nLvl = (uint) mg.num_levels();
		for (uint lvl = nLvl-1; lvl < nLvl; --lvl)
		{
			sh_it_type it = m_sh->begin<Vertex>(i1dSI, lvl);
			sh_it_type it_end = m_sh->end<Vertex>(i1dSI, lvl);

			if (it != it_end)
			{
				if ((dgm.get_status(*it) & (ES_H_MASTER | ES_H_SLAVE | ES_V_SLAVE))
					|| (dgm.get_status(*it) == ES_NONE))
					rankIN1d = locRank;

				vrtIN1d = *it;
				lvl1d = lvl;

				break;
			}
		}
		for (uint lvl = nLvl-1; lvl < nLvl; --lvl)
		{
			sh_it_type it = m_sh->begin<Vertex>(ihdSI, lvl);
			sh_it_type it_end = m_sh->end<Vertex>(ihdSI, lvl);

			if (it != it_end)
			{
				if ((dgm.get_status(*it) & (ES_H_MASTER | ES_H_SLAVE | ES_V_SLAVE))
					|| (dgm.get_status(*it) == ES_NONE))
					rankINhd = locRank;

				vrtINhd = *it;
				lvlhd = lvl;

				break;
			}
		}

		/*
		sv_it_type it = sv.begin<Vertex>(i1dSI, surfWithGhosts, SurfaceView::ALL_BUT_SHADOW_COPY);
		sv_it_type it_end = sv.end<Vertex>(i1dSI, surfWithGhosts, SurfaceView::ALL_BUT_SHADOW_COPY);

		if (it != it_end)
		{
			if (dgm.get_status(*it) & (ES_H_MASTER | ES_H_SLAVE | ES_V_SLAVE))
				rankIN1d = locRank;

			vrtIN1d = *it;
			lvl1d = m_dom->grid()->get_level(vrtIN1d);
		}

		it = sv.begin<Vertex>(ihdSI, surfWithGhosts, SurfaceView::ALL_BUT_SHADOW_COPY);
		it_end = sv.end<Vertex>(ihdSI, surfWithGhosts, SurfaceView::ALL_BUT_SHADOW_COPY);

		if (it != it_end)
		{
			if (dgm.get_status(*it) & (ES_H_MASTER | ES_H_SLAVE | ES_V_SLAVE))
				rankINhd = locRank;

			vrtINhd = *it;
			lvlhd = m_dom->grid()->get_level(vrtINhd);
		}
		*/

		// if surface interface nodes are present, find out whether there are also constrained nodes
		bool constrdExist = false;
		if (vrtIN1d && vrtINhd)
		{
			int constrdSI = intf->constrained_subset_index();
			GridLevel surfNoGhosts = GridLevel(GridLevel::TOP, GridLevel::SURFACE, false);
			sv_it_type it = sv.begin<Vertex>(constrdSI, surfNoGhosts, SurfaceView::ALL_BUT_SHADOW_COPY);
			sv_it_type it_end = sv.end<Vertex>(constrdSI, surfNoGhosts, SurfaceView::ALL_BUT_SHADOW_COPY);

			for (; it != it_end; ++it)
			{
				if (dgm.get_status(*it) & (IS_NORMAL | IS_VSLAVE))
				{
					constrdExist = true;
					break;
				}
			}

//UG_LOGN("intf node 1d state = " << (size_t) dgm.get_status(vrtIN1d));
//UG_LOGN("intf node hd state = " << (size_t) dgm.get_status(vrtINhd));
		}


		// now communicate who is to be h-master
		// (if h-master already exists, they coincide by construction)
		pcl::ProcessCommunicator comm;
		int rankIN1dGlob;
		comm.allreduce(&rankIN1d, &rankIN1dGlob, 1, PCL_RO_MIN);

		int rankINhdGlob;
		comm.allreduce(&rankINhd, &rankINhdGlob, 1, PCL_RO_MIN);

//UG_LOGN("rank_1d: " << rankIN1dGlob);

		// should not happen, but who knows
		UG_COND_THROW(rankIN1dGlob >= nProcs, "No admissible hMaster for 1d interface node.");
		//if (rankIN1dGlob >= nProcs) continue;

		// now send ranks of v-masters with constrained nodes to h-master
		int* recBuff = NULL;
		if (locRank == rankIN1dGlob) recBuff = new int[nProcs];

		if (vrtIN1d && rankIN1d == nProcs && constrdExist) rankIN1d = locRank;
		else rankIN1d = nProcs;
		comm.gather(&rankIN1d, 1, PCL_DT_INT, recBuff, 1, PCL_DT_INT, rankIN1dGlob);

		// on h-master: add new slaves to h-master interface
		if (locRank == rankIN1dGlob)
		{
			for (size_t j = 0; j < (size_t) nProcs; ++j)
			{
				if (recBuff[j] != nProcs && j != (size_t) locRank)
				{
					IntfType& ii = glm.get_layout<Vertex>(INT_H_MASTER).interface(j, lvl1d);

					// push back elem if not already in interface
					IntfType::iterator it = ii.find_insert_pos_sorted(vrtIN1d, gidCmp);
					if (it == ii.end() || gidCmp(vrtIN1d, ii.get_element(it)))
					{
//UG_LOG("Inserting 1d intf vertex " << vrtIN1d << " as master for proc " << j << ".\n");
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
				IntfType& ii = glm.get_layout<Vertex>(INT_H_SLAVE).interface(rankIN1dGlob, lvl1d);

				// push back elem if not already in interface
				IntfType::iterator it = ii.find_insert_pos_sorted(vrtIN1d, gidCmp);
				if (it == ii.end() || gidCmp(vrtIN1d, ii.get_element(it)))
				{
//UG_LOG("Inserting 1d intf vertex " << vrtIN1d << " as slave for proc " << rankIN1dGlob << ".\n");
					ii.insert(vrtIN1d, it);

//for (IntfType::iterator iter = ii.begin(); iter != ii.end(); ++iter)
//UG_LOG(ii.get_element(iter) << "   " << ii.get_local_id(iter) << "\n");
				}
			}
		}

		// delete receive buffer on h-master
		if (locRank == rankIN1dGlob) delete[] recBuff;


		// same again for high-dim intf node

		// should not happen, but who knows
		UG_COND_THROW(rankINhdGlob >= nProcs, "No admissible hMaster for high-dim interface node.");

		// now send ranks of v-masters with constrained nodes to h-master
		recBuff = NULL;
		if (locRank == rankINhdGlob) recBuff = new int[nProcs];
		if (vrtINhd && rankINhd == nProcs && constrdExist) rankINhd = locRank;
		else rankINhd = nProcs;

		comm.gather(&rankINhd, 1, PCL_DT_INT, recBuff, 1, PCL_DT_INT, rankINhdGlob);

		// on h-master: add slaves to h-master interface
		if (locRank == rankINhdGlob)
		{
			for (size_t j = 0; j < (size_t) nProcs; ++j)
			{
				if (recBuff[j] != nProcs && j != (size_t) locRank)
				{
					IntfType& ii = glm.get_layout<Vertex>(INT_H_MASTER).interface(j, lvlhd);

					// push back elem if not already in interface
					IntfType::iterator it = ii.find_insert_pos_sorted(vrtINhd, gidCmp);
					if (it == ii.end() || gidCmp(vrtINhd, ii.get_element(it)))
					{
//UG_LOG("Inserting hd intf vertex " << vrtINhd << " as master for proc " << j << ".\n");
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
				IntfType& ii = glm.get_layout<Vertex>(INT_H_SLAVE).interface(rankINhdGlob, lvlhd);

				// push back elem if not already in interface
				IntfType::iterator it = ii.find_insert_pos_sorted(vrtINhd, gidCmp);
				if (it == ii.end() || gidCmp(vrtINhd, ii.get_element(it)))
				{
//UG_LOG("Inserting hd intf vertex " << vrtINhd << " as slave for proc " << rankINhdGlob << ".\n");
					ii.insert(vrtINhd, it);

//for (IntfType::iterator iter = ii.begin(); iter != ii.end(); ++iter)
//UG_LOG(ii.get_element(iter) << "   " << ii.get_local_id(iter) << "\n");
				}
			}
		}

		// delete receive buffer on hmaster
		if (locRank == rankINhdGlob) delete[] recBuff;
	}

	glm.remove_empty_interfaces();
	dgm.enable_interface_management(true);
	dgm.grid_layouts_changed(false);
}


template <typename TDomain>
void set_distro_adjuster(SmartPtr<TDomain> dom, SmartPtr<InterfaceDistroAdjuster<TDomain> > adj)
{
#ifdef UG_PARALLEL
	DistributedGridManager& dgm = *dom->grid()->distributed_grid_manager();
	dgm.set_distro_adjuster(adj);
#endif
}



// explicit template specializations
#ifdef UG_DIM_1
	template class InterfaceDistroAdjuster<Domain1d>;
	template void set_distro_adjuster<Domain1d>(SmartPtr<Domain1d> dom, SmartPtr<InterfaceDistroAdjuster<Domain1d> > adj);
#endif
#ifdef UG_DIM_2
	template class InterfaceDistroAdjuster<Domain2d>;
	template void set_distro_adjuster<Domain2d>(SmartPtr<Domain2d> dom, SmartPtr<InterfaceDistroAdjuster<Domain2d> > adj);
#endif
#ifdef UG_DIM_3
	template class InterfaceDistroAdjuster<Domain3d>;
	template void set_distro_adjuster<Domain3d>(SmartPtr<Domain3d> dom, SmartPtr<InterfaceDistroAdjuster<Domain3d> > adj);
#endif




} // namespace nernst_planck
} // namespace ug4
