/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Markus Breit
 * Creation date: 2016-05-19
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


#ifndef UG__PLUGINS__NERNST_PLANCK__INTERFACE__INTF_DISTRO_ADJUSTER_H
#define UG__PLUGINS__NERNST_PLANCK__INTERFACE__INTF_DISTRO_ADJUSTER_H

#include <vector>                                              // for vector, allocator

#include "common/util/smart_pointer.h"                         // for SmartPtr, ConstSmartPtr
#include "lib_grid/algorithms/graph/dual_graph.h"              // for DualGraphNeighborCollector
#include "lib_grid/grid_objects/grid_dim_traits.h"             // for grid_dim_traits
#include "lib_grid/parallelization/distro_adjuster.h"          // for DistroAdjuster
#include "lib_grid/tools/subset_handler_multi_grid.h"          // for MGSubsetHandler

#include "../../Parmetis/src/anisotropy_unificator.h"             // for AnisotropyUnificator

#include "interface1d_fv.h"                                    // for IInterface1D

#include "np_config.h"  // contains cmake-configured defines,
                     // esp. whether Parmetis is built along

namespace ug {

// forward declarations
class MGSelector;
template <typename TDomain> class ApproximationSpace;

namespace nernst_planck {
/**
 * @brief This class manages all PNP-relevant distribution issues.
 *
 * (1) Being derived from DistroAdjuster, it is able to manipulate the domain partitioning
 *     after the distribution method has been applied.
 *     This can be used to realize the parallel 1d/3d interfaces where interface nodes
 *     need to be copied to any proc that holds a part of the interface.
 * (2) Being derived from DualGraphNeighborCollector, it is able to influence the way
 *     connections are recorded for the dual graph in the ParMetis distribution method.
 *     This can be used to correctly record the connections from the first element of the
 *     1d extensions to the border elements of the 3d side.
 *     This feature (currently) only works if the whole domain is distributed from one proc.
 * (3) Being derived from AnisotropyProtector, it is able to manipulate the way the
 *     dual graph for ParMetis partitioning is created. This can be used to prevent the
 *     partitioning from having borders along charged membranes.
 *
 * @todo: Make (2) and (3) work in general distribution conditions.
 */
template <typename TDomain>
class PNPDistroManager
: public DistroAdjuster,
#ifdef NPParmetis
  public parmetis::AnisotropyUnificator<TDomain, typename grid_dim_traits<TDomain::dim>::grid_base_object>,
#endif
  public DualGraphNeighborCollector<typename grid_dim_traits<TDomain::dim>::grid_base_object>
{
	public:
		typedef typename grid_dim_traits<TDomain::dim>::grid_base_object elem_type;
		typedef typename grid_dim_traits<TDomain::dim-1>::grid_base_object side_type;
		typedef Attachment<int> AElemIndex;
		typedef Attachment<std::vector<int> > AElemIndices;

	public:
		PNPDistroManager(SmartPtr<ApproximationSpace<TDomain> > approx);
		virtual ~PNPDistroManager() {};

		// inherited from DualGraphNeighborAdjuster
		void collect_neighbors(std::vector<elem_type*>& neighborsOut, elem_type* elem);

		// inherited from DistroAdjuster
		virtual void adjust
		(
			MGSelector& sel,
			bool partitionForLocalProc,
			bool createVerticalInterfaces
		);

#ifdef NPParmetis
		// inherited from AnisotropyUnificator
		void unify
		(
			MultiGrid* mg,
			int lvl,
			int localOffset,
			const Grid::AttachmentAccessor<elem_type, AElemIndex>& aaElemInd,
			const Grid::AttachmentAccessor<side_type, AElemIndices>& aaSideElemInd,
			std::vector<std::pair<int, int> >& unificationPairs
		) const;
#endif

		void adjust_horizontal_interfaces(const GridMessage_Creation& msg);

		void add_interface(SmartPtr<IInterface1D> intf)
		{m_vIntf.push_back(intf);}

	private:
		SmartPtr<ApproximationSpace<TDomain> > m_approx;
		SmartPtr<TDomain> m_dom;
		ConstSmartPtr<MGSubsetHandler> m_sh;
		MessageHub::SPCallbackId m_spGridCreationCallbackID;

		std::vector<SmartPtr<IInterface1D> > m_vIntf;
};


template <typename TDomain>
void set_distro_adjuster(SmartPtr<TDomain> dom, SmartPtr<PNPDistroManager<TDomain> > adj);


} // namespace nernst_planck
} // namespace ug4


#endif // UG__PLUGINS__NERNST_PLANCK__INTERFACE__INTF_DISTRO_ADJUSTER_H
