/*
 * intf_distro_adjuster.h
 *
 *  Created on: 19.05.2016
 *      Author: mbreit
 */

#ifndef UG__PLUGINS__NERNST_PLANCK__INTF_DISTRO_ADJUSTER_H_
#define UG__PLUGINS__NERNST_PLANCK__INTF_DISTRO_ADJUSTER_H_

#include <vector>                                              // for vector, allocator

#include "common/util/smart_pointer.h"                         // for SmartPtr, ConstSmartPtr
#include "lib_grid/algorithms/graph/dual_graph.h"              // for DualGraphNeighborCollector
#include "lib_grid/grid_objects/grid_dim_traits.h"             // for grid_dim_traits
#include "lib_grid/parallelization/distro_adjuster.h"          // for DistroAdjuster
#include "lib_grid/tools/subset_handler_multi_grid.h"          // for MGSubsetHandler

#include "../Parmetis/src/partitioner_parmetis.h"              // for SideSubsetProtector

#include "interface1d_fv.h"                                    // for IInterface1D


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
 *     This feature (currently) only works if the whole domain is distributed from one proc.
 *
 * @todo: Make (2) and (3) work in general distribution conditions.
 */
template <typename TDomain>
class PNPDistroManager
: public DistroAdjuster,
  public parmetis::AnisotropyProtector<TDomain>,
  public DualGraphNeighborCollector<typename grid_dim_traits<TDomain::dim>::grid_base_object>
{
	public:
	typedef typename grid_dim_traits<TDomain::dim>::grid_base_object elem_type;
	typedef typename grid_dim_traits<TDomain::dim-1>::grid_base_object side_type;

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

		//void adjust_horizontal_interfaces(const GridMessage_Creation& msg);

		void add_interface(SmartPtr<IInterface1D> intf)
		{m_vIntf.push_back(intf);}

	private:
		SmartPtr<ApproximationSpace<TDomain> > m_approx;
		SmartPtr<TDomain> m_dom;
		ConstSmartPtr<MGSubsetHandler> m_sh;
		//MessageHub::SPCallbackId m_spGridCreationCallbackID;

		std::vector<SmartPtr<IInterface1D> > m_vIntf;
};


template <typename TDomain>
void set_distro_adjuster(SmartPtr<TDomain> dom, SmartPtr<PNPDistroManager<TDomain> > adj);


} // namespace nernst_planck
} // namespace ug4


#endif // UG__PLUGINS__NERNST_PLANCK__INTF_DISTRO_ADJUSTER_H_
