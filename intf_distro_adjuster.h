/*
 * intf_distro_adjuster.h
 *
 *  Created on: 19.05.2016
 *      Author: mbreit
 */

#ifndef UG__PLUGINS__NERNST_PLANCK__INTF_DISTRO_ADJUSTER_H_
#define UG__PLUGINS__NERNST_PLANCK__INTF_DISTRO_ADJUSTER_H_


#include "lib_grid/parallelization/distro_adjuster.h"
#include "lib_grid/tools/subset_handler_interface.h"
#include "lib_grid/lib_grid_messages.h"
#include "common/util/smart_pointer.h"

#include "lib_grid/algorithms/graph/dual_graph.h"	// DualGraphNeighborCollector
#include "lib_grid/grid_objects/grid_dim_traits.h"

#include "interface1d_fv.h"

#include <vector>


namespace ug {
namespace nernst_planck {

// This class' collect_neighbor method can not be used in Parmetis partitioninf
// for grids without full-dim (w.r.t. TDomain::dim) elements a.t.m.
// To treat that case, you will need to declare this class not only template
// of TDomain, but also of refDim which must then be used in the public derivation
// statement from DualGraphNeighborCollector below.

template <typename TDomain>
class InterfaceDistroAdjuster
	: public DistroAdjuster,
	  public DualGraphNeighborCollector<typename grid_dim_traits<TDomain::dim>::grid_base_object>
{
	public:
	typedef typename grid_dim_traits<TDomain::dim>::grid_base_object elem_type;
	typedef typename grid_dim_traits<TDomain::dim-1>::grid_base_object side_type;

	public:
		InterfaceDistroAdjuster(SmartPtr<ApproximationSpace<TDomain> > approx);
		virtual ~InterfaceDistroAdjuster() {};

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
void set_distro_adjuster(SmartPtr<TDomain> dom, SmartPtr<InterfaceDistroAdjuster<TDomain> > adj);


} // namespace nernst_planck
} // namespace ug4


#endif // UG__PLUGINS__NERNST_PLANCK__INTF_DISTRO_ADJUSTER_H_
