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

#include "interface1d_fv.h"

#include <vector>


namespace ug {
namespace nernst_planck {

template <typename TDomain>
class InterfaceDistroAdjuster
	: public DistroAdjuster
{
	public:
		InterfaceDistroAdjuster(SmartPtr<ApproximationSpace<TDomain> > approx);
		virtual ~InterfaceDistroAdjuster() {};

		// inherited from DistroAdjuster
		virtual void adjust
		(
			MGSelector& sel,
			bool partitionForLocalProc,
			bool createVerticalInterfaces
		);

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
void set_distro_adjuster(SmartPtr<TDomain> dom, SmartPtr<InterfaceDistroAdjuster<TDomain> > adj);


} // namespace nernst_planck
} // namespace ug4


#endif // UG__PLUGINS__NERNST_PLANCK__INTF_DISTRO_ADJUSTER_H_
