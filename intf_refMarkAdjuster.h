/*
 * intf_refMarkAdjuster.h
 *
 *  Created on: 29.04.2016
 *      Author: mbreit
 */

#ifndef UG__PLUGINS__NERNST_PLANCK__INTF_REFMARKADJUSTER_H_
#define UG__PLUGINS__NERNST_PLANCK__INTF_REFMARKADJUSTER_H_

#include <vector>                                             // for vector, allocator

#include "common/util/smart_pointer.h"                        // for SmartPtr, ConstSmartPtr
#include "interface1d_fv.h"                                   // for IInterface1D
#include "lib_grid/refinement/ref_mark_adjuster_interface.h"  // for IRefMarkAdjuster
#include "lib_grid/tools/subset_handler_interface.h"          // for ISubsetHandler


namespace ug {

// forward declarations
class Vertex;
class Edge;
class Face;
class Volume;
class IRefiner;

namespace nernst_planck {


class InterfaceRefMarkAdjuster
	: public IRefMarkAdjuster
{
	public:
		InterfaceRefMarkAdjuster() : IRefMarkAdjuster()	{}

		virtual ~InterfaceRefMarkAdjuster()	{}

		void add_interfaces(const std::vector<SmartPtr<IInterface1D> >& vI)
		{m_vIntf = vI;}

		void set_subset_handler(ConstSmartPtr<ISubsetHandler> ssh)
		{m_ssh = ssh;}

		virtual void ref_marks_changed
		(
			IRefiner& ref,
			const std::vector<Vertex*>& vrts,
			const std::vector<Edge*>& edges,
			const std::vector<Face*>& faces,
			const std::vector<Volume*>& vols
		);

	private:
		std::vector<SmartPtr<IInterface1D> > m_vIntf;

		ConstSmartPtr<ISubsetHandler> m_ssh;
};

void add_interface_ref_mark_adjuster(IRefiner* ref, SmartPtr<InterfaceRefMarkAdjuster> irma);


} // namespace nernst_planck
} // namespace ug

#endif // UG__PLUGINS__NERNST_PLANCK__INTF_REFMARKADJUSTER_H_
