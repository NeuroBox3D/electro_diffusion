/*
 * intf_refMarkAdjuster.h
 *
 *  Created on: 29.04.2016
 *      Author: mbreit
 */

#ifndef UG__PLUGINS__NERNST_PLANCK__INTF_REFMARKADJUSTER_H_
#define UG__PLUGINS__NERNST_PLANCK__INTF_REFMARKADJUSTER_H_

#include "lib_grid/algorithms/refinement/ref_mark_adjuster_interface.h"
#include "interface1d_fv.h"


namespace ug {
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
		//std::vector<int> m_vIntfSI;
		std::vector<SmartPtr<IInterface1D> > m_vIntf;

		ConstSmartPtr<ISubsetHandler> m_ssh;
};


} // namespace nernst_planck
} // namespace ug

#endif // UG__PLUGINS__NERNST_PLANCK__INTF_REFMARKADJUSTER_H_
