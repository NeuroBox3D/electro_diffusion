/*
 * intf_refMarkAdjuster.h
 *
 *  Created on: 29.04.2016
 *      Author: mbreit
 */

#ifndef UG__PLUGINS__NERNST_PLANCK__EXTENSION_REFMARKADJUSTER_H_
#define UG__PLUGINS__NERNST_PLANCK__EXTENSION_REFMARKADJUSTER_H_

#include "lib_grid/refinement/ref_mark_adjuster_interface.h"
#include "lib_grid/refinement/hanging_node_refiner_multi_grid.h"
#include "interface1d_fv.h"


namespace ug {
namespace nernst_planck {


template <typename TDomain>
class ExtensionRefMarkAdjuster
	: public IRefMarkAdjuster
{
	public:
		static const int worldDim = TDomain::dim;

	public:
		/// constructor
		ExtensionRefMarkAdjuster
		(
			SmartPtr<TDomain> dom,
			std::vector<number> dir,
			const std::string useless
		);

		/// destructor
		virtual ~ExtensionRefMarkAdjuster()	{}

		virtual void ref_marks_changed
		(
			IRefiner& ref,
			const std::vector<Vertex*>& vrts,
			const std::vector<Edge*>& edges,
			const std::vector<Face*>& faces,
			const std::vector<Volume*>& vols
		);

	protected:
		void change_edge_marks(IRefiner& ref, Edge* e);
		void change_face_marks(IRefiner& ref, Face* f);
		void change_vol_marks(IRefiner& ref, Volume* v);

	private:
		int m_si;
		ConstSmartPtr<ISubsetHandler> m_ssh;
		typename TDomain::position_accessor_type m_aaPos;
		MathVector<worldDim> m_direction;
};

template <typename TDomain>
void add_extension_ref_mark_adjuster(IRefiner* ref, SmartPtr<ExtensionRefMarkAdjuster<TDomain> > erma);


} // namespace nernst_planck
} // namespace ug

#endif // UG__PLUGINS__NERNST_PLANCK__EXTENSION_REFMARKADJUSTER_H_
