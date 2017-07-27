/*
 * intf_refMarkAdjuster.h
 *
 *  Created on: 29.04.2016
 *      Author: mbreit
 */

#ifndef UG__PLUGINS__NERNST_PLANCK__EXTENSION_REFMARKADJUSTER_H_
#define UG__PLUGINS__NERNST_PLANCK__EXTENSION_REFMARKADJUSTER_H_


#include <vector>                                                 // for vector

#include "common/math/math_vector_matrix/math_vector.h"           // for MathVector
#include "common/types.h"                                         // for number
#include "common/util/smart_pointer.h"                            // for SmartPtr, ConstSmartPtr
#include "lib_grid/refinement/ref_mark_adjuster_interface.h"      // for IRefMarkAdjuster
#include "lib_grid/tools/subset_handler_interface.h"              // for ISubsetHandler


namespace ug {

// forward declarations
class Vertex;
class Edge;
class Face;
class Volume;
class IRefiner;

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
