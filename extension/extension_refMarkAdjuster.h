/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Markus Breit
 * Creation date: 2016-04-29
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

#ifndef UG__PLUGINS__NERNST_PLANCK__EXTENSION__EXTENSION_REFMARKADJUSTER_H
#define UG__PLUGINS__NERNST_PLANCK__EXTENSION__EXTENSION_REFMARKADJUSTER_H


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

#endif // UG__PLUGINS__NERNST_PLANCK__EXTENSION__EXTENSION_REFMARKADJUSTER_H
