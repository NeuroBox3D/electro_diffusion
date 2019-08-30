/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Markus Breit
 * Creation date: 2014-06-06
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


#ifndef UG__PLUGINS__NERNST_PLANCK__EXTENSION__DOMAIN1D_SOLUTION_ADJUSTER_H
#define UG__PLUGINS__NERNST_PLANCK__EXTENSION__DOMAIN1D_SOLUTION_ADJUSTER_H

#include <string>                                        // for string
#include <vector>                                        // for vector, allocator
#include "common/types.h"                                // for number
#include "common/math/math_vector_matrix/math_vector.h"  // for MathVector
#include "common/util/smart_pointer.h"                   // for SmartPtr
#include "lib_disc/function_spaces/grid_function.h"      // for GridFunction


namespace ug {
namespace nernst_planck {


/// Constraint for DoFs that are required to have the same value as their neighbors.
/**
 * This class is intended to make it possible to constrain DoFs on a 1d extension
 * to the same value as their neighbors from a specified subset.
 *
 * These neighbors are found using a nearest neighbor search, where only vector
 * components in a user-defined direction (the direction of the extension) are
 * considered.
 *
 * This class can be used in the PNP problem to constrain the useless vertices
 * in the 1d extensions in such a way that the values on every cross-section
 * through the extension are constant (which is useful for the visualization).
 *
 * \tparam	TDomain				type of Domain
 * \tparam	TAlgebra			type of Algebra
 *
 * \date 07.01.2015
 * \author mbreit
 */
template <typename TDomain, typename TAlgebra>
class Domain1dSolutionAdjuster
{
	public:
		static const int worldDim = TDomain::dim;

	public:
		Domain1dSolutionAdjuster() : m_sortDir(0.0)
		{m_sortDir[0] = 1.0;}

		void add_constrained_subset(const std::string& ss) {m_vConstrdNames.push_back(ss);}
		void add_constrainer_subset(const std::string& ss) {m_vConstrgNames.push_back(ss);}

		void set_sorting_direction(const std::vector<number>& vDir);

		void adjust_solution(SmartPtr<GridFunction<TDomain, TAlgebra> > u);

	private:
		template <typename TBaseElem> void collect_constrainers(SmartPtr<GridFunction<TDomain, TAlgebra> > u);
		template <typename TBaseElem> void adjust_constrained(SmartPtr<GridFunction<TDomain, TAlgebra> > u);

	protected:
		std::vector<std::string> m_vConstrdNames;
		std::vector<std::string> m_vConstrgNames;

		std::vector<int> m_vConstrdSI;
		std::vector<int> m_vConstrgSI;

		MathVector<worldDim> m_sortDir;

		struct DataPoint
		{
			public:
				DataPoint(number coord, number val)
				: m_coord(coord), m_val(val) {}

				struct CompareFunction
				{
					public:
						bool operator() (const DataPoint& a, const DataPoint& b)
						{return a.m_coord < b.m_coord;}
				};
				friend struct CompareFunction;

				number m_coord;
				number m_val;
		};

		std::vector<std::vector<DataPoint> > m_vDataPoints;
};


} // namespace nernst_planck
} // namespace ug


#endif // UG__PLUGINS__NERNST_PLANCK__EXTENSION__DOMAIN1D_SOLUTION_ADJUSTER_H
