/*
 * domain1d_solution_adjuster.cpp
 *
 *  Created on: 07.01.2015
 *      Author: mbreit
 */

#include "domain1d_solution_adjuster.h"

#include <cstddef>                                                 // for size_t
#include <algorithm>                                               // for sort

#include "common/error.h"                                          // for UG_COND_THROW
#include "common/math/math_vector_matrix/math_vector_functions.h"  // for VecNormalize
#include "lib_algebra/cpu_algebra_types.h"                         // for CPUAlgebra
#include "lib_disc/common/multi_index.h"                           // for DoFIndex
#include "lib_disc/dof_manager/dof_distribution.h"                 // for DoFDistribution
#include "lib_disc/domain.h"                                       // for Domain1d, Doma...
#include "lib_disc/function_spaces/dof_position_util.h"            // for InnerDoFPosition
#include "lib_grid/grid/grid_base_objects.h"                       // for Edge, Vertex...
#include "lib_grid/tools/subset_group.h"                           // for SubsetGroup


namespace ug{
namespace nernst_planck{


template <typename TDomain, typename TAlgebra>
void Domain1dSolutionAdjuster<TDomain, TAlgebra>::
set_sorting_direction(const std::vector<number>& vDir)
{
	UG_COND_THROW(vDir.size() < (size_t) TDomain::dim, "Given sorting direction does not have enough components.");

	for (size_t i = 0; i < (size_t) worldDim; ++i)
		m_sortDir[i] = vDir[i];
}


template <typename TDomain, typename TAlgebra>
void Domain1dSolutionAdjuster<TDomain, TAlgebra>::
adjust_solution(SmartPtr<GridFunction<TDomain, TAlgebra> > u)
{
	// translate subset names to indices
	SubsetGroup ssg(u->approx_space()->subset_handler(), m_vConstrdNames);
	m_vConstrdSI = ssg.index_vector();
	ssg.clear();
	ssg.add(m_vConstrgNames);
	m_vConstrgSI = ssg.index_vector();

	// make sure sorting direction is normalized
	VecNormalize(m_sortDir, m_sortDir);

	// resize data points vector
	ConstSmartPtr<DoFDistribution> dd = u->dof_distribution();
	m_vDataPoints.clear();
	m_vDataPoints.resize(dd->num_fct());

	// collect constraining points
	if (dd->max_dofs(VERTEX)) collect_constrainers<Vertex>(u);
	if (dd->max_dofs(EDGE)) collect_constrainers<Edge>(u);

	// sort constrainer points
	typename DataPoint::CompareFunction mcf;
	for (size_t fct = 0; fct < m_vDataPoints.size(); ++fct)
		std::sort(m_vDataPoints[fct].begin(), m_vDataPoints[fct].end(), mcf);

	// query NN  for all constrained DoFs (log search)
	if (dd->max_dofs(VERTEX)) adjust_constrained<Vertex>(u);
	if (dd->max_dofs(EDGE)) adjust_constrained<Edge>(u);
	if (dd->max_dofs(FACE)) adjust_constrained<Face>(u);
	if (dd->max_dofs(VOLUME)) adjust_constrained<Volume>(u);
}


template <typename TDomain, typename TAlgebra>
template <typename TBaseElem>
void Domain1dSolutionAdjuster<TDomain, TAlgebra>::
collect_constrainers(SmartPtr<GridFunction<TDomain, TAlgebra> > u)
{
	ConstSmartPtr<DoFDistribution> dd = u->dof_distribution();

	// loop constraining subsets
	for (size_t i = 0; i < m_vConstrgSI.size(); ++i)
	{
		int si = m_vConstrgSI[i];

		// loop constraining elements
		typename DoFDistribution::traits<TBaseElem>::const_iterator iter, iterEnd;
		iter = dd->begin<TBaseElem>(si);
		iterEnd = dd->end<TBaseElem>(si);

		for (; iter != iterEnd; ++iter)
		{
			// get constrained elem
			TBaseElem* constrg = *iter;

			// loop functions
			size_t numFct = dd->num_fct();
			for (size_t fct = 0; fct < numFct; fct++)
			{
				if (!dd->is_def_in_subset(fct, si)) continue;

				// get inner DoF indices
				std::vector<DoFIndex> constrgInd;
				dd->inner_dof_indices(constrg, fct, constrgInd, false);

				// get corresponding positions
				std::vector<MathVector<worldDim> > globPos;
				InnerDoFPosition<TDomain>(globPos, constrg, *u->domain(), dd->lfeid(fct));

				// loop all DoFs
				size_t nDof = constrgInd.size();
				UG_COND_THROW(globPos.size() != nDof, "#DoF mismatch");

				for (size_t d = 0; d < nDof; ++d)
				{
					// use scalar product with sorting direction as 1d coordinate
					number coord = VecProd(globPos[d], m_sortDir);

					// get solution value at dof
					number val = DoFRef(*u, constrgInd[d]);

					// push back data point
					m_vDataPoints[fct].push_back(DataPoint(coord, val));
				}
			}
		}
	}
}


template <typename TDomain, typename TAlgebra>
template <typename TBaseElem>
void Domain1dSolutionAdjuster<TDomain, TAlgebra>::
adjust_constrained(SmartPtr<GridFunction<TDomain, TAlgebra> > u)
{
	ConstSmartPtr<DoFDistribution> dd = u->dof_distribution();

	// loop constrained subsets
	for (size_t i = 0; i < m_vConstrdSI.size(); ++i)
	{
		int si = m_vConstrdSI[i];

		// loop constrained elements
		typename DoFDistribution::traits<TBaseElem>::const_iterator iter, iterEnd;
		iter = dd->begin<TBaseElem>(si);
		iterEnd = dd->end<TBaseElem>(si);

		for (; iter != iterEnd; ++iter)
		{
			// get constrained elem
			TBaseElem* constrd = *iter;

			// loop functions
			size_t numFct = dd->num_fct();
			for (size_t fct = 0; fct < numFct; fct++)
			{
				if (!dd->is_def_in_subset(fct, si)) continue;

				// get inner DoF indices
				std::vector<DoFIndex> constrdInd;
				dd->inner_dof_indices(constrd, fct, constrdInd, false);

				// get corresponding positions
				std::vector<MathVector<worldDim> > globPos;
				InnerDoFPosition<TDomain>(globPos, constrd, *u->domain(), dd->lfeid(fct));

				// loop all DoFs
				size_t nDof = constrdInd.size();
				UG_COND_THROW(globPos.size() != nDof, "#DoF mismatch");

				for (size_t d = 0; d < nDof; ++d)
				{
					// use scalar product with sorting direction as 1d coordinate
					number coord = VecProd(globPos[d], m_sortDir);

					// get solution value at dof
					number& val = DoFRef(*u, constrdInd[d]);

					// find nearest neighbor
					std::vector<DataPoint>& vDP = m_vDataPoints[fct];
					size_t nDP = vDP.size();

					UG_COND_THROW(!nDP, "No constrainers available.");

					// treat special case with only one constrainer
					if (nDP == 1)
					{
						val = vDP[0].m_val;
						continue;
					}

					// first step: we start with a correctly sized left side of the array
					size_t step = 1;
					while (step < nDP) step = step << 1;
					size_t curr = step >> 1;
					step = step >> 2;

					// if we get into the right side in our search:
					// add as many elements from the left side as necessary to achieve correct size
					if (vDP[curr].m_coord < coord)
						curr = nDP - curr + step;
					else if (vDP[curr-1].m_coord > coord)
						curr -= step;
					else step = 1; // artificially stop as we have found the correct position

					while ((step = step >> 1))
					{
						if (vDP[curr].m_coord < coord)	// enter right part
							curr += step;
						else if (vDP[curr-1].m_coord > coord)	// enter left part
							curr -= step;
						else step = 1;
					}

					// we are now in between the correct elements
					// compare left and right and decide
					if (coord - vDP.at(curr-1).m_coord < vDP.at(curr).m_coord - coord)
						val = vDP.at(curr-1).m_val;
					else
						val = vDP.at(curr).m_val;
				}
			}
		}
	}
}



// explicit template specializations
#ifdef UG_CPU_1
	#ifdef UG_DIM_1
		template class Domain1dSolutionAdjuster<Domain1d, CPUAlgebra>;
	#endif
	#ifdef UG_DIM_2
		template class Domain1dSolutionAdjuster<Domain2d, CPUAlgebra>;
	#endif
	#ifdef UG_DIM_3
		template class Domain1dSolutionAdjuster<Domain3d, CPUAlgebra>;
	#endif
#endif
#ifdef UG_CPU_5
	#ifdef UG_DIM_1
		template class Domain1dSolutionAdjuster<Domain1d, CPUBlockAlgebra<5> >;
	#endif
	#ifdef UG_DIM_2
		template class Domain1dSolutionAdjuster<Domain2d, CPUBlockAlgebra<5> >;
	#endif
	#ifdef UG_DIM_3
		template class Domain1dSolutionAdjuster<Domain3d, CPUBlockAlgebra<5> >;
	#endif
#endif
#ifdef UG_CPU_6
	#ifdef UG_DIM_1
		template class Domain1dSolutionAdjuster<Domain1d, CPUBlockAlgebra<6> >;
	#endif
	#ifdef UG_DIM_2
		template class Domain1dSolutionAdjuster<Domain2d, CPUBlockAlgebra<6> >;
	#endif
	#ifdef UG_DIM_3
		template class Domain1dSolutionAdjuster<Domain3d, CPUBlockAlgebra<6> >;
	#endif
#endif



} // namespace nernst_planck
} // namespace ug
