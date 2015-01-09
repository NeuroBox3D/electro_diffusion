/*
 * copy_neighbor_value_constraint_impl.h
 *
 *  Created on: 07.01.2015
 *      Author: mbreit
 */

#include "copy_neighbor_value_constraint.h"

namespace ug{
namespace nernst_planck{


template <typename TDomain, typename TAlgebra>
CopyNeighborValueConstraint<TDomain, TAlgebra>::CopyNeighborValueConstraint
(
	const char* fcts,
	const char* constrained
)
	: m_siConstr(0)
{
	// store function names
	// transform into vector
	if (fcts == NULL) fcts = "";
	m_vsFct = TokenizeString(fcts);

	// remove white space
	for (size_t i = 0; i < m_vsFct.size(); ++i)
		RemoveWhitespaceFromString(m_vsFct[i]);

	// if no function passed, clear functions
	if (m_vsFct.size() == 1 && m_vsFct[0].empty()) m_vsFct.clear();

	// if functions passed with separator, but not all tokens filled, throw error
	for (size_t i = 0; i < m_vsFct.size(); ++i)
	{
		if (m_vsFct.empty())
		{
			UG_THROW("Error while setting functions: Passed function string lacks a function "
					 "specification at position " << i << " (of " << m_vsFct.size()-1 << ")");
		}
	}

	// store subset name
	m_ssiConstr = std::string(constrained);
}



/// destructor
template <typename TDomain, typename TAlgebra>
CopyNeighborValueConstraint<TDomain, TAlgebra>::~CopyNeighborValueConstraint()
{
	// do nothing
}



///	returns the type of constraints
template <typename TDomain, typename TAlgebra>
int CopyNeighborValueConstraint<TDomain, TAlgebra>::type() const
{
	return CT_CONSTRAINTS;
}



///	sets the approximation space
template <typename TDomain, typename TAlgebra>
void CopyNeighborValueConstraint<TDomain, TAlgebra>::
set_approximation_space(SmartPtr<ApproximationSpace<TDomain> > approxSpace)
{
	// check whether the approximation space has already been set
	bool newApproxSpace = (this->m_spApproxSpace != approxSpace);

	// remember approx space
	this->m_spApproxSpace = approxSpace;

	// invoke callback
	if (newApproxSpace) approximation_space_changed();
}


/// called when the approximation space has changed
template <typename TDomain, typename TAlgebra>
void CopyNeighborValueConstraint<TDomain, TAlgebra>::approximation_space_changed()
{
// get fct indices
	// get dof distribution and domain for later use
	ConstSmartPtr<DoFDistribution> dd = this->approximation_space()->dof_distribution(GridLevel());
	SmartPtr<TDomain> dom = this->approximation_space()->domain();

	// store indices of functions in function group
	ConstSmartPtr<FunctionPattern> spFctPattern = this->approximation_space()->dof_distribution_info();
	FunctionGroup fctGrp(spFctPattern);

	try
	{
		fctGrp.add(m_vsFct);
	}
	UG_CATCH_THROW("Cannot find some symbolic function name.");

	// get those indices
	m_vFct.clear();
	for (size_t fct = 0; fct < fctGrp.size(); ++fct)
	{
		// get unique id of function
		size_t uniqueID = fctGrp[fct];

		m_fctIndexMapper[uniqueID] = m_vFct.size();
		m_vFct.push_back(uniqueID);
	}

// get subset indices
	m_siConstr = dom->subset_handler()->get_subset_index(m_ssiConstr.c_str());

	if (m_siConstr < 0)	// previous call gives -1 if failed
	{
		UG_THROW("The given constrained subset is not contained in the "
				  "geometry handled by this domain's subset handler.");
	}

// find constrainers for every constrained vertex
	fill_constraint_map();
}



/// for every constrained vertex: finds the corresponding constrainer and
/// fills the constrainer map with the pair of corresponding indices
template <typename TDomain, typename TAlgebra>
void CopyNeighborValueConstraint<TDomain, TAlgebra>::fill_constraint_map()
{
	typedef typename MultiGrid::traits<Vertex>::secure_container vrt_list;
	typedef typename MultiGrid::traits<Edge>::secure_container edge_list;
	typedef typename MultiGrid::traits<Face>::secure_container face_list;

	ConstSmartPtr<DoFDistribution> dd = this->m_spApproxSpace->dof_distribution(GridLevel());
	SmartPtr<TDomain> dom = this->approximation_space()->domain();

	// number of functions
	size_t numFct = m_vFct.size();

	// loop constrained vertices
	DoFDistribution::traits<Vertex>::const_iterator iter, iterEnd;
	iter = dd->begin<Vertex>(m_siConstr);
	iterEnd = dd->end<Vertex>(m_siConstr);

	for (; iter != iterEnd; ++iter)
	{
		// get constrained vertex
		Vertex* constrd = *iter;

		// loop associated edges; find the one whose other end is not in the constrained set
		Vertex* constrg = NULL;
		edge_list el;
		dom->grid()->associated_elements(el, constrd);
		for (size_t edge = 0; edge < el.size(); edge++)
		{
			Vertex* other;
			el[edge]->get_opposing_side(constrd, &other);

			if (dom->subset_handler()->get_subset_index(other) != m_siConstr)
			{
				constrg = other;
				break;
			}
		}

		// if not found, we are on the far side of a cuboid
		// in that case, loop associated faces and their vertices:
		// the only vertex not contained in constrained set is the constrainer
		if (!constrg)
		{
			face_list fl;
			dom->grid()->associated_elements(fl, constrd);
			for (size_t face = 0; face < fl.size(); face++)
			{
				vrt_list vl;
				dom->grid()->associated_elements(vl, fl[face]);
				for (size_t vrt = 0; vrt < vl.size(); vrt++)
				{
					if (dom->subset_handler()->get_subset_index(vl[vrt]) != m_siConstr)
					{
						constrg = vl[vrt];
						goto check_constrg;	// yeah! a good use of goto!
					}
				}
			}
		}

		// check that constrainer has now been found
	check_constrg:
		if (!constrg) {UG_THROW("No constrainer found!");}

	// find indices
		// loop functions
		for (size_t fct = 0; fct < numFct; fct++)
		{
			std::vector<size_t> constrdInd;
			std::vector<size_t> constrgInd;

			// constrained index
			if (!dd->is_def_in_subset(m_vFct[fct], m_siConstr))
				{UG_THROW("Function " << m_vFct[fct] << "is not defined on constrained subset " << m_siConstr << ".");}
			dd->inner_algebra_indices_for_fct(constrd, constrdInd, false, m_vFct[fct]);

			// constraining index
			int si = dd->subset_handler()->get_subset_index(constrg);
			if (!dd->is_def_in_subset(m_vFct[fct], si))
				{UG_THROW("Function " << m_vFct[fct] << "is not defined on constraining subset " << si << ".");}
			dd->inner_algebra_indices_for_fct(constrg, constrgInd, false, m_vFct[fct]);

			// should not happen, but for debugging purposes:
			UG_ASSERT(constrgInd.size() == 1, "More (or less) than one function index found on a vertex!");
			UG_ASSERT(constrdInd.size() == 1, "More (or less) than one function index found on a vertex!");

			// fill map with pair of indices
			m_constraintMap[constrdInd[0]] = constrgInd[0];
		}
	}
}


///	adapts jacobian to enforce constraints
template <typename TDomain, typename TAlgebra>
void CopyNeighborValueConstraint<TDomain, TAlgebra>::adjust_jacobian
(	matrix_type& J,
	const vector_type& u,
	ConstSmartPtr<DoFDistribution> dd,
	number time,
	ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
	const number s_a0
)
{
	typedef typename matrix_type::row_iterator row_iterator;

	// loop rows
	size_t nr = J.num_rows();
	for (size_t i = 0; i < nr; i++)
	{
		// write identity row for constrained (first: delete all couplings)
		if (m_constraintMap.find(i) != m_constraintMap.end())
		{
			// we put curly brackets around this so that the row_iterator will be removed afterwards
			// otherwise, we get a negative assert from the matrix implementation
			{
				const row_iterator iterEnd = J.end_row(i);
				for (row_iterator conn = J.begin_row(i); conn != iterEnd; ++conn)
					conn.value() = 0.0;
			}
			J(i, i)	= 1.0;

			// do not alter the columns after that!
			continue;
		}

		// we do not need to adapt the columns for this row, as this constraint
		// is only meant for dofs that are not involved in any elemDisc
	}
}



///	adapts defect to enforce constraints
template <typename TDomain, typename TAlgebra>
void CopyNeighborValueConstraint<TDomain, TAlgebra>::adjust_defect
(	vector_type& d,
	const vector_type& u,
	ConstSmartPtr<DoFDistribution> dd,
	number time,
	ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
	const std::vector<number>* vScaleMass,
	const std::vector<number>* vScaleStiff
)
{
	// loop constrained vertices
	typename std::map<size_t, size_t>::iterator constraintMapEnd = m_constraintMap.end();
	typename std::map<size_t, size_t>::iterator constraintMapIt;
	for (constraintMapIt = m_constraintMap.begin(); constraintMapIt != constraintMapEnd; ++constraintMapIt)
	{
		size_t index = constraintMapIt->first;
		d[index] = 0.0;
	}
}



///	adapts matrix and rhs (linear case) to enforce constraints
template <typename TDomain, typename TAlgebra>
void CopyNeighborValueConstraint<TDomain, TAlgebra>::adjust_linear
(	matrix_type& mat,
	vector_type& rhs,
	ConstSmartPtr<DoFDistribution> dd,
	number time
)
{
	UG_THROW("This feature is not implemented.");
}



///	adapts a rhs to enforce constraints
template <typename TDomain, typename TAlgebra>
void CopyNeighborValueConstraint<TDomain, TAlgebra>::adjust_rhs
(	vector_type& rhs,
	const vector_type& u,
	ConstSmartPtr<DoFDistribution> dd,
	number time
)
{
	UG_THROW("This feature is deactivated in order to test, whether it is needed in the first place.");
}



///	sets the constraints in a solution vector
template <typename TDomain, typename TAlgebra>
void CopyNeighborValueConstraint<TDomain, TAlgebra>::adjust_solution
(	vector_type& u,
	ConstSmartPtr<DoFDistribution> dd,
	number time
)
{
#ifdef UG_PARALLEL
	// check parallel storage type
	if (pcl::NumProcs() > 1 &&
		(!u.has_storage_type(PST_CONSISTENT)))
	{
		UG_THROW("Expected solution to have parallel storage type PST_CONSISTENT,\n"
				 "but it does not. PST is '" << u.get_storage_type() << "'.");
	}
#endif

	// loop constrained vertices
	typename std::map<size_t, size_t>::iterator constraintMapEnd = m_constraintMap.end();
	typename std::map<size_t, size_t>::iterator constraintMapIt;
	for (constraintMapIt = m_constraintMap.begin(); constraintMapIt != constraintMapEnd; ++constraintMapIt)
	{
		size_t constrdInd = constraintMapIt->first;
		size_t constrgInd = constraintMapIt->second;

		u[constrdInd] = u[constrgInd];
	}
}


} // namespace nernst_planck
} // namespace ug
