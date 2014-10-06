/*
 * interface1d_fv1_impl.h
 *
 *  Created on: 06.06.2014
 *      Author: mbreit
 */

#include "interface1d_fv1.h"

namespace ug{
namespace nernst_planck{

template <typename TDomain, typename TAlgebra>
IInterface1DFV1<TDomain, TAlgebra>::IInterface1DFV1
(	const char* fcts,
	const char* high_dim_subset,
	const char* one_dim_subset,
	const char* two_dim_intfNode,
	const char* extension_subset
)
	: m_ssiExt(0), m_ssiIN(0)
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

	// store subset names
	m_sssiExt = std::string(extension_subset);
	m_sssiIN = std::string(two_dim_intfNode);
	m_sssi[0] = std::string(high_dim_subset);
	m_sssi[1] = std::string(one_dim_subset);
}



/// destructor
template <typename TDomain, typename TAlgebra>
IInterface1DFV1<TDomain, TAlgebra>::~IInterface1DFV1()
{
	// do nothing
}



///	returns the type of constraints
template <typename TDomain, typename TAlgebra>
int IInterface1DFV1<TDomain, TAlgebra>::type() const
{
	return CT_CONSTRAINTS;
}



///	sets the approximation space
template <typename TDomain, typename TAlgebra>
void IInterface1DFV1<TDomain, TAlgebra>::
set_approximation_space(SmartPtr<ApproximationSpace<TDomain> > approxSpace)
{
	// check whether the approximation space has already been set
	bool newApproxSpace = (this->m_spApproxSpace != approxSpace);

	// remember approx space
	this->m_spApproxSpace = approxSpace;

	// invoke callback
	if (newApproxSpace) approximation_space_changed();
}

/*
template <typename TDomain, typename TAlgebra>
void IInterface1DFV1<TDomain, TAlgebra>::
check_values_at_interface(const vector_type& sol)
{
	size_t numFct = m_vFct.size();
	for (size_t fct = 0; fct < numFct; fct++)
		std::cout << "u_2d["<<fct<<"] = "<< sol[m_algInd[0][fct]] << "  ";
	std::cout << std::endl;
	for (size_t fct = 0; fct < numFct; fct++)
		std::cout << "u_1d["<<fct<<"] = "<< sol[m_algInd[1][fct]] << "  ";
	std::cout << std::endl << std::endl;
}
*/


/// called when the approximation space has changed
template <typename TDomain, typename TAlgebra>
void IInterface1DFV1<TDomain, TAlgebra>::approximation_space_changed()
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

		m_vFct.push_back(uniqueID);
	}

// get subset indices
	m_ssiExt = dom->subset_handler()->get_subset_index(m_sssiExt.c_str());
	m_ssiIN = dom->subset_handler()->get_subset_index(m_sssiIN.c_str());
	m_ssi[0] = dom->subset_handler()->get_subset_index(m_sssi[0].c_str());
	m_ssi[1] = dom->subset_handler()->get_subset_index(m_sssi[1].c_str());

	if (m_ssiExt < 0 || m_ssiIN < 0 || m_ssi[0] < 0 || m_ssi[1] < 0)	// previous call gives -1 if failed
	{
		UG_THROW("At least one of the given subsets is not contained in the "
				  "geometry handled by this domain's subset handler");
	}

// find algebra indices for interface nodes
// The way this is done here is as follows:
// 1a) Find the constrained vertex on the 1d side uniquely identified by one_dim_subset.
// 2b) Get the unique vertex connected by an edge to it: This is the interface node on the 1d side.
// 2) Find the constrained vertex on the high-dimensional side uniquely identified by two_dim_intfNode.
	typedef typename MultiGrid::traits<Edge>::secure_container edge_list;

	// 1a)
	Vertex* iv1 = NULL;
	Vertex* iv2 = NULL;
	DoFDistribution::traits<Vertex>::const_iterator iter;
	iter = dd->begin<Vertex>(m_ssi[1]);

	if (iter == dd->end<Vertex>(m_ssi[1]))
	{
#ifndef UG_PARALLEL
		UG_THROW("No vertex in constrained subset for 1d side. This is not allowed!");
#else
		if (pcl::NumProcs() <= 1)
		{
			UG_THROW("No vertex in constrained subset for 1d side. This is not allowed!");
		}
		// else do nothing
#endif
	}
	else
	{
		Vertex* constrd = *iter;

		// 1b)
		edge_list el;
		dom->grid()->associated_elements(el, constrd);

		// find edge in the 1d extension domain subset and take other end
		for (size_t i = 0; i < el.size(); i++)
		{
			if (dom->subset_handler()->get_subset_index(el[i]) == m_ssiExt)
			{
				el[i]->get_opposing_side(constrd, &iv1);
				break;
			}
		}
		if (!iv1)
		{
			UG_THROW("Constrained vertex on 1d end is not connected to any edge of the extension."
					 "This case is not permitted.");
		}

		// to be on the safe side
		if (++iter != dd->end<Vertex>(m_ssi[1]))
				{UG_THROW("More than one vertex in constrained subset for 1d side. This is not allowed!");}
	}

	// 2)
	iter = dd->begin<Vertex>(m_ssiIN);
	if (iter == dd->end<Vertex>(m_ssiIN))
	{
#ifndef UG_PARALLEL
		UG_THROW("Vertex for full-dimensional interface node could not be found via its subset index.");
#else
		if (pcl::NumProcs() == 1)
		{
			UG_THROW("Vertex for full-dimensional interface node could not be found via its subset index.");
		}
		// else do nothing
#endif
	}
	else
	{
		iv2 = *iter;

		// to be on the safe side
		if (++iter != dd->end<Vertex>(m_ssiIN))
			{UG_THROW("More than one vertex in full-dim interface node subset. This is not allowed!");}
	}

	// both interface nodes must be on the same processor
	if ((iv1 != NULL) != (iv2 != NULL))
	{
		UG_THROW("The interface nodes of this interface are not both located on the same processor.\n"
				"This is not allowed (as it will typically engender convergence problems).");
	}

	// get algebra indices
	if (iv1 != NULL && iv2 != NULL)
	{
		int ssi1 = dom->subset_handler()->get_subset_index(iv1);
		for (size_t fct = 0; fct < m_vFct.size(); fct++)
		{
			std::vector<size_t> ind1;

			if (!dd->is_def_in_subset(m_vFct[fct], ssi1))
			{
				UG_THROW("Function " << m_vFct[fct] << "is not defined on constraining subset "
						 "for interface node on 1d side!");
			}
			dd->inner_algebra_indices_for_fct(iv1, ind1, false, m_vFct[fct]);

			UG_ASSERT(ind1.size() == 1, "More (or less) than one function index found on a vertex!");

			m_algInd[1].push_back(ind1[0]);
		}

		int ssi2 = dom->subset_handler()->get_subset_index(iv2);
		for (size_t fct = 0; fct < m_vFct.size(); fct++)
		{
			std::vector<size_t> ind2;

			if (!dd->is_def_in_subset(m_vFct[fct], ssi2))
			{
				UG_THROW("Function " << m_vFct[fct] << "is not defined on constraining subset "
						 "for interface node on full-dim side!");
			}
			dd->inner_algebra_indices_for_fct(iv2, ind2, false, m_vFct[fct]);

			UG_ASSERT(ind2.size() == 1, "More (or less) than one function index found on a vertex!");

			m_algInd[0].push_back(ind2[0]);
		}
	}

// find constrainers for every constrained vertex
	fill_constraint_map();

	// disallow constrained nodes without having the interface nodes on the same processor
	if (m_constraintMap.size() && !(iv1 && iv2))
	{
		UG_THROW("Processor has constrained nodes but not the necessary interface nodes. This is not allowed.\n"
				 "Make sure the constrained nodes are not separated from their interface nodes by the partitioning.");
	}
}



/// for every constrained vertex: finds the corresponding constrainer and
/// fills the constrainer map with the pair of corresponding indices
template <typename TDomain, typename TAlgebra>
void IInterface1DFV1<TDomain, TAlgebra>::fill_constraint_map()
{
	typedef typename MultiGrid::traits<Edge>::secure_container edge_list;

	ConstSmartPtr<DoFDistribution> dd = this->m_spApproxSpace->dof_distribution(GridLevel());
	SmartPtr<TDomain> dom = this->approximation_space()->domain();

	// number of functions
	size_t numFct = m_vFct.size();

	// high-dimensional end
	{
		DoFDistribution::traits<Vertex>::const_iterator iter, iterEnd;
		iter = dd->begin<Vertex>(m_ssi[0]);
		iterEnd = dd->end<Vertex>(m_ssi[0]);

		for (; iter != iterEnd; ++iter)
		{
			// get constrained vertex
			Vertex* constrd = *iter;

			// loop edges associated; find the one whose other end is not in the constrained set
			Vertex* constrg = NULL;
			edge_list el;
			dom->grid()->associated_elements(el, constrd);
			for (size_t edge = 0; edge < el.size(); edge++)
			{
				Vertex* other;
				el[edge]->get_opposing_side(constrd, &other);

				if (dom->subset_handler()->get_subset_index(other) != m_ssi[0])
				{
					constrg = other;
					break;
				}
			}

			// ensure that the constrainer has been found
			if (!constrg) {UG_THROW("No constrainer found!");}

		// find indices
			// loop functions
			for (size_t fct = 0; fct < numFct; fct++)
			{
				std::vector<size_t> constrdInd;
				std::vector<size_t> constrgInd;

				// constrained index
				if (!dd->is_def_in_subset(m_vFct[fct], m_ssi[0]))
					{UG_THROW("Function " << m_vFct[fct] << "is not defined on constrained subset " << m_ssi[0] << ".");}
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
				m_constraintMap[constrdInd[0]] = ConstraintInfo(constrgInd[0], fct, 0);
			}
		}
	}

	// 1d end
	{
		DoFDistribution::traits<Vertex>::const_iterator iter, iterEnd;
		iter = dd->begin<Vertex>(m_ssi[1]);
		iterEnd = dd->end<Vertex>(m_ssi[1]);

		for (; iter != iterEnd; ++iter)
		{
			// get constrained vertex
			Vertex* constrd = *iter;

			// get associated edge (there must be only one); take vertex on the other side
			edge_list el;
			dom->grid()->associated_elements(el, constrd);

			// find edge in the 1d extension domain subset and take other end
			Vertex* constrg = NULL;
			for (size_t i = 0; i < el.size(); i++)
			{
				if (dom->subset_handler()->get_subset_index(el[i]) == m_ssiExt)
				{
					el[i]->get_opposing_side(constrd, &constrg);
					break;
				}
			}
			if (!constrg)
			{
				UG_THROW("Constrained vertex on 1d end is not connected to any edge of the extension."
						 "This case is not permitted.");
			}

		// find indices
			// loop functions
			for (size_t fct = 0; fct < numFct; fct++)
			{
				std::vector<size_t> constrdInd;
				std::vector<size_t> constrgInd;

				// constrained index
				if (!dd->is_def_in_subset(m_vFct[fct], m_ssi[1]))
					{UG_THROW("Function " << m_vFct[fct] << "is not defined on constrained subset " << m_ssi[1] << ".");}
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
				m_constraintMap[constrdInd[0]] = ConstraintInfo(constrgInd[0], fct, 1);
			}
		}
	}
/* DEBUGGING
	typedef std::map<Vertex*, Vertex*>::iterator map_it;
	map_it iter = m_constrainerMap.begin();
	map_it iterEnd = m_constrainerMap.end();
	for (; iter != iterEnd; ++iter)
	{
		std::vector<size_t> ind1, ind2;
		dd->inner_algebra_indices_for_fct(iter->first, ind1, false, 0);
		dd->inner_algebra_indices_for_fct(iter->second, ind2, false, 0);
		UG_LOG(ind1[0] << " -> " << ind2[0] << "\n");
	}
*/
}


///	adapts jacobian to enforce constraints
template <typename TDomain, typename TAlgebra>
void IInterface1DFV1<TDomain, TAlgebra>::adjust_jacobian
(	matrix_type& J,
	const vector_type& u,
	ConstSmartPtr<DoFDistribution> dd,
	number time,
	ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
	const number s_a0
)
{
/*
#ifdef UG_PARALLEL
	// check parallel storage type
	if (pcl::NumProcs() > 1 &&
		(!u.has_storage_type(PST_CONSISTENT)))
	{
		UG_THROW("Expected solution to have parallel storage type PST_CONSISTENT,\n"
				 "but it does not. PST is '" << u.get_storage_type() << "'.");
	}

	// J is supposed to have PST "additive" before application of constraints
	// (however, this can not be checked for here since the status is not set but after
	// the application of all constraints in domain_disc_impl.h).
	// It will still have PST "additive" afterwards.
#endif
*/

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

		// adapt current row if it depends on any constrained value
		// including constraints from other equations (if coupled)

		// loop column existing entries of row
		for (row_iterator rit = J.begin_row(i); rit != J.end_row(i); ++rit)
		{
			// if a column index is part of the constrained set, then adjust row
			typename std::map<size_t, ConstraintInfo>::iterator constraintMapIt
				= m_constraintMap.find(rit.index());

			if (constraintMapIt != m_constraintMap.end())
			{
				size_t constrdInd = constraintMapIt->first;
				size_t constrgInd = constraintMapIt->second.constrgInd;
				size_t fct = constraintMapIt->second.fct;
				size_t side = constraintMapIt->second.side;
				size_t intfSideInd = m_algInd[side][fct];
				size_t intfCSideInd = m_algInd[(side+1) % 2][fct];

				// get defect derivatives for this vertex
				typename vector_type::value_type defDeriv[3];
				constraintValueDerivs(defDeriv, u[constrgInd], u[intfSideInd], u[intfCSideInd]);

				// deriv wrt constraining
				J(i, constrgInd) += J(i, constrdInd) * defDeriv[0];

				// deriv wrt interface vertex on constraining side
				J(i, intfSideInd) += J(i, constrdInd) * defDeriv[1];

				// deriv wrt interface vertex on constrained side
				J(i, intfCSideInd) += J(i, constrdInd) * defDeriv[2];

				// deriv wrt constrained vertex
				J(i, constrdInd) = 0.0;
			}
		}
	}
}



///	adapts defect to enforce constraints
template <typename TDomain, typename TAlgebra>
void IInterface1DFV1<TDomain, TAlgebra>::adjust_defect
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
	typename std::map<size_t, ConstraintInfo>::iterator constraintMapEnd = m_constraintMap.end();
	typename std::map<size_t, ConstraintInfo>::iterator constraintMapIt;
	for (constraintMapIt = m_constraintMap.begin(); constraintMapIt != constraintMapEnd; ++constraintMapIt)
	{
		size_t index = constraintMapIt->first;
		d[index] = 0.0;
	}
}



///	adapts matrix and rhs (linear case) to enforce constraints
template <typename TDomain, typename TAlgebra>
void IInterface1DFV1<TDomain, TAlgebra>::adjust_linear
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
void IInterface1DFV1<TDomain, TAlgebra>::adjust_rhs
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
void IInterface1DFV1<TDomain, TAlgebra>::adjust_solution
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
	typename std::map<size_t, ConstraintInfo>::iterator constraintMapEnd = m_constraintMap.end();
	typename std::map<size_t, ConstraintInfo>::iterator constraintMapIt;
	for (constraintMapIt = m_constraintMap.begin(); constraintMapIt != constraintMapEnd; ++constraintMapIt)
	{
		size_t constrdInd = constraintMapIt->first;
		size_t constrgInd = constraintMapIt->second.constrgInd;
		size_t fct = constraintMapIt->second.fct;
		size_t side = constraintMapIt->second.side;
		size_t intfSideInd = m_algInd[side][fct];
		size_t intfCSideInd = m_algInd[(side+1) % 2][fct];

		typename vector_type::value_type def;
		constraintValue(def, u[constrgInd], u[intfSideInd], u[intfCSideInd]);

		u[constrdInd] = def;
	}
}



//////////////////////////////////
//   AdditiveInterface1DFV1   	//
//////////////////////////////////

/// constructor with strings
template <typename TDomain, typename TAlgebra>
AdditiveInterface1DFV1<TDomain, TAlgebra>::AdditiveInterface1DFV1
(	const char* fcts,
	const char* high_dim_subset,
	const char* one_dim_subset,
	const char* two_dim_intfNode,
	const char* extension_subset
) : IInterface1DFV1<TDomain, TAlgebra>(fcts, high_dim_subset, one_dim_subset, two_dim_intfNode, extension_subset) {}


/// function returning the defect value for a constrained node
template <typename TDomain, typename TAlgebra>
void AdditiveInterface1DFV1<TDomain, TAlgebra>::constraintValue
(	typename vector_type::value_type& d,
	const typename vector_type::value_type& u_c,
	const typename vector_type::value_type& u_itf0,
	const typename vector_type::value_type& u_itf1
)
{
	d = u_c + (u_itf1 - u_itf0);
}

/// function returning the defect derivatives for a constrained node
template <typename TDomain, typename TAlgebra>
void AdditiveInterface1DFV1<TDomain, TAlgebra>::constraintValueDerivs
(	typename vector_type::value_type dd[3],
	const typename vector_type::value_type& u_c,
	const typename vector_type::value_type& u_itf0,
	const typename vector_type::value_type& u_itf1
)
{
	dd[0] =  1.0;
	dd[1] = -1.0;
	dd[2] =  1.0;
}




//////////////////////////////////////
//   MultiplicativeInterface1DFV1   //
//////////////////////////////////////

/// constructor with strings
template <typename TDomain, typename TAlgebra>
MultiplicativeInterface1DFV1<TDomain, TAlgebra>::MultiplicativeInterface1DFV1
(	const char* fcts,
	const char* high_dim_subset,
	const char* one_dim_subset,
	const char* two_dim_intfNode,
	const char* extension_subset
) : IInterface1DFV1<TDomain, TAlgebra>(fcts, high_dim_subset, one_dim_subset, two_dim_intfNode, extension_subset) {}


/// function returning the defect value for a constrained node
template <typename TDomain, typename TAlgebra>
void MultiplicativeInterface1DFV1<TDomain, TAlgebra>::constraintValue
(	typename vector_type::value_type& d,
	const typename vector_type::value_type& u_c,
	const typename vector_type::value_type& u_itf0,
	const typename vector_type::value_type& u_itf1
)
{
	if (std::fabs(u_itf0) < 2 * std::numeric_limits<typename vector_type::value_type>::denorm_min())
		{UG_THROW("Denominator practically zero (" << u_itf0 << ").");}

	d = u_c * (u_itf1 / u_itf0);
}

/// function returning the defect derivatives for a constrained node
template <typename TDomain, typename TAlgebra>
void MultiplicativeInterface1DFV1<TDomain, TAlgebra>::constraintValueDerivs
(	typename vector_type::value_type dd[3],
	const typename vector_type::value_type& u_c,
	const typename vector_type::value_type& u_itf0,
	const typename vector_type::value_type& u_itf1
)
{
	if (std::fabs(u_itf0) < 2 * std::numeric_limits<typename vector_type::value_type>::denorm_min())
		{UG_THROW("Denominator practically zero.");}

	dd[0] =  u_itf1 / u_itf0;
	dd[1] = -u_c * u_itf1 / (u_itf0*u_itf0);
	dd[2] =  u_c / u_itf0;
}

} // namespace nernst_planck
} // namespace ug
