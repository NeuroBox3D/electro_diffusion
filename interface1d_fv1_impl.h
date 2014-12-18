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
IInterface1DFV1<TDomain, TAlgebra>::Interface1DMapper::Interface1DMapper(SmartPtr<IAssemble<TAlgebra> > ass)
{
	SmartPtr<AssemblingTuner<TAlgebra> > assTuner = ass->ass_tuner();
	assTuner->set_mapping(this);
}


template <typename TDomain, typename TAlgebra>
void IInterface1DFV1<TDomain, TAlgebra>::Interface1DMapper::add_local_vec_to_global
(
	vector_type& vec,
	const LocalVector& lvec,
	ConstSmartPtr<DoFDistribution> dd
)
{
	const LocalIndices& lind = lvec.get_indices();

	for (size_t fct = 0; fct < lind.num_fct(); fct++)
	{
		for (size_t dof = 0; dof < lind.num_dof(fct); dof++)
		{
			const size_t index = lind.index(fct,dof);
			const size_t comp = lind.comp(fct,dof);

			// find out whether index is constrained in any of the Interfaces associated
			bool constr = false;
			for (size_t intf = 0; intf < m_vspInterface.size(); intf++)
			{
				// is constrained in this interface
				if (m_vspInterface[intf]->m_constraintMap.find(index)
					!= m_vspInterface[intf]->m_constraintMap.end())
				{
					// get index for corresponding 1D interface node
					if (m_vspInterface[intf]->m_fctIndexMapper.find(fct)
						==  m_vspInterface[intf]->m_fctIndexMapper.end())
						UG_THROW("Function for constrained index is unknown in interface.")

					size_t fctIndInIntf = m_vspInterface[intf]->m_fctIndexMapper[fct];
					size_t intf_dofInd = m_vspInterface[intf]->m_algInd[1][fctIndInIntf];

					// add to 1D interface index instead of constrained index
					BlockRef(vec[intf_dofInd], comp) += lvec.value(fct,dof);

					constr = true;
					break;
				}
			}

			// normal mapping and adding in case index is not constrained in any interface
			if (!constr)
				BlockRef(vec[index], comp) += lvec.value(fct,dof);
		}
	}
}


template <typename TDomain, typename TAlgebra>
void IInterface1DFV1<TDomain, TAlgebra>::Interface1DMapper::add_local_mat_to_global
(
	matrix_type& mat,
	const LocalMatrix& lmat,
	ConstSmartPtr<DoFDistribution> dd
)
{
	const LocalIndices& r_lind = lmat.get_row_indices();
	const LocalIndices& c_lind = lmat.get_col_indices();

	for (size_t rfct = 0; rfct < r_lind.num_fct(); rfct++)
	{
		for (size_t rdof = 0; rdof < r_lind.num_dof(rfct); rdof++)
		{
			const size_t r_index = r_lind.index(rfct,rdof);
			const size_t r_comp = r_lind.comp(rfct,rdof);

			// find out whether index is constrained in any of the Interfaces associated
			bool constr = false;
			for (size_t intf = 0; intf < m_vspInterface.size(); intf++)
			{
				// is constrained in this interface
				if (m_vspInterface[intf]->m_constraintMap.find(r_index)
					!= m_vspInterface[intf]->m_constraintMap.end())
				{
					// get index for corresponding 1D interface node
					if (m_vspInterface[intf]->m_fctIndexMapper.find(rfct)
						==  m_vspInterface[intf]->m_fctIndexMapper.end())
						UG_THROW("Function for constrained index is unknown in interface.")

					size_t fctIndInIntf = m_vspInterface[intf]->m_fctIndexMapper[rfct];
					size_t intf_dofInd = m_vspInterface[intf]->m_algInd[1][fctIndInIntf];

					for (size_t cfct = 0; cfct < c_lind.num_fct(); cfct++)
					{
						for (size_t cdof = 0; cdof < c_lind.num_dof(cfct); cdof++)
						{
							const size_t c_index = c_lind.index(cfct,cdof);
							const size_t c_comp = c_lind.comp(cfct,cdof);

							// add to 1D interface index instead of constrained index
							BlockRef(mat(intf_dofInd, c_index), r_comp, c_comp) += lmat.value(rfct, rdof, cfct, cdof);
						}
					}

					constr = true;
					break;
				}
			}

			// normal mapping and adding in case index is not constrained in any interface
			if (!constr)
			{
				for (size_t cfct = 0; cfct < c_lind.num_fct(); cfct++)
				{
					for (size_t cdof = 0; cdof < c_lind.num_dof(cfct); cdof++)
					{
						const size_t c_index = c_lind.index(cfct,cdof);
						const size_t c_comp = c_lind.comp(cfct,cdof);

						// add to 1D interface index instead of constrained index
						BlockRef(mat(r_index, c_index), r_comp, c_comp) += lmat.value(rfct, rdof, cfct, cdof);
					}
				}
			}
		}
	}
}



template <typename TDomain, typename TAlgebra>
void IInterface1DFV1<TDomain, TAlgebra>::Interface1DMapper::add_interface(SmartPtr<interface_type> intf)
{
	m_vspInterface.push_back(intf);
}



template <typename TDomain, typename TAlgebra>
IInterface1DFV1<TDomain, TAlgebra>::IInterface1DFV1
(
	const char* fcts,
	const char* constrained,
	const char* high_dim_intfNode,
	const char* one_dim_intfNode
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

	// store subset names
	m_ssiConstr = std::string(constrained);
	m_ssiIntf[0] = std::string(high_dim_intfNode);
	m_ssiIntf[1] = std::string(one_dim_intfNode);
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

		m_fctIndexMapper[uniqueID] = m_vFct.size();
		m_vFct.push_back(uniqueID);
	}

// get subset indices
	m_siConstr = dom->subset_handler()->get_subset_index(m_ssiConstr.c_str());
	m_siIntf[0] = dom->subset_handler()->get_subset_index(m_ssiIntf[0].c_str());
	m_siIntf[1] = dom->subset_handler()->get_subset_index(m_ssiIntf[1].c_str());

	if (m_siConstr < 0 || m_siIntf[0] < 0 || m_siIntf[1] < 0)	// previous call gives -1 if failed
	{
		UG_THROW("At least one of the given subsets is not contained in the "
				  "geometry handled by this domain's subset handler");
	}

// find algebra indices for interface nodes
	typedef typename MultiGrid::traits<Edge>::secure_container edge_list;

	Vertex* iv1 = NULL;	// high-dim interface node
	Vertex* iv2 = NULL;	// one-dim interface node
	DoFDistribution::traits<Vertex>::const_iterator iter;

	iter = dd->begin<Vertex>(m_siIntf[0]);
	if (iter == dd->end<Vertex>(m_siIntf[0]))
	{
#ifndef UG_PARALLEL
		UG_THROW("No vertex in subset for high-dimensional interface node. This is not allowed!");
#else
		if (pcl::NumProcs() <= 1)
		{
			UG_THROW("No vertex in subset for high-dimensional interface node. This is not allowed!");
		}
		// else do nothing
#endif
	}
	else
	{
		iv1 = *iter;

		// to be on the safe side
		if (++iter != dd->end<Vertex>(m_siIntf[0]))
			{UG_THROW("More than one vertex in subset for high-dimensional interface node. This is not allowed!");}
	}

	iter = dd->begin<Vertex>(m_siIntf[1]);
	if (iter == dd->end<Vertex>(m_siIntf[1]))
	{
#ifndef UG_PARALLEL
		UG_THROW("No vertex in subset for one-dimensional interface node. This is not allowed!");
#else
		if (pcl::NumProcs() <= 1)
		{
			UG_THROW("No vertex in subset for one-dimensional interface node. This is not allowed!");
		}
		// else do nothing
#endif
	}
	else
	{
		iv2 = *iter;

		// to be on the safe side
		if (++iter != dd->end<Vertex>(m_siIntf[1]))
			{UG_THROW("More than one vertex in subset for one-dimensional interface node. This is not allowed!");}
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
				UG_THROW("Function " << m_vFct[fct] << "is not defined on interface node on the high-dimensional side!");
			}
			dd->inner_algebra_indices_for_fct(iv1, ind1, false, m_vFct[fct]);

			UG_ASSERT(ind1.size() == 1, "More (or less) than one function index found on a vertex!");

			m_algInd[0].push_back(ind1[0]);
		}

		int ssi2 = dom->subset_handler()->get_subset_index(iv2);
		for (size_t fct = 0; fct < m_vFct.size(); fct++)
		{
			std::vector<size_t> ind2;

			if (!dd->is_def_in_subset(m_vFct[fct], ssi2))
			{
				UG_THROW("Function " << m_vFct[fct] << "is not defined on interface node on the one-dimensional side!");
			}
			dd->inner_algebra_indices_for_fct(iv2, ind2, false, m_vFct[fct]);

			UG_ASSERT(ind2.size() == 1, "More (or less) than one function index found on a vertex!");

			m_algInd[1].push_back(ind2[0]);
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

		// ensure that the constrainer has been found
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
			m_constraintMap[constrdInd[0]] = ConstraintInfo(constrgInd[0], fct);
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
		std::vector<typename std::map<size_t, ConstraintInfo>::iterator> colIndices;
		for (row_iterator rit = J.begin_row(i); rit != J.end_row(i); ++rit)
		{
			// if a column index is part of the constrained set, then adjust row
			typename std::map<size_t, ConstraintInfo>::iterator constraintMapIt
				= m_constraintMap.find(rit.index());

			if (constraintMapIt != m_constraintMap.end())
			{
				// store constrained col index until row_iterator reaches end
				colIndices.push_back(constraintMapIt);
			}
		}

		// go through row entries again and adapt for constraints
		for (size_t col = 0; col < colIndices.size(); col++)
		{
			size_t constrdInd = colIndices[col]->first;
			size_t constrgInd = colIndices[col]->second.constrgInd;
			size_t fct = colIndices[col]->second.fct;
			size_t IntfNodeHdInd = m_algInd[0][fct];
			size_t IntfNode1dInd = m_algInd[1][fct];

			// get defect derivatives for this vertex
			typename vector_type::value_type defDeriv[3];
			constraintValueDerivs(defDeriv, u[constrgInd], u[IntfNodeHdInd], u[IntfNode1dInd]);

			// deriv wrt constraining
			J(i, constrgInd) += J(i, constrdInd) * defDeriv[0];

			// deriv wrt interface vertex on constraining side
			J(i, IntfNodeHdInd) += J(i, constrdInd) * defDeriv[1];

			// deriv wrt interface vertex on constrained side
			J(i, IntfNode1dInd) += J(i, constrdInd) * defDeriv[2];

			// deriv wrt constrained vertex
			J(i, constrdInd) = 0.0;
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
		size_t intfSideInd = m_algInd[0][fct];
		size_t intfCSideInd = m_algInd[1][fct];

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
	const char* constrained,
	const char* high_dim_intfNode,
	const char* one_dim_intfNode
) : IInterface1DFV1<TDomain, TAlgebra>(fcts, constrained, high_dim_intfNode, one_dim_intfNode) {}


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
	const char* constrained,
	const char* high_dim_intfNode,
	const char* one_dim_intfNode
) : IInterface1DFV1<TDomain, TAlgebra>(fcts, constrained, high_dim_intfNode, one_dim_intfNode) {}


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
