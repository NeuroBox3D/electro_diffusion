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


// fill helper structures
	fill_constrainer_map();
	fill_defect_influence_map();


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

	// the dist for the correct node must be _exactly_ zero
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

	// get algebra indices
	if (iv1)
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
	}

	if (iv2)
	{
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
}



/// for every constrained vertex: finds the corresponding constrainer
template <typename TDomain, typename TAlgebra>
void IInterface1DFV1<TDomain, TAlgebra>::fill_constrainer_map()
{
	typedef typename MultiGrid::traits<Edge>::secure_container edge_list;

	ConstSmartPtr<DoFDistribution> dd = this->m_spApproxSpace->dof_distribution(GridLevel());
	SmartPtr<TDomain> dom = this->approximation_space()->domain();

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

			if (!constrg) {UG_THROW("No constrainer found!");}
			else m_constrainerMap[constrd] = constrg;
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
			Vertex* other = NULL;
			for (size_t i = 0; i < el.size(); i++)
			{
				if (dom->subset_handler()->get_subset_index(el[i]) == m_ssiExt)
				{
					el[i]->get_opposing_side(constrd, &other);
					break;
				}
			}
			if (!other)
			{
				UG_THROW("Constrained vertex on 1d end is not connected to any edge of the extension."
						 "This case is not permitted.");
			}

			if (dom->subset_handler()->get_subset_index(other) != m_ssi[1])
				 m_constrainerMap[constrd] = other;
			else
				{UG_THROW("Connected constrainer is in the constrained subset!");}
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



/// for every constrained vertex: collects the vertices in the constraining set
/// whose defect will be influenced by the constrained one
template <typename TDomain, typename TAlgebra>
void IInterface1DFV1<TDomain, TAlgebra>::fill_defect_influence_map()
{
	typedef typename MultiGrid::traits<elem_type>::secure_container elem_list;
	typedef typename MultiGrid::traits<Vertex>::secure_container vertex_list;

	ConstSmartPtr<DoFDistribution> dd = this->approximation_space()->dof_distribution(GridLevel());
	SmartPtr<TDomain> dom = this->approximation_space()->domain();

	// high-dimensional end
	{
		DoFDistribution::traits<Vertex>::const_iterator iter, iterEnd;
		iter = dd->begin<Vertex>(m_ssi[0]);
		iterEnd = dd->end<Vertex>(m_ssi[0]);

		for (; iter != iterEnd; ++iter)
		{
			// get constrained vertex
			Vertex* constrd = *iter;

			// find the vertices in the constraining subset whose defect depends on this constrained dof
			elem_list el;
			dom->grid()->associated_elements(el, constrd);

			// remember vertices already pushed back to avoid double entries
			std::map<Vertex*, bool> considered;

			// loop all associated elements
			for (size_t elem = 0; elem < el.size(); elem++)
			{
				// get vertices associated to this elem
				vertex_list otherVertices;
				dom->grid()->associated_elements(otherVertices, el[elem]);

				// loop vertices
				for (size_t vrt = 0; vrt < otherVertices.size(); vrt++)
				{
					Vertex* other = otherVertices[vrt];

					// only add if subset index matches constrained set and not yet added
					if (dom->subset_handler()->get_subset_index(other) == m_ssi[0]
						&& !considered[other])
					{
						considered[other] = true;
						UG_ASSERT(m_constrainerMap.find(other) != m_constrainerMap.end(),
								 "Discovered constrained index that is no key in constrainer map!");
						m_defectInfluenceMap[constrd].push_back(m_constrainerMap[other]);
					}
				}
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

			m_defectInfluenceMap[constrd].push_back(m_constrainerMap[constrd]);
		}
	}

/* DEBUGGING
	typedef std::map<Vertex*, std::vector<Vertex*> >::iterator map_it;
	map_it iter = m_defectInfluenceMap.begin();
	map_it iterEnd = m_defectInfluenceMap.end();
	for (; iter != iterEnd; ++iter)
	{
		std::vector<size_t> ind1, ind2;
		dd->inner_algebra_indices_for_fct(iter->first, ind1, false, 0);
		UG_LOG(ind1[0] << " -> ");
		for (size_t j=0; j< iter->second.size(); j++)
		{
			dd->inner_algebra_indices_for_fct(iter->second[j], ind2, true, 0);
			UG_LOG(ind2[0] << " ");
		}
		UG_LOG("\n");
	}
*/
}


template <typename TDomain, typename TAlgebra>
void IInterface1DFV1<TDomain, TAlgebra>::fill_sol_at_intf
(
	std::vector<number> intf_val[2],
	const vector_type& u,
	ConstSmartPtr<DoFDistribution> dd
)
{
	size_t numFct = m_vFct.size();

	// resize
	intf_val[0].resize(numFct);
	intf_val[1].resize(numFct);

#ifdef UG_PARALLEL
if (pcl::NumProcs() > 1)
{
// communicate here: get values for u at interface nodes (as they might not be present on this proc)
	pcl::ProcessCommunicator com;

	// compose one vector of all values u[m_algInd[side][fct]]
	// and add two control values that will indicate how many procs possessed the values
	std::vector<number> comm_vec(2*numFct+2, 0.0);
	if (m_algInd[0].size() == numFct)
	{
		for (size_t fct = 0; fct < numFct; fct++)
			comm_vec[fct] = u[m_algInd[0][fct]];

		comm_vec[2*numFct] = 1.0;
	}
	if (m_algInd[1].size() == numFct)
	{
		for (size_t fct = 0; fct < numFct; fct++)
			comm_vec[numFct+fct] = u[m_algInd[1][fct]];

		comm_vec[2*numFct+1] = 1.0;
	}

	std::vector<number> comm_receive;

	/* TODO: We have a problem here when using GMG:
	 * This method is called by adjust_jacobian() and that one is called after
	 * assembling the jacobian on the base level. If the base solver is gathered
	 * (i.e. executed only on one process), only the executing process will arrive
	 * here and the others will not take part in the communication.
	 * In fact, they are waiting at a defect computation communication point for
	 * the solver which contains the GMG and therefore, this will produce strange
	 * MPI errors, like "MPI ERROR: MPI_ERR_TRUNCATE: message truncated".
	 */

	// communicate: compute the sum over all processes and send it back to all
	com.allreduce(comm_vec, comm_receive, PCL_RO_SUM);

	// now check the result
	if (!comm_receive[2*numFct])
		{UG_THROW("No processor possesses information about the full-dim side interface node!");}
	if (!comm_receive[2*numFct+1])
		{UG_THROW("No processor possesses information about the 1D side interface node!");}

	// in case more than one processor possessed the necessary values: divide by number
	for (size_t fct = 0; fct < numFct; fct++)
	{
		intf_val[0][fct] = comm_receive[fct] / comm_receive[2*numFct];
		intf_val[1][fct] = comm_receive[numFct+fct] / comm_receive[2*numFct+1];
	}
}
else
{
	for (size_t fct = 0; fct < numFct; fct++)
	{
		intf_val[0][fct] = u[m_algInd[0][fct]];
		intf_val[1][fct] = u[m_algInd[1][fct]];
	}
}
#else
	for (size_t fct = 0; fct < numFct; fct++)
	{
		intf_val[0][fct] = u[m_algInd[0][fct]];
		intf_val[1][fct] = u[m_algInd[1][fct]];
	}
#endif
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
	size_t numFct = m_vFct.size();

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
	std::vector<number> intf_val[2];
	fill_sol_at_intf(intf_val, u, dd);

// adjust jacobian parts for one side of the interface,
// then switch roles and do the same thing for the other side
	for (size_t side = 0; side < 2; side++)
	{
		// side is the constraining side; c_side the constrained one
		size_t c_side = (side+1) % 2;

		// loop constrained vertices
		DoFDistribution::traits<Vertex>::const_iterator iter, iterEnd;
		iter = dd->begin<Vertex>(m_ssi[side]);
		iterEnd = dd->end<Vertex>(m_ssi[side]);

		for (; iter != iterEnd; ++iter)
		{
			// get constrained and constraining vertex
			Vertex* constrd = *iter;
			Vertex* constrg = m_constrainerMap[constrd];

			// loop functions
			for (size_t fct = 0; fct < numFct; fct++)
			{
				std::vector<size_t> constrdInd;
				std::vector<size_t> constrgInd;

				// constrained index
				if (!dd->is_def_in_subset(m_vFct[fct], m_ssi[side]))
					{UG_THROW("Function " << m_vFct[fct] << "is not defined on constrained subset " << m_ssi[side] << ".");}
				dd->inner_algebra_indices_for_fct(constrd, constrdInd, false, m_vFct[fct]);

				// constraining index
				int si = dd->subset_handler()->get_subset_index(constrg);
				if (!dd->is_def_in_subset(m_vFct[fct], si))
					{UG_THROW("Function " << m_vFct[fct] << "is not defined on constraining subset " << si << ".");}
				dd->inner_algebra_indices_for_fct(constrg, constrgInd, false, m_vFct[fct]);

				// should not happen, but for debugging purposes:
				UG_ASSERT(constrgInd.size() == 1, "More (or less) than one function index found on a vertex!");
				UG_ASSERT(constrdInd.size() == 1, "More (or less) than one function index found on a vertex!");

				// get defect derivatives for this vertex
				typename vector_type::value_type defDeriv[4];
				constrainedDefectDerivs(defDeriv, u[constrdInd[0]], u[constrgInd[0]],
										intf_val[side][fct], intf_val[c_side][fct]);

				// write Jacobian row for constrained (first: delete all couplings)
				typedef typename matrix_type::row_iterator row_iterator;
				{
					// we put curly brackets around this so that the row_iterator will be removed afterwards
					// otherwise, we get a negative assert from the matrix implementation
					const row_iterator iterEnd = J.end_row(constrdInd[0]);
					for (row_iterator conn = J.begin_row(constrdInd[0]); conn != iterEnd; ++conn)
						conn.value() = 0.0;
				}

				// avoid special cases (where e.g. constrdInd[0] == m_algInd[c_side][fct]))
				// by adding up derivs instead of assigning directly
				J(constrdInd[0], constrdInd[0])	= J(constrdInd[0], constrgInd[0]) = 0.0;
				if (m_algInd[side].size()) J(constrdInd[0], m_algInd[side][fct]) = 0.0;
				if (m_algInd[c_side].size()) J(constrdInd[0], m_algInd[c_side][fct]) = 0.0;

				// as PST is supposedly consistent: set correct values where it is possible
				J(constrdInd[0], constrdInd[0]) += defDeriv[0];
				J(constrdInd[0], constrgInd[0]) += defDeriv[1];
				if (m_algInd[side].size()) J(constrdInd[0], m_algInd[side][fct]) += defDeriv[2];
				if (m_algInd[c_side].size()) J(constrdInd[0], m_algInd[c_side][fct]) += defDeriv[3];

				/*
				// adapt rows for all constraining defects that depend on this constrained
				std::vector<Vertex*>& allConstrainers = m_defectInfluenceMap[constrd];
				for (size_t i = 0; i < allConstrainers.size(); i++)
				{
					Vertex* cv = allConstrainers[i];

					// get index
					std::vector<size_t> cvInd;
					int si = dd->subset_handler()->get_subset_index(cv);
					if (!dd->is_def_in_subset(m_vFct[fct], si))
						{UG_THROW("Function " << m_vFct[fct] << "is not defined on constraining subset " << si << ".");}
					dd->inner_algebra_indices_for_fct(cv, cvInd, false, m_vFct[fct]);
					UG_ASSERT(cvInd.size() == 1, "More (or less) than one function index found on a vertex!");

					// deriv wrt constraining
					J(cvInd[0], constrgInd[0]) += J(cvInd[0], constrdInd[0]) * defDeriv[1];

					// deriv wrt interface vertex on constraining side
					if (m_algInd[side].size()) J(cvInd[0], m_algInd[side][fct]) += J(cvInd[0], constrdInd[0]) * defDeriv[2];

					// deriv wrt interface vertex on constrained side
					if (m_algInd[c_side].size()) J(cvInd[0], m_algInd[c_side][fct]) += J(cvInd[0], constrdInd[0]) * defDeriv[3];

					// deriv wrt constrained vertex
					J(cvInd[0], constrdInd[0]) = 0.0;
				}
				*/
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
	size_t numFct = m_vFct.size();

	// loop constrained vertices on the first side of the interface,
	// then switch roles and do the same thing for the other side
	for (size_t side = 0; side < 2; side++)
	{
		// side is the constraining side; c_side the constrained one
		DoFDistribution::traits<Vertex>::const_iterator iter, iterEnd;
		iter = dd->begin<Vertex>(m_ssi[side]);
		iterEnd = dd->end<Vertex>(m_ssi[side]);

		for (; iter != iterEnd; ++iter)
		{
			// loop functions
			for (size_t fct = 0; fct < numFct; fct++)
			{
				std::vector<size_t> constrdInd;
				//std::vector<size_t> constrgInd;

				// constrained index
				if (!dd->is_def_in_subset(m_vFct[fct], m_ssi[side]))
					{UG_THROW("Function " << m_vFct[fct] << "is not defined on constrained subset " << m_ssi[side] << ".");}
				dd->inner_algebra_indices_for_fct(*iter, constrdInd, false, m_vFct[fct]);

				// better convergence results are yielded
				// by adjusting the solution before computing the defect
				// which automatically reduces the defect to 0
				d[constrdInd[0]] = 0.0;
			}
		}
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
	// in order to activate, you will probably have to comment in the code below
	/*
	size_t numFct = m_vFct.size();

	// loop constrained vertices on the first side of the interface,
	// then switch roles and do the same thing for the other side
	for (size_t side = 0; side < 2; side++)
	{
		// side is the constraining side; c_side the constrained one
		size_t c_side = (side+1) % 2;

		DoFDistribution::traits<Vertex>::const_iterator iter, iterEnd;
		iter = dd->begin<Vertex>(m_ssi[side]);
		iterEnd = dd->end<Vertex>(m_ssi[side]);

		for (; iter != iterEnd; ++iter)
		{
			// get constrained and constraining vertices
			Vertex* constrd = *iter;
			Vertex* constrg = m_constrainerMap[constrd];

			// loop functions
			for (size_t fct = 0; fct < numFct; fct++)
			{
				std::vector<size_t> constrdInd;
				std::vector<size_t> constrgInd;

				// constrained index
				if (!dd->is_def_in_subset(m_vFct[fct], m_ssi[side]))
					{UG_THROW("Function " << m_vFct[fct] << "is not defined on constrained subset " << m_ssi[side] << ".");}
				dd->inner_algebra_indices_for_fct(constrd, constrdInd, false, m_vFct[fct]);

				// constraining index
				int si = dd->subset_handler()->get_subset_index(constrg);
				if (!dd->is_def_in_subset(m_vFct[fct], si))
					{UG_THROW("Function " << m_vFct[fct] << "is not defined on constraining subset " << si << ".");}
				dd->inner_algebra_indices_for_fct(constrg, constrgInd, false, m_vFct[fct]);

				// should not happen, but for debugging purposes:
				UG_ASSERT(constrgInd.size() == 1, "More (or less) than one function index found on a vertex!");
				UG_ASSERT(constrdInd.size() == 1, "More (or less) than one function index found on a vertex!");

				// constrain defect
				constrainedDefect(rhs[constrdIndex[0]], u[constrdInd[0]], u[constrgInd[0]],
								  u[alg_ind[side][fct]], u[alg_ind[c_side][fct]]);
			}
		}
	}
	*/
}



///	sets the constraints in a solution vector
template <typename TDomain, typename TAlgebra>
void IInterface1DFV1<TDomain, TAlgebra>::adjust_solution
(	vector_type& u,
	ConstSmartPtr<DoFDistribution> dd,
	number time
)
{
	size_t numFct = m_vFct.size();

#ifdef UG_PARALLEL
	// check parallel storage type
	if (pcl::NumProcs() > 1 &&
		(!u.has_storage_type(PST_CONSISTENT)))
	{
		UG_THROW("Expected solution to have parallel storage type PST_CONSISTENT,\n"
				 "but it does not. PST is '" << u.get_storage_type() << "'.");
	}
#endif

	std::vector<number> intf_val[2];
	fill_sol_at_intf(intf_val, u, dd);

	// loop constrained vertices on the first side of the interface,
	// then switch roles and do the same thing for the other side
	for (size_t side = 0; side < 2; side++)
	{
		// side is the constraining side; c_side the constrained one
		size_t c_side = (side+1) % 2;

		DoFDistribution::traits<Vertex>::const_iterator iter, iterEnd;
		iter = dd->begin<Vertex>(m_ssi[side]);
		iterEnd = dd->end<Vertex>(m_ssi[side]);

		for (; iter != iterEnd; ++iter)
		{
			// get constrained and constraining vertices
			Vertex* constrd = *iter;
			Vertex* constrg = m_constrainerMap[constrd];

			// loop functions
			for (size_t fct = 0; fct < numFct; fct++)
			{
				std::vector<size_t> constrdInd;
				std::vector<size_t> constrgInd;

				// constrained index
				if (!dd->is_def_in_subset(m_vFct[fct], m_ssi[side]))
					{UG_THROW("Function " << m_vFct[fct] << "is not defined on constrained subset " << m_ssi[side] << ".");}
				dd->inner_algebra_indices_for_fct(constrd, constrdInd, false, m_vFct[fct]);

				// constraining index
				int si = dd->subset_handler()->get_subset_index(constrg);
				if (!dd->is_def_in_subset(m_vFct[fct], si))
					{UG_THROW("Function " << m_vFct[fct] << "is not defined on constraining subset " << si << ".");}
				dd->inner_algebra_indices_for_fct(constrg, constrgInd, false, m_vFct[fct]);

				// should not happen, but for debugging purposes:
				UG_ASSERT(constrgInd.size() == 1, "More (or less) than one function index found on a vertex!");
				UG_ASSERT(constrdInd.size() == 1, "More (or less) than one function index found on a vertex!");

				// constrain defect
				typename vector_type::value_type def;
				constrainedDefect(def, u[constrdInd[0]], u[constrgInd[0]],
								  intf_val[side][fct], intf_val[c_side][fct]);

				u[constrdInd[0]] += def;
			}
		}
	}
}

/// additional linear adjustment for use inside linear solver
template <typename TDomain, typename TAlgebra>
void IInterface1DFV1<TDomain, TAlgebra>::adjust_solution_linear
(	vector_type& u,
	ConstSmartPtr<DoFDistribution> dd,
	SmartPtr<MatrixOperator<matrix_type, vector_type> > J,
	number time
)
{
	size_t numFct = m_vFct.size();

#ifdef UG_PARALLEL
	// check parallel storage type
	if (pcl::NumProcs() > 1 &&
		(!u.has_storage_type(PST_CONSISTENT)))
	{
		UG_THROW("Expected solution to have parallel storage type PST_CONSISTENT,\n"
				 "but it does not. PST is '" << u.get_storage_type() << "'.");
	}
#endif

	std::vector<number> intf_val[2];
	fill_sol_at_intf(intf_val, u, dd);

	// loop constrained vertices on the first side of the interface,
	// then switch roles and do the same thing for the other side
	for (size_t side = 0; side < 2; side++)
	{
		// side is the constraining side; c_side the constrained one
		size_t c_side = (side+1) % 2;

		DoFDistribution::traits<Vertex>::const_iterator iter, iterEnd;
		iter = dd->begin<Vertex>(m_ssi[side]);
		iterEnd = dd->end<Vertex>(m_ssi[side]);

		for (; iter != iterEnd; ++iter)
		{
			// get constrained and constraining vertices
			Vertex* constrd = *iter;
			Vertex* constrg = m_constrainerMap[constrd];

			// loop functions
			for (size_t fct = 0; fct < numFct; fct++)
			{
				std::vector<size_t> constrdInd;
				std::vector<size_t> constrgInd;

				// constrained index
				if (!dd->is_def_in_subset(m_vFct[fct], m_ssi[side]))
					{UG_THROW("Function " << m_vFct[fct] << "is not defined on constrained subset " << m_ssi[side] << ".");}
				dd->inner_algebra_indices_for_fct(constrd, constrdInd, false, m_vFct[fct]);

				// constraining index
				int si = dd->subset_handler()->get_subset_index(constrg);
				if (!dd->is_def_in_subset(m_vFct[fct], si))
					{UG_THROW("Function " << m_vFct[fct] << "is not defined on constraining subset " << si << ".");}
				dd->inner_algebra_indices_for_fct(constrg, constrgInd, false, m_vFct[fct]);

				// should not happen, but for debugging purposes:
				UG_ASSERT(constrgInd.size() == 1, "More (or less) than one function index found on a vertex!");
				UG_ASSERT(constrdInd.size() == 1, "More (or less) than one function index found on a vertex!");

				// constrain defect
				// get defect derivatives for this vertex
				matrix_type& M = J->get_matrix();

				// write Jacobian row for constrained (first: delete all couplings)
				u[constrdInd[0]] = M(constrdInd[0], constrgInd[0])*u[constrgInd[0]];
				if (m_algInd[side].size())
					u[constrdInd[0]] += M(constrdInd[0], m_algInd[side][fct])*intf_val[side][fct];
				if (m_algInd[c_side].size())
					u[constrdInd[0]] += M(constrdInd[0], m_algInd[c_side][fct])*intf_val[c_side][fct];
			}
		}
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
void AdditiveInterface1DFV1<TDomain, TAlgebra>::constrainedDefect
(	typename vector_type::value_type& d,
	const typename vector_type::value_type& u,
	const typename vector_type::value_type& u_c,
	const typename vector_type::value_type& u_itf0,
	const typename vector_type::value_type& u_itf1
)
{
	d = -u + u_c + (u_itf1 - u_itf0);
}

/// function returning the defect derivatives for a constrained node
template <typename TDomain, typename TAlgebra>
void AdditiveInterface1DFV1<TDomain, TAlgebra>::constrainedDefectDerivs
(	typename vector_type::value_type dd[4],
	const typename vector_type::value_type& u,
	const typename vector_type::value_type& u_c,
	const typename vector_type::value_type& u_itf0,
	const typename vector_type::value_type& u_itf1
)
{
	dd[0] = -1.0;
	dd[1] =  1.0;
	dd[2] = -1.0;
	dd[3] =  1.0;
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
void MultiplicativeInterface1DFV1<TDomain, TAlgebra>::constrainedDefect
(	typename vector_type::value_type& d,
	const typename vector_type::value_type& u,
	const typename vector_type::value_type& u_c,
	const typename vector_type::value_type& u_itf0,
	const typename vector_type::value_type& u_itf1
)
{
	if (std::fabs(u_itf0) < 2 * std::numeric_limits<typename vector_type::value_type>::denorm_min())
		{UG_THROW("Denominator practically zero (" << u_itf0 << ").");}

	d = -u + u_c * (u_itf1 / u_itf0);
}

/// function returning the defect derivatives for a constrained node
template <typename TDomain, typename TAlgebra>
void MultiplicativeInterface1DFV1<TDomain, TAlgebra>::constrainedDefectDerivs
(	typename vector_type::value_type dd[4],
	const typename vector_type::value_type& u,
	const typename vector_type::value_type& u_c,
	const typename vector_type::value_type& u_itf0,
	const typename vector_type::value_type& u_itf1
)
{
	if (std::fabs(u_itf0) < 2 * std::numeric_limits<typename vector_type::value_type>::denorm_min())
		{UG_THROW("Denominator practically zero.");}

	dd[0] = -1.0;
	dd[1] =  u_itf1 / u_itf0;
	dd[2] = -u_c * u_itf1 / (u_itf0*u_itf0);
	dd[3] =  u_c / u_itf0;
}

} // namespace nernst_planck
} // namespace ug
