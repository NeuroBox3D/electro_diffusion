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
	const char* one_dim_subset
)
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


/// called when the approximation space has changed
template <typename TDomain, typename TAlgebra>
void IInterface1DFV1<TDomain, TAlgebra>::approximation_space_changed()
{
// get fct indices
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
	m_ssi[0] = this->approximation_space()->domain()->subset_handler()->get_subset_index(m_sssi[0].c_str());
	m_ssi[1] = this->approximation_space()->domain()->subset_handler()->get_subset_index(m_sssi[1].c_str());

	if (m_ssi[0] < 0 || m_ssi[1] < 0)	// previous call gives -1 if failed
	{
		UG_THROW("At least one of the given subsets is not contained in the "
				  "geometry handled by this domain's subset handler");
	}


// fill helper structures
	fill_constrainer_map();
	fill_defect_influence_map();


// find algebra indices for interface nodes
// The way this is done here is as follows:
// 1) Find the constrained vertex on the 1d side uniquely identified by one_dim_subset.
// 2) Get the unique vertex connected by an edge to it: This is the interface node on the 1d side.
// 3) Find the constrained vertex on the high-dimensional side that is nearest to the 1d interface node.
//	  This is the constrained vertex belonging to the high-dimensional interface node
//	  (should have _exactly_ the same coordinates, but you never know...).
	typedef typename MultiGrid::traits<Edge>::secure_container edge_list;
	ConstSmartPtr<DoFDistribution> dd = this->approximation_space()->dof_distribution(GridLevel());

	// 1)
	Vertex* iv1;
	Vertex* iv2 = NULL;
	DoFDistribution::traits<Vertex>::const_iterator iter;
	iter = dd->begin<Vertex>(m_ssi[1]);

	if (iter == dd->end<Vertex>(m_ssi[1]))
		{UG_THROW("No vertex in constrained subset for 1d side. This is not allowed!");}

	Vertex* constrd = *iter;

	// 2)
	edge_list el;
	this->approximation_space()->domain()->grid()->associated_elements(el, constrd);

	if (el.size() != 1)
	{
		UG_THROW("Incorrect number of connected edges for constrained vertex on 1d end."
				 "There should be exactly 1, but " << el.size() << " found.\n"
				 "Correct the error in your geometry!");
	}

	el[0]->get_opposing_side(constrd, &iv1);

	// to be on the safe side
	if (++iter != dd->end<Vertex>(m_ssi[1]))
			{UG_THROW("More than one vertex in constrained subset for 1d side. This is not allowed!");}

	// 3)
	typedef MathVector<TDomain::dim> pos_type;
	typedef Grid::VertexAttachmentAccessor<Attachment<pos_type> > pos_acc_type;
	pos_acc_type pos_acc = this->approximation_space()->domain()->position_accessor();

	pos_type pos_ifn1 = pos_acc[iv1];

	// loop through constrained vertices on high-dim side

	iter = dd->begin<Vertex>(m_ssi[0]);
	DoFDistribution::traits<Vertex>::const_iterator iterEnd = dd->end<Vertex>(m_ssi[0]);

	number dist;
	Vertex* nearest;
	if (iter != iterEnd)
	{
		nearest = *iter;
		pos_type diff;
		VecSubtract(diff, pos_ifn1, pos_acc[nearest]);
		dist = VecTwoNormSq(diff);
		++iter;
	}
	else {UG_THROW("Not even one constrained vertex on the high-dimensional side. This is not admissible.");}

	for (; iter != iterEnd; ++iter)
	{
		pos_type diff;
		VecSubtract(diff, pos_ifn1, pos_acc[*iter]);
		number thisDist = VecTwoNormSq(diff);
		if (thisDist < dist)
		{
			dist = thisDist;
			nearest = *iter;
		}
	}

	iv2 = m_constrainerMap[nearest];

	// get algebra indices
	int ssi1 = this->approximation_space()->domain()->subset_handler()->get_subset_index(iv1);
	int ssi2 = this->approximation_space()->domain()->subset_handler()->get_subset_index(iv2);

	for (size_t fct = 0; fct < m_vFct.size(); fct++)
	{
		std::vector<size_t> ind1;
		std::vector<size_t> ind2;

		if (!dd->is_def_in_subset(m_vFct[fct], ssi1) || !dd->is_def_in_subset(m_vFct[fct], ssi2))
			{UG_THROW("Function " << m_vFct[fct] << "is not defined on constraining subsets!");}
		dd->inner_algebra_indices_for_fct(iv1, ind1, false, m_vFct[fct]);
		dd->inner_algebra_indices_for_fct(iv2, ind2, false, m_vFct[fct]);

		UG_ASSERT(ind1.size() == 1 && ind2.size() == 1,
				  "More (or less) than one function index found on a vertex!");

		m_algInd[0].push_back(ind2[0]);
		m_algInd[1].push_back(ind1[0]);
	}

}



/// for every constrained vertex: finds the corresponding constrainer
template <typename TDomain, typename TAlgebra>
void IInterface1DFV1<TDomain, TAlgebra>::fill_constrainer_map()
{
	typedef typename MultiGrid::traits<Edge>::secure_container edge_list;

	ConstSmartPtr<DoFDistribution> dd = this->m_spApproxSpace->dof_distribution(GridLevel());

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
			this->approximation_space()->domain()->grid()->associated_elements(el, constrd);
			for (size_t edge = 0; edge < el.size(); edge++)
			{
				Vertex* other;
				el[edge]->get_opposing_side(constrd, &other);

				if (this->approximation_space()->domain()->subset_handler()->get_subset_index(other) != m_ssi[0])
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
			this->approximation_space()->domain()->grid()->associated_elements(el, constrd);

			if (el.size() != 1)
			{
				UG_THROW("Incorrect number of connected edges for constrained vertex on 1d end."
					     "There should be exactly 1, but " << el.size() << " found.\n"
					     "Correct the error in your geometry!");
			}

			Vertex* other;
			el[0]->get_opposing_side(constrd, &other);

			if (this->approximation_space()->domain()->subset_handler()->get_subset_index(other) != m_ssi[1])
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

	ConstSmartPtr<DoFDistribution> dd = this->m_spApproxSpace->dof_distribution(GridLevel());

	// high-dimensional end
	{
		DoFDistribution::traits<Vertex>::const_iterator iter, iterEnd;
		iter = dd->begin<Vertex>(m_ssi[0]);
		iterEnd = dd->end<Vertex>(m_ssi[0]);

		for (; iter != iterEnd; ++iter)
		{
			// get constrained and constraining vertices
			Vertex* constrd = *iter;

		// find the vertices in the constraining subset whose defect depends on this constrained dof
			elem_list el;
			this->approximation_space()->domain()->grid()->associated_elements(el, constrd);

			// remember vertices already pushed back to avoid double entries
			std::map<Vertex*, bool> considered;

			// loop all associated elements
			for (size_t elem = 0; elem < el.size(); elem++)
			{
				// get vertices associated to this elem
				vertex_list otherVertices;
				this->approximation_space()->domain()->grid()->associated_elements(otherVertices, el[elem]);

				// loop vertices
				for (size_t vrt = 0; vrt < otherVertices.size(); vrt++)
				{
					Vertex* other = otherVertices[vrt];
					// only add if subset index matches constrained set and not yet added
					if (this->approximation_space()->domain()->subset_handler()->get_subset_index(other) == m_ssi[0]
						&& !considered[other])
					{
						considered[other] = true;
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
			// get constrained and constraining vertices
			Vertex* constrd = *iter;
			Vertex* constrg = m_constrainerMap[constrd];

			m_defectInfluenceMap[constrd].push_back(constrg);
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



/// must be called when the underlying geometry is changed (if this affects the interface)
template <typename TDomain, typename TAlgebra>
void IInterface1DFV1<TDomain, TAlgebra>::geometryChanged()
{
	// TODO: here, at least the constrainerMap will have to be updated!
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
						  	  	  	  	u[m_algInd[side][fct]], u[m_algInd[c_side][fct]]);

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
				J(constrdInd[0], constrdInd[0])
					= J(constrdInd[0], constrgInd[0])
					= J(constrdInd[0], m_algInd[side][fct])
					= J(constrdInd[0], m_algInd[c_side][fct]) = 0.0;

				J(constrdInd[0], constrdInd[0]) += defDeriv[0];
				J(constrdInd[0], constrgInd[0]) += defDeriv[1];
				J(constrdInd[0], m_algInd[side][fct]) += defDeriv[2];
				J(constrdInd[0], m_algInd[c_side][fct]) += defDeriv[3];

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
					J(cvInd[0], m_algInd[side][fct]) += J(cvInd[0], constrdInd[0]) * defDeriv[2];

					// deriv wrt interface vertex on constrained side
					J(cvInd[0], m_algInd[c_side][fct]) += J(cvInd[0], constrdInd[0]) * defDeriv[3];

					// deriv wrt constrained vertex
					J(cvInd[0], constrdInd[0]) = 0.0;
				}
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
				constrainedDefect(d[constrdInd[0]], u[constrdInd[0]], u[constrgInd[0]],
								  u[m_algInd[side][fct]], u[m_algInd[c_side][fct]]);
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
								  u[m_algInd[side][fct]], u[m_algInd[c_side][fct]]);
				u[constrdInd[0]] += def;
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
	const char* one_dim_subset
) : IInterface1DFV1<TDomain, TAlgebra>(fcts, high_dim_subset, one_dim_subset) {}


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
	const char* one_dim_subset
) : IInterface1DFV1<TDomain, TAlgebra>(fcts, high_dim_subset, one_dim_subset) {}


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
		{UG_THROW("Denominator practically zero.");}

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
