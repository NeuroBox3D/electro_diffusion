/*
 * interface1d_fv_impl.h
 *
 *  Created on: 06.06.2014
 *      Author: mbreit
 */

#include "interface1d_fv.h"

namespace ug{
namespace nernst_planck{


template <typename TDomain, typename TAlgebra>
IInterface1D<TDomain, TAlgebra>::Interface1DMapper::Interface1DMapper(SmartPtr<IAssemble<TAlgebra> > ass)
{
	SmartPtr<AssemblingTuner<TAlgebra> > assTuner = ass->ass_tuner();
	assTuner->set_mapping(this);
}


template <typename TDomain, typename TAlgebra>
void IInterface1D<TDomain, TAlgebra>::Interface1DMapper::add_local_vec_to_global
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
					{
						UG_THROW("Function " << fct << " for constrained index " << index
								 << " is unknown in interface " << intf << ".\n");
					}
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
void IInterface1D<TDomain, TAlgebra>::Interface1DMapper::add_local_mat_to_global
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
void IInterface1D<TDomain, TAlgebra>::Interface1DMapper::add_interface(SmartPtr<interface_type> intf)
{
	m_vspInterface.push_back(intf);
}



template <typename TDomain, typename TAlgebra>
IInterface1D<TDomain, TAlgebra>::IInterface1D
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



template <typename TDomain, typename TAlgebra>
IInterface1D<TDomain, TAlgebra>::~IInterface1D()
{
	// do nothing
}



template <typename TDomain, typename TAlgebra>
int IInterface1D<TDomain, TAlgebra>::type() const
{
	return CT_CONSTRAINTS;
}



template <typename TDomain, typename TAlgebra>
void IInterface1D<TDomain, TAlgebra>::
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
void IInterface1D<TDomain, TAlgebra>::
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


template <typename TDomain, typename TAlgebra>
void IInterface1D<TDomain, TAlgebra>::approximation_space_changed()
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

template <typename TDomain, typename TAlgebra>
template <typename TElem, typename TElemDesc, typename TDummy>
IInterface1D<TDomain, TAlgebra>::OrientationOffset<TElem, TElemDesc, TDummy>::
OrientationOffset
(
	std::vector<size_t>& vOrientOffset,
	TElemDesc& target_desc,
	TElem* constrg,
	size_t p
)
{
	UG_THROW("not implemented!")
}

template <typename TDomain, typename TAlgebra>
template <typename TDummy>
IInterface1D<TDomain, TAlgebra>::OrientationOffset<Vertex*, Vertex, TDummy>::
OrientationOffset
(
	std::vector<size_t>& vOrientOffset,
	Vertex& target_desc,
	Vertex* constrg,
	size_t p
)
{
	// do nothing
}

template <typename TDomain, typename TAlgebra>
template <typename TDummy>
IInterface1D<TDomain, TAlgebra>::OrientationOffset<Edge, EdgeDescriptor, TDummy>::
OrientationOffset
(
	std::vector<size_t>& vOrientOffset,
	EdgeDescriptor& target_desc,
	Edge* constrg,
	size_t p
)
{
	ComputeOrientationOffsetLagrange(vOrientOffset, target_desc, constrg, p);
}

template <typename TDomain, typename TAlgebra>
template <typename TDummy>
IInterface1D<TDomain, TAlgebra>::OrientationOffset<Face, FaceDescriptor, TDummy>::
OrientationOffset
(
	std::vector<size_t>& vOrientOffset,
	FaceDescriptor& target_desc,
	Face* constrg,
	size_t p
)
{
	const int n_vrt = target_desc.num_vertices();
	const int id0 = GetVertexIndex(&target_desc, constrg->vertex(0));
	const int id1 = GetVertexIndex(&target_desc, constrg->vertex(1));
	const bool sameOrientation = (id1 == (id0+1)%n_vrt);

	switch (n_vrt)
	{
		case 3:
			MapLagrangeMultiIndexTriangle(vOrientOffset, id0, sameOrientation, p);
			break;
		case 4:
			MapLagrangeMultiIndexQuad(vOrientOffset, id0, sameOrientation, p);
			break;
		default: UG_THROW("No elements with #corners = " << n_vrt << " implemented.");
	}
}


template <typename TDomain, typename TAlgebra>
template <typename TElem, typename TElemDesc, typename TContainingElem, typename TDummy>
IInterface1D<TDomain, TAlgebra>::GetConstrainer<TElem, TElemDesc, TContainingElem, TDummy>::
GetConstrainer(IInterface1D<TDomain, TAlgebra>* const intf, TElem* constrd, TElem** constrg_out)
{
	UG_ASSERT(constrg_out, "Pointer to pointer to out constraining vertex (third argument) must not be NULL.")

	typedef typename MultiGrid::traits<TContainingElem>::secure_container cont_list;
	typedef typename MultiGrid::traits<TElem>::secure_container elem_list;

	SmartPtr<TDomain> dom = intf->approximation_space()->domain();

	// loop elements that have constrd elem as boundary;
	// find the one whose other end is not part of the constrained set
	TElemDesc constrg;
	cont_list cl;
	dom->grid()->associated_elements(cl, constrd);
	size_t cont = 0;
	for (; cont < cl.size(); ++cont)
	{
		if (dom->subset_handler()->get_subset_index(cl[cont]) == intf->m_siConstr)
			continue;

		if (!cl[cont]->get_opposing_side(constrd, constrg))
			{UG_THROW("No opposing side found!");}

		break;
	}

	// ensure that the constrainer has been found
	size_t n_vrt = constrg.num_vertices();
	if (!n_vrt) {UG_THROW("No constrainer found!");}

	// ensure that number of constrainer vertices is the same as number of constrained vertices
	if (n_vrt != constrd->num_vertices())
		UG_THROW("Number of constraining vertices (" << n_vrt << ") does not match number "
				 "of constrained vertices (" << constrd->num_vertices() << ").");

	// get real elem instead of descriptor
	elem_list el;
	dom->grid()->associated_elements(el, cl[cont]);
	for (size_t elem = 0; elem < el.size(); ++elem)
		if (CompareVertices(el[elem], &constrg))
		{
			*constrg_out = el[elem];
			return;
		}

	UG_THROW("No constrainer found!");
}


template <typename TDomain, typename TAlgebra>
template <typename TDummy>
IInterface1D<TDomain, TAlgebra>::GetConstrainer<Vertex, Vertex, Edge, TDummy>::
GetConstrainer(IInterface1D<TDomain, TAlgebra>* const intf, Vertex* constrd, Vertex** constrg_out)
{
	UG_ASSERT(constrg_out, "Pointer to pointer to out constraining vertex (third argument) must not be NULL.")

	typedef typename MultiGrid::traits<Edge>::secure_container edge_list;
	SmartPtr<TDomain> dom = intf->approximation_space()->domain();

	// loop associated edges; find the one not in the constrained set and take other end
	edge_list el;
	dom->grid()->associated_elements(el, constrd);
	for (size_t edge = 0; edge < el.size(); edge++)
	{
		if (dom->subset_handler()->get_subset_index(el[edge]) == intf->m_siConstr)
			continue;

		if (!el[edge]->get_opposing_side(constrd, constrg_out))
			UG_THROW("No opposing side found!");

		return;
	}

	UG_THROW("No constrainer found!");
}


template <typename TDomain, typename TAlgebra>
template <typename TElem, typename TElemDesc, typename TDummy>
IInterface1D<TDomain, TAlgebra>::TargetDescriptor<TElem, TElemDesc, TDummy>::
TargetDescriptor(IInterface1D<TDomain, TAlgebra>* const intf, TElem* constrd, TElemDesc& desc)
{
	size_t n_vrt = constrd->num_vertices();
	for (size_t vrt = 0; vrt < n_vrt; ++vrt)
	{
		// get constrained vertex
		Vertex* constrd_vrt = constrd->vertex(vrt);

		// get associated constrainer vertex
		Vertex* constrg_vrt;
		GetConstrainer<Vertex, Vertex, Edge>(intf, constrd_vrt, &constrg_vrt);

		// set constrg vertex in target descriptor
		desc.set_vertex(vrt, constrg_vrt);
	}
	desc.set_num_vertices(n_vrt);
}


template <typename TDomain, typename TAlgebra>
template <typename TDummy>
IInterface1D<TDomain, TAlgebra>::TargetDescriptor<Vertex, Vertex, TDummy>::
TargetDescriptor(IInterface1D<TDomain, TAlgebra>* const intf, Vertex* constrd, Vertex& desc)
{
	// do nothing
}


template <typename TDomain, typename TAlgebra>
template <typename TDummy>
IInterface1D<TDomain, TAlgebra>::TargetDescriptor<Edge, EdgeDescriptor, TDummy>::
TargetDescriptor(IInterface1D<TDomain, TAlgebra>* const intf, Edge* constrd, EdgeDescriptor& desc)
{
	for (size_t vrt = 0; vrt < 2; ++vrt)
	{
		// get constrained vertex
		Vertex* constrd_vrt = constrd->vertex(vrt);

		// get associated constrainer vertex
		Vertex* constrg_vrt;
		GetConstrainer<Vertex, Vertex, Edge>(intf, constrd_vrt, &constrg_vrt);

		// set constrg vertex in target descriptor
		desc.set_vertex(vrt, constrg_vrt);
	}
}

template <typename TDomain, typename TAlgebra>
template <typename TElem, typename TElemDesc, typename TContainingElem>
void IInterface1D<TDomain, TAlgebra>::fill_constraint_map()
{
	ConstSmartPtr<DoFDistribution> dd = this->m_spApproxSpace->dof_distribution(GridLevel());

	// loop constrained elements
	typename DoFDistribution::traits<TElem>::const_iterator iter, iterEnd;
	iter = dd->begin<TElem>(m_siConstr);
	iterEnd = dd->end<TElem>(m_siConstr);

	for (; iter != iterEnd; ++iter)
	{
		// get constrained elem
		TElem* constrd = *iter;

		// get constrainer
		TElem* constrg;
		GetConstrainer<TElem, TElemDesc, TContainingElem>(this, constrd, &constrg);

		// preparations to find out numbering of constraining descriptor vertices relative
		// to numbering of constrainers of vertices of constrained element
		TElemDesc target_desc;
		TargetDescriptor<TElem, TElemDesc>(this, constrd, target_desc);

		// find indices
		// loop functions
		size_t numFct = m_vFct.size();
		for (size_t fct = 0; fct < numFct; fct++)
		{
			const size_t fct_ind = m_vFct[fct];
			std::vector<size_t> orientationOffsets;
			LFEID lfeid = dd->lfeid(fct_ind);
			switch (lfeid.type())
			{
				case LFEID::LAGRANGE:
					// only calculate orientationOffsets for p > 2, (else only max 1 DoF on sub)
					if (lfeid.order() <= 2) break;

					OrientationOffset<TElem, TElemDesc>(orientationOffsets, target_desc, constrg, lfeid.order());
					break;

				default:
					UG_THROW("This interface is only implemented for Lagrange shape functions.");
			}

			// get inner DoF indices
			std::vector<DoFIndex> constrdInd;
			std::vector<DoFIndex> constrgInd;
			dd->inner_dof_indices(constrd, fct_ind, constrdInd, false);
			dd->inner_dof_indices(constrg, fct_ind, constrgInd, false);

			// should not happen, but for debugging purposes:
			UG_ASSERT(constrdInd.size() ==
						dd->num_fct_dofs(fct_ind, constrd->reference_object_id(), m_siConstr),
					  "Incorrect number of DoFs for fct " << fct_ind << " on subset " << m_siConstr << " with roid "
					  << constrd->reference_object_id() << ": " << constrdInd.size() << " instead of "
					  << dd->num_fct_dofs(fct_ind, constrd->reference_object_id(), m_siConstr));
			UG_ASSERT(constrgInd.size() ==
						dd->num_fct_dofs(fct_ind, constrg->reference_object_id(), dd->subset_handler()->get_subset_index(constrg)),
					  "Incorrect number of DoFs for fct " << fct_ind << " on subset "
					  << dd->subset_handler()->get_subset_index(constrg) << " with roid "
					  << constrg->reference_object_id() << ": " << constrgInd.size() << " instead of "
					  << dd->num_fct_dofs(fct_ind, constrg->reference_object_id(), dd->subset_handler()->get_subset_index(constrg)));

			UG_COND_THROW(constrgInd.size() != constrdInd.size(),
						  "Constraining and constrained elements do not have the same number of DoFs "
						  "(" << constrgInd.size() << "/" << constrdInd.size() << ")");

			// fill map with pairs of indices
			size_t n_ind = constrgInd.size();
			if (orientationOffsets.size())
			{
				UG_COND_THROW(orientationOffsets.size() != n_ind,
							  "Orientation offset vector does not have the required size "
							  "( " << orientationOffsets.size() << " instead of " << n_ind << ")");

				for (size_t i = 0; i < n_ind; ++i)
				{
					UG_ASSERT(!constrdInd[i][1], "Block matrices are not supported (yet)!");
					UG_ASSERT(!constrgInd[orientationOffsets[i]][1], "Block matrices are not supported (yet)!");
					m_constraintMap[constrdInd[i][0]] = ConstraintInfo(constrgInd[orientationOffsets[i]][0], fct);
				}
			}
			else
			{
				for (size_t i = 0; i < n_ind; ++i)
				{
					UG_ASSERT(!constrdInd[i][1], "Block matrices are not supported (yet)!");
					UG_ASSERT(!constrgInd[i][1], "Block matrices are not supported (yet)!");
					m_constraintMap[constrdInd[i][0]] = ConstraintInfo(constrgInd[i][0], fct);
				}
			}
		}
	}
}



template <typename TDomain, typename TAlgebra>
void IInterface1D<TDomain, TAlgebra>::fill_constraint_map()
{
	ConstSmartPtr<DoFDistribution> dd = this->m_spApproxSpace->dof_distribution(GridLevel());

	if (worldDim-1 >= VERTEX && dd->max_dofs(VERTEX) > 0)
		fill_constraint_map<Vertex, Vertex, Edge>();
	if (worldDim-1 >= EDGE && dd->max_dofs(EDGE) > 0)
		fill_constraint_map<Edge, EdgeDescriptor, Face>();
	if (worldDim-1 >= FACE && dd->max_dofs(FACE) > 0)
		fill_constraint_map<Face, FaceDescriptor, Volume>();

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


template <typename TDomain, typename TAlgebra>
void IInterface1D<TDomain, TAlgebra>::adjust_jacobian
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



template <typename TDomain, typename TAlgebra>
void IInterface1D<TDomain, TAlgebra>::adjust_defect
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



template <typename TDomain, typename TAlgebra>
void IInterface1D<TDomain, TAlgebra>::adjust_linear
(	matrix_type& mat,
	vector_type& rhs,
	ConstSmartPtr<DoFDistribution> dd,
	number time
)
{
	UG_THROW("This feature is not implemented.");
}



template <typename TDomain, typename TAlgebra>
void IInterface1D<TDomain, TAlgebra>::adjust_rhs
(	vector_type& rhs,
	const vector_type& u,
	ConstSmartPtr<DoFDistribution> dd,
	number time
)
{
	UG_THROW("This feature is deactivated in order to test, whether it is needed in the first place.");
}



template <typename TDomain, typename TAlgebra>
void IInterface1D<TDomain, TAlgebra>::adjust_solution
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


} // namespace nernst_planck
} // namespace ug
