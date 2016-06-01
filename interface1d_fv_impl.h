/*
 * interface1d_fv_impl.h
 *
 *  Created on: 06.06.2014
 *      Author: mbreit
 */

#include "interface1d_fv.h"

#include "lib_grid/tools/surface_view.h"
#include "lib_grid/algorithms/element_side_util.h"	// GetOpposingSide
#ifdef UG_PARALLEL
	#include "lib_grid/parallelization/distributed_grid.h"
#endif

namespace ug{
namespace nernst_planck{


#if 0
template <typename TDomain, typename TAlgebra>
Interface1D<TDomain, TAlgebra>::Interface1DMapper::Interface1DMapper(SmartPtr<IAssemble<TAlgebra> > ass)
{
	SmartPtr<AssemblingTuner<TAlgebra> > assTuner = ass->ass_tuner();
	assTuner->set_mapping(this);
}


template <typename TDomain, typename TAlgebra>
void Interface1D<TDomain, TAlgebra>::Interface1DMapper::add_local_vec_to_global
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
				if (m_vspInterface[intf]->m_constrainerMap.find(index)
					!= m_vspInterface[intf]->m_constrainerMap.end())
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
void Interface1D<TDomain, TAlgebra>::Interface1DMapper::add_local_mat_to_global
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
				if (m_vspInterface[intf]->m_constrainerMap.find(r_index)
					!= m_vspInterface[intf]->m_constrainerMap.end())
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
void Interface1D<TDomain, TAlgebra>::Interface1DMapper::add_interface(SmartPtr<interface_type> intf)
{
	m_vspInterface.push_back(intf);
}

#endif


template <typename TDomain, typename TAlgebra>
Interface1D<TDomain, TAlgebra>::Interface1D
(
	const char* fcts,
	const char* constrained,
	const char* high_dim_intfNode,
	const char* one_dim_intfNode,
	std::vector<number> dir
)
	: m_siConstr(-1)
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

	// store interface direction
	UG_COND_THROW(dir.size() < worldDim, "Given direction vector does not have enough components "
				  " (" << dir.size() << " instead of " << worldDim <<").");

	for (size_t i = 0; i < worldDim; ++i)
		m_direction[i] = dir[i];

	VecNormalize(m_direction, m_direction);
}



template <typename TDomain, typename TAlgebra>
Interface1D<TDomain, TAlgebra>::~Interface1D()
{
	// do nothing
}



template <typename TDomain, typename TAlgebra>
int Interface1D<TDomain, TAlgebra>::type() const
{
	return CT_MAY_DEPEND_ON_HANGING | CT_CONSTRAINTS;
}



template <typename TDomain, typename TAlgebra>
void Interface1D<TDomain, TAlgebra>::
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
void Interface1D<TDomain, TAlgebra>::approximation_space_changed()
{
// get fct indices
	// get domain for later use
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

	update();
}


template <typename TDomain, typename TAlgebra>
void Interface1D<TDomain, TAlgebra>::update()
{
	std::vector<SmartPtr<DoFDistribution> > vspdd = this->approximation_space()->dof_distributions();
	size_t sz = vspdd.size();

	for (size_t i = 0; i < sz; ++i)
		update(vspdd[i]);
}



template <typename TDomain, typename TAlgebra>
void Interface1D<TDomain, TAlgebra>::update(ConstSmartPtr<DoFDistribution> dd)
{
	bool isSurfDD = dd->grid_level().is_surface();

// find algebra indices for interface nodes
	Vertex* iv1 = NULL;	// high-dim interface node
	Vertex* iv2 = NULL;	// one-dim interface node
	DoFDistribution::traits<Vertex>::const_iterator iter;

	iter = dd->begin<Vertex>(m_siIntf[0]);
	if (iter == dd->end<Vertex>(m_siIntf[0]))
	{
#ifndef UG_PARALLEL
		if (isSurfDD)
			UG_THROW("No vertex in subset for high-dimensional interface node. This is not allowed!");
#else
		if (isSurfDD && pcl::NumProcs() <= 1)
			UG_THROW("No vertex in subset for high-dimensional interface node. This is not allowed!");

		// else do nothing
#endif
	}
	else
	{
		iv1 = *iter;

		// TODO: This might actually happen if interface node is on surface_rim/surface_shadow
		if (++iter != dd->end<Vertex>(m_siIntf[0]))
			{UG_THROW("More than one vertex in subset for high-dimensional interface node. This is not allowed!");}
	}

	iter = dd->begin<Vertex>(m_siIntf[1]);
	if (iter == dd->end<Vertex>(m_siIntf[1]))
	{
#ifndef UG_PARALLEL
		if (isSurfDD)
			UG_THROW("No vertex in subset for one-dimensional interface node. This is not allowed!");
#else
		if (isSurfDD && pcl::NumProcs() <= 1)
		{
			UG_THROW("No vertex in subset for one-dimensional interface node. This is not allowed!");
		}
		// else do nothing
#endif
	}
	else
	{
		iv2 = *iter;

		// TODO: This might actually happen if interface node is on surface_rim/surface_shadow
		if (++iter != dd->end<Vertex>(m_siIntf[1]))
			{UG_THROW("More than one vertex in subset for one-dimensional interface node. This is not allowed!");}
	}

	/*
	// both interface nodes must be on the same processor
	if (iv2 && !iv1)
	{
		UG_THROW("This processor owns the 1d interface node of the interface\n"
				"but not the high-dim interface node.\n"
				"This is not allowed (as it will typically engender convergence problems).");
	}
	*/

	ConstraintInfo& cInfo = m_mConstraints[dd.get()];

	// get algebra indices
	cInfo.algInd[0].clear();
	cInfo.algInd[1].clear();
	if (iv1 != NULL && iv2 != NULL)
	{
		// get subset handler
		ConstSmartPtr<MGSubsetHandler> ssh = this->approximation_space()->subset_handler();

		// check that interface nodes are not hanging
		// (CANNOT happen for 1d end; MUST NOT happen for full-d end)
		if (iv1->is_constrained())
			UG_THROW("Interface nodes must not be hanging, but " << worldDim << "d interface node is.");
		if (iv2->is_constrained())
			UG_THROW("Interface nodes must not be hanging, but 1d interface node is.")

		int ssi1 = ssh->get_subset_index(iv1);
		for (size_t fct = 0; fct < m_vFct.size(); fct++)
		{
			std::vector<size_t> ind1;

			if (!dd->is_def_in_subset(m_vFct[fct], ssi1))
			{
				UG_THROW("Function " << m_vFct[fct] << "is not defined on interface node on the high-dimensional side!");
			}
			dd->inner_algebra_indices_for_fct(iv1, ind1, false, m_vFct[fct]);

			UG_ASSERT(ind1.size() == 1, "More (or less) than one function index found on a vertex!");

			cInfo.algInd[0].push_back(ind1[0]);
		}

		int ssi2 = ssh->get_subset_index(iv2);
		for (size_t fct = 0; fct < m_vFct.size(); fct++)
		{
			std::vector<size_t> ind2;

			if (!dd->is_def_in_subset(m_vFct[fct], ssi2))
			{
				UG_THROW("Function " << m_vFct[fct] << "is not defined on interface node on the one-dimensional side!");
			}
			dd->inner_algebra_indices_for_fct(iv2, ind2, false, m_vFct[fct]);

			UG_ASSERT(ind2.size() == 1, "More (or less) than one function index found on a vertex!");

			cInfo.algInd[1].push_back(ind2[0]);
		}
	}

// find constrainers for every constrained vertex
	fill_constrainer_map(dd);

	// disallow constrained nodes without having the interface nodes on the same processor (for surface dd)
	if (isSurfDD && cInfo.constrainerMap.size() && !(iv1 && iv2))
	{
		UG_THROW("Processor has constrained nodes but not the necessary interface nodes. This is not allowed.\n"
				 "Make sure the constrained nodes are not separated from their interface nodes by the partitioning.");
	}
}


template <typename TDomain, typename TAlgebra>
template <typename TElem, typename TElemDesc, typename TDummy>
Interface1D<TDomain, TAlgebra>::OrientationOffset<TElem, TElemDesc, TDummy>::
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
Interface1D<TDomain, TAlgebra>::OrientationOffset<Vertex, Vertex, TDummy>::
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
Interface1D<TDomain, TAlgebra>::OrientationOffset<Edge, EdgeDescriptor, TDummy>::
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
Interface1D<TDomain, TAlgebra>::OrientationOffset<Face, FaceDescriptor, TDummy>::
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
template <typename TElem, typename TContainingElem>
TElem* Interface1D<TDomain, TAlgebra>::get_constrainer(TElem* constrd)
{
	typedef typename MultiGrid::traits<TContainingElem>::secure_container cont_list;

	SmartPtr<TDomain> dom = this->approximation_space()->domain();

	const typename TDomain::position_accessor_type& aaPos = dom->position_accessor();
	const typename TDomain::position_type& constrdPos = CalculateCenter(constrd, aaPos);

	// loop elements that have constrd elem as boundary;
	// find the one from whose other end the constrained elem is located
	// in the specified direction of the interface
	TElem* constrg = NULL;
	cont_list cl;
	dom->grid()->associated_elements(cl, constrd);
	size_t cont = 0;
	for (; cont < cl.size(); ++cont)
	{
		// continue as long as edge does not point into interface direction
		TElem* opp = GetOpposingSide(*dom->grid(), cl[cont], constrd);
		if (!opp) continue;

		typename TDomain::position_type oppPos = CalculateCenter(opp, aaPos);
		VecSubtract(oppPos, constrdPos, oppPos);
		VecNormalize(oppPos, oppPos);
		if (VecProd(oppPos, m_direction) < 0.95)	// that is about 18 degrees off
			continue;

		//if (dom->subset_handler()->get_subset_index(cl[cont]) == intf->m_siConstr)
		//	continue;

		//if (!cl[cont]->get_opposing_side(constrd, constrg))
		//	{UG_THROW("No opposing side found!");}

		constrg = opp;
		break;
	}

	//UG_COND_THROW(!constrg, "No constrainer found!");

	// ensure that number of constrainer vertices is the same as number of constrained vertices
	//size_t n_vrt = constrg->num_vertices();
	//if (n_vrt != constrd->num_vertices())
	//	UG_THROW("Number of constraining vertices (" << n_vrt << ") does not match number "
	//			 "of constrained vertices (" << constrd->num_vertices() << ").");

	 return constrg;
}

template <typename TDomain, typename TAlgebra>
GridObject* Interface1D<TDomain, TAlgebra>::get_constrainer_object(GridObject* constrd)
{
	switch (constrd->base_object_id())
	{
		case VERTEX:
		{
			Vertex* v = dynamic_cast<Vertex*>(constrd);
			UG_COND_THROW(!v, "Constrained object has base type VERTEX but cannot be cast to it.");
			return get_constrainer<Vertex, Edge>(v);
		}
		case EDGE:
		{
			Edge* e = dynamic_cast<Edge*>(constrd);
			UG_COND_THROW(!e, "Constrained object has base type EDGE but cannot be cast to it.");
			return get_constrainer<Edge, Face>(e);
		}
		case FACE:
		{
			Face* f = dynamic_cast<Face*>(constrd);
			UG_COND_THROW(!f, "Constrained object has base type Face but cannot be cast to it.")
			return get_constrainer<Face, Volume>(f);
		}
		default: UG_THROW("Invalid constrainer base object type " << constrd->base_object_id());
	}
}


#if 0
template <typename TDomain, typename TAlgebra>
template <typename TDummy>
Interface1D<TDomain, TAlgebra>::GetConstrainer<Vertex, Vertex, Edge, TDummy>::
GetConstrainer(Interface1D<TDomain, TAlgebra>* const intf, Vertex* constrd, Vertex** constrg_out)
{
	UG_ASSERT(constrg_out, "Pointer to pointer to out constraining vertex (third argument) must not be NULL.")

	typedef typename MultiGrid::traits<Edge>::secure_container edge_list;
	SmartPtr<TDomain> dom = intf->approximation_space()->domain();

	const typename TDomain::position_accessor_type& aaPos = dom->position_accessor();
	const typename TDomain::position_type& constrdPos = CalculateCenter(constrd, aaPos);

	// loop associated edges; find the one not in the constrained set and take other end
	edge_list el;
	dom->grid()->associated_elements(el, constrd);
	for (size_t edge = 0; edge < el.size(); edge++)
	{
		// continue as long as edge does not point into interface direction
		Vertex* opp;
		if (!el[edge]->get_opposing_side(constrd, &opp)) continue;

		typename TDomain::position_type oppPos = CalculateCenter(opp, aaPos);
		VecSubtract(oppPos, constrdPos, oppPos);
		VecNormalize(oppPos, oppPos);
		if (VecProd(oppPos, m_direction) < 0.95)	// that is about 18 degrees off
			continue;

		/*
		if (dom->subset_handler()->get_subset_index(el[edge]) == intf->m_siConstr
			|| dom->subset_handler()->get_subset_index(el[edge]) == 20		// FIXME: Hack to exclude useless! Do it properly.
			|| dom->subset_handler()->get_subset_index(el[edge]) == intf->m_siIntf[1]) // exclude 1d intf node
			continue;
		*/

		if (!el[edge]->get_opposing_side(constrd, constrg_out))
			UG_THROW("No opposing side found!");

		return;
	}

	UG_THROW("No constrainer found!");
}
#endif

template <typename TDomain, typename TAlgebra>
template <typename TElem, typename TElemDesc, typename TDummy>
Interface1D<TDomain, TAlgebra>::TargetDescriptor<TElem, TElemDesc, TDummy>::
TargetDescriptor(Interface1D<TDomain, TAlgebra>* const intf, TElem* constrd, TElemDesc& desc)
{
	size_t n_vrt = constrd->num_vertices();
	for (size_t vrt = 0; vrt < n_vrt; ++vrt)
	{
		// get constrained vertex
		Vertex* constrd_vrt = constrd->vertex(vrt);

		// get associated constrainer vertex
		Vertex* constrg_vrt = intf->get_constrainer<Vertex, Edge>(constrd_vrt);
		UG_COND_THROW(!constrg_vrt, "Constrainer not found.");

		// set constrg vertex in target descriptor
		desc.set_vertex(vrt, constrg_vrt);
	}
	desc.set_num_vertices(n_vrt);
}


template <typename TDomain, typename TAlgebra>
template <typename TDummy>
Interface1D<TDomain, TAlgebra>::TargetDescriptor<Vertex, Vertex, TDummy>::
TargetDescriptor(Interface1D<TDomain, TAlgebra>* const intf, Vertex* constrd, Vertex& desc)
{
	// do nothing
}


template <typename TDomain, typename TAlgebra>
template <typename TDummy>
Interface1D<TDomain, TAlgebra>::TargetDescriptor<Edge, EdgeDescriptor, TDummy>::
TargetDescriptor(Interface1D<TDomain, TAlgebra>* const intf, Edge* constrd, EdgeDescriptor& desc)
{
	for (size_t vrt = 0; vrt < 2; ++vrt)
	{
		// get constrained vertex
		Vertex* constrd_vrt = constrd->vertex(vrt);

		// get associated constrainer vertex
		Vertex* constrg_vrt = intf->get_constrainer<Vertex, Edge>(constrd_vrt);
		UG_COND_THROW(!constrg_vrt, "Constrainer not found.");

		// set constrg vertex in target descriptor
		desc.set_vertex(vrt, constrg_vrt);
	}
}

template <typename TDomain, typename TAlgebra>
template <typename TElem, typename TElemDesc, typename TContainingElem>
void Interface1D<TDomain, TAlgebra>::fill_constrainer_map(ConstSmartPtr<DoFDistribution> dd)
{
	ConstSmartPtr<SurfaceView> sv = dd->surface_view();
	SmartPtr<MultiGrid> mg = this->m_spApproxSpace->domain()->grid();

	// the constrainer map for this dd
	std::map<size_t, ConstrainerInfo>& constrainerMap = m_mConstraints[dd.get()].constrainerMap;

	// loop constrained elements
	typename DoFDistribution::traits<TElem>::const_iterator iter, iterEnd;
	iter = dd->begin<TElem>(m_siConstr);
	iterEnd = dd->end<TElem>(m_siConstr);

	for (; iter != iterEnd; ++iter)
	{
		// get constrained elem
		TElem* constrd = *iter;

		// check that elem is not hanging
		// (in which case this elem is to be constrained by hanging constraint)
		// TODO: This will only work properly for P1 Lagrange approx spaces.
		if ((*iter)->is_constrained()) continue;


		/* // This is not needed with anisotropic refinement of the interface.
		// on surface rim, this elem might have hanging constrainers;
		// take shadow rim element as constrained instead
		bool isSrfRim = dd->grid_level().is_surface()
						&& sv->surface_state(constrd).contains(SurfaceView::MG_SURFACE_RIM);

		if (isSrfRim)
		{
			constrd = dynamic_cast<TElem*>(mg->get_parent(constrd));
			UG_COND_THROW(!constrd, "Surface_rim constrained does not have a parent (of the correct base type).");
		}
		*/

		// get constrainer
		TElem* constrg = get_constrainer<TElem, TContainingElem>(constrd);

		// might happen in adaptive multi-grid
		if (!constrg)
		{
			// check that no full-dim element is connected to constrained
			// otherwise throw
			typedef typename domain_traits<TDomain::dim>::element_type elem_type;
			typedef typename MultiGrid::traits<elem_type>::secure_container elem_list;
			elem_list el;
			mg->associated_elements(el, constrd);
			UG_COND_THROW(el.size(), "Constrained element without constrainer "
									 "even though full-dim element is connected.");

			// now ignore
			continue;
		}

		// in case of surface rim constrained: use respective children
		/*
		if (isSrfRim)
		{
			UG_COND_THROW(mg->num_children<TElem>(constrd) != 1,
						  "Not exactly one child for shadow rim constrained element.");
			if (!sv->surface_state(constrg).partially_contains(SurfaceView::MG_SHADOW))
			{
#ifdef UG_PARALLEL
				PCL_DEBUG_BARRIER_ALL();
#endif

				typedef typename TDomain::position_accessor_type AAPos;
				AAPos aaPos = this->m_spApproxSpace->domain()->position_accessor();

				UG_THROW("Constrainer of shadow_rim is not shadowed."
						<< "\nPosition constrained: "<<aaPos[(Vertex*)constrd]
						<< "\nPosition constraining: "<<aaPos[(Vertex*)constrg]);
			}
			UG_COND_THROW(mg->num_children<TElem>(constrg) != 1,
						  "Not exactly one child for shadow rim constrainer.");

			constrd = mg->get_child<TElem>(constrd, 0);
			constrg = mg->get_child<TElem>(constrg, 0);

			UG_COND_THROW(!constrd, "Shadow_rim constrained does not have a child.");
			UG_COND_THROW(!constrg, "Shadow_rim constrainer does not have a child.");
		}
		*/

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
					constrainerMap[constrdInd[i][0]] = ConstrainerInfo(constrgInd[orientationOffsets[i]][0], constrg->is_constrained(), fct);
				}
			}
			else
			{
				for (size_t i = 0; i < n_ind; ++i)
				{
					UG_ASSERT(!constrdInd[i][1], "Block matrices are not supported (yet)!");
					UG_ASSERT(!constrgInd[i][1], "Block matrices are not supported (yet)!");
					constrainerMap[constrdInd[i][0]] = ConstrainerInfo(constrgInd[i][0], constrg->is_constrained(), fct);
				}
			}
		}
	}
}



template <typename TDomain, typename TAlgebra>
void Interface1D<TDomain, TAlgebra>::fill_constrainer_map(ConstSmartPtr<DoFDistribution> dd)
{
	// clear current constrainer map
	ConstraintInfo& cInfo = m_mConstraints[dd.get()];
	std::map<size_t, ConstrainerInfo>& constrainerMap = cInfo.constrainerMap;

	constrainerMap.clear();

	if (worldDim-1 >= VERTEX && dd->max_dofs(VERTEX) > 0)
		fill_constrainer_map<Vertex, Vertex, Edge>(dd);
	if (worldDim-1 >= EDGE && dd->max_dofs(EDGE) > 0)
		fill_constrainer_map<Edge, EdgeDescriptor, Face>(dd);
	if (worldDim-1 >= FACE && dd->max_dofs(FACE) > 0)
		fill_constrainer_map<Face, FaceDescriptor, Volume>(dd);

#if 0
// DEBUGGING
	typedef typename std::map<size_t, ConstrainerInfo>::iterator map_it;
	map_it iter = constrainerMap.begin();
	map_it iterEnd = constrainerMap.end();

	UG_LOGN("");
	UG_LOGN("DofDistro " << (dd->grid_level().is_surface() ? "surf " : "level ")
						<< (dd->grid_level().is_surface() ? 0 : dd->grid_level().level()))
	for (; iter != iterEnd; ++iter)
	{
		size_t constrd = iter->first;
		size_t constrg = (iter->second).constrgInd;
		size_t fct = (iter->second).fct;
		bool hanging = (iter->second).constrainerIsHanging;
		UG_LOGN(constrd << " -> " << constrg << "   (" << fct << ")"
				<< (!hanging ? "" : "   hanging constrainer"));
	}
	UG_LOGN("--------------------");
	UG_LOGN("Size: " << constrainerMap.size());
	UG_LOGN("");
#endif
}


template <typename TDomain, typename TAlgebra>
void Interface1D<TDomain, TAlgebra>::adjust_jacobian
(	matrix_type& J,
	const vector_type& u,
	ConstSmartPtr<DoFDistribution> dd,
	int type,
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

/*
// DEBUG (serial!): write constraint map to file
typename std::map<size_t, ConstraintInfo>::iterator it = m_constraintMap.begin();
typename std::map<size_t, ConstraintInfo>::iterator itEnd = m_constraintMap.end();

std::ofstream outFile;
outFile.open("constraintMap.dat", std::ios_base::out);
try
{
	for (; it != itEnd; ++it)
		outFile << it->first << ": " << it->second.constrgInd << " (" << it->second.fct << ")\n";
}
UG_CATCH_THROW("Output file 'debug_ilut.dat' could not be written to.");
outFile.close();
*/

	typedef typename matrix_type::row_iterator row_iterator;
	typename std::map<size_t, ConstrainerInfo>::iterator constrainerMapIt;

	// whether we are dealing with a surface dof distro
	bool isSurf = dd->grid_level().is_surface();

	// get constraint info for dof distro
	typedef typename std::map<const DoFDistribution*, ConstraintInfo>::iterator CIIter;
	CIIter ciit = m_mConstraints.find(dd.get());
	if (ciit == m_mConstraints.end()) update(dd);
	ConstraintInfo& cInfo = m_mConstraints[dd.get()];
	std::map<size_t, ConstrainerInfo>& constrainerMap = cInfo.constrainerMap;

	// loop rows
	size_t nr = J.num_rows();
	for (size_t i = 0; i < nr; i++)
	{
		// adapt current row if it depends on any constrained value

		// loop existing column entries of row
		std::vector<typename std::pair<size_t, ConstrainerInfo> > colIndices;
		for (row_iterator rit = J.begin_row(i); rit != J.end_row(i); ++rit)
		{
			// if a column index is part of the constrained set, then adjust row
			constrainerMapIt = constrainerMap.find(rit.index());

			if (constrainerMapIt != constrainerMap.end())
			{
				// store constrained col index until row_iterator reaches end
				// do that only of constrainer is not hanging and type is CT_CONSTRAINTS
				// OR if constrainer is hanging and type is CT_MAY_DEPEND_ON_HANGING
				if ((!constrainerMapIt->second.constrainerIsHanging && type == CT_CONSTRAINTS)
					|| (constrainerMapIt->second.constrainerIsHanging && type == CT_MAY_DEPEND_ON_HANGING))
					colIndices.push_back(std::make_pair(rit.index(), constrainerMapIt->second));
			}
		}

		// go through row entries again and adapt for constraints
		for (size_t col = 0; col < colIndices.size(); col++)
		{
			size_t constrdInd = colIndices[col].first;
			size_t constrgInd = colIndices[col].second.constrgInd;
			size_t fct = colIndices[col].second.fct;

			// interface nodes might not be present in level dof distro
			if (cInfo.algInd[0].size() < fct+1 || cInfo.algInd[1].size() < fct+1)
			{
				if (isSurf) UG_THROW("Interface node indices not present in surface dof distro.");

				// for level dof distros, only adapt constrained dof and constrainer
				// TODO: only works for linear constraints!
				// get defect derivatives for this vertex
				typename vector_type::value_type defDeriv[3];
				constraintValueDerivs(defDeriv, 1, 1, 1);

				// deriv wrt constraining
				if (J.has_connection(i, constrgInd))
					J(i, constrgInd) += J(i, constrdInd) * defDeriv[0];
				else
				{
					typename matrix_type::value_type help = J(i, constrdInd) * defDeriv[0];
					J(i, constrgInd) = help;
				}

				// deriv wrt constrained vertex
				J(i, constrdInd) = 0.0;
			}
			else
			{
				size_t IntfNodeHdInd = cInfo.algInd[0][fct];
				size_t IntfNode1dInd = cInfo.algInd[1][fct];

				// get defect derivatives for this vertex
				typename vector_type::value_type defDeriv[3];
				constraintValueDerivs(defDeriv, u[constrgInd], u[IntfNodeHdInd], u[IntfNode1dInd]);

				// NOTE: "If a side effect on a scalar object is unsequenced relative to either
				//       another side effect on the same scalar object or a value computation using
				//       the value of the same scalar object, the behavior is undefined."
				//       (http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2013/n3690.pdf)
				//
				//       In both of the following expressions, we might have a side effect in the
				//       evaluation of the lhs (entry might not be present in the CRS matrix and
				//       would need to be created then, possibly changing the value of the memory
				//       the rhs access to the same matrix is associated with)! The result depends
				//       on the order in which lhs and rhs evaluation are carried out, which is
				//       undefined.
				//
				//       J(i, constrgInd) += J(i, constrdInd) * defDeriv[0];
				//		 J(i, constrgInd) = J(i, constrdInd) * defDeriv[0];
				//
				//       We therefore need to take extra caution:

				// deriv wrt constraining
				if (J.has_connection(i, constrgInd))
					J(i, constrgInd) += J(i, constrdInd) * defDeriv[0];
				else
				{
					typename matrix_type::value_type help = J(i, constrdInd) * defDeriv[0];
					J(i, constrgInd) = help;
				}

				// deriv wrt interface vertex on constraining side
				if (J.has_connection(i, IntfNodeHdInd))
					J(i, IntfNodeHdInd) += J(i, constrdInd) * defDeriv[1];
				else
				{
					typename matrix_type::value_type help = J(i, constrdInd) * defDeriv[1];
					J(i, IntfNodeHdInd) = help;
				}

				// deriv wrt interface vertex on constrained side
				if (J.has_connection(i, IntfNode1dInd))
					J(i, IntfNode1dInd) += J(i, constrdInd) * defDeriv[2];
				else
				{
					typename matrix_type::value_type help = J(i, constrdInd) * defDeriv[2];
					J(i, IntfNode1dInd) = help;
				}

				// deriv wrt constrained vertex
				J(i, constrdInd) = 0.0;
			}
		}

		// Now, if the current row belongs to a constrained DoF itself,
		// then we need to add it to its 1d interface DoF and write an identity row instead.
		constrainerMapIt = constrainerMap.find(i);
		if (constrainerMapIt != constrainerMap.end()
			&& ((!constrainerMapIt->second.constrainerIsHanging && type == CT_CONSTRAINTS)
					|| (constrainerMapIt->second.constrainerIsHanging && type == CT_MAY_DEPEND_ON_HANGING)))
		{
			size_t fct = constrainerMapIt->second.fct;

			// on level dof distros, the intf node might not be present
			size_t IntfNode1dInd = (size_t) -1;
			if (cInfo.algInd[1].size() >= fct+1)
				IntfNode1dInd = cInfo.algInd[1][fct];

			// we put curly brackets around this so that the row_iterator will be removed afterwards
			// otherwise, we get a negative assert from the matrix implementation
			{
				const row_iterator iterEnd = J.end_row(i);
				for (row_iterator conn = J.begin_row(i); conn != iterEnd; ++conn)
				{
					// only add to intf node row if it is present
					if (IntfNode1dInd != (size_t) -1)
						J(IntfNode1dInd, conn.index()) += conn.value();

					// in any case, delete constrained row
					conn.value() = 0.0;
				}
			}

			// set unit diagonal
			J(i, i)	= 1.0;

#if 0 // NO! The constrained DoF is no real DoF; its value is determined by adjust_solution()!
			ConstrainerInfo ci = constrainerMapIt->second;
			size_t constrgInd = ci.constrgInd;
			size_t fct = ci.fct;
			size_t IntfNodeHdInd = cInfo.algInd[0][fct];
			size_t IntfNode1dInd = cInfo.algInd[1][fct];

			typename vector_type::value_type defDeriv[3];
			constraintValueDerivs(defDeriv, u[constrgInd], u[IntfNodeHdInd], u[IntfNode1dInd]);

			J(i, constrgInd) = -defDeriv[0];
			J(i, IntfNodeHdInd) = -defDeriv[1];
			J(i, IntfNode1dInd) = -defDeriv[2];
#endif
		}
	}
}



template <typename TDomain, typename TAlgebra>
void Interface1D<TDomain, TAlgebra>::adjust_defect
(	vector_type& d,
	const vector_type& u,
	ConstSmartPtr<DoFDistribution> dd,
	int type,
	number time,
	ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
	const std::vector<number>* vScaleMass,
	const std::vector<number>* vScaleStiff
)
{
	// get constraint info for dof distro
	typedef typename std::map<const DoFDistribution*, ConstraintInfo>::const_iterator CIIter;
	CIIter ciit = m_mConstraints.find(dd.get());
	if (ciit == m_mConstraints.end()) update(dd);
	ConstraintInfo& cInfo = m_mConstraints[dd.get()];
	std::map<size_t, ConstrainerInfo>& constrainerMap = cInfo.constrainerMap;

	// loop constrained vertices
	typename std::map<size_t, ConstrainerInfo>::iterator constrainerMapEnd = constrainerMap.end();
	typename std::map<size_t, ConstrainerInfo>::iterator constrainerMapIt = constrainerMap.begin();
	for (; constrainerMapIt != constrainerMapEnd; ++constrainerMapIt)
	{
		// only once for each constraint type
		if ((constrainerMapIt->second.constrainerIsHanging && type != CT_MAY_DEPEND_ON_HANGING)
			|| (!constrainerMapIt->second.constrainerIsHanging && type != CT_CONSTRAINTS))
			continue;

		// add defect to 1d interface node; then set 0
		size_t index = constrainerMapIt->first;
		size_t fct = constrainerMapIt->second.fct;

		UG_COND_THROW(cInfo.algInd[1].size() < fct+1, "Index for 1d interface node not present.");
		size_t IntfNode1dInd = cInfo.algInd[1][fct];

		d[IntfNode1dInd] += d[index];
		d[index] = 0.0;
	}
}



template <typename TDomain, typename TAlgebra>
void Interface1D<TDomain, TAlgebra>::adjust_linear
(	matrix_type& mat,
	vector_type& rhs,
	ConstSmartPtr<DoFDistribution> dd,
	int type,
	number time
)
{
	UG_THROW("This feature is not implemented.");
}



template <typename TDomain, typename TAlgebra>
void Interface1D<TDomain, TAlgebra>::adjust_rhs
(	vector_type& rhs,
	const vector_type& u,
	ConstSmartPtr<DoFDistribution> dd,
	int type,
	number time
)
{
	UG_THROW("This feature is not implemented.");
}



template <typename TDomain, typename TAlgebra>
void Interface1D<TDomain, TAlgebra>::adjust_solution
(	vector_type& u,
	ConstSmartPtr<DoFDistribution> dd,
	int type,
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

	bool isSurf = dd->grid_level().is_surface();

	// get constraint info for dof distro
	typedef typename std::map<const DoFDistribution*, ConstraintInfo>::const_iterator CIIter;
	CIIter ciit = m_mConstraints.find(dd.get());
	if (ciit == m_mConstraints.end()) update(dd);
	ConstraintInfo& cInfo = m_mConstraints[dd.get()];
	std::map<size_t, ConstrainerInfo>& constrainerMap = cInfo.constrainerMap;

	// loop constrained vertices
	typename std::map<size_t, ConstrainerInfo>::iterator constrainerMapEnd = constrainerMap.end();
	typename std::map<size_t, ConstrainerInfo>::iterator constrainerMapIt = constrainerMap.begin();
	for (; constrainerMapIt != constrainerMapEnd; ++constrainerMapIt)
	{
		// only once for each constraint type
		if ((constrainerMapIt->second.constrainerIsHanging && type != CT_MAY_DEPEND_ON_HANGING)
			|| (!constrainerMapIt->second.constrainerIsHanging && type != CT_CONSTRAINTS))
			continue;

		size_t constrdInd = constrainerMapIt->first;
		size_t constrgInd = constrainerMapIt->second.constrgInd;
		size_t fct = constrainerMapIt->second.fct;

		// interface nodes might not be present in level dof distro
		if (cInfo.algInd[0].size() < fct+1 || cInfo.algInd[1].size() < fct+1)
		{
			if (isSurf)
				{UG_THROW("Interface node indices not present in surface dof distro.");}
			else
				{UG_THROW("Called adjust_solution() for level dof distro; this is not implemented.");}
		}
		else
		{
			size_t intfSideInd = cInfo.algInd[0][fct];
			size_t intfCSideInd = cInfo.algInd[1][fct];

			constraintValue(u[constrdInd], u[constrgInd], u[intfSideInd], u[intfCSideInd]);
		}
	}
}



template <typename TDomain, typename TAlgebra>
void Interface1D<TDomain, TAlgebra>::adjust_correction
(	vector_type& u,
	ConstSmartPtr<DoFDistribution> dd,
	int type,
	number time
)
{
	// get constraint info for dof distro
	typedef typename std::map<const DoFDistribution*, ConstraintInfo>::const_iterator CIIter;
	CIIter ciit = m_mConstraints.find(dd.get());
	if (ciit == m_mConstraints.end()) update(dd);
	ConstraintInfo& cInfo = m_mConstraints[dd.get()];
	std::map<size_t, ConstrainerInfo>& constrainerMap = cInfo.constrainerMap;

	// loop constrained vertices
	typename std::map<size_t, ConstrainerInfo>::iterator constrainerMapEnd = constrainerMap.end();
	typename std::map<size_t, ConstrainerInfo>::iterator constrainerMapIt = constrainerMap.begin();
	for (; constrainerMapIt != constrainerMapEnd; ++constrainerMapIt)
	{
		// only once for each constraint type
		if ((constrainerMapIt->second.constrainerIsHanging && type != CT_MAY_DEPEND_ON_HANGING)
			|| (!constrainerMapIt->second.constrainerIsHanging && type != CT_CONSTRAINTS))
			continue;

		// set to 0
		size_t index = constrainerMapIt->first;
		u[index] = 0.0;
	}
}


/*
template <typename TDomain, typename TAlgebra>
void IInterface1D<TDomain, TAlgebra>::adjust_restriction
(
	vector_type& uCoarse,
	GridLevel coarseLvl,
	const vector_type& uFine,
	GridLevel fineLvl
)
{
	// get dof distro for coarse level
	ConstSmartPtr<DoFDistribution> dd = this->approximation_space()->dof_distribution(coarseLvl);

	// get constraint info for dof distro
	typedef typename std::map<const DoFDistribution*, ConstraintInfo>::const_iterator CIIter;
	CIIter ciit = m_mConstraints.find(dd.get());
	if (ciit == m_mConstraints.end()) update(dd);
	ConstraintInfo& cInfo = m_mConstraints[dd.get()];
	std::map<size_t, ConstrainerInfo>& constrainerMap = cInfo.constrainerMap;

	// loop constrained vertices
	typename std::map<size_t, ConstrainerInfo>::iterator constrainerMapEnd = constrainerMap.end();
	typename std::map<size_t, ConstrainerInfo>::iterator constrainerMapIt = constrainerMap.begin();
	for (; constrainerMapIt != constrainerMapEnd; ++constrainerMapIt)
	{
		// set to 0
		size_t index = constrainerMapIt->first;
		uCoarse[index] = 0.0;
	}
}
*/

template <typename TDomain, typename TAlgebra>
void Interface1D<TDomain, TAlgebra>::adjust_prolongation
(
	matrix_type& P,
	ConstSmartPtr<DoFDistribution> ddFine,
	ConstSmartPtr<DoFDistribution> ddCoarse,
	int type,
	number time
)
{

	// SCREW IT; does not work properly
/*
	// only adapt BEFORE hanging nodes adaption (as coarse grid node can never be hanging)
	if (type != CT_MAY_DEPEND_ON_HANGING)
		return;

	// adjust prolongation in constrainers using correct constraint values (not zero)
	// set prolongation for constrained nodes to 0.0

	typedef typename matrix_type::row_iterator row_iterator;
	typename std::map<size_t, ConstrainerInfo>::iterator constrainerMapIt;

	// get constraint info for dof distros
	typedef typename std::map<const DoFDistribution*, ConstraintInfo>::iterator CIIter;
	CIIter ciit = m_mConstraints.find(ddFine.get());
	if (ciit == m_mConstraints.end()) update(ddFine);
	ConstraintInfo& cInfo_f = m_mConstraints[ddFine.get()];
	std::map<size_t, ConstrainerInfo>& constrainerMap_f = cInfo_f.constrainerMap;

	typedef typename std::map<const DoFDistribution*, ConstraintInfo>::iterator CIIter;
	ciit = m_mConstraints.find(ddCoarse.get());
	if (ciit == m_mConstraints.end()) update(ddCoarse);
	ConstraintInfo& cInfo_c = m_mConstraints[ddCoarse.get()];
	std::map<size_t, ConstrainerInfo>& constrainerMap_c = cInfo_c.constrainerMap;

	// loop rows
	size_t nr = P.num_rows();
	for (size_t i = 0; i < nr; i++)
	{
		// If the current row belongs to a constrained DoF (fine level),
		// then we write a zero row (and do nothing more)
		constrainerMapIt = constrainerMap_f.find(i);
		if (constrainerMapIt != constrainerMap_f.end())
		{
			SetRow(P, i, 0.0);
			continue;
		}


		// adapt current row if it depends on any constrained value

		// loop existing column entries of row
		std::vector<typename std::pair<size_t, ConstrainerInfo> > colIndices;
		for (row_iterator rit = P.begin_row(i); rit != P.end_row(i); ++rit)
		{
			// if a column index is part of the constrained set, then adjust row
			constrainerMapIt = constrainerMap_c.find(rit.index());

			if (constrainerMapIt != constrainerMap_c.end())
			{
				// store constrained col index until row_iterator reaches end
				colIndices.push_back(std::make_pair(rit.index(), constrainerMapIt->second));
			}
		}

		// go through row entries again and adapt for constraints
		for (size_t col = 0; col < colIndices.size(); col++)
		{
			size_t constrdInd = colIndices[col].first;
			size_t constrgInd = colIndices[col].second.constrgInd;
			size_t fct = colIndices[col].second.fct;
			size_t IntfNodeHdInd = (size_t) -1;
			size_t IntfNode1dInd = (size_t) -1;

			// interface nodes might not be present in level dof distro
			if (cInfo_c.algInd[0].size() > fct && cInfo_c.algInd[1].size() > fct)
			{
				IntfNodeHdInd = cInfo_c.algInd[0][fct];
				IntfNode1dInd = cInfo_c.algInd[1][fct];
			}

			// get defect derivatives for this vertex
			typename vector_type::value_type defDeriv[3];
			// FIXME: This will only work if constraints are linear and add up to one.
			constraintValueDerivs(defDeriv, 1, 1, 1);

			// NOTE: "If a side effect on a scalar object is unsequenced relative to either
			//       another side effect on the same scalar object or a value computation using
			//       the value of the same scalar object, the behavior is undefined."
			//       (http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2013/n3690.pdf)
			//
			//       In both of the following expressions, we might have a side effect in the
			//       evaluation of the lhs (entry might not be present in the CRS matrix and
			//       would need to be created then, possibly changing the value of the memory
			//       the rhs access to the same matrix is associated with)! The result depends
			//       on the order in which lhs and rhs evaluation are carried out, which is
			//       undefined.
			//
			//       J(i, constrgInd) += J(i, constrdInd) * defDeriv[0];
			//		 J(i, constrgInd) = J(i, constrdInd) * defDeriv[0];
			//
			//       We therefore need to take extra caution:

			// wrt constraining
			if (P.has_connection(i, constrgInd))
				P(i, constrgInd) += P(i, constrdInd) * defDeriv[0];
			else
			{
				// The following must not be done!
				// We might create the first entry in column constrg,
				// which is the first entry in row constrg of the restriction matrix
				// (as restriction is prolongation transposed by default)
				// and this would mean that the assembled defect for constrg
				// (which might be on the surface in the coarse level)
				// would be overwritten by the restriction!
				// Either make sure that this is not the case or leave it be.
				//typename matrix_type::value_type help = P(i, constrdInd) * defDeriv[0];
				//P(i, constrgInd) = help;
			}

			// wrt interface vertex on constraining side
			if (IntfNodeHdInd != (size_t) -1)
			{
				if (P.has_connection(i, IntfNodeHdInd))
					P(i, IntfNodeHdInd) += P(i, constrdInd) * defDeriv[1];
				else
				{
					//typename matrix_type::value_type help = P(i, constrdInd) * defDeriv[1];
					//P(i, IntfNodeHdInd) = help;
				}
			}

			// wrt interface vertex on constrained side
			if (IntfNode1dInd != (size_t) -1)
			{
				if (P.has_connection(i, IntfNode1dInd))
					P(i, IntfNode1dInd) += P(i, constrdInd) * defDeriv[2];
				else
				{
					//typename matrix_type::value_type help = P(i, constrdInd) * defDeriv[2];
					//P(i, IntfNode1dInd) = help;
				}
			}

			// wrt constrained vertex
			P(i, constrdInd) = 0.0;
		}
	}
*/
}


template <typename TDomain, typename TAlgebra>
void Interface1D<TDomain, TAlgebra>::adjust_restriction
(
	matrix_type& R,
	ConstSmartPtr<DoFDistribution> ddCoarse,
	ConstSmartPtr<DoFDistribution> ddFine,
	int type,
	number time
)
{
	/*
	// adjust restriction like defect,
	// i.e. (only P1): add constrained/constrainer couplings to 1d intf node
	// then set zero row

	SmartPtr<MultiGrid> mg = this->m_spApproxSpace->domain()->grid();

	// get constraint info for coarse dof distro (for 1d intf node indices)
	typedef typename std::map<const DoFDistribution*, ConstraintInfo>::const_iterator CIIter;
	CIIter ciit = m_mConstraints.find(ddCoarse.get());
	if (ciit == m_mConstraints.end()) update(ddCoarse);
	ConstraintInfo& cInfo = m_mConstraints[ddCoarse.get()];


	// loop constrained elements (coarse)
	typename DoFDistribution::traits<Vertex>::const_iterator iter, iterEnd;
	iter = ddCoarse->begin<Vertex>(m_siConstr);
	iterEnd = ddCoarse->end<Vertex>(m_siConstr);

	for (; iter != iterEnd; ++iter)
	{
		// get constrained elem
		Vertex* constrd = *iter;

		// get parent of constrained
		Vertex* constrd_parent = dynamic_cast<Vertex*>(mg->get_parent(constrd));
		if (!constrd_parent) continue;

		// get constrainer of parent
		Vertex* constrg = get_constrainer<Vertex, Edge>(constrd_parent);

		// find indices
		// loop functions
		size_t numFct = m_vFct.size();
		for (size_t fct = 0; fct < numFct; fct++)
		{
			const size_t fct_ind = m_vFct[fct];

			// get inner DoF indices
			std::vector<DoFIndex> constrdInd;
			std::vector<DoFIndex> constrgInd;
			ddCoarse->inner_dof_indices(constrd, fct_ind, constrdInd, false);
			ddFine->inner_dof_indices(constrg, fct_ind, constrgInd, false);

			UG_COND_THROW(constrdInd.size() != 1,
				"Not exactly one index for constrained vertex (" << constrdInd.size() << " instead).");
			UG_COND_THROW(constrgInd.size() != 1,
				"Not exactly one index for constrainer vertex (" << constrgInd.size() << " instead).");

			// add restrictive coupling between coarse constrained and fine constrainer to 1d intf node
			R(cInfo.algInd[1][fct], constrgInd[0][0]) += R(constrdInd[0][0], constrgInd[0][0]);

			// set zero row
			SetRow(R, constrdInd[0], 0.0);
		}
	}
*/
}

} // namespace nernst_planck
} // namespace ug
