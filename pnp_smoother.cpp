/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Markus Breit
 * Creation date: 2017-06-20
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

#include "pnp_smoother.h"

#include <cmath>                                                   // for fabs
#include <algorithm>                                               // for sort
#include <string>                                                  // for operator+, basic_string

#include "common/math/math_vector_matrix/math_vector_functions.h"  // for VecSubtract, VecNormalize, VecCross
#include "common/math/ugmath_types.h"                              // for vector3, vector2
#include "common/util/string_util.h"                               // for TokenizeString
#include "lib_algebra/algebra_common/core_smoothers.h"             // for gs_step_LL
#ifdef UG_PARALLEL
	#include "lib_algebra/parallelization/parallel_index_layout.h"     // for IndexLayout
	#include "lib_algebra/parallelization/parallelization_util.h"      // for MatAddSlaveRowsToMasterRowOverlap0
#endif
#include "lib_algebra/small_algebra/blocks.h"                      // for BlockRef
#include "lib_disc/common/multi_index.h"                           // for DoFIndex
#include "lib_disc/dof_manager/dof_distribution.h"                 // for DoFDistribution, DoFDistribution::traits, DoFDistribution:...
#include "lib_disc/domain.h"                                       // for Domain1d, Domain2d, Domain3d
#include "lib_disc/function_spaces/approximation_space.h"          // for ApproximationSpace
#include "lib_grid/multi_grid.h"                                   // for MultiGrid
#include "lib_grid/algorithms/debug_util.h"                        // for ElementDebugInfo
#include "lib_grid/algorithms/element_side_util.h"                 // for GetOpposingSide
#include "lib_grid/grid/grid_base_objects.h"                       // for Vertex (ptr only), Face, Edge
#include "lib_grid/grid_objects/grid_dim_traits.h"                 // for grid_dim_traits
#include "lib_grid/tools/subset_group.h"                           // for SubsetGroup


namespace ug {
namespace nernst_planck {

template <typename TDomain>
number CalculateAbsNormalProd
(
	const typename TDomain::position_type& vec,
	const typename grid_dim_traits<TDomain::dim>::side_type* side,
	const typename TDomain::position_accessor_type& aaPos
)
{
	return 1.0;
}

template <>
number CalculateAbsNormalProd<Domain2d>
(
	const vector2& vec,
	const Edge* side,
	const Domain2d::position_accessor_type& aaPos
)
{
	vector2 a;
	VecSubtract(a, aaPos[side->vertex(1)], aaPos[side->vertex(0)]);
	std::swap(a[0], a[1]);
	a[1] *= -1.0;
	VecNormalize(a, a);

	return fabs(VecProd(a, vec));
}

template <>
number CalculateAbsNormalProd<Domain3d>
(
	const vector3& vec,
	const Face* side,
	const Domain3d::position_accessor_type& aaPos
)
{
	vector3 a, b;
	VecSubtract(a, aaPos[side->vertex(1)], aaPos[side->vertex(0)]);
	VecSubtract(b, aaPos[side->vertex(2)], aaPos[side->vertex(0)]);
	VecCross(a, a, b);
	VecNormalize(a, a);

	return fabs(VecProd(a, vec));
}



/// sorting utility
struct BlockOrder
{
	BlockOrder(const std::vector<size_t>& _vMinInd)
	: vMinInd(_vMinInd) {};

	bool operator() (size_t i, size_t j)
	{return vMinInd[i] < vMinInd[j];}

	private:
		const std::vector<size_t>& vMinInd;
};

/// sorting utility
struct DoFOrder
{
	bool operator() (const DoFIndex& i, const DoFIndex& j)
	{
		if (i[0] < j[0]) return true;
		if (i[0] > j[0]) return false;
		return i[1] < j[1];
	}
};





template <typename TDomain, typename TAlgebra, template <typename> class TPrecond>
PNPSmoother<TDomain, TAlgebra, TPrecond>::PNPSmoother(SmartPtr<ApproximationSpace<TDomain> > approx)
: m_name(std::string(precond_type::name()) + std::string("_PNP")),
  m_spApprox(approx),
  m_spMG(m_spApprox->domain()->grid()),
  m_spSH(m_spApprox->domain()->subset_handler()),
  m_spDD(m_spApprox->dof_distribution(this->m_gl, true)),
  m_method(0),
  m_ps(0),
  m_bBlockingNeedsReinit(true),
  m_spBM(new MatrixOperator<block_matrix_type, block_vector_type>())
{
	// set this class as listener for adaptation events
	m_spGridAdaptationCallbackID = m_spMG->message_hub()->register_class_callback(this,
		&PNPSmoother<TDomain, TAlgebra, TPrecond>::grid_adaptation_callback);

	// set this class as listener for distribution events
	m_spGridDistributionCallbackID = m_spMG->message_hub()->register_class_callback(this,
		&PNPSmoother<TDomain, TAlgebra, TPrecond>::grid_distribution_callback);
}


template <typename TDomain, typename TAlgebra, template <typename> class TPrecond>
PNPSmoother<TDomain, TAlgebra, TPrecond>::PNPSmoother(const PNPSmoother& parent)
: m_name(parent.m_name),
  m_spApprox(parent.m_spApprox),
  m_spMG(parent.m_spMG),
  m_spSH(parent.m_spSH),
  m_spDD(parent.m_spDD),
  m_vChargeSubsetsPairs(parent.m_vChargeSubsetsPairs),
  m_method(parent.m_method),
  m_ps(parent.m_ps),
  m_bBlockingNeedsReinit(true),
  m_spBM(new MatrixOperator<block_matrix_type, block_vector_type>())
{
	// set this class as listener for adaptation events
	m_spGridAdaptationCallbackID = m_spMG->message_hub()->register_class_callback(this,
		&PNPSmoother<TDomain, TAlgebra, TPrecond>::grid_adaptation_callback);

	// set this class as listener for distribution events
	m_spGridDistributionCallbackID = m_spMG->message_hub()->register_class_callback(this,
		&PNPSmoother<TDomain, TAlgebra, TPrecond>::grid_distribution_callback);
}


template <typename TDomain, typename TAlgebra, template <typename> class TPrecond>
PNPSmoother<TDomain, TAlgebra, TPrecond>::~PNPSmoother()
{}


template <typename TDomain, typename TAlgebra, template <typename> class TPrecond>
const char* PNPSmoother<TDomain, TAlgebra, TPrecond>::name() const
{
	return m_name.c_str();
}


template <typename TDomain, typename TAlgebra, template <typename> class TPrecond>
bool PNPSmoother<TDomain, TAlgebra, TPrecond>::preprocess(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp)
{
	// pointer to the matrix to work with
	matrix_type* pmat;
#ifdef UG_PARALLEL
	matrix_type mat;
	if (pcl::NumProcs() > 1)
	{
		if (m_ps == 0)  // unique/unique
		{
			mat = *pOp;

			// make consistent and set Dirichlet values on slaves
			MatAddSlaveRowsToMasterRowOverlap0(mat);
			std::vector<IndexLayout::Element> vIndex;
			CollectUniqueElements(vIndex, mat.layouts()->slave());
			SetDirichletRow(mat, vIndex);
		}
		else if (m_ps == 1)	// consistent / additive
		{
			mat = *pOp;
			MatMakeConsistentOverlap0(mat);
			//MakeConsistent(*pOp, mat);
		}
		else if (m_ps == 2)	// restricted additive Schwarz (RAS)
		{
			mat = *pOp;
			MatMakeConsistentOverlap0(mat);
			//MakeConsistent(*pOp, mat);
		}
		else UG_THROW("Invalid parallelization strategy " << m_ps << ".");
		pmat = &mat;
	}
	else
#endif

	pmat = &pOp->get_matrix();

	// if blocking needs reinit, do it
	if (m_bBlockingNeedsReinit)
		reinit_blocking();

	// adjust DD to current grid level
	m_spDD = m_spApprox->dof_distribution(this->m_gl, true);

	// write the blocked matrix
	size_t sz = m_vBlocks.size();
	m_spBM->resize_and_clear(sz, sz);

	for (size_t rowInd = 0; rowInd < sz; ++rowInd)
	{
		// loop new blocks
		const std::vector<size_t>& rowBlockInds = m_vBlocks[m_vPerm[rowInd]];
		size_t nRowInd = rowBlockInds.size();
		for (size_t i = 0; i < nRowInd; ++i)
		{
			size_t ri = rowBlockInds[i];

			// loop original matrix row for algebra index
			typename matrix_type::row_iterator rit = pmat->begin_row(ri);
			typename matrix_type::row_iterator ritEnd = pmat->end_row(ri);
			for (; rit != ritEnd; ++rit)
			{
				const size_t ci = rit.index();
				const size_t colInd = m_vInvPerm[m_vIndexToBlock[ci]];

				const std::vector<size_t>& colBlockInds = m_vBlocks[m_vPerm[colInd]];
				const size_t nColInd = colBlockInds.size();

				// copy the whole block //
				const typename matrix_type::value_type& block = rit.value();

				const size_t nr = GetRows(block);
				const size_t nc = GetCols(block);

				bool connCreated = false;
				if (!m_spBM->has_connection(rowInd, colInd))
					connCreated = true;
				typename block_matrix_type::value_type& newBlock = (*m_spBM)(rowInd, colInd);

				// we have to resize and zero-out newly created blocks
				// because it has the size and values of some other entry after initial construction!
				if (connCreated)
					newBlock.resize(nRowInd*nr, nColInd*nc, false);

				// find out which of the indices in new column is ci
				size_t j = (size_t) -1;
				for (size_t k = 0; k < nColInd; ++k)
				{
					if (colBlockInds[k] == ci)
					{
						j = k;
						break;
					}
				}
				UG_COND_THROW(j == (size_t) -1, "Column index not found in new column block indices.");

				// now do the actual copying
				for (size_t r = 0; r < nr; ++r)
					for (size_t c = 0; c < nc; ++c)
						newBlock(i*nr + r, j*nc + c) = BlockRef(block, r, c);
			}
		}
	}

	/*
	// DEBUG: check sizes
	typename matrix_type::value_type dummy;

	size_t nr = m_spBM->num_rows();
	for (size_t r = 0; r < nr; ++r)
	{
		typename block_matrix_type::row_iterator rit = m_spBM->begin_row(r);
		typename block_matrix_type::row_iterator ritEnd = m_spBM->end_row(r);
		for (; rit != ritEnd; ++rit)
		{
			const size_t c = rit.index();
			typename block_matrix_type::value_type& block = rit.value();

			const size_t nbr = block.num_rows();
			const size_t nbc = block.num_cols();

			const size_t nbrExpected = m_vBlocks[m_vPerm[r]].size() * GetRows(dummy);
			const size_t nbcExpected = m_vBlocks[m_vPerm[c]].size() * GetCols(dummy);

			if (nbr != nbrExpected)
				UG_LOGN("Number of rows in block (" << r << "," << c << ") incorrect: "
						<< "should be " << nbrExpected << ", but is " << nbr << ".");

			if (nbc != nbcExpected)
				UG_LOGN("Number of cols in block (" << r << "," << c << ") incorrect: "
						<< "should be " << nbcExpected << ", but is " << nbc << ".");
		}
	}
	*/

	bool success;
	do_preprocess<TPrecond>(success, this, m_spBM);
	return success;
}


template <typename TDomain, typename TAlgebra, template <typename> class TPrecond>
bool PNPSmoother<TDomain, TAlgebra, TPrecond>::step
(
	SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp,
	vector_type& c,
	const vector_type& d
)
{
#ifdef UG_PARALLEL
	// make defect unique
	SmartPtr<vector_type> dTmp = d.clone();

	if (m_ps == 0)
		dTmp->change_storage_type(PST_UNIQUE);
	else if (m_ps == 2)	// restricted additive Schwarz (RAS)
		dTmp->change_storage_type(PST_CONSISTENT);

#else
	const vector_type* dTmp = &d;
#endif


	// copy defect to suitable algebra vector
	// and resize correction along the way
	size_t sz = m_vBlocks.size();
	block_vector_type cb(sz);
	block_vector_type db(sz);
	for (size_t i = 0; i < sz; ++i)
	{
		size_t bsz = m_vBlocks[m_vPerm[i]].size();
		for (size_t j = 0; j < bsz; ++j)
		{
			const typename vector_type::value_type& block = (*dTmp)[m_vBlocks[m_vPerm[i]][j]];

			size_t origBlockSz = GetSize(block);
			if (db[i].size() < bsz * origBlockSz)
			{
				cb[i].resize(bsz * origBlockSz);	// TODO: this is going to be slow due to many small allocations
				db[i].resize(bsz * origBlockSz);
			}

			for (size_t k = 0; k < origBlockSz; ++k)
				db[i][j*origBlockSz + k] = BlockRef(block, k);
		}
	}

	bool success;
	do_step<TPrecond>(success, this, m_spBM, cb, db);

	// copy back correction
	for (size_t i = 0; i < sz; ++i)
	{
		size_t bsz = m_vBlocks[m_vPerm[i]].size();
		for (size_t j = 0; j < bsz; ++j)
		{
			typename vector_type::value_type& block = c[m_vBlocks[m_vPerm[i]][j]];

			size_t origBlockSz = GetSize(block);
			for (size_t k = 0; k < origBlockSz; ++k)
				BlockRef(block, k) = cb[i][j*origBlockSz + k];
		}
	}

#ifdef UG_PARALLEL
	if (m_ps == 0)
		c.set_storage_type(PST_UNIQUE);
	else if (m_ps == 1)
		c.set_storage_type(PST_ADDITIVE);
	else if (m_ps == 2)	// restricted additive Schwarz (RAS)
		c.set_storage_type(PST_UNIQUE); // solution is not really unique, but we ignore non-master corrections
	else UG_THROW("Invalid parallelization strategy " << m_ps << ".");

	c.change_storage_type(PST_CONSISTENT);
#endif

	return success;
}


template <typename TDomain, typename TAlgebra, template <typename> class TPrecond>
bool PNPSmoother<TDomain, TAlgebra, TPrecond>::postprocess()
{
	bool success;
	do_postprocess<TPrecond>(success, this);
	return success;
}


template <typename TDomain, typename TAlgebra, template <typename> class TPrecond>
bool PNPSmoother<TDomain, TAlgebra, TPrecond>::supports_parallel() const
{
	return precond_type::supports_parallel();
}


template <typename TDomain, typename TAlgebra, template <typename> class TPrecond>
void PNPSmoother<TDomain, TAlgebra, TPrecond>::grid_level_has_changed()
{
//std::cout << "Reblocking needed due to changed grid level (" << this->m_gl << ")." << std::endl;
	m_bBlockingNeedsReinit = true;
}



template <typename TDomain, typename TAlgebra, template <typename> class TPrecond>
void PNPSmoother<TDomain, TAlgebra, TPrecond>::add_charge_surface_pair(const std::string& chSsName, const std::string& volSsName)
{
	m_vChargeSubsetsPairs.resize(m_vChargeSubsetsPairs.size() + 1);

	SubsetGroup ssg(m_spSH);
	int chSi, volSi;
	try
	{
		ssg.add(chSsName);
		chSi = ssg[0];
	}
	UG_CATCH_THROW("Something wrong with passed charge subsets. Pair cannot be added.");
	try
	{
		ssg.remove(ssg[0]);
		ssg.add(TokenizeString(volSsName));
		volSi = ssg[0];
	}
	UG_CATCH_THROW("Something wrong with passed volume subsets. Pair cannot be added.");

	m_vChargeSubsetsPairs.back().first = chSi;
	std::vector<int>& vVolSI = m_vChargeSubsetsPairs.back().second;
	vVolSI.resize(ssg.size());
	for (size_t i = 0; i < ssg.size(); ++i)
		vVolSI[i] = ssg[i];
}


template <typename TDomain, typename TAlgebra, template <typename> class TPrecond>
void PNPSmoother<TDomain, TAlgebra, TPrecond>::set_method(int m)
{
	m_method = m;
}


template <typename TDomain, typename TAlgebra, template <typename> class TPrecond>
void PNPSmoother<TDomain, TAlgebra, TPrecond>::reinit_blocking()
{
	// TODO: This is only adequate for P1 approx spaces.
	// find pairs of charged surface / volume vertices;
	// if a volume vertex already belongs to another charged surface vertex,
	// add all together (triplet);
	// throw if a charged surface vertex has more than one volume vertex neighbor
	// (as we only want to do this for hex/prism layers)
	typedef typename DoFDistribution::traits<Vertex>::const_iterator vrt_it_type;
	typedef typename MultiGrid::traits<Edge>::secure_container edge_list;
	typedef typename grid_dim_traits<TDomain::dim>::side_type side_type;
	typedef typename MultiGrid::traits<side_type>::secure_container side_list;

	// preparation: create an attachment that holds the block index for each assigned vertex
	Attachment<size_t> aBlockInd("blockIndex", false);
	if (!m_spMG->has_attachment<Vertex>(aBlockInd))
		m_spMG->attach_to_dv<Vertex>(aBlockInd, (size_t) -1);
	Grid::AttachmentAccessor<Vertex, Attachment<size_t> > aaBlockInd(*m_spMG, aBlockInd);
	const typename TDomain::position_accessor_type& aaPos = m_spApprox->domain()->position_accessor();

	// adjust DD to current grid level
	m_spDD = m_spApprox->dof_distribution(this->m_gl, true);

	// iterate DD vertices (i.e., level vertices for level DDs and surface vertices for surface DDs)
	m_vBlocks.clear();
	m_vIndexToBlock.clear();
	m_vIndexToBlock.resize(m_spDD->num_indices());
	size_t blockInd = 0;
	vrt_it_type vit = m_spDD->begin<Vertex>();
	vrt_it_type vit_end = m_spDD->end<Vertex>();
	for (; vit != vit_end; ++vit)
	{
		Vertex* vrt = *vit;

		// check if already assigned
		if (aaBlockInd[vrt] != (size_t) -1)
			continue;

		// check if belongs to charged surface
		int si = m_spSH->get_subset_index(vrt);
		bool chargeSubset = false;
		std::vector<int>* volSs = NULL;
		size_t nCh = m_vChargeSubsetsPairs.size();
		for (size_t i = 0; i < nCh; ++i)
		{
			if (si == m_vChargeSubsetsPairs[i].first)
			{
				chargeSubset = true;
				volSs = &m_vChargeSubsetsPairs[i].second;
				break;
			}
		}

		// if on charged subset: find corresponding volume vertex:
		// we do not throw if there is none;
		// we do not throw if there is more than one, in that case, we take the one most orthogonal;
		// if that cannot be determined, we just pick the first volume vertex
		Vertex* volVrt = NULL;
		if (chargeSubset)
		{
			// loop associated edges; retain the opposing vertices in the volume subset(s)
			std::vector<Vertex*> vVolVrt;
			edge_list el;
			m_spMG->associated_elements(el, vrt);
			size_t elSz = el.size();
			for (size_t e = 0; e < elSz; ++e)
			{
				Edge* edge = el[e];
				int esi = m_spSH->get_subset_index(edge);

				for (size_t i = 0; i < volSs->size(); ++i)
				{
					if (esi == (*volSs)[i])
					{
						// it is possible the volume vertex is a ghost;
						// then it is not contained in the (ghost-free) level dof distro
						// and we cannot use it
						// (However, this typically means the distribution constraint is violated
						// and the domain is in fact decomposed along the charged subset!)
						Vertex* opp = GetOpposingSide(*m_spMG, edge, vrt);
						if (m_spDD->is_contained(opp))
							vVolVrt.push_back(opp);
					}
				}
			}

			// return if there is exactly one volume neighbor
			size_t nvv = vVolVrt.size();
			if (nvv == 1)
				volVrt = vVolVrt[0];

			// more than one volume neighbor
			else if (nvv >= 2)
			{
				// calculate product of surface normal in vertex with vec to volume vrt
				std::vector<side_type*> vSrfSide;
				side_list sl;
				m_spMG->associated_elements(sl, vrt);

				size_t slSz = sl.size();
				for (size_t s = 0; s < slSz; ++s)
				{
					side_type* side = sl[s];
					int ssi = m_spSH->get_subset_index(side);
					if (ssi == si)
						vSrfSide.push_back(side);
				}

				size_t nss = vSrfSide.size();
				// in the unlikely case where the domain decomposition eliminates all
				// surface neighbors just pick the first volume vertex (should not be too bad!?)
				if (!nss)
					volVrt = vVolVrt[0];

				else
				{
					const typename TDomain::position_type& vrtCoords = aaPos[vrt];

					number maxProd = -1.0;
					size_t vvimax = 0;
					for (size_t v = 0; v < nvv; ++v)
					{
						typename TDomain::position_type vec;
						VecSubtract(vec, aaPos[vVolVrt[v]], vrtCoords);

						number normalProdSum = 0.0;
						for (size_t s = 0; s < nss; ++s)
							normalProdSum += CalculateAbsNormalProd<TDomain>(vec, vSrfSide[s], aaPos);
						normalProdSum /= nss;

						if (normalProdSum > maxProd)
						{
							maxProd = normalProdSum;
							vvimax = v;
						}
					}

					// take the volume vertex with maximal sum
					volVrt = vVolVrt[vvimax];
				}
			}
		}

		if (volVrt)
		{
			// find the algebra indices associated with the surface vertex
			std::vector<size_t> vAlgInd;
			size_t nInd = m_spDD->inner_algebra_indices(vrt, vAlgInd, false);

			// if volume vertex has been assigned before, add surface DoFs to this block
			size_t volBlockInd = aaBlockInd[volVrt];
			if (volBlockInd != (size_t) -1)
			{
				std::vector<size_t>& block = m_vBlocks[volBlockInd];
				size_t blockSz = block.size();
				block.resize(blockSz + nInd);
				for (size_t i = 0; i < nInd; ++i)
				{
					block[blockSz + i] = vAlgInd[i];
					m_vIndexToBlock[vAlgInd[i]] = volBlockInd;
				}

				// retain block index for surface vertex
				aaBlockInd[vrt] = volBlockInd;
			}

			// otherwise also find the volume vertex dofs and combine with surface dofs in new block
			else
			{
				// find the dof indices associated with the volume vertex
				std::vector<size_t> vVolAlgInd;
				size_t nVolInd = m_spDD->inner_algebra_indices(volVrt, vVolAlgInd, false);

				// create new block and add all dofs
				m_vBlocks.resize(blockInd + 1);

				std::vector<size_t>& block = m_vBlocks[blockInd];
				block.resize(nInd + nVolInd);
				for (size_t i = 0; i < nInd; ++i)
				{
					block[i] = vAlgInd[i];
					m_vIndexToBlock[vAlgInd[i]] = blockInd;
				}
				for (size_t i = 0; i < nVolInd; ++i)
				{
					block[nInd + i] = vVolAlgInd[i];
					m_vIndexToBlock[vVolAlgInd[i]] = blockInd;
				}

				// retain block index for both vertices and increment blockInd
				aaBlockInd[vrt] = aaBlockInd[volVrt] = blockInd;
				++blockInd;
			}
		}

		// if not on charge subset, simply put DofIndices to a new block
		else
		{
			// find the dof indices associated with the surface vertex
			std::vector<size_t> vAlgInd;
			size_t nInd = m_spDD->inner_algebra_indices(vrt, vAlgInd, false);

			// create new block and add all dofs
			m_vBlocks.resize(blockInd + 1);

			std::vector<size_t>& block = m_vBlocks[blockInd];
			block.resize(nInd);
			for (size_t i = 0; i < nInd; ++i)
			{
				block[i] = vAlgInd[i];
				m_vIndexToBlock[vAlgInd[i]] = blockInd;
			}

			// retain block index for vertex and increment blockInd
			aaBlockInd[vrt] = blockInd;
			++blockInd;
		}
	}

	// sort each block according to previous order in dof distro
	std::vector<size_t> vMinInd(blockInd);
	for (size_t i = 0; i < blockInd; ++i)
	{
		std::sort(m_vBlocks[i].begin(), m_vBlocks[i].end());
		vMinInd[i] = m_vBlocks[i][0];
	}

	// sort vector of blocks according to previous order in dof distro
	m_vPerm.clear();
	m_vPerm.resize(blockInd);
	for (size_t i = 0; i < blockInd; ++i)
		m_vPerm[i] = i;

	BlockOrder blockOrder(vMinInd);
	std::sort(m_vPerm.begin(), m_vPerm.end(), blockOrder);

	m_vInvPerm.clear();
	m_vInvPerm.resize(blockInd);
	for (size_t i = 0; i < blockInd; ++i)
		m_vInvPerm[m_vPerm[i]] = i;

	/*
	// DEBUG: print block vector
	std::cout << std::endl;
	std::cout << std::endl;
	for (std::size_t i = 0; i < m_vBlocks.size(); ++i)
	{
		for (std::size_t j = 0; j < m_vBlocks[i].size(); ++j)
			std::cout << m_vBlocks[i][j] << " ";
		std::cout << std::endl;
	}

	// DEBUG: check size of blocked matrix
	size_t n = 0;
	for (std::size_t i = 0; i < m_vBlocks.size(); ++i)
		n += m_vBlocks.at(i).size();
	if (n != m_vIndexToBlock.size())
		UG_LOGN("Size mismatch: Number of indices in blockInds is " << n
			<< ", but #indices in DD is " << m_vIndexToBlock.size() << ".");

	// DEBUG: check index to block vector
	for (std::size_t i = 0; i < m_vIndexToBlock.size(); ++i)
	{
		const std::vector<size_t>& blockInds = m_vBlocks[m_vIndexToBlock[i]];
		if (std::find(blockInds.begin(), blockInds.end(), i) == blockInds.end())
			UG_LOGN("Block mapped to index " << i << " does not contain " << i << ".");
	}

	// DEBUG: block size histogram
	std::map<size_t, size_t> blockSizeMap;
	for (std::size_t i = 0; i < m_vBlocks.size(); ++i)
	{
		std::map<size_t, size_t>::iterator it;
		if ((it = blockSizeMap.find(m_vBlocks[i].size())) == blockSizeMap.end())
			blockSizeMap[m_vBlocks[i].size()] = 1;
		else
			++(it->second);
	}

	UG_LOGN("Block size map level " << this->m_gl << ":");
	std::map<size_t, size_t>::const_iterator it = blockSizeMap.begin();
	std::map<size_t, size_t>::const_iterator itEnd = blockSizeMap.end();
	for (; it != itEnd; ++it)
		UG_LOGN("  block size " << it->first << ": " << it->second);
	*/
}


template <typename TDomain, typename TAlgebra, template <typename> class TPrecond>
void PNPSmoother<TDomain, TAlgebra, TPrecond>::grid_adaptation_callback(const GridMessage_Adaption& gma)
{
	// update blocks if grid has been adapted
	if (gma.step_ends())
	{
//std::cout << "Reblocking needed due to grid adaptation." << std::endl;
		m_bBlockingNeedsReinit = true;
	}
}


template <typename TDomain, typename TAlgebra, template <typename> class TPrecond>
void PNPSmoother<TDomain, TAlgebra, TPrecond>::grid_distribution_callback(const GridMessage_Distribution& gmd)
{
	// update blocks if grid has been redistributed
	if (gmd.msg() == GMDT_DISTRIBUTION_STOPS)
	{
//std::cout << "Reblocking needed due to grid distribution." << std::endl;
		m_bBlockingNeedsReinit = true;
	}
}


template <typename TDomain, typename TAlgebra, template <typename> class TPrecond>
template <typename dummy>
PNPSmoother<TDomain, TAlgebra, TPrecond>::do_preprocess<ILU, dummy>::do_preprocess
(
	bool& successOut,
	this_type* pnpSmoother,
	SmartPtr<MatrixOperator<block_matrix_type, block_vector_type> > pOp
)
{
	// do not do a thing if preprocessing disabled
	if (pnpSmoother->m_bDisablePreprocessing)
	{
		successOut = true;
		return;
	}

	block_matrix_type& mat = *pOp;
	pnpSmoother->m_ILU = mat;

	// sort if needed
	if (pnpSmoother->m_bSort)
		pnpSmoother->calc_cuthill_mckee();

	// resize help vector
	pnpSmoother->m_h.resize(mat.num_cols());

	// compute ILU factorization
	if (pnpSmoother->m_beta != 0.0) FactorizeILUBeta(pnpSmoother->m_ILU, pnpSmoother->m_beta);
	else if (block_matrix_type::rows_sorted) FactorizeILUSorted(pnpSmoother->m_ILU, pnpSmoother->m_sortEps);
	else FactorizeILU(pnpSmoother->m_ILU);
	pnpSmoother->m_ILU.defragment();

	successOut = true;
}


template <typename TDomain, typename TAlgebra, template <typename> class TPrecond>
template <typename dummy>
PNPSmoother<TDomain, TAlgebra, TPrecond>::do_step<ILU, dummy>::do_step
(
	bool& successOut,
	this_type* pnpSmoother,
	SmartPtr<MatrixOperator<block_matrix_type, block_vector_type> > pOp,
	block_vector_type& c,
	const block_vector_type& d
)
{
	pnpSmoother->applyLU(c, d, pnpSmoother->m_h);
	successOut = true;
}


template <typename TDomain, typename TAlgebra, template <typename> class TPrecond>
template <typename dummy>
PNPSmoother<TDomain, TAlgebra, TPrecond>::do_step<GaussSeidel, dummy>::do_step
(
	bool& successOut,
	this_type* pnpSmoother,
	SmartPtr<MatrixOperator<block_matrix_type, block_vector_type> > pOp,
	block_vector_type& c,
	const block_vector_type& d
)
{
	gs_step_LL(*pOp, c, d, pnpSmoother->m_relax);
	successOut = true;
}


template <typename TDomain, typename TAlgebra, template <typename> class TPrecond>
void PNPSmoother<TDomain, TAlgebra, TPrecond>::my_clone(SmartPtr<ILinearIterator<vector_type> >& res)
{
	SmartPtr<ILinearIterator<vector_type> > sp(new PNPSmoother(*this));
	res = sp;
}

template <typename TDomain, typename TAlgebra, template <typename> class TPrecond>
void PNPSmoother<TDomain, TAlgebra, TPrecond>::my_clone(SmartPtr<ILinearIterator<block_vector_type> >& res)
{
	SmartPtr<ILinearIterator<block_vector_type> > sp(new PNPSmoother(*this));
	res = sp;
}





// explicit template specialization declaration
#ifdef UG_DIM_1
	#ifdef UG_CPU_1
		template class PNPSmoother<Domain1d, CPUAlgebra, ILU>;
		template class PNPSmoother<Domain1d, CPUAlgebra, GaussSeidel>;
		//template class PNP_ILU<Domain1d, CPUAlgebra>;
	#endif
	#ifdef UG_CPU_5
		template class PNPSmoother<Domain1d, CPUBlockAlgebra<5>, ILU>;
		template class PNPSmoother<Domain1d, CPUBlockAlgebra<5>, GaussSeidel>;
		//template class PNP_ILU<Domain1d, CPUBlockAlgebra<5> >;
	#endif
	#ifdef UG_CPU_6
		template class PNPSmoother<Domain1d, CPUBlockAlgebra<6>, ILU>;
		template class PNPSmoother<Domain1d, CPUBlockAlgebra<6>, GaussSeidel>;
		//template class PNP_ILU<Domain1d, CPUBlockAlgebra<6> >;
	#endif
#endif
#ifdef UG_DIM_2
	#ifdef UG_CPU_1
		template class PNPSmoother<Domain2d, CPUAlgebra, ILU>;
		template class PNPSmoother<Domain2d, CPUAlgebra, GaussSeidel>;
		//template class PNP_ILU<Domain2d, CPUAlgebra>;
	#endif
	#ifdef UG_CPU_5
		template class PNPSmoother<Domain2d, CPUBlockAlgebra<5>, ILU>;
		template class PNPSmoother<Domain2d, CPUBlockAlgebra<5>, GaussSeidel>;
		//template class PNP_ILU<Domain2d, CPUBlockAlgebra<5> >;
	#endif
	#ifdef UG_CPU_6
		template class PNPSmoother<Domain2d, CPUBlockAlgebra<6>, ILU>;
		template class PNPSmoother<Domain2d, CPUBlockAlgebra<6>, GaussSeidel>;
		//template class PNP_ILU<Domain2d, CPUBlockAlgebra<6> >;
	#endif
#endif
#ifdef UG_DIM_3
	#ifdef UG_CPU_1
		template class PNPSmoother<Domain3d, CPUAlgebra, ILU>;
		template class PNPSmoother<Domain3d, CPUAlgebra, GaussSeidel>;
		//template class PNP_ILU<Domain3d, CPUAlgebra>;
	#endif
	#ifdef UG_CPU_5
		template class PNPSmoother<Domain3d, CPUBlockAlgebra<5>, ILU>;
		template class PNPSmoother<Domain3d, CPUBlockAlgebra<5>, GaussSeidel>;
		//template class PNP_ILU<Domain3d, CPUBlockAlgebra<5> >;
	#endif
	#ifdef UG_CPU_6
		template class PNPSmoother<Domain3d, CPUBlockAlgebra<6>, ILU>;
		template class PNPSmoother<Domain3d, CPUBlockAlgebra<6>, GaussSeidel>;
		//template class PNP_ILU<Domain3d, CPUBlockAlgebra<6> >;
	#endif
#endif


} // namespace nernst_planck
} // namespace ug
