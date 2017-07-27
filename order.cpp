/*
 * order.cpp
 *
 *  Created on: 12.08.2015
 *      Author: mbreit
 */

#include "order.h"

#include <cstdlib>                                             // for abs
#include <algorithm>                                           // for min, stable_sort
#include <map>                                                 // for map
#include <set>                                                 // for set

#include "common/assert.h"                                     // for UG_ASSERT
#include "common/error.h"                                      // for UG_COND_THROW, UG_CATCH_THROW, UG_THROW
#include "common/types.h"                                      // for number
#include "common/util/string_util.h"                           // for TokenizeString
#include "lib_disc/dof_manager/dof_distribution.h"             // for DoFDistribution, DoFDistribution::tr...
#include "lib_disc/domain.h"                                   // for Domain1d, Domain2d, Domain3d
#include "lib_disc/function_spaces/approximation_space.h"      // for ApproximationSpace
#include "lib_grid/grid/grid_base_objects.h"                   // for GridBaseObjectId::EDGE, GridBaseObje...


namespace ug {
namespace nernst_planck {


template <typename TBaseElem>
void collect_obj_indices
(
	std::vector<int>& vObjInd,
	SmartPtr<DoFDistribution> dd,
	ConstSmartPtr<MGSubsetHandler> sh,
	const SubsetGroup& csg
)
{
	typename DoFDistribution::traits<TBaseElem>::const_iterator it, itEnd;
	it = dd->begin<TBaseElem>();
	itEnd = dd->end<TBaseElem>();


	std::vector<size_t> vInd;
	for (; it != itEnd; ++it)
	{
		size_t nInd = dd->inner_algebra_indices(*it, vInd, true);

		if (!nInd) continue;

		// obj index is minimum of all indices
		size_t objInd = vInd[0];
		for (size_t i = 1; i < nInd; ++i)
			objInd = std::min(objInd, vInd[i]);

		vObjInd[objInd] = (int) nInd;

		// if subset is constrained: negative value
		int si = sh->get_subset_index(*it);
		for (size_t i = 0; i < csg.size(); ++i)
		{
			if (si == csg[i])
			{
				vObjInd[objInd] *= -1;
				break;
			}
		}
	}
}




template <typename TDomain>
void reorder_dofs(SmartPtr<ApproximationSpace<TDomain> > approxSpace, const char* constrained)
{
	UG_COND_THROW(!approxSpace.valid(), "Approximation space is invalid.");
	std::vector<SmartPtr<DoFDistribution> > vDD = approxSpace->dof_distributions();
	ConstSmartPtr<MGSubsetHandler> sh = approxSpace->subset_handler();

	SubsetGroup sg(sh);
	try	{sg.add(TokenizeString(constrained, ','));}
	UG_CATCH_THROW("Subset group could not be created.");

	// loop dof distributions
	for (size_t i = 0; i < vDD.size(); ++i)
	{
		SmartPtr<DoFDistribution> dd = vDD[i];

		size_t nInd = dd->num_indices();
		std::vector<int> vObjInd(nInd, 0);

		// collect object indices for constrained subsets:
		// the resulting vector will have values != 0 exactly at obj_indices
		// and the absolute value is the number of indices corresponding to the object
		// value is <0 exactly for constrained obj_indices
		if (TDomain::dim >= VERTEX && dd->max_dofs(VERTEX) > 0)
			collect_obj_indices<Vertex>(vObjInd, dd, sh, sg);
		if (TDomain::dim >= EDGE && dd->max_dofs(EDGE) > 0)
			collect_obj_indices<Edge>(vObjInd, dd, sh, sg);
		if (TDomain::dim >= FACE && dd->max_dofs(FACE) > 0)
			collect_obj_indices<Face>(vObjInd, dd, sh, sg);
		if (TDomain::dim >= VOLUME && dd->max_dofs(VOLUME) > 0)
			collect_obj_indices<Volume>(vObjInd, dd, sh, sg);

		// now re-order
		// check integrity
		size_t curr = 0;
		size_t nConstr = 0;
		while (curr < nInd)
		{
			UG_COND_THROW(!vObjInd[curr], "curr index " << curr << " is non-obj, but needs to be obj!");
			UG_COND_THROW(curr+abs(vObjInd[curr]) > nInd, "index vector ends prematurely!");
			if (vObjInd[curr] < 0) nConstr += -vObjInd[curr];
			for (size_t j = 1; j < (size_t) abs(vObjInd[curr]); ++j)
			{
				UG_COND_THROW(vObjInd[curr+j], "curr index " << curr+j << " is obj, but needs to be non-obj!");
			}
			curr += abs(vObjInd[curr]);
		}

		std::vector<size_t> newInd(nInd);
		size_t nUnConstr = 0;
		nConstr = nInd - nConstr;
		curr = 0;
		while (curr < nInd)
		{
			if (vObjInd[curr] > 0)
			{
				for (size_t j = 0; j < (size_t) vObjInd[curr]; ++j)
					newInd[curr+j] = nUnConstr+j;

				nUnConstr += vObjInd[curr];
				curr += vObjInd[curr];
			}
			else
			{
				for (size_t j = 0; j < (size_t) -vObjInd[curr]; ++j)
					newInd[curr+j] = nConstr+j;

				nConstr += -vObjInd[curr];
				curr += -vObjInd[curr];
			}
		}

		// apply new order
		dd->permute_indices(newInd);
	}
}




/*
struct Label
{
	public:
		Label() {};
		bool operator<(const Label& l)
		{
			size_t min = std::min(m_vL.size(), l.m_vL.size());
			for (size_t i = 0; i < min; ++i)
			{
				if (m_vL[i] < l.m_vL[i]) return true;
				if (m_vL[i] > l.m_vL[i]) return false;
			}
			return l.m_vL.size() > min;
		}

	private:
		std::vector<size_t> m_vL;
};
*/



namespace {
struct LabelCmp
{
	LabelCmp(const std::vector<float>& _vL) : vL(_vL) {};
	bool operator()(const size_t& a, const size_t& b)
	{ return vL[a] < vL[b];}

	const std::vector<float>& vL;
};
}


void lex_order
(
	std::vector<size_t>& newIndOut,
	const std::vector<std::vector<size_t> >& vAdj,
	bool preserveConsec
)
{
	size_t nVrt = vAdj.size();

	std::vector<float> vLabel(nVrt, 1.0);
	std::vector<bool> vReached(nVrt, false);
	std::vector<bool> vHandled(nVrt, false);

	// mark all indices without connections as handled
	for (size_t i = 0; i < nVrt; ++i)
		if (!vAdj[i].size())
			vHandled[i] = true;

	// get first unhandled index
	size_t k = 1;
	size_t i_unhandled = nVrt - 1;
	for (; i_unhandled < nVrt; --i_unhandled)
		if (!vHandled[i_unhandled]) break;

	size_t highestLabelVrt = i_unhandled;


	std::vector<size_t> vNewOrder;
	while (true)
	{
		size_t v = highestLabelVrt;

		// assign i to current vertex
		vNewOrder.push_back(v);
		vHandled[v] = true;

		// mark v reached
		vReached[v] = true;
		std::map<number, std::set<size_t> > mReach;

		// mark all unnumbered verts unreached
		for (size_t j = 0; j < nVrt; ++j)
			if (!vHandled[j])
				vReached[j] = false;

		for (size_t j = 0; j < vAdj[v].size(); ++j)
		{
			size_t w = vAdj[v][j];
			if (vHandled[w]) continue;

			mReach[vLabel[w]].insert(w);
			vReached[w] = true;
			vLabel[w] += 0.5;
			// mark (v,w) as an edge of G*a
		}

		// search
		for (size_t j = 1; j <= k; ++j)
		{
			std::set<size_t>& labelSet = mReach[(float) j];
			while (!labelSet.empty())
			{
				size_t w = *labelSet.begin();
				labelSet.erase(w);

				for (size_t k = 0; k < vAdj[w].size(); ++k)
				{
					size_t z = vAdj[w][k];
					if (vReached[z]) continue;

					vReached[z] = true;

					if (vLabel[z] > j)
					{
						mReach[vLabel[z]].insert(z);
						vLabel[z] += 0.5;
						// mark (v,z) as an edge of G*a
					}
					else labelSet.insert(z);
				}
			}
		}

		// sort unnumbered vertices by label value
		std::vector<size_t> vUnhandled;
		for (size_t j = 0; j < nVrt; ++j)
			if (!vHandled[j])
				vUnhandled.push_back(j);

		// end criterion
		if (!vUnhandled.size()) break;

		LabelCmp cmp(vLabel);
		std::stable_sort(vUnhandled.begin(), vUnhandled.end(), cmp);

		// reassign labels to be integers from 0 to k, redefining k appropriately
		float oldLabel = vLabel[vUnhandled[0]];
		k = 1.0;
		vLabel[vUnhandled[0]] = k;
		for (size_t j = 1; j < vUnhandled.size(); ++j)
		{
			float& label = vLabel[vUnhandled[j]];
			if (label - oldLabel > 0.25)
			{
				oldLabel = label;
				k += 1.0;
			}
			label = k;
		}

		highestLabelVrt = vUnhandled.back();
	}

	// now assign new indices
	newIndOut.clear();
	newIndOut.resize(nVrt, (size_t) -1);

	// If we are ordering based on geometrical connectivity,
	// we want to keep together indices of the same geometric object.
	// The input vAdj will therefore only contain adjacency information
	// for the first DoF of every object.
	// During the re-ordering, the other DoFs are ignored and only later
	// inserted into the new order along with their first DoF representative.
	// Attention: The whole thing only works if on every object,
	//            the number of DoFs is the same.
	if (preserveConsec)
	{
		size_t cnt = 0;

		for (size_t oldInd = 0; oldInd < nVrt; ++oldInd)
		{
			// skip non-sorted indices
			if (!vAdj[oldInd].size()) continue;

			// get current entry in vNewOrder
			UG_ASSERT(cnt < vNewOrder.size(), "cnt: " << cnt << ", ordered: " << vNewOrder.size())
			const size_t newInd = vNewOrder[vNewOrder.size() - 1 - cnt];
			++cnt;

			// give the current vNewOrder entry the current index
			UG_ASSERT(newInd < newIndOut.size(), "newInd: " << newInd << ", size: " << newIndOut.size())
			newIndOut[newInd] = oldInd;
		}

		// check that all ordered indices have been written
		if (cnt != vNewOrder.size())
			UG_THROW("Not all indices sorted that must be sorted: "
					<< cnt << " written, but should write: " << vNewOrder.size());

		// fill non-sorted indices (preserving consecutive indexing)
		for (size_t i = 1; i < nVrt; ++i)
			if (newIndOut[i] == (size_t) -1)
				newIndOut[i] = newIndOut[i-1] + 1;
	}

	// If we are ordering based on matrix entries,
	// any unconnected vertex will simply be moved to the bottom.
	else
	{
		size_t newOrdSz = vNewOrder.size();

		for (size_t i = 0; i < newOrdSz; ++i)
		{
			size_t oldInd = vNewOrder[newOrdSz - 1 - i];
			newIndOut[oldInd] = i;
		}

		// move unconnected indices to the bottom of the ordering
		for (size_t i = 0; i < nVrt; ++i)
		{
			if (newIndOut[i] == (size_t) -1)
				newIndOut[i] = newOrdSz++;
		}
	}

	// check that permutation is bijective
	std::vector<size_t> invPerm(nVrt);
	for (size_t i = 0; i < invPerm.size(); ++i)
		invPerm[i] = (size_t) (-1);

	for (size_t i = 0; i < invPerm.size(); ++i)
	{
		UG_COND_THROW(invPerm[newIndOut[i]] != (size_t) (-1), "not a bijective permutation "
			"(double mapping to index " << newIndOut[i] << " by indices " << invPerm[newIndOut[i]] << " and " << i << ")!");
		invPerm[newIndOut[i]] = i;
	}
}


template <typename TDomain>
void reorder_dof_distros_lex(SmartPtr<ApproximationSpace<TDomain> > approx)
{
	UG_COND_THROW(!approx.valid(), "Approximation space is invalid.");
	std::vector<SmartPtr<DoFDistribution> > vDD = approx->dof_distributions();
	ConstSmartPtr<MGSubsetHandler> sh = approx->subset_handler();

	// loop dof distributions
	for (size_t i = 0; i < vDD.size(); ++i)
	{
		SmartPtr<DoFDistribution> dd = vDD[i];

		// get adjacency graph
		std::vector<std::vector<size_t> > vvConnection;
		try {dd->get_connections(vvConnection);}
		UG_CATCH_THROW("No adjacency graph available.");

		// get mapping for LEX M ordering
		std::vector<size_t> vNewIndex;
		lex_order(vNewIndex, vvConnection);

		// reorder indices
		dd->permute_indices(vNewIndex);
	}
}



#ifdef UG_DIM_1
	template void reorder_dofs<Domain1d>(SmartPtr<ApproximationSpace<Domain1d> >, const char*);
	template void reorder_dof_distros_lex<Domain1d>(SmartPtr<ApproximationSpace<Domain1d> >);
#endif
#ifdef UG_DIM_2
	template void reorder_dofs<Domain2d>(SmartPtr<ApproximationSpace<Domain2d> >, const char*);
	template void reorder_dof_distros_lex<Domain2d>(SmartPtr<ApproximationSpace<Domain2d> >);
#endif
#ifdef UG_DIM_3
	template void reorder_dofs<Domain3d>(SmartPtr<ApproximationSpace<Domain3d> >, const char*);
	template void reorder_dof_distros_lex<Domain3d>(SmartPtr<ApproximationSpace<Domain3d> >);
#endif


} // namespace nernst_planck
} // namespace ug


