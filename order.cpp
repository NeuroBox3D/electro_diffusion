/*
 * order.cpp
 *
 *  Created on: 12.08.2015
 *      Author: mbreit
 */

#include "order.h"

#include "common/error.h"
#include "lib_disc/domain.h"


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



#ifdef UG_DIM_1
	template void reorder_dofs<Domain1d>(SmartPtr<ApproximationSpace<Domain1d> >, const char*);
#endif
#ifdef UG_DIM_2
	template void reorder_dofs<Domain2d>(SmartPtr<ApproximationSpace<Domain2d> >, const char*);
#endif
#ifdef UG_DIM_3
	template void reorder_dofs<Domain3d>(SmartPtr<ApproximationSpace<Domain3d> >, const char*);
#endif


} // namespace nernst_planck
} // namespace ug


