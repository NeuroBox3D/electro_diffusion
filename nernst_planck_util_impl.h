/*
 * nernst_planck_util.h
 *
 *  Created on: 17.07.2014
 *      Author: mbreit
 */

#include "nernst_planck_util.h"


namespace ug{
namespace nernst_planck{


template <typename TGridFunction>
number writeResidualsToFile
(
	SmartPtr<TGridFunction> sol1,
	SmartPtr<TGridFunction> sol2,
	const char* cmp,
	const char* fileName
)
{
	FunctionGroup fctGrp1, fctGrp2;
	try
	{
		fctGrp1 = sol1->fct_grp_by_name(cmp);
		fctGrp2 = sol2->fct_grp_by_name(cmp);
	}
	UG_CATCH_THROW("At least one of the functions in '" << cmp
					<< "' is not contained in the approximation space of one of the grid functions"
					" (or something else was wrong).");

	std::ofstream ofs(fileName, std::ios::app);

	if (!ofs.is_open())
	{
		UG_THROW("Could not open output file '" << fileName << "'.");
	}

	//	get vertex iterator for current subset
	typedef typename DoFDistribution::traits<Vertex>::const_iterator itType;

	ConstSmartPtr<DoFDistribution> dd1 = sol1->approx_space()->dof_distribution(GridLevel::TOP);
	ConstSmartPtr<DoFDistribution> dd2 = sol2->approx_space()->dof_distribution(GridLevel::TOP);

	itType it1 = dd1->template begin<Vertex>();
	itType it2 = dd2->template begin<Vertex>();
	itType iterEnd1 = dd1->template end<Vertex>();
	itType iterEnd2 = dd2->template end<Vertex>();

	number l2_sq = 0.0;
	for (; it1 != iterEnd1 && it2 != iterEnd2; ++it1, ++it2)
	{
		for (size_t fi = 0; fi < fctGrp1.size(); fi++)
		{
			std::vector<DoFIndex> ind1, ind2;
			dd1->dof_indices(*it1, fctGrp1[fi], ind1);
			dd2->dof_indices(*it2, fctGrp2[fi], ind2);

			number res = DoFRef(*sol1, ind1[0]) - DoFRef(*sol2, ind2[0]);
			l2_sq += res*res;

			ofs << res << " ";
		}
	}

	if (it1 != iterEnd1 || it2 != iterEnd2)
	{
		UG_THROW("Not the same number of vertices in the grids of the two functions.");
	}

	return l2_sq;
}

template <typename TDomain>
void adjust_geom_after_refinement
(
	SmartPtr<ApproximationSpace<TDomain> > approx,
	const char* innerSubset,
	const char* fullDimIntfNodeSubset,
	const char* lowDimIntfNodeSubset
)
{
	SmartPtr<TDomain> dom = approx->domain();
	SmartPtr<MultiGridSubsetHandler> ssh = dom->subset_handler();
	SmartPtr<MultiGrid> mg = dom->grid();

	// get subset indices for given subsets
	int si_intf = ssh->get_subset_index(fullDimIntfNodeSubset);
	int si_1d_intf = ssh->get_subset_index(lowDimIntfNodeSubset);
	int si_inner = ssh->get_subset_index(innerSubset);

// (A) move interface node of full-dim side if need be (i.e. if it has not already been moved)
	// find full-dim interface node in top_level-1 first
	Vertex* intf = NULL;

	UG_COND_THROW(ssh->num_levels() < 2, "Not enough levels in multigrid. At least two levels must be present "
			"after refinement in order to use this method.");

	uint lvl = ssh->num_levels()-2; // top level - 1
	VertexIterator iter = ssh->begin<Vertex>(si_intf, lvl);
	if (iter == ssh->end<Vertex>(si_intf, lvl))
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
		intf = *iter;

		// to be on the safe side
		if (++iter != ssh->end<Vertex>(si_intf, lvl))
			{UG_THROW("More than one vertex in subset for high-dimensional interface node. This is not allowed!");}
	}

	// stop here if interface node not present on this proc
	if (!intf)
#ifndef UG_PARALLEL
		UG_THROW("No full-dim interface node found on level top-1.");
#else
		return;
	else
	{
		// it may turn out that the level below top is distributed to another
		// process (and therefore exists twice, as vmaster and vslave)
		// in this case, we only need the vslave
		const DistributedGridManager* dgm = mg->distributed_grid_manager();
		if (dgm->contains_status(intf, ES_V_MASTER)) return;
	}
#endif

	// now find new position:
	// find edge connecting full-d intf node to 1d intf node
	// vertex child of this edge is full-d intf node on top level
	typedef typename MultiGrid::traits<Edge>::secure_container edge_list;
	edge_list el;
	dom->grid()->associated_elements(el, intf);
	size_t e = 0;
	for (; e < el.size(); ++e)
	{
		Vertex* other;
		if (!el[e]->get_opposing_side(intf, &other))
			{UG_THROW("No opposing side found!");}

		if (ssh->get_subset_index(other) != si_1d_intf)
			continue;

		break;
	}

	if (e == el.size())
		UG_THROW("New full-dimensional interface node could not be determined.\n"
				 "There might be a processor boundary separating interfaces nodes? "
				 "This is not allowed to happen.");

	Vertex* newIntf = mg->get_child_vertex(el[e]);
	if (!newIntf)
		UG_THROW("New full-dim interface node could not be determined.\n"
				 "The top level interface nodes are involved in a vertical interface. "
				 "This case is not implemented.");

	ssh->assign_subset(newIntf, si_intf);

	// change the subset of the child of top_level-1 full-d intf node to inner subset
	Vertex* prev_intf_top = (*el[e])[0];
	if (!prev_intf_top)
		UG_THROW("New full-dim interface node could not be determined.\n"
				 "The top level interface nodes are involved in a vertical interface. "
				 "This case is not implemented.");
	ssh->assign_subset(mg->get_child_vertex(intf), si_inner);
}


} // namspace calciumDynamics
} // namespace ug
