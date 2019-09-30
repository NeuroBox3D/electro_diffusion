/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Markus Breit
 * Creation date: 2014-07-17
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

#include "nernst_planck_util.h"

#include <cstddef>                                                  // for size_t, NULL
#include <fstream>                                                  // for operator<<, basic_ostream, char...

#include "common/error.h"                                           // for UG_CATCH_THROW, UG_THROW, UG_CO...
#include "lib_algebra/cpu_algebra_types.h"                          // for CPUAlgebra
#include "lib_disc/common/function_group.h"                         // for FunctionGroup
#include "lib_disc/common/multi_index.h"                            // for DoFIndex, DoFRef
#include "lib_disc/dof_manager/dof_distribution.h"                  // for DoFDistribution, DoFDistributio...
#include "lib_disc/domain.h"                                        // for Domain1d, Doma...
#include "lib_disc/function_spaces/approximation_space.h"           // for ApproximationSpace
#include "lib_disc/function_spaces/grid_function.h"                 // for GridFunction
#include "lib_grid/grid/grid_base_object_traits.h"                  // for VertexIterator
#include "lib_grid/grid/grid_base_objects.h"                        // for Vertex (ptr only), GridBaseObje...
#ifdef UG_PARALLEL
	#include "lib_grid/parallelization/distributed_grid.h"          // for DistributedGridManager, Element...
#endif
#include "lib_grid/multi_grid.h"                                    // for MultiGrid
#include "lib_grid/tools/grid_level.h"                              // for GridLevel, GridLevel::::TOP
#include "lib_grid/tools/subset_handler_multi_grid.h"               // for MultiGridSubsetHandler
#ifdef UG_PARALLEL
	#include "pcl/pcl_base.h"                                       // for NumProcs
#endif

namespace ug {
namespace nernst_planck {


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

	ConstSmartPtr<DoFDistribution> dd1 = sol1->approx_space()->dof_distribution(GridLevel());
	ConstSmartPtr<DoFDistribution> dd2 = sol2->approx_space()->dof_distribution(GridLevel());

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

	// find full-dim interface node in level 0
	Vertex* intf = NULL;


	VertexIterator iter = ssh->begin<Vertex>(si_intf, 0);
	if (iter == ssh->end<Vertex>(si_intf, 0))
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
		if (++iter != ssh->end<Vertex>(si_intf, 0))
			{UG_THROW("More than one vertex in subset for high-dimensional interface node. This is not allowed!");}
	}

	// stop here if interface node not present on this proc
	if (!intf)
#ifndef UG_PARALLEL
		UG_THROW("No full-dim interface node found on level top-1.");
#else
		return;

	// it may turn out that the level below top is distributed to another
	// process (and therefore exists twice, as vmaster and vslave)
	// in this case, we only need the vslave
	const DistributedGridManager* dgm = mg->distributed_grid_manager();
	if (dgm->contains_status(intf, ES_V_MASTER)) return;

#endif

	// find edge connecting full-d intf node to 1d intf node
	typedef typename MultiGrid::traits<Edge>::secure_container edge_list;
	edge_list el;
	dom->grid()->associated_elements(el, intf);
	size_t ei = 0;
	for (; ei < el.size(); ++ei)
	{
		Vertex* other;
		if (!el[ei]->get_opposing_side(intf, &other))
			{UG_THROW("No opposing side found!");}

		if (ssh->get_subset_index(other) != si_1d_intf)
			continue;

		break;
	}

	if (ei == el.size())
		UG_THROW("New full-dimensional interface node could not be determined.\n"
				 "There might be a processor boundary separating interfaces nodes? "
				 "This is not allowed to happen.");

	Edge* e = el[ei];


	// now iterate over multigrid levels
	uint nlvl = ssh->num_levels();
	for (uint lvl = 0; lvl < nlvl; ++lvl)
	{
		// if this edge has a vertex child, this is the full-d intf node on the level above
		Vertex* newIntf = mg->get_child_vertex(e);
		if (newIntf)
		{
			ssh->assign_subset(newIntf, si_intf);

			// change the subset of the child of lower-level full-d intf node to inner subset
			// (which is the subset the edge has)
			ssh->assign_subset(mg->get_child_vertex(intf), ssh->get_subset_index(e));
		}
		else break;

		// prepare for next lvl
		intf = newIntf;
		UG_COND_THROW(mg->num_child_edges(e) != 2, "Incorrect number of child edges.");
		Edge* newEdge = mg->get_child_edge(e, 0);
		if (ssh->get_subset_index((*newEdge)[0]) == si_1d_intf || ssh->get_subset_index((*newEdge)[1]) == si_1d_intf)
			e = newEdge;
		else
			e = mg->get_child_edge(e, 1);
	}
}




// template specializations
#ifdef UG_DIM_1
	template void adjust_geom_after_refinement<Domain1d>(SmartPtr<ApproximationSpace<Domain1d> >, const char*, const char*);

	#ifdef UG_CPU_1
		template number writeResidualsToFile<GridFunction<Domain1d, CPUAlgebra> >(SmartPtr<GridFunction<Domain1d, CPUAlgebra> >, SmartPtr<GridFunction<Domain1d, CPUAlgebra> >, const char*, const char*);
	#endif
	#ifdef UG_CPU_5
		template number writeResidualsToFile<GridFunction<Domain1d, CPUBlockAlgebra<5> > >(SmartPtr<GridFunction<Domain1d, CPUBlockAlgebra<5> > >, SmartPtr<GridFunction<Domain1d, CPUBlockAlgebra<5> > >, const char*, const char*);
	#endif
	#ifdef UG_CPU_6
		template number writeResidualsToFile<GridFunction<Domain1d, CPUBlockAlgebra<6> > >(SmartPtr<GridFunction<Domain1d, CPUBlockAlgebra<6> > >, SmartPtr<GridFunction<Domain1d, CPUBlockAlgebra<6> > >, const char*, const char*);
	#endif
#endif
#ifdef UG_DIM_2
	template void adjust_geom_after_refinement<Domain2d>(SmartPtr<ApproximationSpace<Domain2d> >, const char*, const char*);

	#ifdef UG_CPU_1
		template number writeResidualsToFile<GridFunction<Domain2d, CPUAlgebra> >(SmartPtr<GridFunction<Domain2d, CPUAlgebra> >, SmartPtr<GridFunction<Domain2d, CPUAlgebra> >, const char*, const char*);
	#endif
	#ifdef UG_CPU_5
		template number writeResidualsToFile<GridFunction<Domain2d, CPUBlockAlgebra<5> > >(SmartPtr<GridFunction<Domain2d, CPUBlockAlgebra<5> > >, SmartPtr<GridFunction<Domain2d, CPUBlockAlgebra<5> > >, const char*, const char*);
	#endif
	#ifdef UG_CPU_6
		template number writeResidualsToFile<GridFunction<Domain2d, CPUBlockAlgebra<6> > >(SmartPtr<GridFunction<Domain2d, CPUBlockAlgebra<6> > >, SmartPtr<GridFunction<Domain2d, CPUBlockAlgebra<6> > >, const char*, const char*);
	#endif
#endif
#ifdef UG_DIM_3
	template void adjust_geom_after_refinement<Domain3d>(SmartPtr<ApproximationSpace<Domain3d> >, const char*, const char*);

	#ifdef UG_CPU_1
		template number writeResidualsToFile<GridFunction<Domain3d, CPUAlgebra> >(SmartPtr<GridFunction<Domain3d, CPUAlgebra> >, SmartPtr<GridFunction<Domain3d, CPUAlgebra> >, const char*, const char*);
	#endif
	#ifdef UG_CPU_5
		template number writeResidualsToFile<GridFunction<Domain3d, CPUBlockAlgebra<5> > >(SmartPtr<GridFunction<Domain3d, CPUBlockAlgebra<5> > >, SmartPtr<GridFunction<Domain3d, CPUBlockAlgebra<5> > >, const char*, const char*);
	#endif
	#ifdef UG_CPU_6
		template number writeResidualsToFile<GridFunction<Domain3d, CPUBlockAlgebra<6> > >(SmartPtr<GridFunction<Domain3d, CPUBlockAlgebra<6> > >, SmartPtr<GridFunction<Domain3d, CPUBlockAlgebra<6> > >, const char*, const char*);
	#endif
#endif


} // namspace nernst_planck
} // namespace ug
