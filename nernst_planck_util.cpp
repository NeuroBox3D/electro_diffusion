/*
 * nernst_planck_util.cpp
 *
 *  Created on: 17.07.2014
 *      Author: mbreit
 */

#include "nernst_planck_util.h"

#include <cstddef>                                                  // for size_t, NULL
#include <fstream>                                                  // for operator<<, basic_ostream, char...
#include <limits>                                                   // for numeric_limits, numeric_limits<...
#include <sstream>                                                  // ostringstream

#include "common/assert.h"                                          // for UG_ASSERT
#include "common/error.h"                                           // for UG_CATCH_THROW, UG_THROW, UG_CO...
#include "lib_algebra/cpu_algebra_types.h"                          // for CPUAlgebra
#include "lib_disc/common/function_group.h"                         // for FunctionGroup
#include "lib_disc/common/multi_index.h"                            // for DoFIndex, DoFRef
#include "lib_disc/dof_manager/dof_distribution.h"                  // for DoFDistribution, DoFDistributio...
#include "lib_disc/domain.h"                                        // for Domain1d, Doma...
#include "lib_disc/function_spaces/approximation_space.h"           // for ApproximationSpace
#include "lib_disc/function_spaces/dof_position_util.h"             // for InnerDoFPosition
#include "lib_disc/function_spaces/grid_function.h"                 // for GridFunction
#include "lib_disc/local_finite_element/local_finite_element_id.h"  // for LFEID
#include "lib_grid/grid/grid_base_object_traits.h"                  // for VertexIterator
#include "lib_grid/grid/grid_base_objects.h"                        // for Vertex (ptr only), GridBaseObje...
#ifdef UG_PARALLEL
	#include "lib_grid/parallelization/distributed_grid.h"          // for DistributedGridManager, Element...
#endif
#include "lib_grid/multi_grid.h"                                    // for MultiGrid
#include "lib_grid/tools/grid_level.h"                              // for GridLevel, GridLevel::::TOP
#include "lib_grid/tools/subset_group.h"                            // for SubsetGroup
#include "lib_grid/tools/subset_handler_multi_grid.h"               // for MultiGridSubsetHandler
#include "lib_grid/tools/surface_view.h"                            // for SurfaceView::ConstSurfaceViewEl...
#ifdef UG_PARALLEL
	#include "pcl/pcl_base.h"                                       // for NumProcs
#endif

#ifdef NPWithMPM
	#include "../MembranePotentialMapping/vm2ug_rework.h"               // for Mapper
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




// /////////////////////// //
// export solution command //
// /////////////////////// //

// helper function
template <typename TGridFunction, typename TBaseElem>
static void exportSolution
(
	SmartPtr<TGridFunction> solution,
	size_t si,
	size_t fi,
	std::ofstream& ofs
)
{
	typedef typename TGridFunction::domain_type domain_type;

	// retrieve domain and dofDistr from approxSpace
	ConstSmartPtr<domain_type> domain = solution->approx_space()->domain();
	ConstSmartPtr<DoFDistribution> dofDistr = solution->dof_distribution();
	const LFEID lfeid = dofDistr->lfeid(fi);

	//	get elem iterator for current subset and elem type
	typedef typename DoFDistribution::traits<TBaseElem>::const_iterator itType;
	itType iter = dofDistr->template begin<TBaseElem>(si);
	itType iterEnd = dofDistr->template end<TBaseElem>(si);

	// loop over all elems
	for (; iter != iterEnd; ++iter)
	{
		// get current vertex
		TBaseElem* elem = *iter;

		// get coords
		std::vector<typename domain_type::position_type> coords;
		InnerDoFPosition<domain_type>(coords, elem, *domain, lfeid);

		// get multi-indices
		std::vector<DoFIndex> multInd;
		dofDistr->inner_dof_indices(elem, fi, multInd);

		UG_ASSERT(coords.size() == multInd.size(), "#DoF mismatch");

		// get values of DoFs
		number val;
		size_t nDof = multInd.size();
		for (size_t dof = 0; dof < nDof; ++dof)
		{
			val = DoFRef(*solution, multInd[dof]);

			// write solution to file
			for (size_t i = 0; i < coords[dof].size(); ++i)
				ofs << coords[dof][i] << " ";
			ofs << val << std::endl;
		}
	}
}


template <typename TGridFunction>
void exportSolution
(
	SmartPtr<TGridFunction> solution,
	const number time,
	const char* subsetNames,
	const char* functionNames,
	const char* outFileName
)
{
	typedef typename TGridFunction::domain_type domain_type;

	// retrieve dofDistr from solution
	ConstSmartPtr<DoFDistribution> dofDistr = solution->dof_distribution();

	// get subset group to be measured on (if none provided: take all)
	SubsetGroup ssGrp;
	if (!*subsetNames)
	{
		ssGrp.set_subset_handler(dofDistr->subset_handler());
		ssGrp.add_all();
	}
	else
	{
		try {ssGrp = dofDistr->subset_grp_by_name(subsetNames);}
		UG_CATCH_THROW("At least one of the subsets in '" << subsetNames
				<< "' is not contained in the approximation space (or something else was wrong).");
	}

	// get function group to be measured (if none provided: take all)
	FunctionGroup fctGrp;
	if (!*functionNames)
	{
		fctGrp.set_function_pattern(dofDistr->function_pattern());
		fctGrp.add_all();
	}
	else
	{
		try {fctGrp = dofDistr->fct_grp_by_name(functionNames);}
		UG_CATCH_THROW("At least one of the functions in '" << functionNames
						<< "' is not contained in the approximation space (or something else was wrong).");
	}

	// loop functions
	for (size_t fi = 0; fi < fctGrp.size(); fi++)
	{
		// construct outFile name
		std::ostringstream ofnss(outFileName, std::ios_base::app);
		ofnss << "_" << time << "_" << fctGrp.name(fi);

		// create if first time step, append otherwise
		std::ofstream outFile;
		outFile.precision(std::numeric_limits<number>::digits10);
		outFile.open(ofnss.str().c_str(), std::ios_base::out);
		if (!outFile.is_open())
			UG_THROW("Output file '" << ofnss.str() << "' could not be opened.")

		try
		{
			// loop subsets
			for (size_t si = 0; si < ssGrp.size(); si++)
			{
				if (domain_type::dim-1 >= VERTEX && dofDistr->max_fct_dofs(fctGrp[fi], VERTEX, ssGrp[si]) > 0)
					exportSolution<TGridFunction, Vertex>(solution, ssGrp[si], fctGrp[fi], outFile);
				if (domain_type::dim-1 >= EDGE && dofDistr->max_fct_dofs(fctGrp[fi], EDGE, ssGrp[si]) > 0)
					exportSolution<TGridFunction, Edge>(solution, ssGrp[si], fctGrp[fi], outFile);
				if (domain_type::dim-1 >= FACE && dofDistr->max_fct_dofs(fctGrp[fi], FACE, ssGrp[si]) > 0)
					exportSolution<TGridFunction, Face>(solution, ssGrp[si], fctGrp[fi], outFile);
				if (domain_type::dim-1 >= VOLUME && dofDistr->max_fct_dofs(fctGrp[fi], VOLUME, ssGrp[si]) > 0)
					exportSolution<TGridFunction, Volume>(solution, ssGrp[si], fctGrp[fi], outFile);
			}
		}
		UG_CATCH_THROW("Output file '" << ofnss.str() << "' could not be written to.");

		outFile.close();
	}

	return;
}



// /////////////////////// //
// import solution command //
// /////////////////////// //

#ifdef NPWithMPM
// helper function
template <typename TGridFunction, typename TBaseElem>
static void importSolution
(
	SmartPtr<TGridFunction> solution,
	size_t si,
	size_t fi,
	const Mapper<TGridFunction::domain_type::dim, number>& mapper
)
{
	typedef typename TGridFunction::domain_type domain_type;

	// retrieve domain and dofDistr from approxSpace
	ConstSmartPtr<domain_type> domain = solution->approx_space()->domain();
	ConstSmartPtr<DoFDistribution> dofDistr = solution->dof_distribution();
	const LFEID lfeid = dofDistr->lfeid(fi);

	//	get elem iterator for current subset and elem type
	typedef typename DoFDistribution::traits<TBaseElem>::const_iterator itType;
	itType iter = dofDistr->template begin<TBaseElem>(si);
	itType iterEnd = dofDistr->template end<TBaseElem>(si);

	// loop over all elems
	for (; iter != iterEnd; ++iter)
	{
		// get current vertex
		TBaseElem* elem = *iter;

		// get coords
		std::vector<typename domain_type::position_type> coords;
		InnerDoFPosition<domain_type>(coords, elem, *domain, lfeid);

		// get multi-indices
		std::vector<DoFIndex> multInd;
		dofDistr->inner_dof_indices(elem, fi, multInd);

		UG_ASSERT(coords.size() == multInd.size(), "#DoF mismatch");

		// get values of DoFs
		number val;
		size_t nDof = multInd.size();
		for (size_t dof = 0; dof < nDof; ++dof)
		{
			// get value from provider
			try {val = mapper.get_data_from_nearest_neighbor(coords[dof]);}
			UG_CATCH_THROW("No value could be retrieved for DoF at " << coords[dof]);

			DoFRef(*solution, multInd[dof]) = val;
		}
	}
}
#endif

template <typename TGridFunction>
void importSolution
(
	SmartPtr<TGridFunction> solution,
	const char* subsetNames,
	const char* functionNames,
	const char* inFileBaseName
)
{
#ifndef NPWithMPM
	UG_THROW("importSolution uses functionality from the MembranePotentialMapping plugin,"
		"but ewas not compiled with it.");
#else
	typedef typename TGridFunction::domain_type domain_type;
	typedef typename domain_type::position_type pos_type;

	ConstSmartPtr<domain_type> domain = solution->approx_space()->domain();
	ConstSmartPtr<DoFDistribution> dofDistr = solution->dof_distribution();

	// get subset group to be measured on
	SubsetGroup ssGrp;
	if (!*subsetNames)
	{
		ssGrp.set_subset_handler(dofDistr->subset_handler());
		ssGrp.add_all();
	}
	else
	{
		try {ssGrp = dofDistr->subset_grp_by_name(subsetNames);}
		UG_CATCH_THROW("At least one of the subsets in '" << subsetNames
				<< "' is not contained in the approximation space (or something else was wrong).");
	}

	// get function group to be measured (if none provided: take all)
	FunctionGroup fctGrp;
	if (!*functionNames)
	{
		fctGrp.set_function_pattern(dofDistr->function_pattern());
		fctGrp.add_all();
	}
	else
	{
		try {fctGrp = dofDistr->fct_grp_by_name(functionNames);}
		UG_CATCH_THROW("At least one of the functions in '" << functionNames
						<< "' is not contained in the approximation space (or something else was wrong).");
	}


	// loop functions
	for (size_t fi = 0; fi < fctGrp.size(); fi++)
	{
		// construct inFile name
		std::ostringstream ofnss(inFileBaseName, std::ios_base::app);
		ofnss << "_" << fctGrp.name(fi);

		// read values from file and fill mapper structure with it
		Mapper<domain_type::dim, number> valueProvider;
		try {valueProvider.build_tree(ofnss.str(), " ");}
		UG_CATCH_THROW("Underlying mapper object could not build its tree "
					   "on given file (" << ofnss.str() << ").");

		// loop subsets
		for (size_t si = 0; si < ssGrp.size(); si++)
		{
			if (domain_type::dim-1 >= VERTEX && dofDistr->max_fct_dofs(fctGrp[fi], VERTEX, ssGrp[si]) > 0)
				importSolution<TGridFunction, Vertex>(solution, ssGrp[si], fctGrp[fi], valueProvider);
			if (domain_type::dim-1 >= EDGE && dofDistr->max_fct_dofs(fctGrp[fi], EDGE, ssGrp[si]) > 0)
				importSolution<TGridFunction, Edge>(solution, ssGrp[si], fctGrp[fi], valueProvider);
			if (domain_type::dim-1 >= FACE && dofDistr->max_fct_dofs(fctGrp[fi], FACE, ssGrp[si]) > 0)
				importSolution<TGridFunction, Face>(solution, ssGrp[si], fctGrp[fi], valueProvider);
			if (domain_type::dim-1 >= VOLUME && dofDistr->max_fct_dofs(fctGrp[fi], VOLUME, ssGrp[si]) > 0)
				importSolution<TGridFunction, Volume>(solution, ssGrp[si], fctGrp[fi], valueProvider);
		}
	}
#endif
}




// template specializations
#ifdef UG_DIM_1
	template void adjust_geom_after_refinement<Domain1d>(SmartPtr<ApproximationSpace<Domain1d> >, const char*, const char*);

	#ifdef UG_CPU_1
		template number writeResidualsToFile<GridFunction<Domain1d, CPUAlgebra> >(SmartPtr<GridFunction<Domain1d, CPUAlgebra> >, SmartPtr<GridFunction<Domain1d, CPUAlgebra> >, const char*, const char*);
		template void exportSolution<GridFunction<Domain1d, CPUAlgebra> >(SmartPtr<GridFunction<Domain1d, CPUAlgebra> >, const number, const char*, const char*, const char*);
		template void importSolution<GridFunction<Domain1d, CPUAlgebra> >(SmartPtr<GridFunction<Domain1d, CPUAlgebra> >, const char*, const char*, const char*);
	#endif
	#ifdef UG_CPU_5
		template number writeResidualsToFile<GridFunction<Domain1d, CPUBlockAlgebra<5> > >(SmartPtr<GridFunction<Domain1d, CPUBlockAlgebra<5> > >, SmartPtr<GridFunction<Domain1d, CPUBlockAlgebra<5> > >, const char*, const char*);
		template void exportSolution<GridFunction<Domain1d, CPUBlockAlgebra<5> > >(SmartPtr<GridFunction<Domain1d, CPUBlockAlgebra<5> > >, const number, const char*, const char*, const char*);
		template void importSolution<GridFunction<Domain1d, CPUBlockAlgebra<5> > >(SmartPtr<GridFunction<Domain1d, CPUBlockAlgebra<5> > >, const char*, const char*, const char*);
	#endif
	#ifdef UG_CPU_6
		template number writeResidualsToFile<GridFunction<Domain1d, CPUBlockAlgebra<6> > >(SmartPtr<GridFunction<Domain1d, CPUBlockAlgebra<6> > >, SmartPtr<GridFunction<Domain1d, CPUBlockAlgebra<6> > >, const char*, const char*);
		template void exportSolution<GridFunction<Domain1d, CPUBlockAlgebra<6> > >(SmartPtr<GridFunction<Domain1d, CPUBlockAlgebra<6> > >, const number, const char*, const char*, const char*);
		template void importSolution<GridFunction<Domain1d, CPUBlockAlgebra<6> > >(SmartPtr<GridFunction<Domain1d, CPUBlockAlgebra<6> > >, const char*, const char*, const char*);
	#endif
#endif
#ifdef UG_DIM_2
	template void adjust_geom_after_refinement<Domain2d>(SmartPtr<ApproximationSpace<Domain2d> >, const char*, const char*);

	#ifdef UG_CPU_1
		template number writeResidualsToFile<GridFunction<Domain2d, CPUAlgebra> >(SmartPtr<GridFunction<Domain2d, CPUAlgebra> >, SmartPtr<GridFunction<Domain2d, CPUAlgebra> >, const char*, const char*);
		template void exportSolution<GridFunction<Domain2d, CPUAlgebra> >(SmartPtr<GridFunction<Domain2d, CPUAlgebra> >, const number, const char*, const char*, const char*);
		template void importSolution<GridFunction<Domain2d, CPUAlgebra> >(SmartPtr<GridFunction<Domain2d, CPUAlgebra> >, const char*, const char*, const char*);
	#endif
	#ifdef UG_CPU_5
		template number writeResidualsToFile<GridFunction<Domain2d, CPUBlockAlgebra<5> > >(SmartPtr<GridFunction<Domain2d, CPUBlockAlgebra<5> > >, SmartPtr<GridFunction<Domain2d, CPUBlockAlgebra<5> > >, const char*, const char*);
		template void exportSolution<GridFunction<Domain2d, CPUBlockAlgebra<5> > >(SmartPtr<GridFunction<Domain2d, CPUBlockAlgebra<5> > >, const number, const char*, const char*, const char*);
		template void importSolution<GridFunction<Domain2d, CPUBlockAlgebra<5> > >(SmartPtr<GridFunction<Domain2d, CPUBlockAlgebra<5> > >, const char*, const char*, const char*);
	#endif
	#ifdef UG_CPU_6
		template number writeResidualsToFile<GridFunction<Domain2d, CPUBlockAlgebra<6> > >(SmartPtr<GridFunction<Domain2d, CPUBlockAlgebra<6> > >, SmartPtr<GridFunction<Domain2d, CPUBlockAlgebra<6> > >, const char*, const char*);
		template void exportSolution<GridFunction<Domain2d, CPUBlockAlgebra<6> > >(SmartPtr<GridFunction<Domain2d, CPUBlockAlgebra<6> > >, const number, const char*, const char*, const char*);
		template void importSolution<GridFunction<Domain2d, CPUBlockAlgebra<6> > >(SmartPtr<GridFunction<Domain2d, CPUBlockAlgebra<6> > >, const char*, const char*, const char*);
	#endif
#endif
#ifdef UG_DIM_3
	template void adjust_geom_after_refinement<Domain3d>(SmartPtr<ApproximationSpace<Domain3d> >, const char*, const char*);

	#ifdef UG_CPU_1
		template number writeResidualsToFile<GridFunction<Domain3d, CPUAlgebra> >(SmartPtr<GridFunction<Domain3d, CPUAlgebra> >, SmartPtr<GridFunction<Domain3d, CPUAlgebra> >, const char*, const char*);
		template void exportSolution<GridFunction<Domain3d, CPUAlgebra> >(SmartPtr<GridFunction<Domain3d, CPUAlgebra> >, const number, const char*, const char*, const char*);
		template void importSolution<GridFunction<Domain3d, CPUAlgebra> >(SmartPtr<GridFunction<Domain3d, CPUAlgebra> >, const char*, const char*, const char*);
	#endif
	#ifdef UG_CPU_5
		template number writeResidualsToFile<GridFunction<Domain3d, CPUBlockAlgebra<5> > >(SmartPtr<GridFunction<Domain3d, CPUBlockAlgebra<5> > >, SmartPtr<GridFunction<Domain3d, CPUBlockAlgebra<5> > >, const char*, const char*);
		template void exportSolution<GridFunction<Domain3d, CPUBlockAlgebra<5> > >(SmartPtr<GridFunction<Domain3d, CPUBlockAlgebra<5> > >, const number, const char*, const char*, const char*);
		template void importSolution<GridFunction<Domain3d, CPUBlockAlgebra<5> > >(SmartPtr<GridFunction<Domain3d, CPUBlockAlgebra<5> > >, const char*, const char*, const char*);
	#endif
	#ifdef UG_CPU_6
		template number writeResidualsToFile<GridFunction<Domain3d, CPUBlockAlgebra<6> > >(SmartPtr<GridFunction<Domain3d, CPUBlockAlgebra<6> > >, SmartPtr<GridFunction<Domain3d, CPUBlockAlgebra<6> > >, const char*, const char*);
		template void exportSolution<GridFunction<Domain3d, CPUBlockAlgebra<6> > >(SmartPtr<GridFunction<Domain3d, CPUBlockAlgebra<6> > >, const number, const char*, const char*, const char*);
		template void importSolution<GridFunction<Domain3d, CPUBlockAlgebra<6> > >(SmartPtr<GridFunction<Domain3d, CPUBlockAlgebra<6> > >, const char*, const char*, const char*);
	#endif
#endif


} // namspace nernst_planck
} // namespace ug
