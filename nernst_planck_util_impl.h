/*
 * nernst_planck_util.h
 *
 *  Created on: 17.07.2014
 *      Author: mbreit
 */

#include "nernst_planck_util.h"

#include "../MembranePotentialMapping/vm2ug_rework.h"	// Mapper
#include "lib_disc/spatial_disc/disc_util/fv1_geom.h"	// DimFV1Geometry
#include "lib_disc/function_spaces/dof_position_util.h"	// DoFPosition

#include <limits>

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
void exportSolution
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

// helper function
template <typename TGridFunction, typename TBaseElem>
void importSolution
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

template <typename TGridFunction>
void importSolution
(
	SmartPtr<TGridFunction> solution,
	const char* subsetNames,
	const char* functionNames,
	const char* inFileBaseName
)
{
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
}


template <typename TBaseElem, typename TGridFunction>
void scale_dof_indices
(
	ConstSmartPtr<DoFDistribution> dd,
	SmartPtr<TGridFunction> vecOut,
	ConstSmartPtr<TGridFunction> vecIn,
	const std::vector<number>& vScale
)
{
	typename DoFDistribution::traits<TBaseElem>::const_iterator iter, iterEnd;
	std::vector<DoFIndex> vInd;

	try
	{
		// iterate all elements (including SHADOW_RIM_COPY!)
		iter = dd->template begin<TBaseElem>(SurfaceView::ALL);
		iterEnd = dd->template end<TBaseElem>(SurfaceView::ALL);
		for (; iter != iterEnd; ++iter)
		{
			for (size_t fi = 0; fi < dd->num_fct(); ++fi)
			{
				size_t nInd = dd->inner_dof_indices(*iter, fi, vInd);

				// remember multi indices
				for (size_t dof = 0; dof < nInd; ++dof)
					DoFRef(*vecOut, vInd[dof]) = vScale[fi] * DoFRef(*vecIn, vInd[dof]);
			}
		}
	}
	UG_CATCH_THROW("Error while scaling vector.")
}



template <typename TGridFunction>
void scale_dimless_vector
(
	SmartPtr<TGridFunction> scaledVecOut,
	ConstSmartPtr<TGridFunction> dimlessVecIn,
	const std::vector<number>& scalingFactors
)
{
	// check that the correct numbers of scaling factors are given
	size_t n = scalingFactors.size();
	UG_COND_THROW(n != dimlessVecIn->num_fct(), "Number of scaling factors (" << n << ") "
			"does not match number of functions given in dimless vector (" << dimlessVecIn->num_fct() << ").");

	// check that input and output vectors have the same number of components and dofs
	UG_COND_THROW(n != scaledVecOut->num_fct(), "Input and output vectors do not have "
			"the same number of functions (" << n << " vs. " << scaledVecOut->num_fct() << ").");
	for (size_t fct = 0; fct < n; ++fct)
	{
		UG_COND_THROW(dimlessVecIn->num_dofs(fct) != scaledVecOut->num_dofs(fct),
				"Input and output vectors do not have the same number of DoFs for function " << fct
				<< " (" << dimlessVecIn->num_dofs(fct) << " vs. " << scaledVecOut->num_dofs(fct) << ").");
	}

	ConstSmartPtr<DoFDistribution> dd = dimlessVecIn->dof_distribution();

	if (dd->max_dofs(VERTEX))
		scale_dof_indices<Vertex, TGridFunction>(dd, scaledVecOut, dimlessVecIn, scalingFactors);
	if (dd->max_dofs(EDGE))
		scale_dof_indices<Edge, TGridFunction>(dd, scaledVecOut, dimlessVecIn, scalingFactors);
	if (dd->max_dofs(FACE))
		scale_dof_indices<Face, TGridFunction>(dd, scaledVecOut, dimlessVecIn, scalingFactors);
	if (dd->max_dofs(VOLUME))
		scale_dof_indices<Volume, TGridFunction>(dd, scaledVecOut, dimlessVecIn, scalingFactors);
}




} // namspace calciumDynamics
} // namespace ug
