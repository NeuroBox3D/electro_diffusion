/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Markus Breit
 * Creation date: 2018-03-19
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

#include "redistribution_util.h"

#include "lib_disc/domain.h"  // for Domain3d etc.
#include "lib_grid/parallelization/process_hierarchy.h"  // for ProcessHierarchy


namespace ug {
namespace nernst_planck {


template <typename TDomain>
void redistribute
(
	ConstSmartPtr<TDomain> dom,
	SmartPtr<DomainLoadBalancer<TDomain> > loadBalancer,
	const std::string& qualityRecordName,
	size_t firstDistLevel,
	size_t maxDistLevel,
	size_t redistStep,
	size_t minElems,
	size_t nodeSize
)
{
#ifdef UG_PARALLEL
	size_t nProcs = pcl::NumProcs();

	const DomainInfo& domInf = dom->domain_info();
	size_t nLvl = domInf.num_levels();
	std::vector<size_t> numElemsOnLvl(nLvl);
	for (size_t i = 0; i < nLvl; ++i)
		numElemsOnLvl[i] = domInf.num_elements_on_level(i);

	SmartPtr<ProcessHierarchy> ph = make_sp(new ProcessHierarchy());

	UG_LOG("Process hierarchy: ");
	nLvl = std::min(nLvl, (size_t) maxDistLevel);
	size_t lastProcs = 1;
	for (size_t i = firstDistLevel; i < nLvl; i += redistStep)
	{
		size_t nProcOnDistroLevel = (numElemsOnLvl[i] / minElems / nodeSize) * nodeSize;
		if (nProcOnDistroLevel > nProcs)
			nProcOnDistroLevel = nProcs;
		size_t redistProcs = nProcOnDistroLevel /= lastProcs;
		redistProcs = std::max((size_t) 1, redistProcs);
		ph->add_hierarchy_level(i, redistProcs);
		lastProcs *= redistProcs;

		UG_LOG("(l" << i << ", " << nProcOnDistroLevel << ") ");
	}
	UG_LOGN("");

	loadBalancer->set_next_process_hierarchy(ph);
	loadBalancer->rebalance();
	loadBalancer->create_quality_record(qualityRecordName.c_str());
#endif
}


#ifdef UG_DIM_1
	template void redistribute<Domain1d>(ConstSmartPtr<Domain1d>, SmartPtr<DomainLoadBalancer<Domain1d> >,
		const std::string&, size_t, size_t, size_t, size_t, size_t);
#endif
#ifdef UG_DIM_2
	template void redistribute<Domain2d>(ConstSmartPtr<Domain2d>, SmartPtr<DomainLoadBalancer<Domain2d> >,
		const std::string&, size_t, size_t, size_t, size_t, size_t);
#endif
#ifdef UG_DIM_3
	template void redistribute<Domain3d>(ConstSmartPtr<Domain3d>, SmartPtr<DomainLoadBalancer<Domain3d> >,
		const std::string&, size_t, size_t, size_t, size_t, size_t);
#endif


} // namespace nernst_planck
} // namespace ug

