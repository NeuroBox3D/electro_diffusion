/*
 * redistribution_util.cpp
 *
 *  Created on: 2018-03-19
 *      Author: mbreit
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

