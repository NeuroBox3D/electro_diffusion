/*
 * redistribution_util.h
 *
 *  Created on: 2018-03-19
 *      Author: mbreit
 */

#ifndef REDISTRIBUTION_UTIL_H_
#define REDISTRIBUTION_UTIL_H_

#include <cstddef>  // for size_t
#include "common/util/smart_pointer.h"  // for SmartPtr
#include "lib_disc/parallelization/domain_load_balancer.h"  // for DomainLoadBalancer


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
);

} // namespace nernst_planck
} // namespace ug

#endif // REDISTRIBUTION_UTIL_H_
