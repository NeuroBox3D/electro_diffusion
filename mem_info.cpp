/*
 * mem_info.cpp
 *
 *  Created on: 2019-06-28
 *      Author: mbreit
 */

#include "mem_info.h"

#include "pcl/pcl_base.h"  // for NumProcs
#include "pcl/pcl_process_communicator.h"  // for NumProcs


namespace ug {
namespace nernst_planck {


void MemInfo::communicate_process_values()
{
#ifdef UG_PARALLEL
	if (pcl::NumProcs() > 1)
	{
		pcl::ProcessCommunicator pc;
		pc.allreduce(&m_locRes, &m_gloRes, 1, PCL_RO_SUM);
		pc.allreduce(&m_locVirt, &m_gloVirt, 1, PCL_RO_SUM);
		pc.allreduce(&m_locRes, &m_maxRes, 1, PCL_RO_MAX);
		pc.allreduce(&m_locVirt, &m_maxVirt, 1, PCL_RO_MAX);
	}
#endif
}

number MemInfo::local_resident_memory() const
{
	return m_locRes;
}

number MemInfo::local_virtual_memory() const
{
	return m_locVirt;
}

number MemInfo::global_resident_memory() const
{
	return m_gloRes;
}

number MemInfo::global_virtual_memory() const
{
	return m_gloVirt;
}

number MemInfo::max_resident_memory() const
{
	return m_maxRes;
}

number MemInfo::max_virtual_memory() const
{
	return m_maxVirt;
}


} // namespace nernst_planck
} // namespace ug
