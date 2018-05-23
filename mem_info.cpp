/*
 * mem_info.cpp
 *
 *  Created on: 15.02.2018
 *      Author: mbreit
 */

#include "mem_info.h"

#include <mach/mach.h>  // for task_info etc.

#include "common/error.h"  // for UG_COND_THROW
#include "pcl/pcl_base.h"  // for NumProcs
#include "pcl/pcl_process_communicator.h"  // for NumProcs

namespace ug {
namespace nernst_planck {


void MemInfo::memory_consumption()
{
	task_basic_info t_info;
	mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;
	UG_COND_THROW(task_info(mach_task_self(), TASK_BASIC_INFO,
			                (task_info_t) &t_info, &t_info_count) != KERN_SUCCESS,
		"Task info could not be obtained.");

	m_locRes = t_info.resident_size;
	m_locVirt = t_info.virtual_size;

	// global
#ifdef UG_PARALLEL
	if (pcl::NumProcs() > 1)
	{
		pcl::ProcessCommunicator pc;
		pc.allreduce(&m_locRes, &m_gloRes, 1, PCL_RO_SUM);
		pc.allreduce(&m_locVirt, &m_gloVirt, 1, PCL_RO_SUM);
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


} // namespace nernst_planck
} // namespace ug
