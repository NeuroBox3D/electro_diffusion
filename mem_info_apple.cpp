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

	communicate_process_values();
}


} // namespace nernst_planck
} // namespace ug
