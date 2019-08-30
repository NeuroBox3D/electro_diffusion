/*
 * mem_info.cpp
 *
 *  Created on: 2019-06-28
 *      Author: mbreit
 *      from: Don Wakefield's answer to
 *      https://stackoverflow.com/questions/669438/how-to-get-memory-usage-at-runtime-using-c
 */

#include "mem_info.h"

#include <unistd.h>
#include <ios>
#include <fstream>
#include <string>


namespace ug {
namespace nernst_planck {


void MemInfo::memory_consumption()
{
	// 'file' stat seems to give the most reliable results
	std::ifstream stat_stream("/proc/self/stat", std::ios_base::in);

	// dummy vars for leading entries in stat that we don't care about
	std::string dummy;
	long rss;

	stat_stream >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy
			   >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy
			   >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy
			   >> dummy >> dummy >> dummy >> m_locVirt >> rss;

	stat_stream.close();

	long page_size_kb = sysconf(_SC_PAGE_SIZE);
	m_locRes = rss * page_size_kb;

	communicate_process_values();
}


} // namespace nernst_planck
} // namespace ug
