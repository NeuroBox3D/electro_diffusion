/*
 * mem_info.h
 *
 *  Created on: 15.02.2018
 *      Author: mbreit
 */

#ifndef UG__PLUGINS__NERNST_PLANCK__MEM_INFO_H_
#define UG__PLUGINS__NERNST_PLANCK__MEM_INFO_H_

#include "common/types.h"  // for number

namespace ug {
namespace nernst_planck {

class MemInfo
{
	public:
		void memory_consumption();
		number local_resident_memory() const;
		number local_virtual_memory() const;
		number global_resident_memory() const;
		number global_virtual_memory() const;

	private:
		size_t m_locRes;
		size_t m_locVirt;
		size_t m_gloRes;
		size_t m_gloVirt;
};

} // namespace nernst_planck
} // namespace ug


#endif // UG__PLUGINS__NERNST_PLANCK__MEM_INFO_H_
