/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Markus Breit
 * Creation date: 2019-06-28
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
