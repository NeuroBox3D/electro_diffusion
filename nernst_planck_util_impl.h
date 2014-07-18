/*
 * nernst_planck_util.h
 *
 *  Created on: 17.07.2014
 *      Author: mbreit
 */

#include "nernst_planck_util.h"


namespace ug{
namespace nernst_planck{


template <typename TGridFunction>
number writeResidualsToFile(SmartPtr<TGridFunction> sol1, SmartPtr<TGridFunction> sol2, const char* fileName)
{
	size_t nDof = sol1->num_dofs();

	if (sol2->num_dofs() != nDof)
	{
		UG_THROW("Number of DoFs is not the same for both grid functions.");
	}

	std::ofstream ofs(fileName, std::ios::app);

	if (!ofs.is_open())
	{
		UG_THROW("Could not open output file '" << fileName << "'.");
	}


	number l2_sq = 0.0;
	for (size_t i = 0; i < nDof; i++)
	{
		number res = (*sol1)[i] - (*sol2)[i];
		l2_sq += res*res;

		ofs << res << " ";
	}

	return l2_sq;
}



} // namspace calciumDynamics
} // namespace ug
