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
number writeResidualsToFile
(
	SmartPtr<TGridFunction> sol1,
	SmartPtr<TGridFunction> sol2,
	const char* cmp,
	const char* fileName
)
{
	FunctionGroup fctGrp1, fctGrp2;
	try
	{
		fctGrp1 = sol1->fct_grp_by_name(cmp);
		fctGrp2 = sol2->fct_grp_by_name(cmp);
	}
	UG_CATCH_THROW("At least one of the functions in '" << cmp
					<< "' is not contained in the approximation space of one of the grid functions"
					" (or something else was wrong).");

	std::ofstream ofs(fileName, std::ios::app);

	if (!ofs.is_open())
	{
		UG_THROW("Could not open output file '" << fileName << "'.");
	}

	//	get vertex iterator for current subset
	typedef typename DoFDistribution::traits<Vertex>::const_iterator itType;

	ConstSmartPtr<DoFDistribution> dd1 = sol1->approx_space()->dof_distribution(GridLevel::TOP);
	ConstSmartPtr<DoFDistribution> dd2 = sol2->approx_space()->dof_distribution(GridLevel::TOP);

	itType it1 = dd1->template begin<Vertex>();
	itType it2 = dd2->template begin<Vertex>();
	itType iterEnd1 = dd1->template end<Vertex>();
	itType iterEnd2 = dd2->template end<Vertex>();

	number l2_sq = 0.0;
	for (; it1 != iterEnd1 && it2 != iterEnd2; ++it1, ++it2)
	{
		for (size_t fi = 0; fi < fctGrp1.size(); fi++)
		{
			std::vector<DoFIndex> ind1, ind2;
			dd1->dof_indices(*it1, fctGrp1[fi], ind1);
			dd2->dof_indices(*it2, fctGrp2[fi], ind2);

			number res = DoFRef(*sol1, ind1[0]) - DoFRef(*sol2, ind2[0]);
			l2_sq += res*res;

			ofs << res << " ";
		}
	}

	if (it1 != iterEnd1 || it2 != iterEnd2)
	{
		UG_THROW("Not the same number of vertices in the grids of the two functions.");
	}

	return l2_sq;
}



} // namspace calciumDynamics
} // namespace ug
