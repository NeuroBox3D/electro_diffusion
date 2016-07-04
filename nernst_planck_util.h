/*
 * nernst_planck_util.h
 *
 *  Created on: 17.07.2014
 *      Author: mbreit
 */

#ifndef UG__PLUGINS__EXPERIMENTAL__NERNST_PLANCK__NERNST_PLANCK_UTIL_H
#define UG__PLUGINS__EXPERIMENTAL__NERNST_PLANCK__NERNST_PLANCK_UTIL_H


#include "common/common.h"
#include "lib_disc/function_spaces/approximation_space.h"
#include "lib_disc/function_spaces/grid_function.h"

#include <iostream>
#include <fstream>


namespace ug {
namespace nernst_planck {


template <typename TGridFunction>
number writeResidualsToFile
(
	SmartPtr<TGridFunction> sol1,
	SmartPtr<TGridFunction> sol2,
	const char* cmp,
	const char* fileName
);


/// adjusts interface after a refinement of the geometry
/**
 * When a geometry with an interface is refined, the interface node on the full-dimensional
 * side will no longer be correctly located on the top level (one layer too far from the interface).
 * We therefore need to re-locate this interface node every time we refine the interface (globally).
 */
template <typename TDomain>
void adjust_geom_after_refinement
(
	SmartPtr<ApproximationSpace<TDomain> > approx,
	const char* fullDimIntfNodeSubset,
	const char* lowDimIntfNodeSubset
);



/**
 * \brief outputs solution for specified functions on specified subsets to file
 *
 *	The solutions (separately for each function) are written to (a) continuing
 *	file(s) specified by the user. A new filename is created for each time step.
 *	At the moment, this is only functional, if UG is not used in parallel mode
 *	- for my lack of parallelizing knowledge.
 *	The function could later serve for saving results of a simulation that can
 *	later be used to "feed" the unknown functions with start values in a
 *	consecutive simulation.
 *
 * \param solution		the vector containing the unknowns of the problem
 * \param approx		the underlying approximation space
 * \param time			the simulation time which the averaged value is taken at
 * 						and shall be accorded to in the output file
 * \param subsetNames	contains the names of the subsets that the averaging
 * 						is to be performed on, separated by commas;
 * 						empty string for all subsets
 * \param functionNames	contains the names of the functions that the averaging
 * 						is to be performed for, separated by commas;
 * 						for each function, a separate file with the function name
 * 						as a suffix will be created;
 * 						empty string for all functions
 * \param outFileName	the name of the output file(s), i.e. their prefix
 *
 * \warning This function is very old and will probably not work properly.
 */
template <typename TVector, typename TDomain>
void exportSolution
(
	const TVector& solution,
	const ApproximationSpace<TDomain>& approx,
	const number time,
	const char* subsetNames,
	const char* functionNames,
	const char* outFileName
);


/**
 * \brief Inputs solution for a specified function on specified subsets from file
 * 		  and sets them as a constraint to the solution.
 *
 *	The solutions (separately for each function) are read from a file specified
 *	by the user.
 *	At the moment, this is only functional, if UG is not used in parallel mode.
 *
 * \param solution		the vector the solution is to be written to
 * \param subsetNames	subsets the solution is to be specified on
 * \param functionName	function the solution is to be specified for
 * \param inFileName	the name of the input file
 */
template <typename TGridFunction>
void importSolution
(
	SmartPtr<TGridFunction> solution,
	const char* subsetNames,
	const char* functionName,
	const char* inFileName
);



template <typename TGridFunction>
void scale_dimless_vector
(
	SmartPtr<TGridFunction> scaledVecOut,
	ConstSmartPtr<TGridFunction> dimlessVecIn,
	const std::vector<number>& scalingFactors
);



#if 0
// does not work: would (at least) need implementation of TestLayout for multi_level_layout_tag
#ifdef UG_PARALLEL
#include "lib_algebra/parallelization/parallelization_util.h"
#include "pcl/pcl_layout_tests.h"

struct VertexCoordinate
{
	VertexCoordinate(SmartPtr<MultiGrid> mg, size_t dim)
	{
		m_mg = mg;
		m_dim = dim;
	}

	number operator() (const Vertex& v)
	{
		Grid::VertexAttachmentAccessor<APosition> aaPos(*m_mg, aPosition);
		return aaPos[&v][m_dim];
	}

	SmartPtr<MultiGrid> m_mg;
	size_t m_dim;
};


template <typename TDomain>
void TestPositions(SmartPtr<TDomain> dom)
{
	SmartPtr<MultiGrid> mg = dom->grid();
	pcl::InterfaceCommunicator<VertexLayout> intfCom;
	DistributedGridManager& dgm = *mg->distributed_grid_manager();
	GridLayoutMap& glm = dgm.grid_layout_map();
	VertexLayout* master_layout;
	VertexLayout* slave_layout;

	if (glm.has_layout<Vertex>(INT_H_SLAVE))
		master_layout = &glm.get_layout<Vertex>(INT_H_SLAVE);
	if (glm.has_layout<Vertex>(INT_H_MASTER))
		slave_layout = &glm.get_layout<Vertex>(INT_H_MASTER);


	for (size_t cmp = 0; cmp < TDomain::dim; ++cmp)
	{
		if (!pcl::TestLayout<VertexLayout, number>(pcl::ProcessCommunicator(), intfCom,
                *master_layout, *slave_layout, false, VertexCoordinate(mg, cmp), true))
		{
			UG_LOGN("Something went wrong.");
		}
	}
}

#endif
#endif

} // namspace calciumDynamics
} // namespace ug

#include "nernst_planck_util_impl.h"

#endif // UG__PLUGINS__EXPERIMENTAL__NERNST_PLANCK__NERNST_PLANCK_UTIL_H
