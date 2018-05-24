/*
 * NeckRecorder.h
 *
 *  Created on: 2017-11-29
 *      Author: mbreit
 */

#ifndef UG__PLUGINS__NERNST_PLANCK__NECK_RECORDER_H
#define UG__PLUGINS__NERNST_PLANCK__NECK_RECORDER_H

#include "common/types.h"  // number
#include "common/util/smart_pointer.h"  // SmartPtr

#include "lib_disc/function_spaces/approximation_space.h" // ApproximationSpace
#include "lib_disc/function_spaces/grid_function.h" // ApproximationSpace

#include <string>
#include <vector>


namespace ug {
namespace nernst_planck {


struct NeckRecorderBase
{
	const number F; // 96485 C/mol
	const number R; // 8.31451 J/(mol*K)

	enum {_K_ = 0, _NA_, _CL_, _A_, _PHI_};

	struct IPData
	{
		number weight;
		number detJ;
		std::vector<number> vShapes; // vector traverses corners of volume
		std::vector<number> vGradZ;  // vector traverses corners of volume
	};

	struct IntegrationVolume
	{
		std::vector<std::vector<DoFIndex> > vDofIndex;	// vector traverses corners of volume and then functions
		std::vector<IPData> vIPData;					// vector traverses quadrature IPs of volume
	};

	NeckRecorderBase();
};


/**
 * NeckRecorder class
 *
 * This class can perform measurements of current across and potential on
 * defined cutting planes through a spine neck.
 * Cutting planes are always parallel to the xy plane in 3d and to the x axis in 2d.
 *
 * @note Only implemented for FV1 (hanging nodes allowed) and NoUpwind.
 *
 * @todo implement for arbitrary upwinding (not trivial)!
 * @todo implement for arbitrary FV order
 */
template <typename TDomain, typename TAlgebra>
class NeckRecorder
: public NeckRecorderBase
{
	public:
		static const int dim = TDomain::dim;

	public:
		/// constructor
		NeckRecorder(SmartPtr<ApproximationSpace<TDomain> > approx);

		/// add measurement zone
		void add_measurement_zone(number zCoord, const std::string& name);

		/// set cytosolic subset name
		void set_cytosolic_subset(const std::string& ss);

		/**
		 * @brief set diffusion constants
		 *
		 * Order must be: K, Na, Cl, A.
		 * Default values: K: 1.96e-9, Na: 1.33e-9, Cl: 2.03e-9, A: 2.0e-9 (m^2/s)
		**/
		void set_diffusion_constants(const std::vector<number>& diffConsts);

		/// set temperature (default: 298.15K)
		void set_temperature(number t);

		/// measure current
		void record_current
		(
			const std::string& fileName,
			number time,
			ConstSmartPtr<GridFunction<TDomain, TAlgebra> > u,
			number scale
		);

		/// measure average potential
		void record_potential
		(
			const std::string& fileName,
			number time,
			ConstSmartPtr<GridFunction<TDomain, TAlgebra> > u
		);

	protected:
		/// find intersecting volumes, prepare integration
		void prepare
		(
			std::vector<std::vector<IntegrationVolume> >& integData,
			bool useSCVFMode
		);

	protected:
		/// approx space
		SmartPtr<ApproximationSpace<TDomain> > m_spApprox;

		/// cytosolic subset index
		int m_siCyt;

		/// measurement zones (identified by z coordinate)
		std::vector<number> m_vMeasCoords;
		std::vector<std::string> m_vMeasZoneNames;

		/// diffusion constants (order: K, Na, Cl, A)
		std::vector<number> m_vDiffConst;

		/// temperature in K
		number m_temp;

		/// data for integration
		std::vector<std::vector<IntegrationVolume> > m_vIntegrationDataCurrent; // first index for measurement zone
		std::vector<std::vector<IntegrationVolume> > m_vIntegrationDataPot; // first index for measurement zone

		/// integration data has been prepared
		bool bPreparedCurrent;
		bool bPreparedPot;
};


} // namespace nernst_planck
} // namespace ug

#endif // UG__PLUGINS__NERNST_PLANCK__NECK_RECORDER_H
