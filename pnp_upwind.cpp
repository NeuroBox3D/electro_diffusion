/*
 * pnp_upwind.cpp
 *
 *  Created on: 11.07.2016
 *      Author: mbreit
 */

#include "pnp_upwind.h"

#include <cmath>                                                   // for exp
#include <cstddef>                                                 // for size_t, NULL
#include <vector>                                                  // for vector, allocator

#include "common/assert.h"                                         // for UG_ASSERT
#include "common/error.h"                                          // for UG_COND_THROW
#include "common/math/math_vector_matrix/math_vector_functions.h"  // for VecDot, VecNormalize, VecScale
#include "common/util/metaprogramming_util.h"                      // for Int2Type
#include "lib_disc/spatial_disc/disc_util/fv1_geom.h"              // for DimFV1Geometry, FV1Geometry
#include "lib_disc/spatial_disc/disc_util/hfv1_geom.h"             // for HFV1Geometry


//#define DEBUG_PNP_UPWIND

namespace ug {

// forward declarations
class RegularEdge;
class Triangle;
class Quadrilateral;
class Tetrahedron;
class Pyramid;
class Prism;
class Octahedron;
class Hexahedron;

namespace nernst_planck {


template <int worldDim>
PNPUpwind<worldDim>::PNPUpwind()
: m_alpha(1e-3)
{
	// The shapes do not depend on the diffusion. Thus, we can set the
	// derivative to be always zero w.r.t. the diffusion for all shapes.
	set_non_zero_deriv_diffusion_flag(false);

	// register evaluation functions
	register_func(Int2Type<worldDim>());
}

template <int worldDim>
PNPUpwind<worldDim>::~PNPUpwind()
{}



template <int worldDim>
void PNPUpwind<worldDim>::set_alpha(number a)
{
	m_alpha = a;
}



#ifdef DEBUG_PNP_UPWIND
template<int worldDim>
struct DebugHelper
{
	public:
		static DebugHelper& instance()
		{
			static DebugHelper instance;
			return instance;
		}

		void write_to_file(GridObject* elem, const std::string& s)
		{
			// are we in the next iteration?
			bool nextIt = false;
			int elemIt = 0;
			std::map<GridObject*, int>::iterator it = m_elemMap.find(elem);
			if (it == m_elemMap.end())
			{
				elemIt = m_elemMap[elem] = 0;
				if (m_currIt == -1)
				{
					nextIt = true;
					++m_currIt;
				}
			}
			else
			{
				if (it->second == 5*(m_currIt+1)-1)
				{
					nextIt = true;
					++m_currIt;
				}
				elemIt = ++it->second;
			}

			if (nextIt)
			{
				// construct new outFile name
				std::ostringstream ofnss("pnpUpwind", std::ios_base::app);
				ofnss << "_" << m_currIt << ".csv";
				outFileName = ofnss.str();
			}

			std::ofstream outFile;
			outFile.open(outFileName.c_str(), std::ios_base::app);

			if (nextIt)
			{
				std::string header = std::string("coordX");
				if (worldDim >= 2) header.append(std::string(",coordY"));
				if (worldDim >= 3) header.append(std::string(",coordZ"));
				header.append(std::string(",fluxX"));
				if (worldDim >= 2) header.append(std::string(",fluxY"));
				if (worldDim >= 3) header.append(std::string(",fluxZ"));
				outFile << header << std::endl;
			}

			if (elemIt % 5 == 0)
				outFile << s;

			outFile.close();
		}

	private:
		DebugHelper() : m_currIt(-1) {}

		DebugHelper(DebugHelper const&);    // do not implement
		void operator=(DebugHelper const&); // do not implement

	private:
		std::map<GridObject*, int> m_elemMap;
		int m_currIt;
		std::string outFileName;
};
#endif


template <int worldDim>
template <typename TFVGeom>
bool PNPUpwind<worldDim>::
update(const TFVGeom* geo,
       const MathVector<worldDim>* vel,
       const MathMatrix<worldDim, worldDim>* diff,
       bool computeDeriv)
{
#ifdef DEBUG_PNP_UPWIND
	// DEBUG: write upwind directions to file
	static void* onlyConsiderMe = (void*) this;
	std::ostringstream oss;
#endif

    PROFILE_BEGIN_GROUP(PNPUpwind_update, "Upwind");

	UG_ASSERT(geo != NULL, "Null pointer passed as FV geometry.");
	UG_ASSERT(vel != NULL, "Null pointer passed as velocity.");

	const size_t numSh = geo->num_sh();

	const number beta = 20.0;

	// loop sub-control volume faces
	const size_t numSCVF = geo->num_scvf();
	for (size_t s = 0; s < numSCVF; ++s)
	{
		const typename TFVGeom::SCVF& scvf = geo->scvf(s);

		UG_COND_THROW(numSh != scvf.num_sh(),
			"Number of geometry shapes does not equal number of SCVF shapes.");

		// compute flux
		const number flux = VecDot(scvf.normal(), vel[s]);

		// calculate exponentially weighed convection shapes
		const number velNorm = VecTwoNorm(vel[s]);

		// for performance, only resize if arrays have to grow;
		// keep maximal size, simply ignore excess entries
		if (m_vDir.size() < numSh)
		{
			m_vDir.resize(numSh);
			m_vDirVelProdExp.resize(numSh);
		}

		number dirVelProdExpSum = 0.0;
		for (size_t sh = 0; sh < numSh; ++sh)
		{
			VecScaleAdd(m_vDir[sh], 1.0, scvf.global_ip(), -1.0, geo->global_node_position(sh));
			VecNormalize(m_vDir[sh], m_vDir[sh]);
			number dirVelProd = VecDot(m_vDir[sh], vel[s]);
			m_vDirVelProdExp[sh] = exp(beta*dirVelProd/(velNorm + m_alpha));
			dirVelProdExpSum += m_vDirVelProdExp[sh];
		}

		for (size_t sh = 0; sh < numSh; ++sh)
			conv_shape(s, sh) = m_vDirVelProdExp[sh] / dirVelProdExpSum * flux;


		// compute derivatives (w.r.t. vel) if needed
		if (computeDeriv)
		{
			for (size_t sh = 0; sh < numSh; ++sh)
				VecSet(D_vel(s, sh), 0.0);

			MathVector<worldDim> dirDiff;
			MathVector<worldDim> innerDeriv;
			for (size_t sh = 0; sh < numSh; ++sh)
			{
				for (size_t sh2 = 0; sh2 < numSh; ++sh2)
				{
					VecScaleAdd(dirDiff, 1.0, m_vDir[sh], -1.0, m_vDir[sh2]);
					if (velNorm > 1e-8)
						VecScaleAdd(innerDeriv, 1.0/(velNorm + m_alpha), dirDiff,
							-VecDot(vel[s], dirDiff)/(velNorm*(velNorm + m_alpha)*(velNorm + m_alpha)), vel[s]);
					else
						VecScale(innerDeriv, dirDiff, 1.0/(velNorm + m_alpha));
					VecScaleAppend(D_vel(s, sh), beta*flux*m_vDirVelProdExp[sh2], innerDeriv,
							m_vDirVelProdExp[sh2], scvf.normal());
				}
				VecScale(D_vel(s, sh), D_vel(s, sh), m_vDirVelProdExp[sh] / (dirVelProdExpSum*dirVelProdExpSum));
			}
		}

#ifdef DEBUG_PNP_UPWIND
		if (onlyConsiderMe == (void*) this)
		{
			MathVector<worldDim> temp;
			MathVector<worldDim> dir = 0.0;

			for (size_t sh = 0; sh < numSh; ++sh)
			{
				VecScaleAdd(temp, -1.0, scvf.global_ip(), 1.0, geo->global_node_position(sh));
				VecScaleAdd(dir, 1.0, dir, dirVelProdExp[sh] / dirVelProdExpSum, temp);
			}

			for (size_t i = 0; i < worldDim; ++i)
				oss << scvf.global_ip()[i] << ", ";
			for (size_t i = 0; i < worldDim-1; ++i)
				oss << dir[i] << ", ";
			oss << dir[worldDim-1] << "\n";
		}
#endif
	}

#ifdef DEBUG_PNP_UPWIND
	if (onlyConsiderMe == (void*) this)
		DebugHelper<worldDim>::instance().write_to_file(geo->elem(), oss.str());
#endif

	return true;
}



template <int worldDim>
void PNPUpwind<worldDim>::register_func(Int2Type<1>)
{
	register_func_for_refDim<1>();

	register_func_for_elem<RegularEdge>();
}

template <int worldDim>
void PNPUpwind<worldDim>::register_func(Int2Type<2>)
{
	register_func(Int2Type<1>());

	register_func_for_refDim<2>();

	register_func_for_elem<Triangle>();
	register_func_for_elem<Quadrilateral>();
}
template <> void PNPUpwind<1>::register_func(Int2Type<2>) {}

template <int worldDim>
void PNPUpwind<worldDim>::register_func(Int2Type<3>)
{
	register_func(Int2Type<2>());

	register_func_for_refDim<3>();

	register_func_for_elem<Tetrahedron>();
	register_func_for_elem<Pyramid>();
	register_func_for_elem<Prism>();
	register_func_for_elem<Hexahedron>();
	register_func_for_elem<Octahedron>();
}
template <> void PNPUpwind<1>::register_func(Int2Type<3>) {}
template <> void PNPUpwind<2>::register_func(Int2Type<3>) {}


template <int worldDim>
template <typename TElem>
void PNPUpwind<worldDim>::register_func_for_elem()
{
	typedef FV1Geometry<TElem, worldDim> TGeom;
	typedef bool (this_type::*TFunc)
			(const TGeom*, const MathVector<worldDim>*, const MathMatrix<worldDim, worldDim>*, bool);

	base_type::template register_update_func<TGeom, TFunc>(&this_type::template update<TGeom>);

	typedef HFV1Geometry<TElem, worldDim> THGeom;
	typedef bool (this_type::*THFunc)
			(const THGeom*, const MathVector<worldDim>*, const MathMatrix<worldDim, worldDim>*, bool);

	base_type::template register_update_func<THGeom, THFunc>(&this_type::template update<THGeom>);
}

template <int worldDim>
template <int refDim>
void PNPUpwind<worldDim>::register_func_for_refDim()
{
	typedef DimFV1Geometry<refDim, worldDim> TGeom;
	typedef bool (this_type::*TFunc)
			(const TGeom*, const MathVector<worldDim>*, const MathMatrix<worldDim, worldDim>*, bool);

	base_type::template register_update_func<TGeom, TFunc>(&this_type::template update<TGeom>);
}


// explicit template specializations
#ifdef UG_DIM_1
	template class PNPUpwind<1>;
#endif
#ifdef UG_DIM_2
	template class PNPUpwind<2>;
#endif
#ifdef UG_DIM_3
	template class PNPUpwind<3>;
#endif


} // namespace nernst_planck
} // namespace ug
