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
: m_expFact(1e2)
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
void PNPUpwind<worldDim>::set_exp_factor(number a)
{
	m_expFact = a;
}



template <int worldDim>
template <typename TFVGeom>
bool PNPUpwind<worldDim>::
update(const TFVGeom* geo,
       const MathVector<worldDim>* vel,
       const MathMatrix<worldDim, worldDim>* diff,
       bool computeDeriv)
{
/* DEBUG: write upwind directions to file
	static void* onlyConsiderMe = (void*) this;
*/

	UG_ASSERT(geo != NULL, "Null pointer passed as FV geometry.");
	UG_ASSERT(vel != NULL, "Null pointer passed as velocity.");

	const size_t numSH = geo->num_sh();

/* DEBUG: write upwind directions to file
// construct outFile name
std::ostringstream ofnss("pnpUpwind.csv", std::ios_base::app);
std::ofstream outFile;

if (onlyConsiderMe == (void*) this)
	outFile.open(ofnss.str().c_str(), std::ios_base::app);
*/

	// loop sub-control volume faces
	for (size_t ip = 0; ip < geo->num_scvf(); ++ip)
	{
		const typename TFVGeom::SCVF& scvf = geo->scvf(ip);

		// compute flux
		const number flux = VecDot(scvf.normal(), vel[ip]);

		// set all shapes to 0.0
		for (size_t sh = 0; sh < numSH; ++sh)
			conv_shape(ip, sh) = 0.0;

		// calculate exponentially weighed convection shapes
		size_t num_scvf_sh = scvf.num_sh();
		UG_COND_THROW(geo->num_scv_ips() < num_scvf_sh,
					  "Less SCVs present than SCVF shape functions.");
		std::vector<MathVector<worldDim> > dir(num_scvf_sh);
		std::vector<number> dirVelProd(num_scvf_sh);
		std::vector<number> dirVelProdExp(num_scvf_sh);
		number dirVelProdExpSum = 0.0;
		for (size_t sh = 0; sh < num_scvf_sh; ++sh)
		{
			VecScaleAdd(dir[sh], 1.0, scvf.global_ip(), -1.0, geo->scv_global_ips()[sh]);
			VecNormalize(dir[sh], dir[sh]);
			dirVelProd[sh] = VecDot(dir[sh], vel[ip]);
			dirVelProdExp[sh] = exp(m_expFact*dirVelProd[sh]);
			dirVelProdExpSum += dirVelProdExp[sh];
		}

		for (size_t sh = 0; sh < num_scvf_sh; ++sh)
			conv_shape(ip, sh) = dirVelProdExp[sh] / dirVelProdExpSum * flux;


		// compute derivatives (w.r.t. vel) if needed
		if (computeDeriv)
		{
			for (size_t sh = 0; sh < numSH; ++sh)
				VecSet(D_vel(ip, sh), 0.0);

			for (size_t sh = 0; sh < num_scvf_sh; ++sh)
			{
				for (size_t sh2 = 0; sh2 < num_scvf_sh; ++sh2)
				{
					MathVector<worldDim> dirDiff;
					VecScaleAdd(dirDiff, m_expFact*flux, dir[sh], -m_expFact*flux, dir[sh2]);
					VecScaleAdd(D_vel(ip, sh), 1.0, D_vel(ip, sh), dirVelProdExp[sh2], dirDiff, dirVelProdExp[sh2], scvf.normal());
				}
				VecScale(D_vel(ip, sh), D_vel(ip, sh), dirVelProdExp[sh] / (dirVelProdExpSum*dirVelProdExpSum));
			}
		}

/* DEBUG: write upwind directions to file
if (onlyConsiderMe == (void*) this)
{
	MathVector<worldDim> temp;
	MathVector<worldDim> dir = 0.0;

	for (size_t sh = 0; sh < num_scvf_sh; ++sh)
	{
		VecScaleAdd(temp, -1.0, scvf.global_ip(), 1.0, geo->scv_global_ips()[sh]);
		if (flux)
			VecScaleAdd(dir, 1.0, dir, conv_shape(ip, sh)/flux, temp);
	}
	try
	{
		for (size_t i = 0; i < worldDim; ++i)
			outFile << scvf.global_ip()[i] << ", ";
		for (size_t i = 0; i < worldDim; ++i)
			outFile << dir[i] << ", ";
		outFile << "0\n";
	}
	UG_CATCH_THROW("Output file" << ofnss.str() << "could not be written to.");
}
*/
	}

/* DEBUG: write upwind directions to file
if (onlyConsiderMe == (void*) this)
	outFile.close();
*/

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
