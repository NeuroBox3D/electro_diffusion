/*
 * charge_marking.h
 *
 *  Created on: 30.03.2016
 *      Author: mbreit
 */

#ifndef UG__PLUGINS__NERNST_PLANCK__CHARGE_MARKING_H_
#define UG__PLUGINS__NERNST_PLANCK__CHARGE_MARKING_H_

#include "lib_disc/function_spaces/error_elem_marking_strategy.h"
#include "lib_disc/function_spaces/approximation_space.h"
#include "common/util/smart_pointer.h"

#include "interface1d_fv.h"

namespace ug {
namespace nernst_planck {


/// mark elements adjacent to given charged surfaces
template <typename TDomain>
class ChargeMarking : public IElementMarkingStrategy<TDomain>
{

public:
	typedef IElementMarkingStrategy<TDomain> base_type;

	ChargeMarking(number tol, size_t max_level)
	: m_tol(tol), m_max_level(max_level) {};

	void set_tolerance(number tol) {m_tol = tol;}
	void set_max_level(size_t max_level) {m_max_level = max_level;}

	//void add_interface(SmartPtr<IInterface1D> intf) {m_vIntf.push_back(intf);}

	void add_surface(int surf_si, int adj_vol_si);
	void remove_surface(int surf_si, int adj_vol_si);

	void mark
	(
		typename base_type::elem_accessor_type& aaError,
		IRefiner& refiner,
		ConstSmartPtr<DoFDistribution> dd
	);

	void mark_without_error(SmartPtr<IRefiner> refiner, SmartPtr<ApproximationSpace<TDomain> > approx);

protected:
	void mark(IRefiner& refiner, ConstSmartPtr<DoFDistribution> dd);

protected:
	number m_tol;
	size_t m_max_level;
	//std::vector<SmartPtr<IInterface1D> > m_vIntf;

	/// vector holding pairs of charged surface and adjacent
	/// element subset indices
	std::vector<std::pair<int, int> > m_vpSurfaces;
};


} // namespace nernst_planck
} // namespace ug

#endif // UG__PLUGINS__NERNST_PLANCK__CHARGE_MARKING_H_
