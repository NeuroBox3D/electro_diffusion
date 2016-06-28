/*
 * charge_marking.cpp
 *
 *  Created on: 30.03.2016
 *      Author: mbreit
 */

#include "charge_marking.h"
#include "intf_refMarkAdjuster.h"

#include "lib_disc/domain.h"
#include "lib_grid/refinement/hanging_node_refiner_base.h"

#include <algorithm>

namespace ug {
namespace nernst_planck {


template <typename TDomain>
void ChargeMarking<TDomain>::add_surface(int surf_si, int adj_vol_si)
{
	m_vpSurfaces.push_back(std::make_pair(surf_si, adj_vol_si));
}

template <typename TDomain>
void ChargeMarking<TDomain>::remove_surface(int surf_si, int adj_vol_si)
{
	std::pair<int,int> pair = std::make_pair(surf_si, adj_vol_si);
	std::vector<std::pair<int,int> >::iterator match =
		std::find(m_vpSurfaces.begin(), m_vpSurfaces.end(), pair);

	if (match != m_vpSurfaces.end())
		m_vpSurfaces.erase(match);
}



template <typename TDomain>
void ChargeMarking<TDomain>::mark
(
	typename base_type::elem_accessor_type& aaError,
	IRefiner& refiner,
	ConstSmartPtr<DoFDistribution> dd
)
{
	number minElemErr;
	number maxElemErr;
	number errTotal;
	size_t numElem;

	ComputeMinMaxTotal(aaError, dd, minElemErr, maxElemErr, errTotal, numElem);
	if (errTotal <= m_tol || dd->multi_grid()->num_levels() > m_max_level)
		return;

	mark(refiner, dd);
}


template <typename TDomain>
void ChargeMarking<TDomain>::mark_without_error
(
	SmartPtr<IRefiner> refiner,
	SmartPtr<ApproximationSpace<TDomain> > approx
)
{
	// get surface dof distribution
	ConstSmartPtr<DoFDistribution> dd = approx->dof_distribution(GridLevel(), false);

	mark(*refiner, dd);
}



template <typename TDomain>
void ChargeMarking<TDomain>::mark(IRefiner& refiner, ConstSmartPtr<DoFDistribution> dd)
{
	typedef typename MultiGrid::traits<typename domain_traits<TDomain::dim>::element_type>::secure_container elem_list_type;
	typedef typename domain_traits<TDomain::dim>::side_type side_type;
	typedef typename MultiGrid::traits<side_type>::secure_container side_list_type;
	typedef typename DoFDistribution::traits<Vertex>::const_iterator iter_type;

	Grid* grid = refiner.grid();
	ConstSmartPtr<ISubsetHandler> sh = dd->subset_handler();

	size_t sz = m_vpSurfaces.size();
	for (size_t i = 0; i < sz; ++i)
	{
		int surf_si = m_vpSurfaces[i].first;
		int vol_si = m_vpSurfaces[i].second;
		iter_type iter = dd->template begin<Vertex>(surf_si);
		iter_type iterEnd = dd->template end<Vertex>(surf_si);

		//	loop elements for marking
		for (; iter != iterEnd; ++iter)
		{
			Vertex* v = *iter;
			elem_list_type el;
			grid->associated_elements(el, v);
			size_t el_sz = el.size();
			for (size_t j = 0; j < el_sz; ++j)
			{
				int elem_si = sh->get_subset_index(el[j]);
				if (elem_si == vol_si)
				{
					// mark the element for anisotropic refinement
					refiner.mark(el[j], RM_REFINE);
				}
			}
		}
	}
/*
	if (m_vIntf.size())
	{
		HangingNodeRefinerBase<MGSelector>* hnr = dynamic_cast<HangingNodeRefinerBase<MGSelector>*>(&refiner);
		if (hnr)
		{
			SmartPtr<InterfaceRefMarkAdjuster> intfRMA = make_sp(new InterfaceRefMarkAdjuster());
			intfRMA->set_subset_handler(sh);
			intfRMA->add_interfaces(m_vIntf);
			hnr->add_ref_mark_adjuster(intfRMA);
		}
		else
		{
			UG_LOGN("Warning: ChargeMarking is supposed to work with an instance of HangingNodeRefiner\n"
					"         and might not work properly with the given refiner type.");
		}
	}
*/
}



// explicit template specializations
#ifdef UG_DIM_1
	template class ChargeMarking<Domain1d>;
#endif
#ifdef UG_DIM_2
	template class ChargeMarking<Domain2d>;
#endif
#ifdef UG_DIM_3
	template class ChargeMarking<Domain3d>;
#endif



} // namespace nernst_planck
} // namespace ug