/*
 * pnp1d_fv.cpp
 *
 *  Created on: 13.08.2015
 *      Author: mbreit
 */

#include "pnp1d_fv.h"

namespace ug{
namespace nernst_planck{

//////////////////////////
//   Constructors   	//
//////////////////////////

// constructor with strings
template <typename TDomain>
PNP1D_FV<TDomain>::PNP1D_FV
(
	SmartPtr<ApproximationSpace<TDomain> > approx,
	const char* functions,
	const char* subsets
)
: IElemDisc<TDomain>(functions, subsets),
  m_spApproxSpace(approx),
  m_bNonRegularGrid(false),
  m_aRadius("radius"),
  m_nIon(0), _PHI_(0),
  m_R(8.31451), m_T(298.15), m_F(96485.0),
  m_reprDim(size_t(3))
{
	// throw: Poisson equation is not instationary!
	UG_THROW("Do not use the PNP1D classes. They are faulty and need repair.");

	// check dimensionality of subsets
	SubsetGroup ssg(m_spApproxSpace->subset_handler(), this->m_vSubset);
	int d = ssg.get_highest_subset_dimension();
	if (d != 1)
	{
		UG_THROW("Error in PNP1D_FV1: This elem disc is only supposed to work on 1D subsets.\n"
				 "Yet it has been initialized with a subset of dimension " << d << ".");
	}

	this->clear_add_fct();
}

// constructor with vectors
template <typename TDomain>
PNP1D_FV<TDomain>::PNP1D_FV
(
	SmartPtr<ApproximationSpace<TDomain> > approx,
	const std::vector<std::string>& vFct,
	const std::vector<std::string>& vSubset
)
: IElemDisc<TDomain>(vFct, vSubset),
  m_spApproxSpace(approx),
  m_bNonRegularGrid(false),
  m_aRadius("radius"),
  m_nIon(0), _PHI_(0),
  m_R(8.31451), m_T(298.15), m_F(96485.0),
  m_reprDim(size_t(3))
{
	// throw: Poisson equation is not instationary!
	UG_THROW("Do not use the PNP1D classes. They are faulty and need repair.");

	// check dimensionality of subsets
	SubsetGroup ssg(m_spApproxSpace->subset_handler(), this->m_vSubset);
	int d = ssg.get_highest_subset_dimension();
	if (d != 1)
	{
		UG_THROW("Error in PNP1D_FV1: This elem disc is only supposed to work on 1D subsets.\n"
				"Yet it has been initialized with a subset of dimension " << d << ".");
	}

	this->clear_add_fct();
}


//////////////////////////////////////////
//   Setter functions for parameters   	//
//////////////////////////////////////////

// set the names of the ion species present
template <typename TDomain>
void PNP1D_FV<TDomain>::set_ions(const std::vector<std::string>& vIons)
{
	 m_nIon = vIons.size();
	 _PHI_ = m_nIon;

	 m_vIonNames = vIons;
}

// set ion valencies
template <typename TDomain>
void PNP1D_FV<TDomain>::set_valencies(const std::vector<int>& vValencies)
{
	if (vValencies.size() != m_nIon)
	{
		UG_THROW("Number of valency values given in set_valencies() does not match the number of ions.\n"
				 "Maybe you forgot to define the ions present by a call to set_ions()?");
	}

	m_vValency = vValencies;
}

// set ion reversal potentials
template <typename TDomain>
void PNP1D_FV<TDomain>::set_reversal_potentials(const std::vector<number>& vRevPot)
{
	if (vRevPot.size() != m_nIon)
	{
		UG_THROW("Number of valency values given in set_reversal_potentials() does not match the number of ions.\n"
				 "Maybe you forgot to define the ions present by a call to set_ions()?");
	}

	m_vRevPot = vRevPot;
}

// set ion specific conductances of the plasma membrane
template <typename TDomain>
void PNP1D_FV<TDomain>::set_specific_conductances(const std::vector<number>& vSpecCond)
{
	if (vSpecCond.size() != m_nIon)
	{
		UG_THROW("Number of valency values given in set_specific_conductances() does not match the number of ions.\n"
				 "Maybe you forgot to define the ions present by a call to set_ions()?");
	}

	m_vSpecConductance = vSpecCond;
}

// set ion specific capacities for influency effects at the plasma membrane
template <typename TDomain>
void PNP1D_FV<TDomain>::set_specific_capacities(const std::vector<number>& vSpecCap)
{
	if (vSpecCap.size() != m_nIon)
	{
		UG_THROW("Number of valency values given in set_specific_capacities() does not match the number of ions.\n"
				 "Maybe you forgot to define the ions present by a call to set_ions()?");
	}

	m_vSpecCapacity = vSpecCap;
}

// set ion diffusion constants
template <typename TDomain>
void PNP1D_FV<TDomain>::set_diffusion_constants(const std::vector<number>& vDiff)
{
	if (vDiff.size() != m_nIon)
	{
		UG_THROW("Number of valency values given in set_diffusion_constants() does not match the number of ions.\n"
				 "Maybe you forgot to define the ions present by a call to set_ions()?");
	}

	m_vDiffusionTensor.resize(m_nIon);
	for (size_t i = 0; i < m_nIon; i++)
	{
		// make sure, all entries are 0.0
		m_vDiffusionTensor[i] *= 0.0;

		// create isotrope diagonal matrix
		for (int j = 0; j < worldDim; j++)
			m_vDiffusionTensor[i](j,j) = vDiff[i];
	}
}

// set effective electric permettivity in the dendrite (eps_0 * eps_r)
template <typename TDomain>
void PNP1D_FV<TDomain>::set_permettivities(const number eps_dend, const number eps_mem)
{
	// make sure all entries are 0.0
	m_permittivity_dend = 0.0;

	// create isotrope diagonal matrix
	for (int j = 0; j < worldDim; j++)
	{
		m_permittivity_dend(j,j) = eps_dend;
	}

	m_permittivity_mem = eps_mem;
}

// set constant dendritic radius
template <typename TDomain>
void PNP1D_FV<TDomain>::set_membrane_thickness(const number d)
{
	m_mem_thickness = d;
}

// set constant dendritic radius
template <typename TDomain>
void PNP1D_FV<TDomain>::set_dendritic_radius(const number r)
{
	// handle the attachment
	if (m_spApproxSpace->domain()->grid()->has_vertex_attachment(m_aRadius))
		UG_THROW("Radius attachment necessary for PNP1D_FV elem disc "
				 "could not be made, since it already exists.");
	m_spApproxSpace->domain()->grid()->attach_to_vertices_dv(m_aRadius, r);

	m_aaRadius = Grid::AttachmentAccessor<Vertex, ANumber>(*m_spApproxSpace->domain()->grid(), m_aRadius);
}

// set constant dendritic radius
template <typename TDomain>
void PNP1D_FV<TDomain>::set_rtf(number R, number T, number F)
{
	m_R = R;
	m_T = T;
	m_F = F;
}

template <typename TDomain>
void PNP1D_FV<TDomain>::set_represented_dimension(size_t dim)
{
	switch (dim)
	{
		case 2: m_reprDim = 2; break;
		case 3: m_reprDim = 3; break;
		default: UG_THROW("Can only represent either dimension 2 or 3.");
	}
}


//////////////////////////////////////////////
//   Assembling functions from IElemDisc   	//
//////////////////////////////////////////////

template<typename TDomain>
void PNP1D_FV<TDomain>::
prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid)
{
	//	check grid
	if (bNonRegularGrid)
		UG_THROW("PNP1D_FV: Only regular grid implemented.");

	// check number
	if (vLfeID.size() < 1)
		UG_THROW("PNP1D_FV: Wrong number of functions given. "
				 "This discretization needs at least 1.");

	for (std::vector<LFEID>::size_type i = 0; i < vLfeID.size(); ++i)
	{
		// check that Lagrange
		if (vLfeID[i].type() != LFEID::LAGRANGE)
			UG_THROW("PNP1D_FV: Only Lagrange type supported.");

		// check that not adaptive
		if (vLfeID[i].order() < 1)
			UG_THROW("PNP1D_FV: Adaptive order not implemented.");

		// check that orders are all the same
		if (i > 0 && vLfeID[i].order() != vLfeID[i-1].order())
			UG_THROW("PNP1D_FV: FV order needs to be the same for all unknowns involved.");

		// check that dims are all the same
		if (i > 0 && vLfeID[i].dim() != vLfeID[i-1].dim())
			UG_THROW("PNP1D_FV: FV dim needs to be the same for all unknowns involved.");
	}

	m_lfeID = vLfeID[0];
	m_quadOrder = m_lfeID.order()+1;

	// update assemble functions
	register_all_funcs();
}

template<typename TDomain>
bool PNP1D_FV<TDomain>::use_hanging() const
{
	return true;
}


//	virtual prepares the loop over all elements of one type
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void PNP1D_FV<TDomain>::prep_elem_loop(const ReferenceObjectID roid, const int si)
{
	//	request geometry
	TFVGeom& geo = GeomProvider<TFVGeom>::get(m_lfeID, m_quadOrder);
	try	{geo.update_local(roid, m_lfeID, m_quadOrder);}
	UG_CATCH_THROW("PNP1D_FV::prep_elem_loop: Cannot update local finite volume geometry.");
}


//	prepares one element for assembling
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void PNP1D_FV<TDomain>::prep_elem(const LocalVector& u, GridObject* elem, const ReferenceObjectID roid, const MathVector<worldDim> vCornerCoords[])
{
	// update geometry for this element
	TFVGeom& geo = GeomProvider<TFVGeom>::get(m_lfeID, m_quadOrder);

	try {geo.update(elem, vCornerCoords, &(this->subset_handler()));}
	UG_CATCH_THROW("PNP1D_FV::prep_elem: Cannot update finite volume geometry.");
}


//	postprocesses the loop over all elements of one type
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void PNP1D_FV<TDomain>::fsh_elem_loop()
{}


// virtual Assembling of Defect (Stiffness part)
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void PNP1D_FV<TDomain>::add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<worldDim> vCornerCoords[])
{
	// get finite volume geometry
	static const TFVGeom& fvgeom = GeomProvider<TFVGeom>::get(m_lfeID, m_quadOrder);

	// cast elem to appropriate type
	TElem* pElem = dynamic_cast<TElem*>(elem);
	if (!pElem) {UG_THROW("Wrong element type.");}

	// radii at elem end points
	number r0 = m_aaRadius[pElem->vertex(0)];
	number r1 = m_aaRadius[pElem->vertex(1)];

// divergence terms (flux through the dendrite)
	// loop SCVF
	for (size_t s = 0; s < fvgeom.num_scvf(); s++)
	{
		// get current SCVF
		const typename TFVGeom::SCVF& scvf = fvgeom.scvf(s);

		// get radius at IP (linear interpolation) - there is only one IP!
		number loc_coord = scvf.local_ip(0)[0];
		number radiusAtIP = r0 + loc_coord*(r1-r0);

		// scvf area
		number scvfArea = (m_reprDim == 3 ? PI * radiusAtIP * radiusAtIP : 2*radiusAtIP);

	// potential
		// to compute eps \nabla phi
		MathVector<worldDim> eps_grad_phi, grad_phi;

		// compute gradient and shape at ip
		VecSet(grad_phi, 0.0);
		for (size_t sh = 0; sh < scvf.num_sh(); sh++)
			VecScaleAppend(grad_phi, u(_PHI_,sh), scvf.global_grad(0,sh));

		// scale by permittivity tensor
		MatVecMult(eps_grad_phi, m_permittivity_dend, grad_phi);

		// compute flux density
		number flux = VecDot(eps_grad_phi, scvf.normal());

		// scale with hidden area of dendritic crosssection
		flux *= scvfArea;

		// add to local defect
		d(_PHI_, scvf.from()) -= flux;
		d(_PHI_, scvf.to()  ) += flux;

	// ion concentrations
		for (size_t i = 0; i < m_nIon; i++)
		{
		// diffusion
			// to compute D \nabla c_i
			MathVector<worldDim> D_grad, grad;

			// compute gradient and shape at ip
			VecSet(grad, 0.0);
			for (size_t sh = 0; sh < scvf.num_sh(); sh++)
				VecScaleAppend(grad, u(i,sh), scvf.global_grad(0,sh));

			// scale by diffusion tensor
			MatVecMult(D_grad, m_vDiffusionTensor[i], grad);

			// compute flux density
			number diff_flux = VecDot(D_grad, scvf.normal());

			// scale with hidden area of dendritic cross-section
			diff_flux *= scvfArea;

			// add to local defect
			d(i, scvf.from()) -= diff_flux;
			d(i, scvf.to()  ) += diff_flux;

		// electricity
			// get concentration at ip
			number c_i = 0.0;
			for (size_t sh = 0; sh < scvf.num_sh(); sh++)
				c_i += u(i,sh) * scvf.shape(0,sh);

			// scale electric field by diffusion tensor
			MathVector<worldDim> D_grad_phi;
			MatVecMult(D_grad_phi, m_vDiffusionTensor[i], grad_phi);

			// compute inner product
			number el_flux = VecDot(D_grad_phi, scvf.normal());

			// scale by concentration and constants
			el_flux *= scvfArea * m_vValency[i] * m_F/(m_R*m_T) * c_i;

			// add to local defect
			d(i, scvf.from()) -= el_flux;
			d(i, scvf.to()  ) += el_flux;
		}
	}

// volume parts
	// loop SCV
	for (size_t s = 0; s < fvgeom.num_scv(); ++s)
	{
		// get current SCV
		const typename TFVGeom::SCV& scv = fvgeom.scv(s);

		// get associated node
		const int co = scv.node_id();

		// calculate integral (no mass lumping here)
		for (size_t ip = 0; ip < scv.num_ip(); ++ip)
		{
			// get radius at IP (linear interpolation)
			number loc_coord = scv.local_ip(ip)[0];
			number radiusAtIP = r0 + loc_coord*(r1-r0);

			// area/volume factor for dendritic scv (TODO: probably not completely correct)
			number scvAreaWeight = scv.weight(ip) * (m_reprDim == 3 ? 2*PI * radiusAtIP : 2);
			number scvVolWeight = scv.weight(ip) * (m_reprDim == 3 ? PI * radiusAtIP * radiusAtIP : 2*radiusAtIP);

			number phiAtIP = 0.0;
			for (size_t sh = 0; sh < scv.num_sh(); ++sh)
				phiAtIP += u(_PHI_, sh) * scv.shape(ip, sh);

			// potential equation: charge density term
			for (size_t i = 0; i < m_nIon; ++i)
			{
				// solution at ip
				number ionAtIP = 0.0;
				for (size_t sh = 0; sh < scv.num_sh(); ++sh)
					ionAtIP += u(i, sh) * scv.shape(ip, sh);

				d(_PHI_, co) += ionAtIP * m_vValency[i] * m_F * scvVolWeight;

				// ion concentrations: passive channel flux
				number flux = m_vSpecConductance[i] * (phiAtIP - m_vRevPot[i]);
				d(i, co) += flux / (m_vValency[i] * m_F) * scvAreaWeight;
			}
		}
	}
}

// virtual Assembling of Defect (Mass part)
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void PNP1D_FV<TDomain>::add_def_M_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<worldDim> vCornerCoords[])
{
	// get finite volume geometry
	static const TFVGeom& fvgeom = GeomProvider<TFVGeom>::get(m_lfeID, m_quadOrder);

	// cast elem to appropriate type
	TElem* pElem = dynamic_cast<TElem*>(elem);
	if (!pElem) {UG_THROW("Wrong element type.");}

	// radii at elem end points
	number r0 = m_aaRadius[pElem->vertex(0)];
	number r1 = m_aaRadius[pElem->vertex(1)];

// defects wrt time (not using lumping)
	// loop SCV
	for (size_t s = 0; s < fvgeom.num_scv(); ++s)
	{
		// get current SCV
		const typename TFVGeom::SCV& scv = fvgeom.scv(s);

		// get associated node
		const int co = scv.node_id();

		// calculate integral (no mass lumping here)
		for (size_t ip = 0; ip < scv.num_ip(); ++ip)
		{
			// get radius at IP (linear interpolation)
			number loc_coord = scv.local_ip(ip)[0];
			number radiusAtIP = r0 + loc_coord*(r1-r0);

			// area/volume factor for dendritic scv (TODO: probably not completely correct)
			number scvAreaWeight = scv.weight(ip) * (m_reprDim == 3 ? 2*PI * radiusAtIP : 2);
			number scvVolWeight = scv.weight(ip) * (m_reprDim == 3 ? PI * radiusAtIP * radiusAtIP : 2*radiusAtIP);

			number phiAtIP = 0.0;
			for (size_t sh = 0; sh < scv.num_sh(); ++sh)
				phiAtIP += u(_PHI_, sh) * scv.shape(ip, sh);

			for (size_t i = 0; i < m_nIon; ++i)
			{
				// potential in ion equations
				d(i, co) += phiAtIP * m_vSpecCapacity[i] / (m_vValency[i] * m_F) * scvAreaWeight;

				// ion concentrations
				number ionAtIP = 0.0;
				for (size_t sh = 0; sh < scv.num_sh(); ++sh)
					ionAtIP += u(i, sh) * scv.shape(ip, sh);
				d(i, co) += ionAtIP * scvVolWeight;
			}
		}
	}
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void PNP1D_FV<TDomain>::add_rhs_elem(LocalVector& rhs, GridObject* elem, const MathVector<worldDim> vCornerCoords[])
{}


// Assembling of Jacobian (Stiffness part)
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void PNP1D_FV<TDomain>::add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<worldDim> vCornerCoords[])
{
	// get finite volume geometry
	static const TFVGeom& fvgeom = GeomProvider<TFVGeom>::get(m_lfeID, m_quadOrder);

	// cast elem to appropriate type
	TElem* pElem = dynamic_cast<TElem*>(elem);
	if (!pElem) {UG_THROW("Wrong element type.");}

// divergence terms (flux through the dendrite)
	// loop SCVF
	number r0 = m_aaRadius[pElem->vertex(0)];
	number r1 = m_aaRadius[pElem->vertex(1)];
	for (size_t s = 0; s < fvgeom.num_scvf(); ++s)
	{
		// get current SCVF
		const typename TFVGeom::SCVF& scvf = fvgeom.scvf(s);

		// get radius at IP (linear interpolation)
		number loc_coord = scvf.local_ip(0)[0];
		number radiusAtIP = r0 + loc_coord*(r1-r0);

		// scvf area
		number scvfArea = (m_reprDim == 3 ? PI * radiusAtIP * radiusAtIP : 2*radiusAtIP);

		// compute grad(phi)
		MathVector<worldDim> D_grad_phi, grad_phi;
		VecSet(grad_phi, 0.0);
		for (size_t sh = 0; sh < scvf.num_sh(); sh++)
			VecScaleAppend(grad_phi, u(_PHI_,sh), scvf.global_grad(0, sh));

	// potential
		// to compute eps \nabla phi
		MathVector<worldDim> eps_grad;

		// compute gradient and shape at ip
		for (size_t sh = 0; sh < scvf.num_sh(); sh++)
		{
			// scale by permittivity tensor
			MatVecMult(eps_grad, m_permittivity_dend, scvf.global_grad(0, sh));

			// compute flux density
			number flux = VecDot(eps_grad, scvf.normal());

			// scale with hidden area of dendritic crosssection
			flux *= scvfArea;

			// add to local Jacobian
			J(_PHI_, scvf.from(), _PHI_, sh) -= flux;
			J(_PHI_, scvf.to(),   _PHI_, sh) += flux;
		}

	// ion concentrations
		for (size_t i = 0; i < m_nIon; i++)
		{
		// diffusion
			// to compute D \nabla c_i
			MathVector<worldDim> D_grad;

			// compute gradient and shape at ip
			for (size_t sh = 0; sh < scvf.num_sh(); sh++)
			{
				// scale by diffusion tensor
				MatVecMult(D_grad, m_vDiffusionTensor[i], scvf.global_grad(0, sh));

				// compute flux density
				number diff_flux = VecDot(D_grad, scvf.normal());

				// scale with hidden area of dendritic crosssection
				diff_flux *= scvfArea;

				// add to local Jacobian
				J(i, scvf.from(), i, sh) -= diff_flux;
				J(i, scvf.to(),   i, sh) += diff_flux;
			}

		// electricity
			// get concentration at ip
			number c_i = 0.0;
			for (size_t sh = 0; sh < scvf.num_sh(); sh++)
				c_i += u(i,sh) * scvf.shape(0, sh);

			for (size_t sh = 0; sh < scvf.num_sh(); sh++)
			{
			// deriv wrt concs
				// scale electric field by diffusion tensor
				MatVecMult(D_grad_phi, m_vDiffusionTensor[i], grad_phi);

				// compute inner product
				number el_flux = VecDot(D_grad_phi, scvf.normal());

				// scale by concentration and constants
				el_flux *= scvfArea * m_vValency[i] * m_F/(m_R*m_T) * scvf.shape(0, sh);

				// add to local Jacobian
				J(i, scvf.from(), i, sh) -= el_flux;
				J(i, scvf.to(),   i, sh) += el_flux;

			// deriv wrt potential
				// scale electric field by diffusion tensor
				MatVecMult(D_grad, m_vDiffusionTensor[i], scvf.global_grad(0, sh));

				// compute inner product
				el_flux = VecDot(D_grad, scvf.normal());

				// scale by concentration and constants
				el_flux *= scvfArea * m_vValency[i] * m_F/(m_R*m_T) * c_i;

				// add to local defect
				J(i, scvf.from(), _PHI_, sh) -= el_flux;
				J(i, scvf.to(),   _PHI_, sh) += el_flux;
			}
		}
	}

// volume parts (using lumping)
	// loop SCV
	for (size_t s = 0; s < fvgeom.num_scv(); ++s)
	{
		// get current SCV
		const typename TFVGeom::SCV& scv = fvgeom.scv(s);

		// get associated node
		const int co = scv.node_id();

		// calculate integral (no mass lumping here)
		for (size_t ip = 0; ip < scv.num_ip(); ++ip)
		{
			// get radius at IP (linear interpolation)
			number loc_coord = scv.local_ip(ip)[0];
			number radiusAtIP = r0 + loc_coord*(r1-r0);

			// area/volume factor for dendritic scv (TODO: probably not completely correct)
			number scvAreaWeight = scv.weight(ip) * (m_reprDim == 3 ? 2*PI * radiusAtIP : 2);
			number scvVolWeight = scv.weight(ip) * (m_reprDim == 3 ? PI * radiusAtIP * radiusAtIP : 2*radiusAtIP);

			for (size_t i = 0; i < m_nIon; ++i)
			{
				for (size_t sh = 0; sh < scv.num_sh(); ++sh)
				{
					// potential equation: charge density term
					J(_PHI_, co, i, sh) += scv.shape(ip, sh) * m_vValency[i] * m_F * scvVolWeight;

					// ion concentrations: passive channel flux
					number flux = m_vSpecConductance[i] * scv.shape(ip, sh);
					J(i, co, _PHI_, sh) += flux / (m_vValency[i] * m_F) * scvAreaWeight;
				}
			}
		}
	}
}

// Assembling of Jacobian (Mass part)
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void PNP1D_FV<TDomain>::add_jac_M_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<worldDim> vCornerCoords[])
{
	// get finite volume geometry
	static const TFVGeom& fvgeom = GeomProvider<TFVGeom>::get(m_lfeID, m_quadOrder);

	// cast elem to appropriate type
	TElem* pElem = dynamic_cast<TElem*>(elem);
	if (!pElem) {UG_THROW("Wrong element type.");}

	// radii at elem end points
	number r0 = m_aaRadius[pElem->vertex(0)];
	number r1 = m_aaRadius[pElem->vertex(1)];

// derivs wrt time (not using lumping)
	// loop SCV
	for (size_t s = 0; s < fvgeom.num_scv(); s++)
	{
		// get current SCV
		const typename TFVGeom::SCV& scv = fvgeom.scv(s);

		// get associated node
		const int co = scv.node_id();

		// calculate integral (no mass lumping here)
		for (size_t ip = 0; ip < scv.num_ip(); ++ip)
		{
			// get radius at IP (linear interpolation)
			number loc_coord = scv.local_ip(ip)[0];
			number radiusAtIP = r0 + loc_coord*(r1-r0);

			// area/volume factor for dendritic scv (TODO: probably not completely correct)
			number scvAreaWeight = scv.weight(ip) * (m_reprDim == 3 ? 2*PI * radiusAtIP : 2);
			number scvVolWeight = scv.weight(ip) * (m_reprDim == 3 ? PI * radiusAtIP * radiusAtIP : 2*radiusAtIP);

			for (size_t i = 0; i < m_nIon; ++i)
			{
				for (size_t sh = 0; sh < scv.num_sh(); ++sh)
				{
					// potential in ion equations
					J(i, co, _PHI_, sh) += scv.shape(ip, sh) * m_vSpecCapacity[i] / (m_vValency[i] * m_F) * scvAreaWeight;

					// ion concentrations
					J(i, co, i, sh) += scv.shape(ip, sh) * scvVolWeight;
				}
			}
		}
	}
}

//////////////////////////////////////////////
//   Error estimation assembling methods	//
//////////////////////////////////////////////

//	virtual prepares the loop over all elements of one type for the computation of the error estimator
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void PNP1D_FV<TDomain>::prep_err_est_elem_loop(const ReferenceObjectID roid, const int si)
{
	UG_THROW("PNP1D_FV::prep_err_est_elem_loop not yet implemented.");
}

//	virtual prepares the loop over all elements of one type for the computation of the error estimator
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void PNP1D_FV<TDomain>::prep_err_est_elem(const LocalVector& u, GridObject* elem, const MathVector<worldDim> vCornerCoords[])
{
	UG_THROW("PNP1D_FV::prep_err_est_elem not yet implemented.");
}

//	virtual compute the error estimator (stiffness part) contribution for one element
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void PNP1D_FV<TDomain>::compute_err_est_A_elem(const LocalVector& u, GridObject* elem, const MathVector<worldDim> vCornerCoords[], const number& scale)
{
	UG_THROW("PNP1D_FV::compute_err_est_A_elem not yet implemented.");
}

//	virtual compute the error estimator (mass part) contribution for one element
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void PNP1D_FV<TDomain>::compute_err_est_M_elem(const LocalVector& u, GridObject* elem, const MathVector<worldDim> vCornerCoords[], const number& scale)
{
	UG_THROW("PNP1D_FV::compute_err_est_M_elem not yet implemented.");
}

//	virtual postprocesses the loop over all elements of one type in the computation of the error estimator
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void PNP1D_FV<TDomain>::fsh_err_est_elem_loop()
{
	UG_THROW("PNP1D_FV::fsh_err_est_elem_loop not yet implemented.");
}


//////////////////////////////////////////
//   registering assembling functions   //
//////////////////////////////////////////

template<typename TDomain>
void PNP1D_FV<TDomain>::
register_all_funcs()
{
//	switch assemble functions
	register_func<RegularEdge, DimFVGeometry<worldDim, 1> >();
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void PNP1D_FV<TDomain>::
register_func()
{
	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;

	this->clear_add_fct(id);
	this->set_prep_elem_loop_fct(	id, &this_type::template prep_elem_loop<TElem, TFVGeom>);
	this->set_prep_elem_fct(	 	id, &this_type::template prep_elem<TElem, TFVGeom>);
	this->set_fsh_elem_loop_fct( 	id, &this_type::template fsh_elem_loop<TElem, TFVGeom>);
	this->set_add_jac_A_elem_fct(	id, &this_type::template add_jac_A_elem<TElem, TFVGeom>);
	this->set_add_jac_M_elem_fct(	id, &this_type::template add_jac_M_elem<TElem, TFVGeom>);
	this->set_add_def_A_elem_fct(	id, &this_type::template add_def_A_elem<TElem, TFVGeom>);
	this->set_add_def_M_elem_fct(	id, &this_type::template add_def_M_elem<TElem, TFVGeom>);
	this->set_add_rhs_elem_fct(	 	id, &this_type::template add_rhs_elem<TElem, TFVGeom>);

	// error estimator parts
	this->set_prep_err_est_elem_loop(	id, &this_type::template prep_err_est_elem_loop<TElem, TFVGeom>);
	this->set_prep_err_est_elem(		id, &this_type::template prep_err_est_elem<TElem, TFVGeom>);
	this->set_compute_err_est_A_elem(	id, &this_type::template compute_err_est_A_elem<TElem, TFVGeom>);
	this->set_compute_err_est_M_elem(	id, &this_type::template compute_err_est_M_elem<TElem, TFVGeom>);
	this->set_fsh_err_est_elem_loop(	id, &this_type::template fsh_err_est_elem_loop<TElem, TFVGeom>);
}


//////////////////////////////////
//   explicit specializations   //
//////////////////////////////////

#ifdef UG_DIM_1
	template class PNP1D_FV<Domain1d>;
#endif
#ifdef UG_DIM_2
	template class PNP1D_FV<Domain2d>;
#endif
#ifdef UG_DIM_3
	template class PNP1D_FV<Domain3d>;
#endif

} // namespace nernst_planck
} // namespace ug
