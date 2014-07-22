/*
 * PNP_1D.cpp
 *
 *  Created on: 08.07.2014
 *      Author: mbreit
 */

#include "PNP_1D.h"

namespace ug{
namespace nernst_planck{

//////////////////////////
//   Constructors   	//
//////////////////////////

// constructor with strings
template <typename TDomain>
PNP_1D<TDomain>::PNP_1D
(
	SmartPtr<ApproximationSpace<TDomain> > approx,
	const char* functions,
	const char* subsets
)
: IElemDisc<TDomain>(functions, subsets),
  R(8.31451), T(298.15), F(96485.0),
  m_spApproxSpace(approx),
  m_bNonRegularGrid(false),
  m_aRadius("radius"),
  m_nIon(0), _PHI_(0)
{
	SubsetGroup ssg(m_spApproxSpace->subset_handler(), this->m_vSubset);
	int d = ssg.get_highest_subset_dimension();
	if (d != 1)
	{
		UG_THROW("Error in PNP_1D: This elem disc is only supposed to work on 1D subsets.\n"
				 "Yet it has been initialized with a subset of dimension " << d << ".");
	}

	register_all_funcs();
}

// constructor with vectors
template <typename TDomain>
PNP_1D<TDomain>::PNP_1D
(
	SmartPtr<ApproximationSpace<TDomain> > approx,
	const std::vector<std::string>& vFct,
	const std::vector<std::string>& vSubset
)
: IElemDisc<TDomain>(vFct, vSubset),
  R(8.314), T(310.0), F(96485.0),
  m_spApproxSpace(approx),
  m_bNonRegularGrid(false),
  m_aRadius("radius"),
  m_nIon(0), _PHI_(0)
{
	SubsetGroup ssg(m_spApproxSpace->subset_handler(), this->m_vSubset);
	int d = ssg.get_highest_subset_dimension();
	if (d != 1)
	{
		UG_THROW("Error in PNP_1D: This elem disc is only supposed to work on 1D subsets.\n"
				"Yet it has been initialized with a subset of dimension " << d << ".");
	}

	register_all_funcs();
}


//////////////////////////////////////////
//   Setter functions for parameters   	//
//////////////////////////////////////////

// set the names of the ion species present
template <typename TDomain>
void PNP_1D<TDomain>::set_ions(const std::vector<std::string>& vIons)
{
	 m_nIon = vIons.size();
	 _PHI_ = m_nIon;

	 m_vIonNames = vIons;
}

// set ion valencies
template <typename TDomain>
void PNP_1D<TDomain>::set_valencies(const std::vector<int>& vValencies)
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
void PNP_1D<TDomain>::set_reversal_potentials(const std::vector<number>& vRevPot)
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
void PNP_1D<TDomain>::set_specific_conductances(const std::vector<number>& vSpecCond)
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
void PNP_1D<TDomain>::set_specific_capacities(const std::vector<number>& vSpecCap)
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
void PNP_1D<TDomain>::set_diffusion_constants(const std::vector<number>& vDiff)
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
		for (size_t j = 0; j < worldDim; j++)
			m_vDiffusionTensor[i](j,j) = vDiff[i];
	}
}

// set effective electric permettivity in the dendrite (eps_0 * eps_r)
template <typename TDomain>
void PNP_1D<TDomain>::set_permettivities(const number eps_dend, const number eps_mem)
{
	// make sure, all entries are 0.0
	m_permittivity_dend *= 0.0;

	// create isotrope diagonal matrix
	for (size_t j = 0; j < worldDim; j++)
	{
		m_permittivity_dend(j,j) = eps_dend;
	}

	m_permittivity_mem = eps_mem;
}

// set constant dendritic radius
template <typename TDomain>
void PNP_1D<TDomain>::set_membrane_thickness(const number d)
{
	m_mem_thickness = d;
}

// set constant dendritic radius
template <typename TDomain>
void PNP_1D<TDomain>::set_dendritic_radius(const number r)
{
	// handle the attachment
	if (m_spApproxSpace->domain()->grid()->has_vertex_attachment(m_aRadius))
		UG_THROW("Radius attachment necessary for PNP_1D elem disc "
				 "could not be made, since it already exists.");
	m_spApproxSpace->domain()->grid()->attach_to_vertices_dv(m_aRadius, r);

	m_aaRadius = Grid::AttachmentAccessor<Vertex, ANumber>(*m_spApproxSpace->domain()->grid(), m_aRadius);
}


//////////////////////////////////////////////
//   Assembling functions from IElemDisc   	//
//////////////////////////////////////////////


template<typename TDomain>
void PNP_1D<TDomain>::prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid)
{
	// check that Lagrange 1st order
	for (std::size_t i = 0; i < vLfeID.size(); i++)
		if (vLfeID[i].type() != LFEID::LAGRANGE || vLfeID[i].order() != 1)
			{UG_THROW("PNP_1D: 1st order Lagrange expected.");}

	// update assemble functions
	m_bNonRegularGrid = bNonRegularGrid;
	register_all_funcs();
}

template<typename TDomain>
bool PNP_1D<TDomain>::use_hanging() const
{
	return true;
}



//	virtual prepares the loop over all elements of one type
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void PNP_1D<TDomain>::prep_elem_loop(const ReferenceObjectID roid, const int si)
{
	/*
	// set local positions
	if (!TFVGeom::usesHangingNodes)
	{

		static const int refDim = TElem::dim;
		TFVGeom& geo = GeomProvider<TFVGeom>::get();
		const MathVector<refDim>* vSCVFip = geo.scvf_local_ips();
		const size_t numSCVFip = geo.num_scvf_ips();
		const MathVector<refDim>* vSCVip = geo.scv_local_ips();
		const size_t numSCVip = geo.num_scv_ips();
		m_vimDiffusion.template 		set_local_ips<refDim>(vSCVFip,numSCVFip, false);
		m_imVelocity.template 		set_local_ips<refDim>(vSCVFip,numSCVFip, false);
		m_imFlux.template 			set_local_ips<refDim>(vSCVFip,numSCVFip, false);
		m_imSource.template 		set_local_ips<refDim>(vSCVip,numSCVip, false);
		m_imVectorSource.template 	set_local_ips<refDim>(vSCVFip,numSCVFip, false);
		m_imReactionRate.template 	set_local_ips<refDim>(vSCVip,numSCVip, false);
		m_imReaction.template 		set_local_ips<refDim>(vSCVip,numSCVip, false);
		m_imReactionRateExpl.template set_local_ips<refDim>(vSCVip,numSCVip, false);
		m_imReactionExpl.template 	set_local_ips<refDim>(vSCVip,numSCVip, false);
		m_imSourceExpl.template 	set_local_ips<refDim>(vSCVip,numSCVip, false);
		m_imMassScale.template 		set_local_ips<refDim>(vSCVip,numSCVip, false);
		m_imMass.template 			set_local_ips<refDim>(vSCVip,numSCVip, false);
	}
	*/
}

//	virtual prepare one elements for assembling
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void PNP_1D<TDomain>::prep_elem(const LocalVector& u, GridObject* elem, const MathVector<worldDim> vCornerCoords[])
{
	// update geometry for this element
	static TFVGeom& geo = GeomProvider<TFVGeom>::get();

	try {geo.update(elem, vCornerCoords, &(this->subset_handler()));}
	UG_CATCH_THROW("PNP_1D::prep_elem: Cannot update finite volume geometry.");

	/*
	// set local positions
	if (TFVGeom::usesHangingNodes)
	{
		const int refDim = TElem::dim;
		const MathVector<refDim>* vSCVFip = geo.scvf_local_ips();
		const size_t numSCVFip = geo.num_scvf_ips();
		const MathVector<refDim>* vSCVip = geo.scv_local_ips();
		const size_t numSCVip = geo.num_scv_ips();
		m_vimDiffusion.template 		set_local_ips<refDim>(vSCVFip,numSCVFip);
		m_imVelocity.template 		set_local_ips<refDim>(vSCVFip,numSCVFip);
		m_imFlux.template 			set_local_ips<refDim>(vSCVFip,numSCVFip);
		m_imSource.template 		set_local_ips<refDim>(vSCVip,numSCVip);
		m_imVectorSource.template 	set_local_ips<refDim>(vSCVFip,numSCVFip);
		m_imReactionRate.template 	set_local_ips<refDim>(vSCVip,numSCVip);
		m_imReaction.template 		set_local_ips<refDim>(vSCVip,numSCVip);
		m_imReactionRateExpl.template 	set_local_ips<refDim>(vSCVip,numSCVip);
		m_imReactionExpl.template 	set_local_ips<refDim>(vSCVip,numSCVip);
		m_imSourceExpl.template		set_local_ips<refDim>(vSCVip,numSCVip);
		m_imMassScale.template 		set_local_ips<refDim>(vSCVip,numSCVip);
		m_imMass.template 			set_local_ips<refDim>(vSCVip,numSCVip);
	}

	//	set global positions
	const MathVector<dim>* vSCVFip = geo.scvf_global_ips();
	const size_t numSCVFip = geo.num_scvf_ips();
	const MathVector<dim>* vSCVip = geo.scv_global_ips();
	const size_t numSCVip = geo.num_scv_ips();
	m_vimDiffusion.			set_global_ips(vSCVFip, numSCVFip);
	m_imVelocity.			set_global_ips(vSCVFip, numSCVFip);
	m_imFlux.				set_global_ips(vSCVFip, numSCVFip);
	m_imSource.				set_global_ips(vSCVip, numSCVip);
	m_imVectorSource.		set_global_ips(vSCVFip, numSCVFip);
	m_imReactionRate.		set_global_ips(vSCVip, numSCVip);
	m_imReactionRateExpl.	set_global_ips(vSCVip, numSCVip);
	m_imReactionExpl.		set_global_ips(vSCVip, numSCVip);
	m_imSourceExpl.			set_global_ips(vSCVip, numSCVip);
	m_imReaction.			set_global_ips(vSCVip, numSCVip);
	m_imMassScale.			set_global_ips(vSCVip, numSCVip);
	m_imMass.				set_global_ips(vSCVip, numSCVip);
	*/
}

//	virtual postprocesses the loop over all elements of one type
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void PNP_1D<TDomain>::fsh_elem_loop()
{}

// virtual Assembling of Defect (Stiffness part)
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void PNP_1D<TDomain>::add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<worldDim> vCornerCoords[])
{
	// get finite volume geometry
	static const TFVGeom& fvgeom = GeomProvider<TFVGeom>::get();

	// cast elem to appropriate type
	TElem* pElem = dynamic_cast<TElem*>(elem);
	if (!pElem) {UG_THROW("Wrong element type.");}

// divergence terms (flux through the dendrite)
	// loop SCVF
	for (size_t s = 0; s < fvgeom.num_scvf(); s++)
	{
		// get current SCVF
		const typename TFVGeom::SCVF& scvf = fvgeom.scvf(s);

		// get radius as mean value of the neighbouring vertex radii
		number radiusAtIP = 0.5 * (  m_aaRadius[pElem->vertex(scvf.from())]
								   + m_aaRadius[pElem->vertex(scvf.to())]);

		// scvf area
		number scvfArea = PI * radiusAtIP * radiusAtIP;

	// potential
		// to compute eps \nabla phi
		MathVector<worldDim> eps_grad_phi, grad_phi;

		// compute gradient and shape at ip
		VecSet(grad_phi, 0.0);
		for (size_t sh = 0; sh < scvf.num_sh(); sh++)
			VecScaleAppend(grad_phi, u(_PHI_,sh), scvf.global_grad(sh));

		// scale by permittivity tensor
		MatVecMult(eps_grad_phi, m_permittivity_dend, grad_phi);

		// compute flux density
		number flux = VecDot(eps_grad_phi, scvf.normal());

		// scale with hidden area of dendritic crosssection and Faraday const
		flux *= scvfArea / F;

		// add to local defect
		d(_PHI_, scvf.from()) -= flux;
		d(_PHI_, scvf.to()  ) += flux;

//std::cout << "potential div: " << flux << std::endl;

	// ion concentrations
		for (size_t i = 0; i < m_nIon; i++)
		{
		// diffusion
			// to compute D \nabla c_i
			MathVector<worldDim> D_grad, grad;

			// compute gradient and shape at ip
			VecSet(grad, 0.0);
			for (size_t sh = 0; sh < scvf.num_sh(); sh++)
				VecScaleAppend(grad, u(i,sh), scvf.global_grad(sh));

			// scale by diffusion tensor
			MatVecMult(D_grad, m_vDiffusionTensor[i], grad);

			// compute flux density
			number diff_flux = VecDot(D_grad, scvf.normal());

			// scale with hidden area of dendritic crosssection
			diff_flux *= scvfArea;

			// add to local defect
			d(i, scvf.from()) -= diff_flux;
			d(i, scvf.to()  ) += diff_flux;

//std::cout << "\t\t\t\t\t\t\tion diffusion " << i << ": " << diff_flux << std::endl;

		// electricity
			// get concentration at ip
			number c_i = 0.0;
			for (size_t sh = 0; sh < scvf.num_sh(); sh++)
				c_i += u(i,sh) * scvf.shape(sh);

			// scale electric field by diffusion tensor
			MathVector<worldDim> D_grad_phi;
			MatVecMult(D_grad_phi, m_vDiffusionTensor[i], grad_phi);

			// compute inner product
			number el_flux = VecDot(D_grad_phi, scvf.normal());

			// scale by concentration and constants
			el_flux *= scvfArea * m_vValency[i] * F/(R*T) * c_i;

			// add to local defect
			d(i, scvf.from()) -= el_flux;
			d(i, scvf.to()  ) += el_flux;

//std::cout << "\t\t\t\t\t\t\tion electric " << i << ": " << el_flux << std::endl;
		}
	}

// volume parts (using lumping)
	// loop SCV
	for (size_t s = 0; s < fvgeom.num_scv(); s++)
	{
		// get current SCV
		const typename TFVGeom::SCV& scv = fvgeom.scv(s);

		// get associated node
		const int co = scv.node_id();

		// radius at integration point
		number radiusAtIP = m_aaRadius[pElem->vertex(co)];

		// surface area of dendritic scv
		number scvArea = 2*PI * radiusAtIP * scv.volume();

		// volume of dendritic scv
		number scvVolume = PI * radiusAtIP * radiusAtIP * scv.volume();

	// potential equation

		// E field over cylinder surface area (supposedly constant)
		number el_surf_flux = m_permittivity_mem * u(_PHI_,co) / m_mem_thickness;

		// scale by surface area and Faraday constant
		el_surf_flux *= scvArea / F;

		d(_PHI_, co) -= el_surf_flux;

		// charge density term
		for (size_t i = 0; i < m_nIon; i++)
		{
			d(_PHI_, co) += u(i, co) * m_vValency[i] * scvVolume;
//std::cout << "charge density " << i << ": " << u(i, co) * m_vValency[i] * scvVolume << std::endl;
		}
	// ion concentrations
		for (size_t i = 0; i < m_nIon; i++)
		{
		// passive channel flux
			// calculate current
			number flux = m_vSpecConductance[i] * (u(_PHI_,co) - m_vRevPot[i]);

			// scale by constants and volume
			flux *= scvArea / (m_vValency[i] * F);

			// add to local defect
			d(i, co) += flux;
//std::cout << "\t\t\t\t\t\t\tpassive leak " << i << ": " << flux << std::endl;
		}
	}
//std::cout << "-----------------------------------------------------------------------------" << std::endl;
}

// virtual Assembling of Defect (Mass part)
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void PNP_1D<TDomain>::add_def_M_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<worldDim> vCornerCoords[])
{
	// get finite volume geometry
	static const TFVGeom& fvgeom = GeomProvider<TFVGeom>::get();

	// cast elem to appropriate type
	TElem* pElem = dynamic_cast<TElem*>(elem);
	if (!pElem) {UG_THROW("Wrong element type.");}

// defects wrt time (using lumping)
	// loop SCV
	for (size_t s = 0; s < fvgeom.num_scv(); s++)
	{
		// get current SCV
		const typename TFVGeom::SCV& scv = fvgeom.scv(s);

		// get associated node
		const int co = scv.node_id();

		// radius at integration point
		number radiusAtIP = m_aaRadius[pElem->vertex(co)];

		// surface area of dendritic scv
		number scvArea = 2*PI * radiusAtIP * scv.volume();

		// volume of dendritic scv
		number scvVolume = PI * radiusAtIP * radiusAtIP * scv.volume();

		for (size_t i = 0; i < m_nIon; i++)
		{
		// potential in ion equations
			// calculate current
			number flux = u(_PHI_, co) * m_vSpecCapacity[i];

			// scale by constants and volume
			flux *= scvArea / (m_vValency[i] * F);

			d(i, co) += flux;
//std::cout << "%\t\t\t\t\t\t\tcapacitary " << i << ": " << flux << std::endl;

		// ion concentrations
			// add to local matrix
			d(i, co) += u(i,co) * scvVolume;
//std::cout << "%\t\t\t\t\t\t\tion mass " << i << ": " << u(i,co) * scvVolume << std::endl;
		}
	}
//std::cout << "-----------------------------------------------------------------------------" << std::endl;
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void PNP_1D<TDomain>::add_rhs_elem(LocalVector& rhs, GridObject* elem, const MathVector<worldDim> vCornerCoords[])
{}

// Assembling of Jacobian (Stiffness part)
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void PNP_1D<TDomain>::add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<worldDim> vCornerCoords[])
{
	// get finite volume geometry
	static const TFVGeom& fvgeom = GeomProvider<TFVGeom>::get();

	// cast elem to appropriate type
	TElem* pElem = dynamic_cast<TElem*>(elem);
	if (!pElem) {UG_THROW("Wrong element type.");}

// divergence terms (flux through the dendrite)
	// loop SCVF
	for (size_t s = 0; s < fvgeom.num_scvf(); s++)
	{
		// get current SCVF
		const typename TFVGeom::SCVF& scvf = fvgeom.scvf(s);

		// get radius as mean value of the neighbouring vertex radii
		number radiusAtIP = 0.5 * (  m_aaRadius[pElem->vertex(scvf.from())]
								   + m_aaRadius[pElem->vertex(scvf.to())]);

		// scvf area
		number scvfArea = PI * radiusAtIP * radiusAtIP;

		// compute grad(phi)
		MathVector<worldDim> D_grad_phi, grad_phi;
		VecSet(grad_phi, 0.0);
		for (size_t sh = 0; sh < scvf.num_sh(); sh++)
			VecScaleAppend(grad_phi, u(_PHI_,sh), scvf.global_grad(sh));

	// potential
		// to compute eps \nabla phi
		MathVector<worldDim> eps_grad;

		// compute gradient and shape at ip
		for (size_t sh = 0; sh < scvf.num_sh(); sh++)
		{
			// scale by permittivity tensor
			MatVecMult(eps_grad, m_permittivity_dend, scvf.global_grad(sh));

			// compute flux density
			number flux = VecDot(eps_grad, scvf.normal());

			// scale with hidden area of dendritic crosssection and Faraday const
			flux *= scvfArea / F;

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
				MatVecMult(D_grad, m_vDiffusionTensor[i], scvf.global_grad(sh));

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
				c_i += u(i,sh) * scvf.shape(sh);

			for (size_t sh = 0; sh < scvf.num_sh(); sh++)
			{
			// deriv wrt concs
				// scale electric field by diffusion tensor
				MatVecMult(D_grad_phi, m_vDiffusionTensor[i], grad_phi);

				// compute inner product
				number el_flux = VecDot(D_grad_phi, scvf.normal());

				// scale by concentration and constants
				el_flux *= scvfArea * m_vValency[i] * F/(R*T) * scvf.shape(sh);

				// add to local Jacobian
				J(i, scvf.from(), i, sh) -= el_flux;
				J(i, scvf.to(),   i, sh) += el_flux;

			// deriv wrt potential
				// scale electric field by diffusion tensor
				MatVecMult(D_grad, m_vDiffusionTensor[i], scvf.global_grad(sh));

				// compute inner product
				el_flux = VecDot(D_grad, scvf.normal());

				// scale by concentration and constants
				el_flux *= scvfArea * m_vValency[i] * F/(R*T) * c_i;

				// add to local defect
				J(i, scvf.from(), _PHI_, sh) -= el_flux;
				J(i, scvf.to(),   _PHI_, sh) += el_flux;
			}
		}
	}

// volume parts (using lumping)
	// loop SCV
	for (size_t s = 0; s < fvgeom.num_scv(); s++)
	{
		// get current SCV
		const typename TFVGeom::SCV& scv = fvgeom.scv(s);

		// get associated node
		const int co = scv.node_id();

		// radius at integration point
		number radiusAtIP = m_aaRadius[pElem->vertex(co)];

		// surface area of dendritic scv
		number scvArea = 2*PI * radiusAtIP * scv.volume();

		// volume of dendritic scv
		number scvVolume = PI * radiusAtIP * radiusAtIP * scv.volume();

	// potential equation
		// E field over cylinder surface area
		J(_PHI_, co, _PHI_, co) -= m_permittivity_mem / m_mem_thickness * scvArea / F;

		// charge density term
		for (size_t i = 0; i < m_nIon; i++)
			J(_PHI_, co, i, co) += m_vValency[i] * scvVolume;

	// ion concentrations
		for (size_t i = 0; i < m_nIon; i++)
		{
		// passive channel flux
			// calculate current
			number flux = m_vSpecConductance[i];

			// scale by constants and volume
			flux *= scvArea / (m_vValency[i] * F);

			// add to local Jacobian
			J(i, co, _PHI_, co) += flux;
		}
	}
}

// Assembling of Jacobian (Mass part)
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void PNP_1D<TDomain>::add_jac_M_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<worldDim> vCornerCoords[])
{
	// get finite volume geometry
	static const TFVGeom& fvgeom = GeomProvider<TFVGeom>::get();

	// cast elem to appropriate type
	TElem* pElem = dynamic_cast<TElem*>(elem);
	if (!pElem) {UG_THROW("Wrong element type.");}

// derivs wrt time (using lumping)
	// loop SCV
	for (size_t s = 0; s < fvgeom.num_scv(); s++)
	{
		// get current SCV
		const typename TFVGeom::SCV& scv = fvgeom.scv(s);

		// get associated node
		const int co = scv.node_id();

		// radius at integration point
		number radiusAtIP = m_aaRadius[pElem->vertex(co)];

		// surface area of dendritic scv
		number scvArea = 2*PI * radiusAtIP * scv.volume();

		// volume of dendritic scv
		number scvVolume = PI * radiusAtIP * radiusAtIP * scv.volume();

		for (size_t i = 0; i < m_nIon; i++)
		{
		// potential in ion equations
			// calculate current
			number flux = m_vSpecCapacity[i];

			// scale by constants and volume
			flux *= scvArea / (m_vValency[i] * F);

			// add to Jacobian
			J(i, co, _PHI_, co) += flux;

			// ion concentrations
			// add to Jacobian
			J(i, co, i, co) += scvVolume;
		}
	}
}

//////////////////////////////////////////////
//   Error estimation assembling methods	//
//////////////////////////////////////////////

//	virtual prepares the loop over all elements of one type for the computation of the error estimator
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void PNP_1D<TDomain>::prep_err_est_elem_loop(const ReferenceObjectID roid, const int si)
{
	UG_THROW("PNP_1D::prep_err_est_elem_loop not yet implemented.");
}

//	virtual prepares the loop over all elements of one type for the computation of the error estimator
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void PNP_1D<TDomain>::prep_err_est_elem(const LocalVector& u, GridObject* elem, const MathVector<worldDim> vCornerCoords[])
{
	UG_THROW("PNP_1D::prep_err_est_elem not yet implemented.");
}

//	virtual compute the error estimator (stiffness part) contribution for one element
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void PNP_1D<TDomain>::compute_err_est_A_elem(const LocalVector& u, GridObject* elem, const MathVector<worldDim> vCornerCoords[], const number& scale)
{
	UG_THROW("PNP_1D::compute_err_est_A_elem not yet implemented.");
}

//	virtual compute the error estimator (mass part) contribution for one element
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void PNP_1D<TDomain>::compute_err_est_M_elem(const LocalVector& u, GridObject* elem, const MathVector<worldDim> vCornerCoords[], const number& scale)
{
	UG_THROW("PNP_1D::compute_err_est_M_elem not yet implemented.");
}

//	virtual postprocesses the loop over all elements of one type in the computation of the error estimator
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void PNP_1D<TDomain>::fsh_err_est_elem_loop()
{
	UG_THROW("PNP_1D::fsh_err_est_elem_loop not yet implemented.");
}


//////////////////////////////////////////
//   registering assembling functions   //
//////////////////////////////////////////

template<typename TDomain>
void PNP_1D<TDomain>::
register_all_funcs()
{
//	get all grid element types in this
	typedef typename domain_traits<1>::DimElemList ElemList;

//	switch assemble functions
	boost::mpl::for_each<ElemList>(Register(this));
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void PNP_1D<TDomain>::
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
	template class PNP_1D<Domain1d>;
#endif
#ifdef UG_DIM_2
	template class PNP_1D<Domain2d>;
#endif
#ifdef UG_DIM_3
	template class PNP_1D<Domain3d>;
#endif

} // namespace nernst_planck
} // namespace ug
