/*
 * pnp1d_fv.h
 *
 *  Created on: 13.08.2015
 *      Author: mbreit
 */

#ifndef UG__PLUGINS__EXPERIMENTAL__NERNST_PLANCK__PNP1D_FV_H_
#define UG__PLUGINS__EXPERIMENTAL__NERNST_PLANCK__PNP1D_FV_H_

#include "common/common.h"
#include "common/math/ugmath.h"
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"
#include "lib_disc/spatial_disc/disc_util/fvho_geom.h"
#include "lib_disc/spatial_disc/disc_util/geom_provider.h"


namespace ug{
namespace nernst_planck{

/// \addtogroup plugin_nernst_planck
/// \{

/// 1D PNP element discretization for higher-order FV schemes
/**	This class implements the IElemDisc interface to provide a 1D simplified
 *	FV discretization for the Poisson-Nernst-Planck (PNP) problem on branched
 *	tubular objects of variable radius (intended usage: dendrites), incorporating
 *	a leakage flux across the tubule's surface boundary as well as a capacitive
 *	flux across the same boundary simulating the movement of ions into / out of
 *	the ionic layer that exists there because of influency effects (the membrane
 *	is supposed to be positively charged).
 *
 *	The simplification to 1D assumes that the geometry is perfectly rotationally
 *	symmetric and that the unknowns only vary in value axially. It is then the
 *	axial dimension that solely remains to be resolved in space.
 *
 *	Thus, the following set of equations is discretized by the FV scheme:
 *
 *	1) Equation for ion species i:
 *  \f[
 * 		\int\limits_{B} \frac{\partial c_i}{\partial t} d\mu =
 * 		\int\limits_{{\partial B}_{ax}} D_i\left( \vec{\nabla} c_i +
 * 										 \frac{z_iF}{RT} c_i \vec{\nabla} \phi \right)
 * 			\cdot \vec{n}_{{\partial B}_{ax}} d\nu
 * 		-\int\limits_{{\partial B}_{r}} \frac{1}{z_iF}
 * 			\left( g_i \left( \phi-E_i \right) + C_i \frac{\partial \phi}{\partial t} \right) d\nu
 *	\f]
 *
 *	2) Equation for the potential:
 *	\f[
 * 		\int\limits_B \sum_i z_i F c_i d\mu =
 * 		\int\limits_{{\partial B}_{ax}} \varepsilon_0 \varepsilon_r \vec{\nabla}\phi \cdot \vec{n} d\nu
 *	\f]
 *
 *	Both of the equations are supposed to hold for any arbitrarily thin
 *	cylindrical piece \f$ B \f$ of the tubular object whose boundary is decomposed
 *	into the two caps of the cylinder (\f$ {\partial B}_{ax} \f$) and the surface
 *	area (\f$ {\partial B}_r \f$).
 *
 *	Meaning of the symbols:
 * <ul>
 * <li>	\f$ c_i \f$ the i-th unknown ion concentration [mol/m^3],
 * <li>	\f$ \phi \f$ the unknown potential [V],
 * <li>	\f$ D_i \f$ diffusion tensor for the i-th ion species [m^2/s],
 * <li>	\f$ z_i \f$ valency of the i-th ion species [1],
 * <li>	\f$ g_i \f$ specific membrane conductance for the i-th ion species [C/Vsm^2],
 * <li>	\f$ C_i \f$ specific membrane layer capacity for the i-th ion species [C/(Vm^2)],
 * <li>	\f$ E_i \f$ membrane reversal potential for the i-th ion species [V],
 * <li>	\f$ \varepsilon_r \f$ relative permettivity of the tubule [1],
 * <li>	\f$ \varepsilon_0 \f$ vacuum permettivity (physical constant) [C/Vm],
 * <li>	\f$ R \f$ universal gas constant (physical constant) [VC/molK],
 * <li>	\f$ T \f$ temperature [K],
 * <li>	\f$ F \f$ Faraday constant [C/mol],
 * </ul>
 *
 * \note As the potential equation is stationary in nature but the whole
 * 		 discretization needs to be treated as time-dependent, it is a
 * 		 requirement to use the backward Euler method for time discretization.
 * 		 Otherwise, the potential equation will be assembled incorrectly.
 *
 * \todo The current implementation calculates with more volume and surface
 * 		 at branching points than do exist. This may lead to slower dynamics there.
 * 		 However, as the assumption of piecewise cylindrical compartments is inexact
 *		 anyway, this might not be all too bad.
 *
 * \date 13.08.2015
 * \author mbreit
*/

template <typename TDomain>
class PNP1D_FV: public IElemDisc<TDomain>
{
	public:
		typedef PNP1D_FV<TDomain> this_type;

		typedef typename TDomain::position_attachment_type pos_att_type;
		typedef typename TDomain::position_accessor_type pos_acc_type;

		/// error estimator type
		typedef MultipleSideAndElemErrEstData<TDomain> err_est_type;

		static const int worldDim = TDomain::dim;

		///	constructor
		PNP1D_FV(SmartPtr<ApproximationSpace<TDomain> > approx,
			   const char* functions = "", const char* subsets = "");

		///	constructor
		PNP1D_FV(SmartPtr<ApproximationSpace<TDomain> > approx,
			   const std::vector<std::string>& vFct, const std::vector<std::string>& vSubset);

		/// destructor
		virtual ~PNP1D_FV() {};

	// inherited methods
		/// requests assembling for trial spaces and grid type
		/**
		 * This function is called before the assembling starts. The
		 * IElemDisc-Implementation is supposed to checks if it can assemble the set
		 * of LFEID and the grid type. It may register corresponding assembling
		 * functions or perform other initialization.
		 * If the ElemDisc does not support the setting it should throw an exception.
		 *
		 * \param[in] vLfeID			vector of Local Finite Element IDs
		 * \param[in] bNonRegularGrid	regular grid type
		 */
			virtual void prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid);

		///	returns if discretization acts on hanging nodes if present
		/**
		 * This function returns if a discretization really needs the hanging nodes
		 * in case of non-regular grid. This may not be the case for e.g. finite
		 * element assemblings but is needed for finite volumes
		 */
			virtual bool use_hanging() const;


	////////////////////////////
	// assembling functions
	////////////////////////////

		/// prepare the timestep
		//virtual void prep_timestep_elem(number future_time, const number time, const LocalVector& u, GridObject* elem, const MathVector<worldDim> vCornerCoords[]);

		///	virtual prepares the loop over all elements of one type
		template <typename TElem, typename TFVGeom>
		void prep_elem_loop(const ReferenceObjectID roid, const int si);

		///	virtual prepare one elements for assembling
		template <typename TElem, typename TFVGeom>
		void prep_elem(const LocalVector& u, GridObject* elem, const ReferenceObjectID roid, const MathVector<worldDim> vCornerCoords[]);

		///	virtual postprocesses the loop over all elements of one type
		template <typename TElem, typename TFVGeom>
		void fsh_elem_loop();

		/// virtual Assembling of Defect (Stiffness part)
		template <typename TElem, typename TFVGeom>
		void add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<worldDim> vCornerCoords[]);

		/// virtual Assembling of Defect (Mass part)
		template <typename TElem, typename TFVGeom>
		void add_def_M_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<worldDim> vCornerCoords[]);

		///	assembles the local right hand side
		template<typename TElem, typename TFVGeom>
		void add_rhs_elem(LocalVector& rhs, GridObject* elem, const MathVector<worldDim> vCornerCoords[]);

		/// Assembling of Jacobian (Stiffness part)
		template <typename TElem, typename TFVGeom>
		void add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<worldDim> vCornerCoords[]);

		/// Assembling of Jacobian (Mass part)
		template <typename TElem, typename TFVGeom>
		void add_jac_M_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<worldDim> vCornerCoords[]);

		///	virtual prepares the loop over all elements of one type for the computation of the error estimator
		template <typename TElem, typename TFVGeom>
		void prep_err_est_elem_loop(const ReferenceObjectID roid, const int si);

		///	virtual prepares the loop over all elements of one type for the computation of the error estimator
		template <typename TElem, typename TFVGeom>
		void prep_err_est_elem(const LocalVector& u, GridObject* elem, const MathVector<worldDim> vCornerCoords[]);

		///	virtual compute the error estimator (stiffness part) contribution for one element
		template <typename TElem, typename TFVGeom>
		void compute_err_est_A_elem(const LocalVector& u, GridObject* elem, const MathVector<worldDim> vCornerCoords[], const number& scale);

		///	virtual compute the error estimator (mass part) contribution for one element
		template <typename TElem, typename TFVGeom>
		void compute_err_est_M_elem(const LocalVector& u, GridObject* elem, const MathVector<worldDim> vCornerCoords[], const number& scale);

		///	virtual postprocesses the loop over all elements of one type in the computation of the error estimator
		template <typename TElem, typename TFVGeom>
		void fsh_err_est_elem_loop();

	public:
	// setters for PNP parameters
		/// set the names of the ion species present
		void set_ions(const std::vector<std::string>& vIons);

		/// set ion valencies
		void set_valencies(const std::vector<int>& vValencies);

		/// set ion reversal potentials
		void set_reversal_potentials(const std::vector<number>& vRevPot);

		/// set ion specific conductances of the plasma membrane
		void set_specific_conductances(const std::vector<number>& vSpecCond);

		/// set ion specific capacities for influency effects at the plasma membrane
		void set_specific_capacities(const std::vector<number>& vSpecCap);

		/// set ion diffusion constants
		void set_diffusion_constants(const std::vector<number>& vDiff);

		/// set effective electric permettivity in the dendrite (eps_0 * eps_r)
		void set_permettivities(const number eps_dend, const number eps_mem);

		/// set constant dendritic radius
		void set_membrane_thickness(const number d);

		/// set constant dendritic radius
		void set_dendritic_radius(const number r);

		/// set constants R, T, F
		void set_rtf(number R, number T, number F);

		/// set represented dimensionality
		void set_represented_dimension(size_t dim);


	protected:
		/// approximation space (needed before given by domainDisc)
		SmartPtr<ApproximationSpace<TDomain> > m_spApproxSpace;

		///	current regular grid flag
		bool m_bNonRegularGrid;

		/// dendritic radius attachment and accessor
		ANumber m_aRadius;
		Grid::AttachmentAccessor<Vertex, ANumber> m_aaRadius;

		/// number of ions == index of potential unknown
		size_t m_nIon;
		size_t _PHI_;

		/// ion names
		std::vector<std::string> m_vIonNames;

		/// ion valencies
		std::vector<int> m_vValency;

		/// ion reversal potentials
		std::vector<number> m_vRevPot;

		/// ion specific conductances of the plasma membrane
		std::vector<number> m_vSpecConductance;

		/// ion specific "capacities" for influence effects at the plasma membrane
		std::vector<number> m_vSpecCapacity;

		/// ion diffusion tensor
		std::vector<MathMatrix<worldDim,worldDim> > m_vDiffusionTensor;

		/// electric permettivity (eps_r * eps_0)
		MathMatrix<worldDim,worldDim> m_permittivity_dend;
		number m_permittivity_mem;

		/// membrane thickness
		number m_mem_thickness;

		/// universal gas constant
    	number m_R;

    	/// temperature
    	number m_T;

    	/// Faraday constant
		number m_F;

		/// represented dimension (default: 3)
		size_t m_reprDim;

	protected:
		///	shape function set
		LFEID m_lfeID;

		///	integration order
		int m_quadOrder;


	private:
		///	register utils
		///	@{
		void register_all_funcs();
		template <typename TElem, typename TFVGeom> void register_func();
		/// @}
};

/// \}

} // namespace nernst_planck
} // namespace ug

#endif // UG__PLUGINS__EXPERIMENTAL__NERNST_PLANCK__PNP1D_FV_H_
