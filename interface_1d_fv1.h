/*
 * interface_1d_fv1.h
 *
 *  Created on: 27.05.2014
 *      Author: mbreit
 */

#ifndef INTERFACE_1D_FV1_H_
#define INTERFACE_1D_FV1_H_

#include "common/common.h"
#include "common/math/ugmath.h"
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"


namespace ug{

template <typename TDomain>
class IInterface1DFV1: public IElemDisc<TDomain>
{
	public:
		static const int worldDim = TDomain::dim;

		struct InterfaceNodes
		{
			friend class IInterface1DFV1<TDomain>;
			private:
				MathVector<worldDim> coords[2];
		};

		///	constructor
		IInterface1DFV1(const char* functions = "", const char* subsets = "");

		///	constructor
		IInterface1DFV1(const std::vector<std::string>& vFct, const std::vector<std::string>& vSubset);

		/// destructor
		virtual ~IInterface1DFV1() {};

	// inherited methods
		///	 returns the type of elem disc
		//virtual int type() const {return EDT_ELEM | EDT_SIDE;}

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
			virtual void prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid) = 0;

		///	returns if discretization acts on hanging nodes if present
		/**
		 * This function returns if a discretization really needs the hanging nodes
		 * in case of non-regular grid. This may not be the case for e.g. finite
		 * element assemblings but is needed for finite volumes
		 */
			virtual bool use_hanging() const {return false;}

		////////////////////////////
		// assembling functions
		////////////////////////////
		public:
		/// prepare the timestep
			virtual void prep_timestep_elem(const number time, const LocalVector& u, GridObject* elem, const MathVector<worldDim> vCornerCoords[]);

		///	virtual prepares the loop over all elements of one type
			virtual void prep_elem_loop(const ReferenceObjectID roid, const int si);

		///	virtual prepare one elements for assembling
			virtual void prep_elem(const LocalVector& u, GridObject* elem, const MathVector<worldDim> vCornerCoords[]);

		///	virtual postprocesses the loop over all elements of one type
			virtual void fsh_elem_loop();

		/// virtual finish the timestep
			virtual void fsh_timestep_elem(const number time, const LocalVector& u, GridObject* elem, const MathVector<worldDim> vCornerCoords[]);

		/// Assembling of Jacobian (Stiffness part)
			virtual void add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<worldDim> vCornerCoords[]);

		/// Assembling of Jacobian (Mass part)
			virtual void add_jac_M_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<worldDim> vCornerCoords[]);

		/// virtual Assembling of Defect (Stiffness part)
			virtual void add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<worldDim> vCornerCoords[]);

		/// defect for explicit terms
			virtual void add_def_A_expl_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<worldDim> vCornerCoords[]);

		/// virtual Assembling of Defect (Mass part)
			virtual void add_def_M_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<worldDim> vCornerCoords[]);

		/// virtual Assembling of Right-Hand Side
			virtual void add_rhs_elem(LocalVector& rhs, GridObject* elem, const MathVector<worldDim> vCornerCoords[]);

		///	virtual prepares the loop over all elements of one type for the computation of the error estimator
			virtual void prep_err_est_elem_loop(const ReferenceObjectID roid, const int si);

		///	virtual prepares the loop over all elements of one type for the computation of the error estimator
			virtual void prep_err_est_elem(const LocalVector& u, GridObject* elem, const MathVector<worldDim> vCornerCoords[]);

		///	virtual compute the error estimator (stiffness part) contribution for one element
			virtual void compute_err_est_A_elem(const LocalVector& u, GridObject* elem, const MathVector<worldDim> vCornerCoords[], const number& scale);

		///	virtual compute the error estimator (mass part) contribution for one element
			virtual void compute_err_est_M_elem(const LocalVector& u, GridObject* elem, const MathVector<worldDim> vCornerCoords[], const number& scale);

		///	virtual compute the error estimator (rhs part) contribution for one element
			virtual void compute_err_est_rhs_elem(GridObject* elem, const MathVector<worldDim> vCornerCoords[], const number& scale);

		///	virtual postprocesses the loop over all elements of one type in the computation of the error estimator
			virtual void fsh_err_est_elem_loop();



	protected:
		InterfaceNodes m_in;

};


} // namespace ug

#endif /* INTERFACE_1D_FV1_H_ */
