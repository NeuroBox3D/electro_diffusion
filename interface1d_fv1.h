/*
 * interface1d_fv1.h
 *
 *  Created on: 06.06.2014
 *      Author: mbreit
 */

#ifndef INTERFACE1D_FV1_H_
#define INTERFACE1D_FV1_H_


#include <map>
#include <limits>	// for numeric_limits<>::max()

#include "common/common.h"
#include "common/math/ugmath.h"
#include "lib_algebra/operator/interface/matrix_operator.h"
#include "lib_disc/spatial_disc/constraints/constraint_interface.h"
#include "lib_disc/spatial_disc/disc_util/fv1_geom.h"
#include "lib_disc/spatial_disc/disc_util/geom_provider.h"


namespace ug{
namespace nernst_planck{

/// Interface class for coupling high-dimensional (2D/3D) discretizations with 1D simplifications.
/**
 * This class is intended to be the base class for implementations of interfaces
 * linking a spatially fully detailed 2D/3D model, typically a "window" of a much
 * larger geometry where only local dynamics are of interest, but where processes
 * in remote locations still matter (e.g. for realistic boundary conditions in the
 * window), with a 1D simplified model for the rest of the geometry.
 *
 * This class provides an implementation for FV discretizations of both sides of
 * the interface.
 *
 * ------------------------
 * Geometry specifications
 * ------------------------
 * Suppose you have the following situation: left: 2D, right: 1D
 *
 *     --- o --- o
 *         |  /  |
 *     --- o --- O       O --- o ---
 *         |  /  |
 *     --- o --- o
 *
 * The two nodes marked by a capital O will be referred to as "interface" nodes/vertices.
 *
 * Then you have to introduce constrained nodes (x) on both sides that will simulate
 * the behaviour of a non-existing extension of the respective subdomain:
 *
 *     --- o --- o  ---  x
 *         |  /  |       |
 *     --- o  -  O/x === x/O --- o ---
 *         |  /  |       |
 *     --- o --- o  ---  x
 *
 * There are two vertices located at the same coordinates in the location of the
 * interface nodes. One is the interface vertex, the other is a constrained vertex
 * connected to the other interface vertex.
 * Each constrained vertex is connected to exactly one non-constrained vertex, which
 * is referred to as its "constrainer" or "constraining vertex".
 * The two subdomains remain disconnected!
 *
 * The constrained vertices on the high-dimensional side, in conjunction with their
 * respective constrainers, form new high-dimensional elements (quadrilaterals in 2D,
 * prisms or hexahedra in 3D). A discretization on the high-dimensional side with
 * an equivalent on the 1D side that is to be coupled by means of this interface
 * will have to be extended to those new elements.
 * This holds for the 1D side analogously.
 *
 * Subsets:
 * ---------
 * The constrained manifold (at least the vertices) of the high-dimensional side
 * must be (exclusively) assigned to a subset that needs to be specified in the
 * interface constructor. The same is necessary for the (single) 1D constrained
 * vertex.
 *
 * -------------
 * How it works
 * -------------
 * The coupling is realized by constraints (this is why this base class is derived
 * from IDomainConstraint). As the respective discretizations on both sides of the
 * interface which are to be coupled are extended to the respective new elements
 * nothing has to be assembled. The one thing that has to be done is assigning
 * meaningful values to the constrained vertices.
 *
 * The constraint takes the general form:   u* = f(u,u1,u2)
 * where u* is the value of a function at a constrained vertex, u the respective
 * value at its constrainer, u1 is the value on the constraining side, u2 that on
 * the constrained one.
 *
 * The constraint function f is problem-specific and needs to be carefully chosen.
 * In order to define such a function, you will have to derive your own interface
 * class from this base class and implement the two virtual methods
 *
 * - constraintValue
 * - constraintValueDerivs
 *
 * constraintValue must calculate the value of the constraints, i.e.
 * \f[
 *      f(u,u_1,u_2);
 * \f]
 * constraintValueDerivs must calculate the partial derivations wrt
 * \f$ u \f$, \f$ u_1 \f$, \f$ u_2 \f$ thereof.
 * Cf. documentation of the virtual methods.
 *
 * Two default implementations are realized in the derived classes
 * AdditiveInterface1DFV1 and MultiplicativeInterface1DFV1 with:
 *
 * \f[
 *		f(u,u_1,u_2) = u + (u_2-u_1);
 * \f]
 * \f[
 * 		f(u,u_1,u_2) = u \cdot (u_2/u_1),
 * \f]
 * respectively.
 *
 *
 * \note At the moment, if you want to use parallel GMG with a gathered base solver
 * with this implementation, you have to set gmg:set_rap(true); otherwise you will
 * experience MPI failures, since the gathering process will try to communicate
 * with all others in order to assemble the constraints given here, whereas the
 * other processes will not reach this point of communication.
 *
 * \tparam	TDomain				type of Domain
 * \tparam	TAlgebra			type of Algebra
 *
 * \date 08.07.2014
 * \author mbreit
 */


template <typename TDomain, typename TAlgebra>
class IInterface1DFV1: public IDomainConstraint<TDomain, TAlgebra>
{
	public:
		/// own type
		typedef IInterface1DFV1<TDomain, TAlgebra> this_type;

		/// world dimension
		static const int worldDim = TDomain::dim;

		///	domain type
		typedef TDomain domain_type;

		///	algebra type
		typedef TAlgebra algebra_type;

		///	type of algebra matrix
		typedef typename algebra_type::matrix_type matrix_type;

		///	type of algebra vector
		typedef typename algebra_type::vector_type vector_type;

		// element type
		typedef typename domain_traits<worldDim>::element_type elem_type;

		///	constructor
		/**
		 * @param fcts				the functions to which the interface applies
		 * 							(must be defined on both sides of the interface!)
		 * @param high_dim_subset	the name of the constrained subset on the high-dimensional side
		 * @param one_dim_subset	the name of the constrained subset on the one-dimensional side
		 * @param extension_subset	the name of the extension domain subset (1D side)
		 */
		IInterface1DFV1(const char* fcts, const char* high_dim_subset, const char* one_dim_subset,
						const char* two_dim_intfNode, const char* extension_subset);

		/// destructor
		virtual ~IInterface1DFV1();


	// inherited methods (from IConstraint)

		///	returns the type of constraints
		virtual int type() const;

		///	sets the approximation space
		/**
		 * Is called by the domain discretization the constraint is added to
		 * before constraint application.
		 *
		 * @param approxSpace	the approximation space to be set
		 */
		virtual void set_approximation_space(SmartPtr<ApproximationSpace<TDomain> > approxSpace);

		//void check_values_at_interface(const vector_type& sol);

		///	adapts jacobian to enforce constraints
		virtual void adjust_jacobian
		(	matrix_type& J,
			const vector_type& u,
			ConstSmartPtr<DoFDistribution> dd,
			number time = 0.0,
			ConstSmartPtr<VectorTimeSeries<vector_type> > vSol = SPNULL,
			const number s_a0 = 1.0
		);

		///	adapts defect to enforce constraints
		virtual void adjust_defect
		(	vector_type& d,
			const vector_type& u,
			ConstSmartPtr<DoFDistribution> dd,
			number time = 0.0,
			ConstSmartPtr<VectorTimeSeries<vector_type> > vSol = SPNULL,
			const std::vector<number>* vScaleMass = NULL,
			const std::vector<number>* vScaleStiff = NULL
		);

		///	adapts matrix and rhs (linear case) to enforce constraints
		virtual void adjust_linear
		(	matrix_type& mat,
			vector_type& rhs,
			ConstSmartPtr<DoFDistribution> dd,
			number time = 0.0
		);

		///	adapts a rhs to enforce constraints
		virtual void adjust_rhs
		(	vector_type& rhs,
			const vector_type& u,
			ConstSmartPtr<DoFDistribution> dd,
			number time = 0.0
		);

		///	sets the constraints in a solution vector
		virtual void adjust_solution
		(	vector_type& u,
			ConstSmartPtr<DoFDistribution> dd,
			number time = 0.0
		);

	protected:
		/// called when the approximation space has changed
		void approximation_space_changed();

	public:
		// methods to be implemented by a concretization of this interface

		/// function returning the defect value for a constrained node
		/**
		 * The defect is: f(constrainingDoF, interfaceDoF0, interfaceDoF1),
		 * where f is any inter-/extrapolating function that is suited to calculate the value
		 * for the constrainedDoF from the values of the constraining DoF and the two interface
		 * DoFs. Of course, f should be differentiable.
		 * interfaceDoF0 is on the same side as the constrainingDoF,
		 * interfaceDoF1 is on the side of the constrainedDoF.
		 *
		 * @param d			defect value at constrained vertex; this is the output
		 * @param u_c		solution at corresponding constrainer vertex
		 * @param u_itf0	solution at interface vertex on constraining side
		 * @param u_itf1	solution at interface vertex on constrained side
		 */
		virtual void constraintValue
		(	typename vector_type::value_type& d,
			const typename vector_type::value_type& u_c,
			const typename vector_type::value_type& u_itf0,
			const typename vector_type::value_type& u_itf1
		) = 0;

		/// function returning the defect derivatives for a constrained node
		/**
		 * The defect is: f(constrainingDoF, interfaceDoF0, interfaceDoF1),
		 * where f is any inter-/extrapolating function that is suited to calculate the value
		 * for the constrainedDoF from the values of the constraining DoF and the two interface
		 * DoFs. Of course, f should be differentiable.
		 * interfaceDoF0 is on the same side as the constrainingDoF,
		 * interfaceDoF1 is on the side of the constrainedDoF.
		 *
		 * @param dd		defect derivs at constrained vertex (in the same order as input variables);
		 * 					this is the output
		 * @param u_c		solution at corresponding constrainer vertex
		 * @param u_itf0	solution at interface vertex on constraining side
		 * @param u_itf1	solution at interface vertex on constrained side
		 */
		virtual void constraintValueDerivs
		(	typename vector_type::value_type dd[3],
			const typename vector_type::value_type& u_c,
			const typename vector_type::value_type& u_itf0,
			const typename vector_type::value_type& u_itf1
		) = 0;

	public:
		/// struct containing information about the dependencies of a constrained index
		struct ConstraintInfo
		{
			public:

			/// default constructor
			ConstraintInfo()
				: constrgInd(0), fct(0), side(0) {};

			/// custom constructor
			ConstraintInfo(size_t _constrgInd, size_t _fct, size_t _side)
				: constrgInd(_constrgInd), fct(_fct), side(_side) {};

			/// the constraining index
			size_t constrgInd;

			/// the function this constrained index belongs to (relative to fcts given to this class)
			size_t fct;

			/// which side of the interface the constrained index belongs to
			/// (0 for high-dim side, 1 for 1d side)
			size_t side;
		};

	private:
		/// for every constrained vertex: finds the corresponding constrainer and
		/// fills the constrainer map with the pair of corresponding indices
		void fill_constraint_map();

	protected:
		/// constrained functions
		std::vector<size_t> m_vFct;
		std::vector<std::string> m_vsFct;

		/// subset index of extension domain
		int m_ssiExt;
		std::string m_sssiExt;

		/// subset index of 2d interface node
		int m_ssiIN;
		std::string m_sssiIN;

		/// subset indices of constrained vertices
		// m_ssi[0] MUST contain the ssi of the constrained vertices for the high-dim. end;
		// m_ssi[1] the ssi of the constrained vertex for the 1d end
		int m_ssi[2];
		std::string m_sssi[2];

		/// algebraic indices for the interface nodes (and all functions)
		std::vector<size_t> m_algInd[2];

		/// solution values at the interface nodes (for all functions)
		std::vector<number> m_intf_val[2];

		/// algebraic indices of constrained nodes and their respective constrainers
		std::map<size_t, ConstraintInfo> m_constraintMap;
};


template <typename TDomain, typename TAlgebra>
class AdditiveInterface1DFV1 : public IInterface1DFV1<TDomain, TAlgebra>
{
	public:
		///	type of algebra vector
		typedef typename TAlgebra::vector_type vector_type;

	public:
		///	constructor
		AdditiveInterface1DFV1(const char* fcts, const char* high_dim_subset,
							   const char* one_dim_subset, const char* one_dim_intfNode,
							   const char* extension_subset);

		/// destructor
		virtual ~AdditiveInterface1DFV1() {};

		// inherited from IInterface1DFV1

		/// \copydoc IInterface1DFV1::constraintValue()
		virtual void constraintValue
		(	typename vector_type::value_type& d,
			const typename vector_type::value_type& u_c,
			const typename vector_type::value_type& u_itf0,
			const typename vector_type::value_type& u_itf1
		);

		/// \copydoc IInterface1DFV1::constraintValueDerivs()
		virtual void constraintValueDerivs
		(	typename vector_type::value_type dd[3],
			const typename vector_type::value_type& u_c,
			const typename vector_type::value_type& u_itf0,
			const typename vector_type::value_type& u_itf1
		);
};


template <typename TDomain, typename TAlgebra>
class MultiplicativeInterface1DFV1: public IInterface1DFV1<TDomain, TAlgebra>
{
	public:
		///	type of algebra vector
		typedef typename TAlgebra::vector_type vector_type;

	public:
		///	constructor
		MultiplicativeInterface1DFV1(const char* fcts, const char* high_dim_subset,
				   const char* two_dim_subset, const char* one_dim_intfNode,
				   const char* extension_subset);

		/// destructor
		virtual ~MultiplicativeInterface1DFV1() {};

		// inherited from IInterface1DFV1

		/// \copydoc IInterface1DFV1::constraintValue()
		virtual void constraintValue
		(	typename vector_type::value_type& d,
			const typename vector_type::value_type& u_c,
			const typename vector_type::value_type& u_itf0,
			const typename vector_type::value_type& u_itf1
		);

		/// \copydoc IInterface1DFV1::constraintValueDerivs()
		virtual void constraintValueDerivs
		(	typename vector_type::value_type dd[3],
			const typename vector_type::value_type& u_c,
			const typename vector_type::value_type& u_itf0,
			const typename vector_type::value_type& u_itf1
		);
};

} // namespace nernst_planck
} // namespace ug


#include "interface1d_fv1_impl.h"

#endif /* INTERFACE1D_FV1_H_ */
