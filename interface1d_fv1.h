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
#include "lib_disc/dof_manager/orientation.h"	// for MapLagrangeMultiIndexTriangle etc.


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
 * This class provides an implementation for FV (of arbitrary order) discretizations
 * of both sides of the interface.
 *
 * ------------------------
 * Geometry specifications
 * ------------------------
 * Suppose you have the following situation: left: 2D, right: 1D
 *
 *     --- o --- o
 *         |  /  |
 *     --- o --- I       I --- o ---
 *         |  /  |
 *     --- o --- o
 *
 * The two nodes marked by a capital I will be referred to as "interface" nodes/vertices.
 *
 * Then you have to introduce constrained nodes (x) on both sides that will simulate
 * the behavior of a non-existing extension of the respective subdomain:
 *
 *     --- o --- o --- x
 *         |  /  |     |
 *     --- o  -  I --- I --- o ---
 *         |  /  |     |
 *     --- o --- o --- x
 *
 * Each constrained vertex is connected to exactly one unconstrained ("normal") vertex
 * in the high-dimensional subdomain which is referred to as its "constrainer" or
 * "constraining vertex".
 *
 * The constrained vertices, in conjunction with their respective constrainers,
 * form new high-dimensional elements (quadrilaterals in 2D, prisms or hexahedra in 3D).
 * A discretization on the high-dimensional side with an equivalent on the 1D side that
 * is to be coupled by means of this interface will have to be extended to those new elements.
 *
 * A practical point: As, for the time being, there is no handling of lower-dimensional
 * elements not part of a full-dimensional one in the process of distributing the geometry
 * when simulating in parallel, the one-dimensional elements will typically be extended to
 * form full-dimensional elements by some nodes of a subset "useless" which will then be
 * set to an arbitrary value using a constraint (e.g. Dirichlet boundary), like this:
 *
 *	   --- o --- o --- x  u --- u --- u ---
 *         |  /  |     | /    /	    /
 *     --- o  -  I --- I --- o --- o ---
 *         |  /  |     |
 *     --- o --- o --- x
 *
 * Subsets:
 * ---------
 * The constrained manifold (vertices and edges (and faces in 3D)) of the high-dimensional
 * side must be (exclusively) assigned to a subset that needs to be specified in the
 * interface constructor. The same is necessary for the interface node on both sides of
 * the interface which both need their own subset.
 *
 * -------------
 * How it works
 * -------------
 *
 * (1)
 * The coupling is realized by constraints (this is why this base class is derived
 * from IDomainConstraint). As the respective discretizations on both sides of the
 * interface which are to be coupled are extended to the respective new elements
 * nothing has to be assembled. What has to be done is assigning meaningful values
 * to the constrained vertices.
 *
 * The constraint takes the general form:   u* = f(u,u1,u2)
 * where u* is the value of a function at a constrained vertex, u the respective
 * value at its constrainer, u1 is the value on the full-D interface node,
 * u2 the value on the 1D interface node.
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
 * constraintValueDerivs must calculate the partial derivatives wrt.
 * \f$ u \f$, \f$ u_1 \f$, \f$ u_2 \f$ thereof.
 * Cf. documentation of the virtual methods.
 *
 * Two default implementations are realized in the derived classes
 * AdditiveIInterface1D and MultiplicativeIInterface1D with:
 *
 * \f[
 *		f(u,u_1,u_2) = u + (u_2-u_1);
 * \f]
 * \f[
 * 		f(u,u_1,u_2) = u \cdot (u_2/u_1),
 * \f]
 * respectively.
 *
 * (2)
 * The Interface base class implements a nested class Interface1DMapper.
 * In order for the interface to work properly, an instance of this class will have
 * to be created in the lua script via "interfaceMapper = Interface1DMapper(domainDisc)"
 * and all the interfaces of the simulation need to be passed to it using its method
 * add_interface().
 * This will make sure the local assemblings in the constrained nodes will be added
 * to the Jacobian and defect components of the 1D interface node. This is necessary for
 * conservation of the flowing quantity.
 *
 * (3)
 * The interface must not be divided by parallel distribution!
 * In order to make sure this does not happen, you will have to do two things:
 *
 *   (a) Use "metisReweigh" as distribution method and use ProtectSubsetPartitionWeighting
 *   	 as weighting function. Set high weights on the constrained subset as well as the
 *   	 two interface node subsets.
 *   (b) For this to work, there needs to be a proper full-D element connecting the full-D
 *   	 subdomain with the 1D subdomain. This can be achieved by creating an element that
 *   	 consists of the two interface nodes and the nearest node of the "useless" subset
 *   	 (see graphical representation above).
 *
 *
 * @note This implementation is only designed to work for Finite Volume
 * 		 discretizations and will probably fail miserably if used in
 * 		 another setting.
 *
 * \tparam	TDomain				type of Domain
 * \tparam	TAlgebra			type of Algebra
 *
 * \date 08.07.2014
 * \author mbreit
 */


template <typename TDomain, typename TAlgebra>
class IInterface1D: public IDomainConstraint<TDomain, TAlgebra>
{
	public:
		// nested class for mapping local "shadow" DoFs to global 1D interface DoF
		class Interface1DMapper : public ILocalToGlobalMapper<TAlgebra>
		{
			public:
				///	algebra type
				typedef TAlgebra algebra_type;

				///	type of algebra matrix
				typedef typename algebra_type::matrix_type matrix_type;

				///	type of algebra vector
				typedef typename algebra_type::vector_type vector_type;

				/// type of the interface
				typedef IInterface1D<TDomain,TAlgebra> interface_type;

			public:
				///	default constructor
				Interface1DMapper(SmartPtr<IAssemble<TAlgebra> > ass);

				///	destructor
				virtual ~Interface1DMapper() {};

				///	adds a local vector to the global one
				virtual void add_local_vec_to_global(vector_type& vec, const LocalVector& lvec, ConstSmartPtr<DoFDistribution> dd);

				///	adds a local matrix to the global one
				virtual void add_local_mat_to_global(matrix_type& mat, const LocalMatrix& lmat, ConstSmartPtr<DoFDistribution> dd);

				///	modifies local solution vector for adapted defect computation (we do not require this here)
				virtual void modify_LocalSol(LocalVector& vecMod, const LocalVector& lvec, ConstSmartPtr<DoFDistribution> dd) {};

				/// adds an interface to the mapper
				void add_interface(SmartPtr<interface_type> intf);
			private:
				std::vector<SmartPtr<interface_type> > m_vspInterface;
		};


	// make private members accessible for mapper class (esp. for accessing index values)
	friend class Interface1DMapper;

	public:
		/// own type
		typedef IInterface1D<TDomain, TAlgebra> this_type;

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
		 * @param constrained		the name of the constrained subset
		 * @param high_dim_intfNode	the name of the subset of the high-dim interface node
		 * @param one_dim_intfNode	the name of the subset of the one-dim interface node
		 */
		IInterface1D(const char* fcts, const char* constrained,
						const char* high_dim_intfNode, const char* one_dim_intfNode);

		/// destructor
		virtual ~IInterface1D();


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
				: constrgInd(0), fct(0){};

			/// custom constructor
			ConstraintInfo(size_t _constrgInd, size_t _fct)
				: constrgInd(_constrgInd), fct(_fct) {};

			/// the constraining index
			size_t constrgInd;

			/// the function this constrained index belongs to (relative to fcts given to this class)
			size_t fct;
		};

	private:

		/// computes offsets for algebra indices
		/**
		 *  Computation gives offsets with respect to an order implicitly defined by the
		 *  target element descriptor (consult in-code documentation in orientation.cpp).
		 *
		 *  This functionality is needed when the constraintMap is being filled in the case
		 *  of FV of order p > 2.
		 *
		 *  The element passed as third argument might be rotated and mirrored with respect
		 *  to the target descriptor. Offsets are computed such that elem's i-th (inner) algebra
		 *  index w.r.t. the order prescribed by the target descriptor can be obtained by
		 *
		 *      std::vector<DoFIndex> ind;
		 *      dd->inner_dof_indices(elem, fct, ind, false);
		 *      DoFIndex dof_i = ind[offsets[i]];
		 *
		 * @note Why usage of structs and Dummy template?
		 *		 Implementation is realized using that partial specialization of nested template
		 *		 structs is allowed (in contrast to specialization of templated methods) in not
		 *		 fully specialized template classes.
		 */
		/// @{
		template <typename TElem, typename TElemDesc, typename TDummy = void>
		struct OrientationOffset
		{
			OrientationOffset
			(
				std::vector<size_t>& vOrientOffset,
				TElemDesc& target_desc,
				TElem* constrg,
				size_t p
			);
		};
		template <typename TDummy>
		struct OrientationOffset<Vertex*, Vertex, TDummy>
		{
			OrientationOffset
			(
				std::vector<size_t>& vOrientOffset,
				Vertex& target_desc,
				Vertex* constrg,
				size_t p
			);
		};
		template <typename TDummy>
		struct OrientationOffset<Edge, EdgeDescriptor, TDummy>
		{
			OrientationOffset
			(
				std::vector<size_t>& vOrientOffset,
				EdgeDescriptor& target_desc,
				Edge* constrg,
				size_t p
			);
		};
		template <typename TDummy>
		struct OrientationOffset<Face, FaceDescriptor, TDummy>
		{
			OrientationOffset
			(
				std::vector<size_t>& vOrientOffset,
				FaceDescriptor& target_desc,
				Face* constrg,
				size_t p
			);
		};
		/// @}

		/// computes the constrainer of a constrained element
		/**
		 *  The element can be any type of element contained in the constrained subset.
		 *  Its constrainer is defined as the unique element of the same type and same
		 *  dimension (dim) such that there is a base element of dimension (dim+1) not
		 *  contained in the constrained subsets which connects constrained and constraining
		 *  element.
		 *
		 *  This method is used when the constraintMap is being filled.
		 *
		 * @note Why usage of structs and Dummy template?
		 *		 Implementation is realized using that partial specialization of nested template
		 *		 structs is allowed (in contrast to specialization of templated methods) in not
		 *		 fully specialized template classes.
		 */
		/// @{
		template <typename TElem, typename TElemDesc, typename TContainingElem, typename TDummy = void>
		struct GetConstrainer
		{
			GetConstrainer(IInterface1D<TDomain, TAlgebra>* const intf, TElem* constrd, TElem** constrg_out);
		};
		template <typename TDummy>
		struct GetConstrainer<Vertex, Vertex, Edge, TDummy>
		{
			GetConstrainer(IInterface1D<TDomain, TAlgebra>* const intf, Vertex* constrd, Vertex** constrg_out);
		};
		/// @}
		template <typename TElem, typename TElemDesc, typename TContainingElem, typename TDummy>
		friend struct GetConstrainer;

		/// computes the target descriptor for a constrained element
		/**
		 * The target descriptor is a descriptor for the constrainer (as defined for
		 * GetConstrainer) which has the same orientation and whose 0-th vertex is the
		 * constrainer of the constrained element's 0-th vertex.
		 * By this construction, its inner DoF indices have the same numbering as those
		 * of the constrained element and it can thus be used to compute offsets using
		 * OrientationOffset.
		 *
		 * As this is never necessary for vertices, the specialized template for vertices
		 * does not do anything.
		 *
		 * @note Why usage of structs and Dummy template?
		 *		 Implementation is realized using that partial specialization of nested template
		 *		 structs is allowed (in contrast to specialization of templated methods) in not
		 *		 fully specialized template classes.
		 */
		/// @{
		template <typename TElem, typename TElemDesc, typename TDummy = void>
		struct TargetDescriptor
		{
			TargetDescriptor(IInterface1D<TDomain, TAlgebra>* const intf, TElem* constrd, TElemDesc& desc);
		};
		template <typename TDummy>
		struct TargetDescriptor<Vertex, Vertex, TDummy>
		{
			TargetDescriptor(IInterface1D<TDomain, TAlgebra>* const intf, Vertex* constrd, Vertex& desc);
		};
		template <typename TDummy>
		struct TargetDescriptor<Edge, EdgeDescriptor, TDummy>
		{
			TargetDescriptor(IInterface1D<TDomain, TAlgebra>* const intf, Edge* constrd, EdgeDescriptor& desc);
		};
		/// @}
		template <typename TElem, typename TElemDesc, typename TDummy>
		friend struct TargetDescriptor;

		/// fill the constraintMap with inner indices of a specific element type
		template <typename TElem, typename TElemDesc, typename TContainingElem>
		void fill_constraint_map();

		/// fill the constraintMap
		/**
		 *  After a call to this method, the constraintMap will contain any algebra index
		 *  that belongs to the constrained subset as a key to an info struct containing
		 *  (a) its constraining algebra index,
		 *  (b) the function index this index belongs to.
		 *
		 *  The map is filled per element type containg indices using
		 *  template <typename TElem, typename TElemDesc, typename TContainingElem> void fill_constraint_map().
		 */
		void fill_constraint_map();

	protected:
		/// constrained functions
		std::vector<size_t> m_vFct;
		std::vector<std::string> m_vsFct;
		std::map<size_t, size_t> m_fctIndexMapper;

		/*
		/// subset index of extension domain
		int m_ssiExt;
		std::string m_sssiExt;
		*/

		/*
		/// subset index of 2d interface node
		int m_ssiIN;
		std::string m_sssiIN;
		*/

		/// subset indices of constrained vertices
		int m_siConstr;
		std::string m_ssiConstr;

		// [0] MUST contain the subset index of the interface node for the high-dim. end;
		// [1] the subset index of the interface node for the 1d end
		/// subset indices of interface nodes
		int m_siIntf[2];
		std::string m_ssiIntf[2];

		/// algebraic indices for the interface nodes (and all functions)
		// [0] MUST contain the index of the interface node for the high-dim. end;
		// [1] the index of the interface node for the 1d end
		std::vector<size_t> m_algInd[2];

		/// solution values at the interface nodes (for all functions)
		// [0] MUST contain the values at the interface node for the high-dim. end;
		// [1] the values at the interface node for the 1d end
		std::vector<number> m_intf_val[2];

		/// algebraic indices of constrained nodes and their respective constrainers
		std::map<size_t, ConstraintInfo> m_constraintMap;
};


template <typename TDomain, typename TAlgebra>
class AdditiveInterface1D : public IInterface1D<TDomain, TAlgebra>
{
	public:
		///	type of algebra vector
		typedef typename TAlgebra::vector_type vector_type;

	public:
		///	constructor
		AdditiveInterface1D
		(
			const char* fcts,
			const char* constrained,
			const char* high_dim_intfNode,
			const char* one_dim_intfNode
		)
	 	: IInterface1D<TDomain, TAlgebra>
		  (fcts, constrained, high_dim_intfNode, one_dim_intfNode)
		{}

		/// destructor
		virtual ~AdditiveInterface1D() {};

		// inherited from IInterface1D

		/// \copydoc IInterface1D::constraintValue()
		virtual void constraintValue
		(	typename vector_type::value_type& d,
			const typename vector_type::value_type& u_c,
			const typename vector_type::value_type& u_itf0,
			const typename vector_type::value_type& u_itf1
		)
		{
			Proxy<typename vector_type::value_type>::constraintValue(d, u_c, u_itf0, u_itf1);
		}

		/// \copydoc IInterface1D::constraintValueDerivs()
		virtual void constraintValueDerivs
		(	typename vector_type::value_type dd[3],
			const typename vector_type::value_type& u_c,
			const typename vector_type::value_type& u_itf0,
			const typename vector_type::value_type& u_itf1
		)
		{
			Proxy<typename vector_type::value_type>::constraintValueDerivs(dd, u_c, u_itf0, u_itf1);
		}


		// The following proxy implementation is used to distinguish the cases:
		// vector_type::value_type = number (CPU1)
		// vector_type::value_type = ...	(other CPUs)
		// As the standard explicitly forbids the specialization of templated class methods
		// in a not fully specialized template class, it is not possible to use
		// template <typename TValueType> virtual void constraintValue(...)
		// and specialize it for TValueType=number.
		// However, it is allowed to _partially_ specialize template classes, this
		// is where the proxy comes in.
	private:
		// _partial_ specialization is allowed (for that: existence of dummy)
		template <typename TValue, typename DUMMY = void>
		struct Proxy
		{
			static void constraintValue
			(
				TValue& d,
				const TValue& u_c,
				const TValue& u_itf0,
				const TValue& u_itf1
			)
			{
				UG_THROW("This class is not implemented for a different algebra type than CPU1.");
			}

			static void constraintValueDerivs
			(
				TValue dd[3],
				const TValue& u_c,
				const TValue& u_itf0,
				const TValue& u_itf1
			)
			{
				UG_THROW("This class is not implemented for a different algebra type than CPU1.");
			}
		};

		template <typename DUMMY>
		struct Proxy<number, DUMMY>
		{
			static void constraintValue
			(
				number& d,
				const number& u_c,
				const number& u_itf0,
				const number& u_itf1
			)
			{
				d = u_c + (u_itf1 - u_itf0);
			}

			static void constraintValueDerivs
			(
				number dd[3],
				const number& u_c,
				const number& u_itf0,
				const number& u_itf1
			)
			{
				dd[0] =  1.0;
				dd[1] = -1.0;
				dd[2] =  1.0;
			}
		};
};


template <typename TDomain, typename TAlgebra>
class MultiplicativeInterface1D: public IInterface1D<TDomain, TAlgebra>
{
	public:
		///	type of algebra vector
		typedef typename TAlgebra::vector_type vector_type;

	public:
		///	constructor
		MultiplicativeInterface1D
		(
			const char* fcts,
			const char* constrained,
			const char* high_dim_intfNode,
			const char* one_dim_intfNode
		)
		: IInterface1D<TDomain, TAlgebra>
		  (fcts, constrained, high_dim_intfNode, one_dim_intfNode)
		{};

		/// destructor
		virtual ~MultiplicativeInterface1D() {};

		// inherited from IInterface1D

		/// \copydoc IInterface1D::constraintValue()
		virtual void constraintValue
		(
			typename vector_type::value_type& d,
			const typename vector_type::value_type& u_c,
			const typename vector_type::value_type& u_itf0,
			const typename vector_type::value_type& u_itf1
		)
		{
			Proxy<typename vector_type::value_type>::constraintValue(d, u_c, u_itf0, u_itf1);
		}

		/// \copydoc IInterface1D::constraintValueDerivs()
		virtual void constraintValueDerivs
		(
			typename vector_type::value_type dd[3],
			const typename vector_type::value_type& u_c,
			const typename vector_type::value_type& u_itf0,
			const typename vector_type::value_type& u_itf1
		)
		{
			Proxy<typename vector_type::value_type>::constraintValueDerivs(dd, u_c, u_itf0, u_itf1);
		}

		// The following proxy implementation is used to distinguish the cases:
		// vector_type::value_type = number (CPU1)
		// vector_type::value_type = ...	(other CPUs)
		// As the standard explicitly forbids the specialization of templated class methods
		// in a not fully specialized template class, it is not possible to use
		// template <typename TValueType> virtual void constraintValue(...)
		// and specialize it for TValueType=number.
		// However, it is allowed to _partially_ specialize template classes, this
		// is where the proxy comes in.
	private:
		// _partial_ specialization is allowed (for that: existence of dummy)
		template <typename TValue, typename DUMMY = void>
		struct Proxy
		{
			static void constraintValue
			(
				TValue& d,
				const TValue& u_c,
				const TValue& u_itf0,
				const TValue& u_itf1
			)
			{
				UG_THROW("This class is not implemented for a different algebra type than CPU1.");
			}

			static void constraintValueDerivs
			(
				TValue dd[3],
				const TValue& u_c,
				const TValue& u_itf0,
				const TValue& u_itf1
			)
			{
				UG_THROW("This class is not implemented for a different algebra type than CPU1.");
			}
		};

		template <typename DUMMY>
		struct Proxy<number, DUMMY>
		{
			static void constraintValue
			(
				number& d,
				const number& u_c,
				const number& u_itf0,
				const number& u_itf1
			)
			{
				if (std::fabs(u_itf0) < 2 * std::numeric_limits<number>::denorm_min())
					{UG_THROW("Denominator practically zero (" << u_itf0 << ").");}

				d = u_c * (u_itf1 / u_itf0);
			}

			static void constraintValueDerivs
			(
				number dd[3],
				const number& u_c,
				const number& u_itf0,
				const number& u_itf1
			)
			{
				if (std::fabs(u_itf0) < 2 * std::numeric_limits<typename vector_type::value_type>::denorm_min())
					{UG_THROW("Denominator practically zero.");}

				dd[0] =  u_itf1 / u_itf0;
				dd[1] = -u_c * u_itf1 / (u_itf0*u_itf0);
				dd[2] =  u_c / u_itf0;
			}
		};
};

} // namespace nernst_planck
} // namespace ug


#include "interface1d_fv1_impl.h"

#endif /* INTERFACE1D_FV1_H_ */
