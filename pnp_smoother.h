/*
 * pnp_smoother.h
 *
 *  Created on: 20.06.2017
 *      Author: mbreit
 */

#ifndef UG__PLUGINS__NERNST_PLANCK__PNP_SMOOTHER_H_
#define UG__PLUGINS__NERNST_PLANCK__PNP_SMOOTHER_H_

#include <cstddef>                                                             // for size_t
#include <utility>                                                             // for pair
#include <vector>                                                              // for vector

#include "common/types.h"                                                      // for number
#include "common/util/message_hub.h"                                           // for MessageHub, MessageHub::SPCallbackId
#include "common/util/smart_pointer.h"                                         // for SmartPtr, SPNULL
#include "lib_algebra/cpu_algebra/sparsematrix.h"                              // for SparseMatrix, SparseMatrix<>::value_type
#include "lib_algebra/cpu_algebra/vector.h"                                    // for Vector
#include "lib_algebra/cpu_algebra_types.h"                                     // for CPUAlgebra, CPUVariableBlockAlgebra::...
#include "lib_algebra/operator/interface/linear_iterator.h"                    // for ILinearIterator
#include "lib_algebra/operator/interface/matrix_operator.h"                    // for MatrixOperator
#include "lib_algebra/operator/preconditioner/gauss_seidel.h"                  // for GaussSeidel
#include "lib_algebra/operator/preconditioner/ilu.h"                           // for ILU
#ifdef UG_PARALLEL
	#include "lib_algebra/parallelization/parallel_matrix.h"                   // for ParallelMatrix
	#include "lib_algebra/parallelization/parallel_vector.h"                   // for ParallelVector<>::this_type, Parallel...
#endif
#include "lib_algebra/small_algebra/small_matrix/densematrix.h"                // for DenseMatrix, operator*
#include "lib_algebra/small_algebra/small_matrix/densevector.h"                // for DenseVector
#include "lib_algebra/small_algebra/storage/variable_array.h"                  // for VariableArray2, VariableArray1
#include "lib_disc/operator/linear_operator/level_preconditioner_interface.h"  // for ILevelPreconditioner
#include "lib_grid/lib_grid_messages.h"                                        // for GridMessage_Adaption, GridMessage_Dis...
#include "lib_grid/tools/subset_handler_multi_grid.h"                          // for MGSubsetHandler


namespace ug {

// forward declarations
class DoFDistribution;
template <typename TDomain> class ApproximationSpace;

namespace nernst_planck {

// We need these two helper classes to treat the clone() method in both interfaces,
// otherwise we would have to implement two clone() methods with differing return types.
template <template <typename> class TPrecond>
class PNPSmootherInnerHelper : public TPrecond<CPUVariableBlockAlgebra>
{
	public:
		typedef ILinearIterator<CPUVariableBlockAlgebra::vector_type> clone_result_type;

	public:
		virtual ~PNPSmootherInnerHelper() {}

		virtual void my_clone(SmartPtr<clone_result_type>& res) = 0;

		virtual SmartPtr<clone_result_type> clone()
		{
			SmartPtr<clone_result_type> res(SPNULL);
			my_clone(res);
			return res;
		}
};

template <typename TAlgebra>
class PNPSmootherHelper : public ILevelPreconditioner<TAlgebra>
{
	public:
		typedef ILinearIterator<typename TAlgebra::vector_type> clone_result_type;

	public:
		virtual ~PNPSmootherHelper() {}

		virtual void my_clone(SmartPtr<clone_result_type>& res) = 0;

		virtual SmartPtr<clone_result_type> clone()
		{
			SmartPtr<clone_result_type> res(SPNULL);
			my_clone(res);
			return res;
		}
};


template <typename TDomain, typename TAlgebra, template <typename> class TPrecond>
class PNPSmoother
: public PNPSmootherHelper<TAlgebra>,
  protected PNPSmootherInnerHelper<TPrecond> // in order to use the preconditioner's protected methods
{
	public:
		typedef PNPSmoother<TDomain, TAlgebra, TPrecond> this_type;
		typedef TAlgebra algebra_type;
		typedef typename TAlgebra::vector_type vector_type;
		typedef typename TAlgebra::matrix_type matrix_type;
		typedef TPrecond<CPUVariableBlockAlgebra> precond_type;

#ifdef UG_PARALLEL
		// TODO: There is room for optimization here:
		// - Maybe do not copy matrix entries, but work on a wrapper matrix class;
		// - maybe do copy, but into another structure, using a sparse matrix
		//   of 5x5 PNP block matrices that can be stored and inverted very efficiently.
		typedef ParallelMatrix<SparseMatrix<DenseMatrix<VariableArray2<number> > > > block_matrix_type;
		typedef ParallelVector<Vector<DenseVector<VariableArray1<number> > > > block_vector_type;
#else
		typedef SparseMatrix<DenseMatrix<VariableArray2<number> > > block_matrix_type;
		typedef Vector<DenseVector<VariableArray1<number> > > block_vector_type;
#endif

	public:
		/// constructor
		PNPSmoother(SmartPtr<ApproximationSpace<TDomain> > approx);

		/// clone constructor
		PNPSmoother(const PNPSmoother& parent);

		/// destructor
		virtual ~PNPSmoother();

		// just to silence shadowing virtual warning
		using precond_type::preprocess;
		using precond_type::step;
		using precond_type::postprocess;

		/// @copydoc IPreconditioner::name
		virtual const char* name() const;

		/// @copydoc IPreconditioner::preprocess
		virtual bool preprocess(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp);

		/// @copydoc IPreconditioner::step
		virtual bool step
		(
			SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp,
			vector_type& c,
			const vector_type& d
		);

		/// @copydoc IPreconditioner::postprocess
		virtual bool postprocess();

		virtual bool supports_parallel() const;

		virtual void grid_level_has_changed();

	protected:
		/// helper struct to enable specialization-specific implementation of preprocess
		template <template <typename> class TPC, typename dummy = void>
		struct do_preprocess
		{
			/// default implementation just returns true and does nothing
			do_preprocess
			(
				bool& successOut,
				this_type* pnpSmoother,
				SmartPtr<MatrixOperator<block_matrix_type, block_vector_type> > pOp
			)
			{
				successOut = true;
			}
		};
		template <typename dummy>
		struct do_preprocess<ILU, dummy>
		{
			/// ILU implementation needs to perform the LU decomposition
			do_preprocess
			(
				bool& successOut,
				this_type* pnpSmoother,
				SmartPtr<MatrixOperator<block_matrix_type, block_vector_type> > pOp
			);
		};

		/// helper struct to enable specialization-specific implementation of step
		template <template <typename> class TPC, typename dummy = void>
		struct do_step
		{
			/// default implementation just returns true and does nothing
			do_step
			(
				bool& successOut,
				this_type* pnpSmoother,
				SmartPtr<MatrixOperator<block_matrix_type, block_vector_type> > pOp,
				block_vector_type& c,
				const block_vector_type& d
			)
			{
				successOut = true;
			}
		};
		template <typename dummy>
		struct do_step<ILU, dummy>
		{
			/// ILU implementation needs to perform the LU decomposition
			do_step
			(
				bool& successOut,
				this_type* pnpSmoother,
				SmartPtr<MatrixOperator<block_matrix_type, block_vector_type> > pOp,
				block_vector_type& c,
				const block_vector_type& d
			);
		};
		template <typename dummy>
		struct do_step<GaussSeidel, dummy>
		{
			/// ILU implementation needs to perform the LU decomposition
			do_step
			(
				bool& successOut,
				this_type* pnpSmoother,
				SmartPtr<MatrixOperator<block_matrix_type, block_vector_type> > pOp,
				block_vector_type& c,
				const block_vector_type& d
			);
		};

		/// helper struct to enable specialization-specific implementation of postprocess
		template <template <typename> class TPC, typename dummy = void>
		struct do_postprocess
		{
			/// default implementation just returns true and does nothing
			do_postprocess
			(
				bool& successOut,
				this_type* pnpSmoother
			)
			{
				successOut = true;
			}
		};

	public:
		/// add charged surface and corresponding volume subsets
		void add_charge_surface_pair(const std::string& chSsName, const std::string& volSsName);

		/// set modus operandi
		void set_method(int m);

		/**
		 * @brief Set parallelization strategy
		 * Possible values:
		 *   0   unique matrix, unique defect
		 *   1   consistent matrix, additive defect
		**/
		void set_parallelization_strategy(int ps) {m_ps = ps;}

	protected:
		/// (re-)calculate the matrix blocks
		void reinit_blocking();

		/// handling grid adaptation events
		void grid_adaptation_callback(const GridMessage_Adaption& gma);

		/// handling grid redistribution events
		void grid_distribution_callback(const GridMessage_Distribution& gmd);

	protected:
		virtual void my_clone(SmartPtr<ILinearIterator<vector_type> >& res);
		virtual void my_clone(SmartPtr<ILinearIterator<block_vector_type> >& res);

	protected:
		std::string m_name;

		MessageHub::SPCallbackId m_spGridAdaptationCallbackID;
		MessageHub::SPCallbackId m_spGridDistributionCallbackID;

		/// underlying approx space
		SmartPtr<ApproximationSpace<TDomain> > m_spApprox;

		/// multi-grid
		SmartPtr<MultiGrid> m_spMG;

		/// subset handler
		SmartPtr<MGSubsetHandler> m_spSH;

		/// dof distro
		SmartPtr<DoFDistribution> m_spDD;

		/// pairs of charge subsets (charge carrying surface subset, corresponding volume subset)
		std::vector<std::pair<int, std::vector<int> > > m_vChargeSubsetsPairs;

		/// modus operandi: 0 = non-overlapping 2- or 3-blocks orthogonal to charged surfaces
		int m_method;

		int m_ps;


		/// matrix blocks: each block holds the algebra indices given here
		std::vector<std::vector<size_t> > m_vBlocks;

		/// status of blocking info
		bool m_bBlockingNeedsReinit;

		/// matrix blocks: each index is held by the block given here
		std::vector<size_t> m_vIndexToBlock;


		/**
		 * @brief block sorting info
		 * - access new block i by m_vBlocks[m_vPerm[i]],
		 * - access new index of old block j by m_vInvPerm[j].
		 **/
		std::vector<size_t> m_vPerm;
		std::vector<size_t> m_vInvPerm;

		/// the blocked matrix
		SmartPtr<MatrixOperator<block_matrix_type, block_vector_type> > m_spBM;
};


} // namespace nernst_planck
} // namespace ug

#endif // UG__PLUGINS__NERNST_PLANCK__PNP_SMOOTHER_H_
