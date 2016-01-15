/*
 * constrained_ilu.h
 *
 *  Created on: 23.09.2014
 *      Author: mbreit
 */

#ifndef UG__PLUGINS__EXPERIMENTAL__NERNST_PLANCK__CONSTRAINED_ILU_H
#define UG__PLUGINS__EXPERIMENTAL__NERNST_PLANCK__CONSTRAINED_ILU_H

#include "common/util/smart_pointer.h"
#include "lib_algebra/operator/interface/preconditioner.h"
#include "lib_disc/spatial_disc/constraints/constraint_interface.h"

#ifdef UG_PARALLEL
	#include "pcl/pcl_util.h"
	#include "lib_algebra/parallelization/parallelization_util.h"
	#include "lib_algebra/parallelization/parallel_matrix_overlap_impl.h"
#endif
#include "lib_algebra/algebra_common/permutation_util.h"
#include "interface1d_fv.h"

#include <vector>

namespace ug{
namespace nernst_planck{


///	constrained preconditioners
template <typename TDomain, typename TAlgebra, template<typename> class TPrecond>
class ConstrainedPreconditioner : public TPrecond<TAlgebra>
{
	public:
		typedef typename TAlgebra::vector_type vector_type; ///< vector type
		typedef typename TAlgebra::matrix_type matrix_type; ///< matrix type
		typedef typename IPreconditioner<TAlgebra>::matrix_operator_type matrix_operator_type; ///< matrix operator type
		typedef GridFunction<TDomain, TAlgebra> gf_type;
		typedef TPrecond<TAlgebra> base_type; ///< base type

	public:
		///	constructor
		ConstrainedPreconditioner(ConstSmartPtr<ApproximationSpace<TDomain> > approx)
		: base_type(), m_spDD(approx->dof_distribution(GridLevel())), m_time(0.0)
		{}

		/// clone constructor
		ConstrainedPreconditioner(const ConstrainedPreconditioner<TDomain, TAlgebra, TPrecond> &parent)
		: base_type(parent), m_spDD(parent.m_spDD), m_time(parent.m_time)
		{}

		///	clone
		virtual SmartPtr<ILinearIterator<vector_type> > clone()
		{
			return make_sp(new ConstrainedPreconditioner<TDomain, TAlgebra, TPrecond>(*this));
		}

		///	destructor
		virtual ~ConstrainedPreconditioner(){}

	protected:
		/// @copydoc ILinearIterator::name()
		virtual const char* name() const
		{
			return base_type::name();
		}

		/// @copydoc IPreconditioner::preprocess()
		virtual bool preprocess(SmartPtr<matrix_operator_type> pOp)
		{
			return base_type::preprocess(pOp);
		}

		/// @copydoc IPreconditioner::step()
		virtual bool step(SmartPtr<matrix_operator_type> pOp, vector_type& c, const vector_type& d)
		{
			// perform normal step
			base_type::step(pOp, c, d);

			// apply constraints
			try
			{
				gf_type* gf = dynamic_cast<gf_type*>(&c);
				if (gf && gf->dof_distribution().valid())
				{
					for (size_t i = 0; i < m_vConstraint.size(); i++)
						m_vConstraint[i]->adjust_solution(c, gf->dof_distribution(), m_time);
				}
				else if (m_spDD.valid())
				{
					UG_LOGN("No DoF distribution found. Using surface DoF distribution.");
					for (size_t i = 0; i < m_vConstraint.size(); i++)
						m_vConstraint[i]->adjust_solution(c, m_spDD, m_time);
				}
				else
					UG_THROW("DoF distribution is not valid.")
			} UG_CATCH_THROW(" Cannot adjust solution.");


		//	we're done
			return true;
		}

		/// @copydoc IPreconditioner::postprocess()
		virtual bool postprocess()
		{
			return base_type::postprocess();
		}

	public:
		/// @copydoc ILinearIterator::supports_parallel()
		virtual bool supports_parallel() const
		{
			return base_type::supports_parallel();
		}

		/// setter for time
		void set_time(number time)
		{
			m_time = time;
		};

		/// adding a constraint
		void add_constraint(SmartPtr<IDomainConstraint<TDomain, TAlgebra> > constraint)
		{
			m_vConstraint.push_back(constraint);
		}

	protected:
		///	storage for factorization
		std::vector<SmartPtr<IDomainConstraint<TDomain, TAlgebra> > > m_vConstraint;
		ConstSmartPtr<DoFDistribution> m_spDD;
		number m_time;
};

} // end namespace nernst_planck
} // end namespace ug


#endif // UG__PLUGINS__EXPERIMENTAL__NERNST_PLANCK__CONSTRAINED_ILU_H
