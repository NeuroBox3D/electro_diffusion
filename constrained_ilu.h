/*
 * constrained_ilu.h
 *
 *  Created on: 23.09.2014
 *      Author: mbreit
 */

#ifndef CONSTRAINED_ILU_H_
#define CONSTRAINED_ILU_H_

#include "common/util/smart_pointer.h"
#include "lib_algebra/operator/interface/preconditioner.h"
#include "lib_disc/spatial_disc/constraints/constraint_interface.h"

#ifdef UG_PARALLEL
	#include "pcl/pcl_util.h"
	#include "lib_algebra/parallelization/parallelization_util.h"
	#include "lib_algebra/parallelization/parallel_matrix_overlap_impl.h"
#endif
#include "lib_algebra/algebra_common/permutation_util.h"
#include "../plugins/experimental/nernst_planck/interface1d_fv1.h"

#include <vector>

namespace ug{
namespace nernst_planck{

// ILU(0) solver, i.e. static pattern ILU w/ P=P(A)
// (cf. Y Saad, Iterative methods for Sparse Linear Systems, p. 270)
template<typename Matrix_type>
bool FactorizeILU_c(Matrix_type &A)
{
	PROFILE_FUNC_GROUP("algebra ILU");
	typedef typename Matrix_type::row_iterator row_iterator;
	typedef typename Matrix_type::value_type block_type;

	// for all rows
	for(size_t i=1; i < A.num_rows(); i++)
	{
		// eliminate all entries A(i, k) with k<i with rows A(k, .) and k<i
		const row_iterator rowEnd = A.end_row(i);
		for(row_iterator it_k = A.begin_row(i);
								it_k != rowEnd && (it_k.index() < i); ++it_k)
		{
			const size_t k = it_k.index();
			block_type &a_ik = it_k.value();

			// add row k to row i by A(i, .) -= A(k,.)  A(i,k) / A(k,k)
			// so that A(i,k) is zero.
			// safe A(i,k)/A(k,k) in A(i,k)
			if(fabs(BlockNorm(A(k,k))) < 1e-15*BlockNorm(A(i,k)))
				UG_THROW("Diag is Zero for k="<<k<<", cannot factorize ILU.");

			a_ik /= A(k,k);

			row_iterator it_j = it_k;
			for(++it_j; it_j != rowEnd; ++it_j)
			{
				const size_t j = it_j.index();
				block_type& a_ij = it_j.value();
				bool bFound;
				row_iterator p = A.get_connection(k,j, bFound);
				if(bFound)
				{
					const block_type a_kj = p.value();
					a_ij -= a_ik *a_kj;
				}
			}
		}
	}

	return true;
}

// ILU(0)-beta solver, i.e.
// -> static pattern ILU w/ P=P(A)
// -> Fill-in is computed and lumped onto the diagonal
template<typename Matrix_type>
bool FactorizeILUBeta_c(Matrix_type &A, number beta)
{
	PROFILE_FUNC_GROUP("algebra ILUBeta");
	typedef typename Matrix_type::row_iterator row_iterator;
	typedef typename Matrix_type::value_type block_type;

	// for all rows i=1:n do
	for(size_t i=1; i < A.num_rows(); i++)
	{
		block_type &Aii = A(i,i);
		block_type Nii(Aii); Nii*=0.0;

		// for k=1:(i-1) do
		// eliminate all entries A(i, k) with k<i with rows A(k, .) and k<i
		const row_iterator it_iEnd = A.end_row(i);
		for(row_iterator it_ik = A.begin_row(i);
							it_ik != it_iEnd && (it_ik.index() < i); ++it_ik)
		{

			// add row k to row i by A(i, .) -=  [A(i,k) / A(k,k)] A(k,.)
			// such that A(i,k) is zero.

			// 1) Contribution to L part:
			// store A(i,k)/A(k,k) in A(i,k)
			const size_t k = it_ik.index();
			block_type &a_ik = it_ik.value();
			a_ik /= A(k,k);

			// 2) Contribution to U part:
			// compute contributions from row k for j=k:N
			const row_iterator it_kEnd = A.end_row(k);
			for (row_iterator it_kj=A.begin_row(k); it_kj != it_kEnd ;++it_kj)
			{
				const size_t j = it_kj.index();
				if (j<i) continue;  // index j belongs L part?

				// this index j belongs U part
				const block_type& a_kj = it_kj.value();

				bool aijFound;
				row_iterator pij = A.get_connection(i,j, aijFound);
				if(aijFound) {
					// entry belongs to pattern
					// -> proceed with standard elimination
					block_type a_ij = pij.value();
					a_ij -= a_ik *a_kj ;

				} else {
					// entry DOES NOT belong to pattern
					// -> we lump it onto the diagonal
					// TODO : non square matrices!!!
					Nii -=  a_ik * a_kj;
				}

			}
		}

		// add fill-in to diagonal
		AddMult(Aii, beta, Nii);
	}

	return true;
}

template<typename Matrix_type>
bool FactorizeILUSorted_c(Matrix_type &A)
{
	PROFILE_FUNC_GROUP("algebra ILU");
	typedef typename Matrix_type::row_iterator row_iterator;
	typedef typename Matrix_type::value_type block_type;

	static const number __eps = 1e-50;

	// for all rows
	for(size_t i=1; i < A.num_rows(); i++)
	{

		// eliminate all entries A(i, k) with k<i with rows A(k, .) and k<i
		for(row_iterator it_k = A.begin_row(i);
							it_k != A.end_row(i) && (it_k.index() < i); ++it_k)
		{
			const size_t k = it_k.index();
			block_type &a_ik = it_k.value();
			block_type &a_kk = A(k,k);

			// add row k to row i by A(i, .) -= A(k,.)  A(i,k) / A(k,k)
			// so that A(i,k) is zero.
			// safe A(i,k)/A(k,k) in A(i,k)
			a_ik /= a_kk;

			if(fabs(BlockNorm(A(k,k))) < __eps * BlockNorm(A(i,k)))
				UG_THROW("ILU: Blocknorm of diagonal is near-zero for k="<<k<<
				         " with eps: "<<__eps<<", ||A_kk||="<<fabs(BlockNorm(A(k,k)))
				         <<", ||A_ik||="<<BlockNorm(A(i,k)));

			typename Matrix_type::row_iterator it_ij = it_k; // of row i
			++it_ij; // skip a_ik
			typename Matrix_type::row_iterator it_kj = A.begin_row(k); // of row k

			while(it_ij != A.end_row(i) && it_kj != A.end_row(k))
			{
				if(it_ij.index() > it_kj.index())
					++it_kj;
				else if(it_ij.index() < it_kj.index())
					++it_ij;
				else
				{
					block_type &a_ij = it_ij.value();
					const block_type &a_kj = it_kj.value();
					a_ij -= a_ik * a_kj;
					++it_kj; ++it_ij;
				}
			}
		}
	}

	return true;
}


// solve x = L^-1 b
template<typename Matrix_type, typename Vector_type>
bool invert_L_c(const Matrix_type &A, Vector_type &x, const Vector_type &b)
{
	PROFILE_FUNC_GROUP("algebra ILU");
	typedef typename Matrix_type::const_row_iterator const_row_iterator;

	typename Vector_type::value_type s;
	for(size_t i=0; i < x.size(); i++)
	{
		s = b[i];
		for(const_row_iterator it = A.begin_row(i); it != A.end_row(i); ++it)
		{
			if(it.index() >= i) continue;
			MatMultAdd(s, 1.0, s, -1.0, it.value(), x[it.index()]);
		}
		x[i] = s;
	}

	return true;
}

// solve x = U^-1 * b
template<typename Matrix_type, typename Vector_type>
bool invert_U_c(const Matrix_type &A, Vector_type &x, const Vector_type &b)
{
	PROFILE_FUNC_GROUP("algebra ILU");
	typedef typename Matrix_type::const_row_iterator const_row_iterator;

	typename Vector_type::value_type s;

	static const number __eps = 1e-8;
	size_t numNearZero=0;

	// last row diagonal U entry might be close to zero with corresponding close to zero rhs
	// when solving Navier Stokes system, therefore handle separately
	if(x.size() > 0)
	{
		size_t i=x.size()-1;
		s = b[i];

		// check if diag part is significantly smaller than rhs
		// This may happen when matrix is indefinite with one eigenvalue
		// zero. In that case, the factorization on the last row is
		// nearly zero due to round-off errors. In order to allow ill-
		// scaled matrices (i.e. small matrix entries row-wise) this
		// is compared to the rhs, that is small in this case as well.
		if (BlockNorm(A(i,i)) <= __eps * BlockNorm(s))
		{
			if(numNearZero++<5)
			{	UG_LOG("ILU Warning: Near-zero diagonal entry "
					"with norm "<<BlockNorm(A(i,i))<<" in last row of U "
					" with corresponding non-near-zero rhs with norm "
					<< BlockNorm(s) << ". Setting rhs to zero.\n");
			}
			// set correction to zero
			x[i] = 0;
		} else {
			// c[i] = s/uii;
			InverseMatMult(x[i], 1.0, A(i,i), s);
		}
	}
	if(x.size() <= 1) return true;

	// handle all other rows
	for(size_t i = x.size()-2; ; --i)
	{
		s = b[i];
		for(const_row_iterator it = A.begin_row(i); it != A.end_row(i); ++it)
		{
			if(it.index() <= i) continue;
			// s -= it.value() * x[it.index()];
			MatMultAdd(s, 1.0, s, -1.0, it.value(), x[it.index()]);

		}
		// x[i] = s/A(i,i);
		InverseMatMult(x[i], 1.0, A(i,i), s);
		if(i == 0) break;
	}

	if(numNearZero>=5)
	{	UG_LOG("...\nILU Warning: " << numNearZero << " ( out of " << x.size() << ") near-zero diagonal entries in last row of U.\n");	}

	return true;
}

///	ILU / ILU(beta) preconditioner
template <typename TDomain, typename TAlgebra>
class ILUC : public IPreconditioner<TAlgebra>
{
	public:
	///	Domain type
		typedef TDomain domain_type;

	///	Algebra type
		typedef TAlgebra algebra_type;

	///	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	///	Matrix type
		typedef typename TAlgebra::matrix_type matrix_type;

	///	Matrix Operator type
		typedef typename IPreconditioner<TAlgebra>::matrix_operator_type matrix_operator_type;

	///	Base type
		typedef IPreconditioner<TAlgebra> base_type;

	protected:
		using base_type::set_debug;
		using base_type::debug_writer;
		using base_type::write_debug;

	public:
	//	Constructor
		ILUC(ConstSmartPtr<ApproximationSpace<TDomain> > approx, double beta=0.0)
			: m_beta(beta), m_bSort(false), m_spDD(approx->dof_distribution(GridLevel())), m_time(0.0) {};

	/// clone constructor
		ILUC(const ILUC<TDomain, TAlgebra> &parent)
			: base_type(parent),
			  m_beta(parent.m_beta), m_bSort(parent.m_bSort), m_spDD(parent.m_spDD), m_time(parent.m_time)
		{}

	///	Clone
		virtual SmartPtr<ILinearIterator<vector_type> > clone()
		{
			return make_sp(new ILUC<domain_type, algebra_type>(*this));
		}

	///	Destructor
		virtual ~ILUC(){}

	///	returns if parallel solving is supported
		virtual bool supports_parallel() const {return true;}

	///	set factor for \f$ ILU_{\beta} \f$
		void set_beta(double beta) {m_beta = beta;}

	/// set cuthill-mckee sort on/off
		void set_sort(bool b)
		{
			m_bSort = b;
		}


	protected:
	//	Name of preconditioner
		virtual const char* name() const {return "ILUC";}

	private:
		// cuthill-mckee sorting
		void calc_cuthill_mckee()
		{
			PROFILE_BEGIN_GROUP(ILU_ReorderCuthillMcKey, "ilu algebra");
			GetCuthillMcKeeOrder(m_ILU, m_newIndex);
			m_bSortIsIdentity = GetInversePermutation(m_newIndex, m_oldIndex);

			if(!m_bSortIsIdentity)
			{
				matrix_type mat;
				mat = m_ILU;
				SetMatrixAsPermutation(m_ILU, mat, m_newIndex);
			}
		}

	protected:

	//	Preprocess routine
		virtual bool preprocess(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp)
		{
			matrix_type &mat = *pOp;
			PROFILE_BEGIN_GROUP(ILU_preprocess, "algebra ILU");
		//	Debug output of matrices
			write_debug(mat, "ILU_BeforeMakeConsistent");

			m_ILU = mat;
#ifdef 	UG_PARALLEL
			MatAddSlaveRowsToMasterRowOverlap0(m_ILU);

		//	set zero on slaves
			std::vector<IndexLayout::Element> vIndex;
			CollectUniqueElements(vIndex,  m_ILU.layouts()->slave());
			SetDirichletRow(m_ILU, vIndex);
#endif


			if(m_bSort)
				calc_cuthill_mckee();

		//	Debug output of matrices
			write_debug(m_ILU, "ILU_BeforeFactorize");

		//	resize help vector
			m_h.resize(mat.num_cols());

		// 	Compute ILU Factorization
			if (m_beta!=0.0) FactorizeILUBeta_c(m_ILU, m_beta);
			else if(matrix_type::rows_sorted) FactorizeILUSorted_c(m_ILU);
			else FactorizeILU_c(m_ILU);
			m_ILU.defragment();

		//	Debug output of matrices
			write_debug(m_ILU, "ILU_AfterFactorize");

		//	we're done
			return true;
		}


		void applyLU(vector_type &c, const vector_type &d, vector_type &tmp)
		{
			if(!m_bSort || m_bSortIsIdentity)
			{
				// 	apply iterator: c = LU^{-1}*d
				invert_L_c(m_ILU, tmp, d); // h := L^-1 d
				invert_U_c(m_ILU, c, tmp); // c := U^-1 h = (LU)^-1 d
			}
			else
			{
				// we save one vector here by renaming
				SetVectorAsPermutation(tmp, d, m_newIndex);
				invert_L_c(m_ILU, c, tmp); // c = L^{-1} d
				invert_U_c(m_ILU, tmp, c); // tmp = (LU)^{-1} d
				SetVectorAsPermutation(c, tmp, m_oldIndex);
			}
		}

	//	Stepping routine
		virtual bool step(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp, vector_type& c, const vector_type& d)
		{
			PROFILE_BEGIN_GROUP(ILU_step, "algebra ILU");
		//	\todo: introduce damping
#ifdef UG_PARALLEL
		//	for debug output (only for application is written)
			static bool first = true;

		//	make defect unique
			SmartPtr<vector_type> spDtmp = d.clone();
			spDtmp->change_storage_type(PST_UNIQUE);

			applyLU(c, *spDtmp, m_h);

		//	Correction is always consistent
			c.set_storage_type(PST_ADDITIVE);

		//	write debug
			if(first) write_debug(c, "ILU_c");

			c.change_storage_type(PST_CONSISTENT);

		//	write debug
			if(first) {write_debug(c, "ILU_cConsistent"); first = false;}

#else
			applyLU(c, d, m_h);
#endif
		// apply constraints
			try
			{
				if (m_spDD.valid())
					for (size_t i = 0; i < m_vConstraint.size(); i++)
						m_vConstraint[i]->adjust_solution_linear(c, m_spDD, pOp, m_time);
				else
					UG_THROW("DoF distribution is not valid.")
			} UG_CATCH_THROW(" Cannot adjust solution.");


		//	we're done
			return true;
		}

	///	Postprocess routine
		virtual bool postprocess() {return true;}

	public:
	/// setter for time
		void set_time(number time) {m_time = time;};

	/// adding a constraint
		void add_constraint(SmartPtr<IInterface1DFV1<TDomain, TAlgebra> > constraint) {m_vConstraint.push_back(constraint);}

	protected:
	///	storage for factorization
		matrix_type m_ILU;

	///	help vector
		vector_type m_h;

	/// Factor for ILU-beta
		number m_beta;

	/// for cuthill-mckee reordering
		std::vector<size_t> m_newIndex, m_oldIndex;
		bool m_bSortIsIdentity;
		bool m_bSort;

	/// constraint members
		std::vector<SmartPtr<IInterface1DFV1<TDomain, TAlgebra> > > m_vConstraint;
		ConstSmartPtr<DoFDistribution> m_spDD;
		number m_time;
};

} // end namespace nernst_planck
} // end namespace ug


#endif /* CONSTRAINED_ILU_H_ */
