/*
 * copyNeighborValueConstraint.h
 *
 *  Created on: 06.06.2014
 *      Author: mbreit
 */

#ifndef UG__PLUGINS__EXPERIMENTAL__NERNST_PLANCK__COPY_NEIGHBOR_VALUE_CONSTRAINT_H
#define UG__PLUGINS__EXPERIMENTAL__NERNST_PLANCK__COPY_NEIGHBOR_VALUE_CONSTRAINT_H


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

#if 0

/// Constraint for DoFs that are required to have the same value as their neighbors.
/**
 * This class is intended to make it possible to constrain DoFs of a specified
 * subset to the same value as their neighbors from a different subset.
 *
 * These neighbors need to be uniquely defined as the vertices connected by
 * an edge and belonging to a different subset than the constrained vertices.
 *
 * This class can be used in the PNP problem to constrain the useless vertices
 * in the 1d extensions in such a way that the values on every cross-section
 * through the extension are constant (which is useful for the visualization).
 *
 * \tparam	TDomain				type of Domain
 * \tparam	TAlgebra			type of Algebra
 *
 * \date 07.01.2015
 * \author mbreit
 */


template <typename TDomain, typename TAlgebra>
class CopyNeighborValueConstraint: public IDomainConstraint<TDomain, TAlgebra>
{
	public:
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
		 */
		CopyNeighborValueConstraint(const char* fcts, const char* constrained);

		/// destructor
		virtual ~CopyNeighborValueConstraint();


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

		///	adapts jacobian to enforce constraints
		virtual void adjust_jacobian
		(	matrix_type& J,
			const vector_type& u,
			ConstSmartPtr<DoFDistribution> dd,
			int type,
			number time = 0.0,
			ConstSmartPtr<VectorTimeSeries<vector_type> > vSol = SPNULL,
			const number s_a0 = 1.0
		);

		///	adapts defect to enforce constraints
		virtual void adjust_defect
		(	vector_type& d,
			const vector_type& u,
			ConstSmartPtr<DoFDistribution> dd,
			int type,
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
			int type,
			number time = 0.0
		);

		///	adapts a rhs to enforce constraints
		virtual void adjust_rhs
		(	vector_type& rhs,
			const vector_type& u,
			ConstSmartPtr<DoFDistribution> dd,
			int type,
			number time = 0.0
		);

		///	sets the constraints in a solution vector
		virtual void adjust_solution
		(	vector_type& u,
			ConstSmartPtr<DoFDistribution> dd,
			int type,
			number time = 0.0
		);

	protected:
		/// called when the approximation space has changed
		void approximation_space_changed();

	private:
		/// for every constrained vertex: finds the corresponding constrainer and
		/// fills the constrainer map with the pair of corresponding indices
		void fill_constraint_map();

	protected:
		/// constrained functions
		std::vector<size_t> m_vFct;
		std::vector<std::string> m_vsFct;
		std::map<size_t, size_t> m_fctIndexMapper;

		/// subset indices of constrained vertices
		int m_siConstr;
		std::string m_ssiConstr;

		/// algebraic indices of constrained nodes and their respective constrainers
		std::map<size_t, size_t> m_constraintMap;
};

#endif

template <typename TDomain, typename TAlgebra>
class Domain1dSolutionAdjuster
{
	public:
		static const int worldDim = TDomain::dim;

	public:
		Domain1dSolutionAdjuster() : m_sortDir(0.0)
		{m_sortDir[0] = 1.0;}

		void add_constrained_subset(const char* ss) {m_vConstrdNames.push_back(ss);}
		void add_constrainer_subset(const char* ss) {m_vConstrgNames.push_back(ss);}

		void set_sorting_direction(std::vector<number> vDir);

		void adjust_solution(SmartPtr<GridFunction<TDomain, TAlgebra> > u);

	private:
		template <typename TBaseElem> void collect_constrainers(SmartPtr<GridFunction<TDomain, TAlgebra> > u);
		template <typename TBaseElem> void adjust_constrained(SmartPtr<GridFunction<TDomain, TAlgebra> > u);

	protected:
		std::vector<std::string> m_vConstrdNames;
		std::vector<std::string> m_vConstrgNames;

		std::vector<int> m_vConstrdSI;
		std::vector<int> m_vConstrgSI;

		MathVector<worldDim> m_sortDir;

		struct DataPoint
		{
			public:
				DataPoint(number coord, number val)
				: m_coord(coord), m_val(val) {}

				struct CompareFunction
				{
					public:
						bool operator() (const DataPoint& a, const DataPoint& b)
						{return a.m_coord < b.m_coord;}
				};
				friend struct CompareFunction;

				number m_coord;
				number m_val;
		};

		std::vector<std::vector<DataPoint> > m_vDataPoints;
};


} // namespace nernst_planck
} // namespace ug


#include "copy_neighbor_value_constraint_impl.h"

#endif // UG__PLUGINS__EXPERIMENTAL__NERNST_PLANCK__COPY_NEIGHBOR_VALUE_CONSTRAINT_H
