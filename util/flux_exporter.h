/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Markus Breit
 * Creation date: 2016-07-04
 *
 * This file is part of NeuroBox, which is based on UG4.
 *
 * NeuroBox and UG4 are free software: You can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 3
 * (as published by the Free Software Foundation) with the following additional
 * attribution requirements (according to LGPL/GPL v3 §7):
 *
 * (1) The following notice must be displayed in the appropriate legal notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 *
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 *
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating PDE based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * "Stepniewski, M., Breit, M., Hoffer, M. and Queisser, G.
 *   NeuroBox: computational mathematics in multiscale neuroscience.
 *   Computing and visualization in science (2019).
 * "Breit, M. et al. Anatomically detailed and large-scale simulations studying
 *   synapse loss and synchrony using NeuroBox. Front. Neuroanat. 10 (2016), 8"
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#ifndef UG__PLUGINS__NERNST_PLANCK__UTIL__FLUX_EXPORTER_H
#define UG__PLUGINS__NERNST_PLANCK__UTIL__FLUX_EXPORTER_H

#include <cmath>                                                           // for fabs
#include <cstddef>                                                         // for size_t
#include <algorithm>                                                       // for max
#include <map>                                                             // for operator!=
#include <set>                                                             // for set
#include <utility>                                                         // for pair
#include <vector>                                                          // for vector

#include "common/types.h"                                                  // for number
#include "common/math/math_vector_matrix/math_vector.h"                    // for MathVector
#include "common/util/smart_pointer.h"                                     // for SmartPtr, Cons...
#include "lib_algebra/cpu_algebra_types.h"                                 // for CPUAlgebra
#include "lib_disc/common/function_group.h"                                // for FunctionGroup
#include "lib_disc/common/local_algebra.h"                                 // for LocalVector
#include "lib_disc/function_spaces/grid_function.h"                        // for GridFunction
#include "lib_disc/io/vtkoutput.h"                                         // for VTKOutput
#include "lib_disc/local_finite_element/local_finite_element_id.h"         // for LFEID
#include "lib_disc/spatial_disc/constraints/continuity_constraints/p1_continuity_constraints.h"  // for SymP1Constraint...
#include "lib_grid/grid/grid_base_objects.h"                               // for GridObject (pt...
#include "lib_grid/grid_objects/grid_dim_traits.h"                         // for grid_dim_traits
#include "lib_grid/tools/subset_group.h"                                   // for SubsetGroup
#include "pcl/pcl_base.h"                                                  // for NumProcs


namespace ug {

// forward declarations
class MultiGrid;
template <int dim> class IConvectionShapes;
template <typename TAlgebra> class IConstraint;

namespace nernst_planck {


template <typename TGridFunction>
class FluxExporter
{
	public:
		typedef typename TGridFunction::algebra_type algebra_type;
		typedef typename TGridFunction::domain_type dom_type;
		static const int dim = dom_type::dim;
		typedef typename TGridFunction::domain_type::subset_handler_type sh_type;
		typedef typename grid_dim_traits<dim>::element_type elem_type;
		typedef typename TGridFunction::const_element_iterator it_type;

	public:
		FluxExporter(SmartPtr<TGridFunction> u, const std::string& cmp_name_spec, const std::string& cmp_name_pot);
		virtual ~FluxExporter() {};

		void set_diff_const(number diff_const) {m_diffConst = diff_const;}
		void set_conv_const(number conv_const) {m_convConst = conv_const;}
		void set_quad_order(number quadOrder) {m_quadOrder = quadOrder;}
		void set_upwind(SmartPtr<IConvectionShapes<dim> > shapes) {m_spConvShape = shapes;}
		void set_hanging_constraint(SmartPtr<IDomainConstraint<dom_type, algebra_type> > constr);

		void set_subsets(const std::vector<std::string>& vSubsets);
		void set_subsets(const char* cSubsets);

		void write_flux
		(
			SmartPtr<VTKOutput<dim> > vtkOutput,
			std::string filename,
			size_t step,
			number time,
			std::string fluxName,
			number scale_factor
		);

		void write_flux
		(
			SmartPtr<VTKOutput<dim> > vtkOutput,
			std::string filename,
			size_t step,
			number time
		) {write_flux(vtkOutput, filename, step, time, std::string("flux"), 1.0);}

		void write_flux
		(
			SmartPtr<VTKOutput<dim> > vtkOutput,
			std::string filename,
			size_t step,
			number time,
			number scale_factor
		) {write_flux(vtkOutput, filename, step, time, std::string("flux"), scale_factor);}

		void write_flux
		(
			SmartPtr<VTKOutput<dim> > vtkOutput,
			std::string filename,
			size_t step,
			number time,
			std::string fluxName
		) {write_flux(vtkOutput, filename, step, time, fluxName, 1.0);}

		// only FV1
		void write_box_fluxes
		(
			std::string filename,
			size_t step,
			number time,
			std::string fluxName,
			number scale_factor
		);

	protected:
		template <typename TElem>
		void add_side_subsets(MultiGrid& mg, GridObject* elem, std::set<int>& sss);

		template <bool hanging, int order>
		struct AssembleWrapper
		{
			AssembleWrapper
			(
				FluxExporter<TGridFunction>* _flEx,
				SmartPtr<GridFunction<dom_type, CPUAlgebra> > _flux,
				SmartPtr<GridFunction<dom_type, CPUAlgebra> > _vol,
				int _si
			)
			: flEx(_flEx), flux(_flux), vol(_vol), si(_si) {}

			template <typename TElem>
			void operator() (TElem&)
			{
				flEx->assemble<hanging, order, TElem>(flux, vol, si);
			}

			FluxExporter<TGridFunction>* flEx;
			SmartPtr<GridFunction<dom_type, CPUAlgebra> > flux;
			SmartPtr<GridFunction<dom_type, CPUAlgebra> > vol;
			int si;
		};
		template <bool hanging, int order>
		friend struct AssembleWrapper;

		template <bool hanging, int order, typename TElem>
		void assemble
		(
			SmartPtr<GridFunction<dom_type, CPUAlgebra> > flux,
			SmartPtr<GridFunction<dom_type, CPUAlgebra> > vol,
			int si
		);

		template <int order>
		void assemble
		(
			SmartPtr<GridFunction<dom_type, CPUAlgebra> > flux,
			SmartPtr<GridFunction<dom_type, CPUAlgebra> > vol
		);


		// PREP ELEM LOOP
		template <typename TFVGeom, typename Dummy = void>
		struct prep_elem_loop
		{
			prep_elem_loop(FluxExporter<TGridFunction>* flEx, const ReferenceObjectID roid);
		};
		template <typename TFVGeom, typename Dummy>
		friend struct prep_elem_loop;

		// specialize for FV1Geometry and HFV1Geometry
		// (little bit tricky, only these two have the signature template <class, int> class)
		template <typename TElem, int dim, template <class, int> class TFV1Geom, typename Dummy>
		struct prep_elem_loop<TFV1Geom<TElem, dim>, Dummy>
		{
			prep_elem_loop(FluxExporter<TGridFunction>* flEx, const ReferenceObjectID roid);
		};


		// PREP ELEM
		template <typename TFVGeom, typename Dummy = void>
		struct prep_elem
		{
			prep_elem
			(
				FluxExporter<TGridFunction>* flEx,
				GridObject* elem,
				const std::vector<MathVector<FluxExporter<TGridFunction>::dim> >& vCornerCoords
			);
		};
		template <typename TFVGeom, typename Dummy>
		friend struct prep_elem;

		// specialize for FV1Geometry and HFV1Geometry
		// (little bit tricky, only these two have the signature template <class, int> class)
		template <typename TElem, int dim, template <class, int> class TFV1Geom, typename Dummy>
		struct prep_elem<TFV1Geom<TElem, dim>, Dummy>
		{
			prep_elem
			(
				FluxExporter<TGridFunction>* flEx,
				GridObject* elem,
				const std::vector<MathVector<FluxExporter<TGridFunction>::dim> >& vCornerCoords
			);
		};


		// ASSEMBLE FLUX
		template <typename TFVGeom, typename Dummy = void>
		struct assemble_flux_elem
		{
			assemble_flux_elem
			(
				FluxExporter<TGridFunction>* flEx,
				LocalVector& f,
				const LocalVector& u,
				GridObject* elem,
				const std::vector<MathVector<FluxExporter<TGridFunction>::dim> >& vCornerCoords
			);
		};
		template <typename TFVGeom, typename Dummy>
		friend struct assemble_flux_elem;

		// specialize for FV1Geometry and HFV1Geometry
		// (little bit tricky, only these two have the signature template <class, int> class)
		template<typename TElem, int dim, template <class, int> class TFV1Geom, typename Dummy>
		struct assemble_flux_elem<TFV1Geom<TElem, dim>, Dummy>
		{
			assemble_flux_elem
			(
				FluxExporter<TGridFunction>* flEx,
				LocalVector& f,
				const LocalVector& u,
				GridObject* elem,
				const std::vector<MathVector<FluxExporter<TGridFunction>::dim> >& vCornerCoords
			);
		};


		// ASSEMBLE VOLUME
		template <typename TFVGeom, typename Dummy = void>
		struct assemble_vol_elem
		{
			assemble_vol_elem(const FluxExporter<TGridFunction>* flEx, LocalVector& vol);
		};
		template <typename TFVGeom, typename Dummy>
		friend struct assemble_vol_elem;

		// specialize for FV1Geometry and HFV1Geometry
		// (little bit tricky, only these two have the signature template <class, int> class)
		template<typename TElem, int dim, template <class, int> class TFV1Geom, typename Dummy>
		struct assemble_vol_elem<TFV1Geom<TElem, dim>, Dummy>
		{
			assemble_vol_elem(const FluxExporter<TGridFunction>* flEx, LocalVector& vol);
		};


		template <typename TBaseElem>
		void div_flux_by_vol_and_scale
		(
			SmartPtr<GridFunction<dom_type, CPUAlgebra> > flux,
			SmartPtr<GridFunction<dom_type, CPUAlgebra> > vol,
			number scale_factor
		);


		SmartPtr<GridFunction<typename TGridFunction::domain_type, CPUAlgebra> >
		calc_flux(number scale_factor);

	protected:
		struct CmpCoords
		{
		    bool operator() (const MathVector<dim>& a, const MathVector<dim>& b) const
		    {
		    	for (size_t d = 0; d < (size_t) dim; ++d)
		    	{
		    		const number a_d = a[d];
		    		const number b_d = b[d];
		    		const bool equal = fabs(a_d-b_d) < 1e-8*(1e-8+std::max(fabs(a_d),fabs(b_d)));
					if (!equal) return a_d < b_d;
		    	}

		    	return false;
		    }
		};

	protected:
		SmartPtr<TGridFunction> m_u;
		ConstSmartPtr<sh_type> m_sh;
		FunctionGroup m_fg;
		SubsetGroup m_sg;
		std::vector<std::string> m_vSubset;

		size_t m_cmpSpec;
		size_t m_cmpPot;

		LFEID m_lfeid;

		number m_diffConst;
		number m_convConst;

		int m_quadOrder;

		bool m_bWriteFluxMatrix;
		typedef typename CPUAlgebra::matrix_type matrix_type;
		matrix_type m_fluxMatrix;
		SmartPtr<DoFDistribution> m_ddBoxFluxes;
		//typedef std::vector<std::pair<MathVector<dim>, MathVector<dim> > > FluxMap;
		//FluxMap m_fluxMap;

		SmartPtr<IDomainConstraint<dom_type, CPUAlgebra> > m_hangingConstraintFlux;
		SmartPtr<IDomainConstraint<dom_type, CPUAlgebra> > m_hangingConstraintVol;
		SmartPtr<IDomainConstraint<dom_type, CPUAlgebra> > m_hangingConstraintBoxFlux;

		SmartPtr<IConvectionShapes<dim> > m_spConvShape;
};

} // namespace nernst_planck
} // namespace ug

#endif // UG__PLUGINS__NERNST_PLANCK__UTIL__FLUX_EXPORTER_H
