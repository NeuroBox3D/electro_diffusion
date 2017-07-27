/*
 * flux_exporter.h
 *
 *  Created on: 04.07.2016
 *      Author: mbreit
 */

#ifndef UG__PLUGINS__NERNST_PLANCK__FLUX_EXPORTER_H
#define UG__PLUGINS__NERNST_PLANCK__FLUX_EXPORTER_H

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
#include "lib_disc/common/function_group.h"                                // for FunctionGroup
#include "lib_disc/common/local_algebra.h"                                 // for LocalVector
#include "lib_disc/local_finite_element/local_finite_element_id.h"         // for LFEID
#include "lib_grid/grid/grid_base_objects.h"                               // for GridObject (pt...
#include "lib_grid/grid_objects/grid_dim_traits.h"                         // for grid_dim_traits
#include "lib_grid/tools/subset_group.h"                                   // for SubsetGroup
#include "pcl/pcl_base.h"                                                  // for NumProcs


namespace ug {

// forward declarations
class MultiGrid;
template <int TDim> class VTKOutput;
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
		FluxExporter(SmartPtr<TGridFunction> u, std::string cmp_name_spec, std::string cmp_name_pot);
		virtual ~FluxExporter() {};

		void set_diff_const(number diff_const) {m_diffConst = diff_const;}
		void set_conv_const(number conv_const) {m_convConst = conv_const;}
		void set_quad_order(number quadOrder) {m_quadOrder = quadOrder;}
		void set_hanging_constraint(SmartPtr<IConstraint<algebra_type> > constr) {m_hangingConstraint = constr;}
		void set_upwind(SmartPtr<IConvectionShapes<dim> > shapes) {m_spConvShape = shapes;}

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
				SmartPtr<TGridFunction> _flux,
				SmartPtr<TGridFunction> _vol,
				int _si
			)
			: flEx(_flEx), flux(_flux), vol(_vol), si(_si) {}

			template <typename TElem>
			void operator() (TElem&)
			{
				flEx->assemble<hanging, order, TElem>(flux, vol, si);
			}

			FluxExporter<TGridFunction>* flEx;
			SmartPtr<TGridFunction> flux;
			SmartPtr<TGridFunction> vol;
			int si;
		};
		template <bool hanging, int order>
		friend struct AssembleWrapper;

		template <bool hanging, int order, typename TElem>
		void assemble(SmartPtr<TGridFunction> flux, SmartPtr<TGridFunction> vol, int si);

		template <int order>
		void assemble(SmartPtr<TGridFunction> flux, SmartPtr<TGridFunction> vol);


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
		void div_flux_by_vol_and_scale(SmartPtr<TGridFunction> flux, SmartPtr<TGridFunction> vol, number scale_factor);

		SmartPtr<TGridFunction> calc_flux(number scale_factor);

	protected:
		struct CmpCoords
		{
		    bool operator() (const MathVector<dim>& a, const MathVector<dim>& b) const
		    {
		    	for (size_t d = 0; d < dim; ++d)
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

		bool m_bWriteFluxMap;
		typedef std::map<MathVector<dim>, std::pair<MathVector<dim>, number>, CmpCoords> FluxMap;
		FluxMap m_fluxMap;

		SmartPtr<IConstraint<algebra_type> > m_hangingConstraint;
		SmartPtr<IConvectionShapes<dim> > m_spConvShape;
};

} // namespace nernst_planck
} // namespace ug

#endif // UG__PLUGINS__NERNST_PLANCK__FLUX_EXPORTER_H
