/*
 * flux_exporter.h
 *
 *  Created on: 04.07.2016
 *      Author: mbreit
 */

#ifndef UG__PLUGINS__NERNST_PLANCK__FLUX_EXPORTER_H
#define UG__PLUGINS__NERNST_PLANCK__FLUX_EXPORTER_H


#include "lib_disc/spatial_disc/constraints/constraint_interface.h"	// IConstraint
#include "lib_disc/io/vtkoutput.h"									// VTKOutput

#include <vector>
#include <string>

namespace ug {
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

	protected:
		template <typename TElem>
		void add_side_subsets(MultiGrid& mg, GridObject* elem, std::set<int>& sss);

		template <bool hanging, int order>
		struct AssembleWrapper
		{
			AssembleWrapper
			(
				const FluxExporter<TGridFunction>* _flEx,
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

			const FluxExporter<TGridFunction>* flEx;
			SmartPtr<TGridFunction> flux;
			SmartPtr<TGridFunction> vol;
			int si;
		};
		template <bool hanging, int order>
		friend struct AssembleWrapper;

		template <bool hanging, int order, typename TElem>
		void assemble(SmartPtr<TGridFunction> flux, SmartPtr<TGridFunction> vol, int si) const;

		template <int order>
		void assemble(SmartPtr<TGridFunction> flux, SmartPtr<TGridFunction> vol);


		// PREP ELEM LOOP
		template <typename TFVGeom, typename Dummy = void>
		struct prep_elem_loop
		{
			prep_elem_loop(const FluxExporter<TGridFunction>* flEx, const ReferenceObjectID roid);
		};
		template <typename TFVGeom, typename Dummy>
		friend struct prep_elem_loop;

		// specialize for FV1Geometry and HFV1Geometry
		// (little bit tricky, only these two have the signature template <class, int> class)
		template <typename TElem, int dim, template <class, int> class TFV1Geom, typename Dummy>
		struct prep_elem_loop<TFV1Geom<TElem, dim>, Dummy>
		{
			prep_elem_loop(const FluxExporter<TGridFunction>* flEx, const ReferenceObjectID roid);
		};


		// PREP ELEM
		template <typename TFVGeom, typename Dummy = void>
		struct prep_elem
		{
			prep_elem
			(
				const FluxExporter<TGridFunction>* flEx,
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
				const FluxExporter<TGridFunction>* flEx,
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
				const FluxExporter<TGridFunction>* flEx,
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
				const FluxExporter<TGridFunction>* flEx,
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

		SmartPtr<IConstraint<algebra_type> > m_hangingConstraint;
};

} // namespace nernst_planck
} // namespace ug

#endif // UG__PLUGINS__NERNST_PLANCK__FLUX_EXPORTER_H