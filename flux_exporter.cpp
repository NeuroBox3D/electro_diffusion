/*
 * flux_exporter.cpp
 *
 *  Created on: 04.07.2016
 *      Author: mbreit
 */

#include "flux_exporter.h"
#include "vtk_export_ho.h"

#include "lib_disc/function_spaces/grid_function.h"				// GridFunction
#include "lib_disc/spatial_disc/disc_util/geom_provider.h"		// GeomProvider
#include "lib_disc/spatial_disc/disc_util/fvho_geom.h"			// FVGeom, DimFVGeom
#include "lib_disc/spatial_disc/disc_util/fv1_geom.h"			// FV1Geometry
#include "lib_disc/spatial_disc/disc_util/hfv1_geom.h"			// HFV1Geometry
#include "lib_disc/common/groups_util.h"						// CreateFunctionIndexMapping
#include "lib_grid/grid_objects/grid_dim_traits.h"				// grid_dim_traits
#include "lib_disc/function_spaces/dof_position_util.h"			// DoFPosition
#include "lib_disc/function_spaces/grid_function_user_data.h"	// GridFunctionVectorData
#include "lib_grid/grid_objects/grid_objects.h"					// geometry_traits<TElem>::REFERENCE_OBJECT_ID

#include <set>

namespace ug {
namespace nernst_planck {



template <typename TGridFunction>
FluxExporter<TGridFunction>::FluxExporter
(
	SmartPtr<TGridFunction> u,
	std::string cmp_name_spec,
	std::string cmp_name_pot

)
: m_u(u), m_diffConst(0.0), m_convConst(0.0), m_quadOrder(-1)
{
	// set subset handler
	m_sh = m_u->domain()->subset_handler();

	// find out indices of species and potential in grid function & create fctGrp
	m_fg.set_function_pattern(m_u->dof_distribution_info());
	try
	{
		m_fg.add(cmp_name_spec);
		m_fg.add(cmp_name_pot);
	}
	UG_CATCH_THROW("Function names are faulty.");

	try
	{
		m_cmpSpec = m_u->fct_id_by_name(cmp_name_spec.c_str());
		m_cmpPot = m_u->fct_id_by_name(cmp_name_pot.c_str());
	}
	UG_CATCH_THROW("Species and/or potential names could not be found in passed grid function.");

	// find out order and lfeid of species discretization
	m_lfeid = m_u->lfeid(m_cmpSpec);

	// check that species and potential functions have the same shape functions
	UG_COND_THROW(m_u->lfeid(m_cmpPot).order() != m_lfeid.order(),
			"Potential function does not have the same order as species function ("
			<< m_u->lfeid(m_cmpPot).order() << " / " << m_lfeid.order() << ")");
	UG_COND_THROW(m_u->lfeid(m_cmpPot).type() != m_lfeid.type(),
			"Potential function does not have the same type as species function ("
			<< m_u->lfeid(m_cmpPot).type() << " / " << m_lfeid.type() << ")");
}


template <typename TGridFunction>
void FluxExporter<TGridFunction>::set_subsets(const std::vector<std::string>& vSubsets)
{
	// check valid subset handler
	UG_COND_THROW(!m_sh.valid(), "Subset handler of passed grid function is invalid.")

	// check subsets: only allow full-dim subsets
	m_sg.set_subset_handler(m_sh);
	try
	{
		if (vSubsets.empty())
			m_sg.add_all();
		else
			m_sg.add(vSubsets);
	}
	UG_CATCH_THROW("Subsets are faulty.");

	for (size_t s = 0; s < m_sg.size(); ++s)
	{
		const int si = m_sg[s];
		const int ssdim = DimensionOfSubset(*m_sh, si);
		if (ssdim != dim)
		{
			// remove lower-dim subsets introduced by automated addition of all subsets
			if (vSubsets.empty())
				m_sg.remove(si);

			// if chosen by user, throw
			else
				UG_COND_THROW(ssdim != dim , "Only subsets of full domain dim allowed.")
		}
	}

	// construct vector of subsets on which to distribute unknowns for flux and vol;
	// add chosen element subsets and all subsets that are part of their borders
	std::set<int> reqSubsetsForDD;
	MultiGrid& mg = *m_u->domain()->grid();
	for (size_t s = 0; s < m_sg.size(); ++s)
	{
		const int si = m_sg[s];

		reqSubsetsForDD.insert(si);

		it_type it = m_u->template begin<elem_type>(si);
		it_type it_end = m_u->template end<elem_type>(si);
		for (; it != it_end; ++it)
		{
			switch (dim)
			{
				case 3:
					add_side_subsets<Face>(mg, *it, reqSubsetsForDD);
					// no break; also enter lower-dim cases
				case 2:
					add_side_subsets<Edge>(mg, *it, reqSubsetsForDD);
					// no break; also enter lower-dim cases
				case 1:
					add_side_subsets<Vertex>(mg, *it, reqSubsetsForDD);
					break;
				default: UG_THROW("Dimension other than 1, 2, 3 not implemented.")
			}
		}
	}

	std::set<int>::const_iterator it = reqSubsetsForDD.begin();
	std::set<int>::const_iterator it_end = reqSubsetsForDD.end();

#ifdef UG_PARALLEL	// communicate subset ids
	size_t nSI = reqSubsetsForDD.size();
	std::vector<int> vSI(nSI);
	std::vector<int> vSIglob;
	for (; it != it_end; ++it)
		vSI.push_back(*it);
	pcl::ProcessCommunicator pc;
	pc.allgatherv(vSIglob, vSI);
	reqSubsetsForDD.insert(vSIglob.begin(), vSIglob.end());
	it = reqSubsetsForDD.begin();	// reset iterator!
#endif

	m_vSubset.clear();
	for (; it != it_end; ++it)
		m_vSubset.push_back(std::string(m_sh->get_subset_name(*it)));
}

template <typename TGridFunction>
void FluxExporter<TGridFunction>::set_subsets(const char* cSubsets)
{
	set_subsets(TokenizeTrimString(cSubsets));
}


template <typename TGridFunction>
void FluxExporter<TGridFunction>::write_flux
(
	SmartPtr<VTKOutput<dim> > vtkOutput,
	std::string filename,
	size_t step,
	number time,
	std::string fluxName,
	number scale_factor
)
{
	// some checks first
	UG_COND_THROW(!m_u.valid(), "Grid function passed in constructor is not (or no longer) valid.");
	UG_COND_THROW(!m_sh.valid(), "Subset handler of passed grid function is invalid.");
	UG_COND_THROW(!m_fg.size(), "No functions selected for output.");
	UG_COND_THROW(!m_sg.size(), "No subsets selected for output.");
	UG_COND_THROW(!m_vSubset.size(), "No subsets selected for output.");
	UG_COND_THROW(!m_diffConst, "No diffusion constant (!= 0) specified.");
	UG_COND_THROW(!m_convConst, "No convection constant (!= 0) specified.");

	// calculate flux as vector-valued grid function
	SmartPtr<TGridFunction> flux = calc_flux(scale_factor);

	// export to vtk
	std::vector<std::string> flux_cmp_names(dim);
	if (dim >= 1) flux_cmp_names[0] = std::string("flux_x");
	if (dim >= 2) flux_cmp_names[1] = std::string("flux_y");
	if (dim >= 3) flux_cmp_names[2] = std::string("flux_z");
	vtkOutput->clear_selection();
	vtkOutput->select(flux_cmp_names, fluxName.c_str());

	vtk_export_ho<TGridFunction, dim>(flux, flux_cmp_names, m_lfeid.order(), vtkOutput, filename.c_str(), step, time, m_sg);

	//size_t sgSz = m_sg.size();
	//for (size_t s = 0; s < sgSz; ++s)
	//	vtkOutput->print_subset(filename.c_str(), *flux, m_sg[s], step, time);

	// make sure that flux grid function is deleted!
	//vtkOutput->clear_data_selection();
}




template <typename TGridFunction>
template <typename TFVGeom, typename Dummy>
FluxExporter<TGridFunction>::prep_elem_loop<TFVGeom, Dummy>::prep_elem_loop
(
	const FluxExporter<TGridFunction>* flEx,
	const ReferenceObjectID roid
)
{
	TFVGeom& geo = GeomProvider<TFVGeom>::get(flEx->m_lfeid, flEx->m_quadOrder);
	try {geo.update_local(roid, flEx->m_lfeid, flEx->m_quadOrder);}
	UG_CATCH_THROW("Failed updating Finite Volume Geometry for elem type.");
}

template <typename TGridFunction>
template <typename TElem, int dim, template <class, int> class TFV1Geom, typename Dummy>
FluxExporter<TGridFunction>::prep_elem_loop<TFV1Geom<TElem, dim>, Dummy>::prep_elem_loop
(
	const FluxExporter<TGridFunction>* flEx,
	const ReferenceObjectID roid
)
{
	// convection shapes?
}



template <typename TGridFunction>
template <typename TFVGeom, typename Dummy>
FluxExporter<TGridFunction>::prep_elem<TFVGeom, Dummy>::prep_elem
(
	const FluxExporter<TGridFunction>* flEx,
	GridObject* elem,
	const std::vector<MathVector<FluxExporter<TGridFunction>::dim> >& vCornerCoords
)
{
	TFVGeom& geo = GeomProvider<TFVGeom>::get(flEx->m_lfeid, flEx->m_quadOrder);
	try {geo.update(elem, &vCornerCoords[0], flEx->m_sh.get());}
	UG_CATCH_THROW("Failed updating Finite Volume Geometry for elem.");
}



// specialize for FV1Geometry and HFV1Geometry
// (little bit tricky, only these two have the signature template <class, int> class)
template <typename TGridFunction>
template<typename TElem, int dim, template <class, int> class TFV1Geom, typename Dummy>
FluxExporter<TGridFunction>::prep_elem<TFV1Geom<TElem, dim>, Dummy>::prep_elem
(
	const FluxExporter<TGridFunction>* flEx,
	GridObject* elem,
	const std::vector<MathVector<FluxExporter<TGridFunction>::dim> >& vCornerCoords
)
{
	TFV1Geom<TElem, dim>& geo = GeomProvider<TFV1Geom<TElem, dim> >::get();
	try {geo.update(elem, &vCornerCoords[0], flEx->m_sh.get());}
	UG_CATCH_THROW("Failed updating Finite Volume Geometry for elem.");
}



template <typename TGridFunction>
template<typename TFVGeom, typename Dummy>
FluxExporter<TGridFunction>::assemble_flux_elem<TFVGeom, Dummy>::assemble_flux_elem
(
	const FluxExporter<TGridFunction>* flEx,
	LocalVector& f,
	const LocalVector& u,
	GridObject* elem,
	const std::vector<MathVector<FluxExporter<TGridFunction>::dim> >& vCornerCoords
)
{
	// request FV geometry
	const TFVGeom& geo = GeomProvider<TFVGeom>::get(flEx->m_lfeid, flEx->m_quadOrder);
	const size_t dim = FluxExporter<TGridFunction>::dim;
	const size_t _C_ = 0;
	const size_t _P_ = 1;


	// get dof positions
	std::vector<MathVector<dim> > vDoFPos;
	if (!DoFPosition<dim>(vDoFPos, elem->reference_object_id(), vCornerCoords, flEx->m_lfeid))
		UG_THROW("Failed getting exact DoF positions.");


	// loop SCVFs
	for (size_t s = 0; s < geo.num_scvf(); ++s)
	{
		const typename TFVGeom::SCVF& scvf = geo.scvf(s);

		MathVector<dim> flux(0.0);

		// loop integration points
		for (size_t ip = 0; ip < scvf.num_ip(); ++ip)
		{
			number scalarFluxIP;

			// diffusive term
			MathVector<dim> grad(0.0);
			for (size_t sh = 0; sh < scvf.num_sh(); ++sh)
				VecScaleAppend(grad, u(_C_,sh), scvf.global_grad(ip, sh));
			VecScale(grad, grad, flEx->m_diffConst);
			scalarFluxIP = -VecDot(grad, scvf.normal());

			// convective term
			number cAtIP = 0.0;
			for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
				cAtIP += u(_C_, sh) * scvf.shape(ip, sh);
			VecSet(grad, 0.0);
			for (size_t sh = 0; sh < scvf.num_sh(); ++sh)
				VecScaleAppend(grad, u(_P_,sh), scvf.global_grad(ip, sh));
			scalarFluxIP -= flEx->m_convConst * cAtIP * VecDot(grad, scvf.normal());

			// multiply by vector from box center to ip and add to flux
			for (size_t i = 0; i < dim; ++i)
			{
				size_t from = scvf.from();
				VecScaleAdd(grad, 1.0, scvf.global_ip(ip), -1.0, vDoFPos[from]);
				f(i, from) += scvf.weight(ip) * scalarFluxIP * grad[i];

				size_t to = scvf.to();
				VecScaleAdd(grad, 1.0, scvf.global_ip(ip), -1.0, vDoFPos[to]);
				f(i, to) -= scvf.weight(ip) * scalarFluxIP * grad[i];
			}
		}
	}
}

template <typename TGridFunction>
template<typename TElem, int dim, template <class, int> class TFV1Geom, typename Dummy>
FluxExporter<TGridFunction>::assemble_flux_elem<TFV1Geom<TElem, dim>, Dummy>::assemble_flux_elem
(
	const FluxExporter<TGridFunction>* flEx,
	LocalVector& f,
	const LocalVector& u,
	GridObject* elem,
	const std::vector<MathVector<FluxExporter<TGridFunction>::dim> >& vCornerCoords
)
{
	// request FV geometry
	const TFV1Geom<TElem, dim>& geo = GeomProvider<TFV1Geom<TElem, dim> >::get();

	const size_t _C_ = 0;
	const size_t _P_ = 1;

	// loop SCVFs
	for (size_t s = 0; s < geo.num_scvf(); ++s)
	{
		const typename TFV1Geom<TElem, dim>::SCVF& scvf = geo.scvf(s);

		number scalarFluxIP;

		// diffusive term
		MathVector<dim> grad(0.0);
		for (size_t sh = 0; sh < scvf.num_sh(); ++sh)
			VecScaleAppend(grad, u(_C_,sh), scvf.global_grad(sh));
		VecScale(grad, grad, flEx->m_diffConst);
		scalarFluxIP = -VecDot(grad, scvf.normal());

		// convective term
		number cAtIP = 0.0;
		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
			cAtIP += u(_C_, sh) * scvf.shape(sh);
		VecSet(grad, 0.0);
		for (size_t sh = 0; sh < scvf.num_sh(); ++sh)
			VecScaleAppend(grad, u(_P_,sh), scvf.global_grad(sh));
		scalarFluxIP -= flEx->m_convConst * cAtIP * VecDot(grad, scvf.normal());

		// multiply by vector from box center to ip and add to flux
		size_t from = scvf.from();
		const typename TFV1Geom<TElem, dim>::SCV& scvFrom = geo.scv(from);
		size_t to = scvf.to();
		const typename TFV1Geom<TElem, dim>::SCV& scvTo = geo.scv(to);

		VecScaleAdd(grad, 1.0, scvf.global_ip(), -1.0, scvFrom.global_ip());
		for (size_t i = 0; i < dim; ++i)
			f(i, from) += scalarFluxIP * grad[i];

		VecScaleAdd(grad, 1.0, scvf.global_ip(), -1.0, scvTo.global_ip());
		for (size_t i = 0; i < dim; ++i)
			f(i, to) -= scalarFluxIP * grad[i];
	}
}


template <typename TGridFunction>
template<typename TFVGeom, typename Dummy>
FluxExporter<TGridFunction>::assemble_vol_elem<TFVGeom, Dummy>::assemble_vol_elem
(
	const FluxExporter<TGridFunction>* flEx,
	LocalVector& vol
)
{
	// request FV geometry
	const TFVGeom& geo = GeomProvider<TFVGeom>::get(flEx->m_lfeid, flEx->m_quadOrder);

	// loop SCVFs
	for (size_t s = 0; s < geo.num_scv(); ++s)
	{
		const typename TFVGeom::SCV& scv = geo.scv(s);

		// loop integration points
		for (size_t ip = 0; ip < scv.num_ip(); ++ip)
			vol(0, scv.node_id()) += scv.weight(ip);
	}
}

template <typename TGridFunction>
template<typename TElem, int dim, template <class, int> class TFV1Geom, typename Dummy>
FluxExporter<TGridFunction>::assemble_vol_elem<TFV1Geom<TElem, dim>, Dummy>::assemble_vol_elem
(
	const FluxExporter<TGridFunction>* flEx,
	LocalVector& vol
)
{
	// request FV geometry
	const TFV1Geom<TElem, dim>& geo = GeomProvider<TFV1Geom<TElem, dim> >::get();

	// loop SCVFs
	for (size_t s = 0; s < geo.num_scv(); ++s)
	{
		const typename TFV1Geom<TElem, dim>::SCV& scv = geo.scv(s);

		vol(0, scv.node_id()) += scv.volume();
	}
}




// Choose FVGeom as in ConvectionDiffusionFV.

// standard FVGeom
template <int dim, typename TElem, bool hanging, int order, int quadOrder = order+1>
struct ChooseProperFVGeom
{
	typedef DimFVGeometry<dim> FVGeom;
};

// for order 1 and !hanging: use FV1Geom
template <int dim, typename TElem, int quadOrder>
struct ChooseProperFVGeom<dim, TElem, false, 1, quadOrder>
{
	typedef FV1Geometry<TElem, dim> FVGeom;
};

// for order 1 and hanging: use HFV1Geom
template <int dim, typename TElem, int quadOrder>
struct ChooseProperFVGeom<dim, TElem, true, 1, quadOrder>
{
	typedef HFV1Geometry<TElem, dim> FVGeom;
};

// take FVGeom instead if dim > 1 and 1 < order < 4 and quadOrder == order+1
template <typename TElem, bool hanging>
struct ChooseProperFVGeom<2, TElem, hanging, 2>
{
	typedef FVGeometry<2, TElem, 2> FVGeom;
};
template <typename TElem, bool hanging>
struct ChooseProperFVGeom<2, TElem, hanging, 3>
{
	typedef FVGeometry<3, TElem, 2> FVGeom;
};
template <typename TElem, bool hanging>
struct ChooseProperFVGeom<3, TElem, hanging, 2>
{
	typedef FVGeometry<2, TElem, 3> FVGeom;
};
template <typename TElem, bool hanging>
struct ChooseProperFVGeom<3, TElem, hanging, 3>
{
	typedef FVGeometry<3, TElem, 3> FVGeom;
};

// take DimFVGeometry again instead if dim == 3 and TElem == Pyramid, Octahedron
template <bool hanging>
struct ChooseProperFVGeom<3, Pyramid, hanging, 2>
{
	typedef DimFVGeometry<3> FVGeom;
};
template <bool hanging>
struct ChooseProperFVGeom<3, Pyramid, hanging, 3>
{
	typedef DimFVGeometry<3> FVGeom;
};
template <bool hanging>
struct ChooseProperFVGeom<3, Octahedron, hanging, 2>
{
	typedef DimFVGeometry<3> FVGeom;
};
template <bool hanging>
struct ChooseProperFVGeom<3, Octahedron, hanging, 3>
{
	typedef DimFVGeometry<3> FVGeom;
};




template <typename TGridFunction>
template <bool hanging, int order, typename TElem>
void FluxExporter<TGridFunction>::assemble
(
	SmartPtr<TGridFunction> flux,
	SmartPtr<TGridFunction> vol,
	int si
) const
{
	typedef typename TGridFunction::domain_type dom_type;
	typedef typename DoFDistribution::traits<TElem>::const_iterator it_type;
	static const int dim = dom_type::dim;

	ConstSmartPtr<DoFDistribution> ddU = m_u->dd();
	ConstSmartPtr<DoFDistribution> ddFlux = flux->dd();
	ConstSmartPtr<DoFDistribution> ddVol = vol->dd();

	// get elem iterators
	it_type it = ddU->template begin<TElem>(si);
	it_type it_end = ddU->template end<TElem>(si);

	if (it == it_end) return;

	// reference object id
	static const ReferenceObjectID roid = geometry_traits<TElem>::REFERENCE_OBJECT_ID;

	// prep elem loop
	if (m_quadOrder == order+1)
	{
		typedef typename ChooseProperFVGeom<dim, TElem, hanging, order>::FVGeom FVGeom;
		prep_elem_loop<FVGeom>(this, roid);
	}
	else
	{
		typedef typename ChooseProperFVGeom<dim, TElem, hanging, order, order>::FVGeom FVGeom;
		prep_elem_loop<FVGeom>(this, roid);
	}

	// local indices and local algebra
	LocalIndices indU;
	LocalIndices indFlux;
	LocalIndices indVol;
	LocalVector locU;
	LocalVector locFlux;
	LocalVector locVol;

	std::vector<MathVector<dim> > vCornerCoords(TElem::NUM_VERTICES);

	// loop elements
	for (; it != it_end; ++it)
	{
		TElem* elem = *it;

		// get corner coords
		FillCornerCoordinates(&vCornerCoords[0], *elem, *m_u->approx_space()->domain());

		// get global indices
		ddU->indices(elem, indU, true);
		ddFlux->indices(elem, indFlux, true);
		ddVol->indices(elem, indVol, true);

		// resize local algebra
		locU.resize(indU);
		locFlux.resize(indFlux);
		locVol.resize(indVol);

		// read local values for solution
		GetLocalVector(locU, *m_u);

		// use only indices 0 and 1 for spec and pot
		FunctionIndexMapping fim;
		try {CreateFunctionIndexMapping(fim, m_fg, m_fg.function_pattern());}
		UG_CATCH_THROW("Failed creating function index mapping.");
		locU.access_by_map(fim);

		// reset local vectors
		locFlux = 0.0;
		locVol = 0.0;

		// prepare elem; assemble flux; assemble volume
		ConstSmartPtr<dom_type> dom = m_u->approx_space()->domain();
		if (m_quadOrder == order+1)
		{
			typedef typename ChooseProperFVGeom<dim, TElem, hanging, order>::FVGeom FVGeom;
			prep_elem<FVGeom>(this, elem, vCornerCoords);
			assemble_flux_elem<FVGeom>(this, locFlux, locU, elem, vCornerCoords);
			assemble_vol_elem<FVGeom>(this, locVol);
		}
		else
		{
			typedef typename ChooseProperFVGeom<dim, TElem, hanging, order, order>::FVGeom FVGeom;
			prep_elem<FVGeom>(this, elem, vCornerCoords);
			assemble_flux_elem<FVGeom>(this, locFlux, locU, elem, vCornerCoords);
			assemble_vol_elem<FVGeom>(this, locVol);
		}

		// add local vecs to global
		AddLocalVector(*flux, locFlux);
		AddLocalVector(*vol, locVol);
	}
}


template <typename TGridFunction>
template <int order>
void FluxExporter<TGridFunction>::assemble
(
	SmartPtr<TGridFunction> flux,
	SmartPtr<TGridFunction> vol
)
{
	typedef typename TGridFunction::domain_type dom_type;
	static const int dim = dom_type::dim;

	// set quad order if not set by user
	if (m_quadOrder == -1)
		m_quadOrder = m_lfeid.order() + 1;

	// loop subsets
	for (size_t s = 0; s < m_sg.size(); ++s)
	{
		const int si = m_sg[s];
		bool bNonRegularGrid = !m_sg.regular_grid(s);

		try
		{
			// get all grid element types in this dim
			typedef typename domain_traits<dim>::DimElemList ElemList;

			// switch assemble functions
			if (bNonRegularGrid)
			{
				UG_COND_THROW(order > 1, "Function order is >1, but hanging nodes detected.\n"
					"This should not happen as hanging nodes are only supported with order 1.");

				boost::mpl::for_each<ElemList>(AssembleWrapper<true, order>(this, flux, vol, si));
			}
			else
			{
				boost::mpl::for_each<ElemList>(AssembleWrapper<true, order>(this, flux, vol, si));
			}
		}
		UG_CATCH_THROW("assemble: Assembling of elements of dimension "
						<< dim << " in subset " << si << " failed.");
	}
}


template <typename TGridFunction>
template <typename TBaseElem>
void FluxExporter<TGridFunction>::div_flux_by_vol_and_scale
(
	SmartPtr<TGridFunction> flux,
	SmartPtr<TGridFunction> vol,
	number scale_factor
)
{
	typedef typename DoFDistribution::traits<TBaseElem>::const_iterator it_type;

	ConstSmartPtr<DoFDistribution> dd = vol->dd();

	// in the parallel case: make flux and volume consistent (is additive)
#ifdef UG_PARALLEL
	if (pcl::NumProcs() > 1)
	{
		flux->set_storage_type(PST_ADDITIVE);
		vol->set_storage_type(PST_ADDITIVE);

		flux->change_storage_type(PST_CONSISTENT);
		vol->change_storage_type(PST_CONSISTENT);
	}
#endif

	// get elem iterators
	it_type it = dd->template begin<TBaseElem>();
	it_type it_end = dd->template end<TBaseElem>();

	for (; it != it_end; ++it)
	{
		TBaseElem* elem = *it;

		std::vector<DoFIndex> volInd;
		vol->inner_dof_indices(elem, 0, volInd, false);
		size_t volIndSz = volInd.size();

		std::vector<DoFIndex> fluxInd;
		for (size_t fct = 0; fct < flux->num_fct(); ++fct)
		{
			flux->inner_dof_indices(elem, fct, fluxInd, true);

			UG_COND_THROW(fluxInd.size() != volIndSz, "Index numbers do not match (flux: "
				<< fluxInd.size() << ", vol: " << volIndSz << ").");

			for (size_t i = 0; i < volIndSz; ++i)
			{
				number& fluxVal = DoFRef(*flux, fluxInd[i]);
				number& volVal = DoFRef(*vol, volInd[i]);
				if (volVal == 0.0)
				{
					if (fluxVal != 0.0)
						UG_THROW("Zero volume.")
				}
				else fluxVal *= scale_factor / volVal;
			}
		}
	}
}


template <typename TGridFunction>
template <typename TElem>
void FluxExporter<TGridFunction>::add_side_subsets
(
	MultiGrid& mg,
	GridObject* elem,
	std::set<int>& sss
)
{
	typedef typename MultiGrid::traits<TElem>::secure_container elem_list;
	elem_list el;
	mg.associated_elements(el, elem);
	size_t elsz = el.size();
	for (size_t e = 0; e < elsz; ++e)
	{
		int esi = m_sh->get_subset_index(el[e]);
		sss.insert(esi);
	}
}

template <typename TGridFunction>
SmartPtr<TGridFunction> FluxExporter<TGridFunction>::calc_flux(number scale_factor)
{
	// set up approx space with flux function of the same order as solution grid function
	SmartPtr<ApproximationSpace<dom_type> > approxFlux =
		make_sp(new ApproximationSpace<dom_type>(m_u->approx_space()->domain()));

	std::vector<std::string> flux_cmp_names(dim);
	if (dim >= 1) flux_cmp_names[0] = std::string("flux_x");
	if (dim >= 2) flux_cmp_names[1] = std::string("flux_y");
	if (dim >= 3) flux_cmp_names[2] = std::string("flux_z");

	try {approxFlux->add(flux_cmp_names, m_lfeid, m_vSubset);}
	UG_CATCH_THROW("Failed to add functions to approximation space on specified subsets.");
	approxFlux->init_top_surface();

	// set up approx space for box volumes
	SmartPtr<ApproximationSpace<dom_type> > approxVol =
			make_sp(new ApproximationSpace<dom_type>(m_u->approx_space()->domain()));
	std::vector<std::string> vol_cmp_names(1);
	vol_cmp_names[0] = std::string("vol");
	try {approxVol->add(vol_cmp_names, m_lfeid, m_vSubset);}
		UG_CATCH_THROW("Failed to add functions to approximation space on specified subsets.");
	approxVol->init_top_surface();


	// create flux and volume grid functions from their respective approx spaces
	SmartPtr<TGridFunction> flux = make_sp(new TGridFunction(approxFlux));
	SmartPtr<TGridFunction> vol = make_sp(new TGridFunction(approxVol));
	flux->set(0.0);	// in order to init CONSISTENT
	vol->set(0.0);	// in order to init CONSISTENT


	// calculate flux and volumes
	switch (m_lfeid.order())
	{
		case 1:	assemble<1>(flux, vol);
				break;
		case 2:	assemble<2>(flux, vol);
				break;
		case 3:	assemble<3>(flux, vol);
				break;
		default: assemble<4>(flux, vol);
	}


	// apply constraints (notably for hanging nodes)
	if (m_hangingConstraint.valid())
	{
		m_hangingConstraint->adjust_defect(*flux, *m_u, flux->dd(), CT_HANGING);
		m_hangingConstraint->adjust_defect(*vol, *m_u, vol->dd(), CT_HANGING);
	}

	// box-wise division of flux value by box volume to get true flux density; scaling
	if (dim >= VOLUME && vol->max_dofs(VOLUME))
			div_flux_by_vol_and_scale<Volume>(flux, vol, scale_factor);
	if (dim >= FACE && vol->max_dofs(FACE))
			div_flux_by_vol_and_scale<Face>(flux, vol, scale_factor);
	if (dim >= EDGE && vol->max_dofs(EDGE))
			div_flux_by_vol_and_scale<Edge>(flux, vol, scale_factor);
	if (dim >= VERTEX && vol->max_dofs(VERTEX))
			div_flux_by_vol_and_scale<Vertex>(flux, vol, scale_factor);

	return flux;
}


// explicit template specializations
#ifdef UG_CPU_1
	#ifdef UG_DIM_1
		template class FluxExporter<GridFunction<Domain1d, CPUAlgebra> >;
	#endif
	#ifdef UG_DIM_2
		template class FluxExporter<GridFunction<Domain2d, CPUAlgebra> >;
	#endif
	#ifdef UG_DIM_3
		template class FluxExporter<GridFunction<Domain3d, CPUAlgebra> >;
	#endif
#endif


} // namespace nernst_planck
} // namespace ug
