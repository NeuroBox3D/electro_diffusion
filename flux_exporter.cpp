/*
 * flux_exporter.cpp
 *
 *  Created on: 04.07.2016
 *      Author: mbreit
 */

#include "flux_exporter.h"

#include <boost/mpl/for_each.hpp>                                    // for for_each
#include <limits>                                                    // for numeric_...
#include <fstream>                                                   // for ofstream
#include <mpi.h>                                                     // for MPI_File_write ...

#include "common/error.h"                                            // for UG_COND_...
#include "common/math/math_vector_matrix/math_matrix.h"              // for MathMatrix
#include "common/math/math_vector_matrix/math_vector_functions.h"    // for VecScale...
#ifdef UG_PARALLEL
	#include "lib_algebra/parallelization/parallel_storage_type.h"       // for Parallel...
	#include "lib_algebra/parallelization/parallelization_util.h"        // for MatMakeConsistentOverlap0
#endif
#include "lib_disc/common/groups_util.h"                             // for CreateFunctionIndexMapping...
#include "lib_disc/common/multi_index.h"                             // for DoFIndex
#include "lib_disc/dof_manager/dof_distribution.h"                   // for DoFDistr...
#include "lib_disc/domain.h"                                         // for Domain1d, Doma...
#include "lib_disc/domain_traits.h"                                  // for domain_traits
#include "lib_disc/domain_util.h"                                    // for FillCornerCoordinates
#include "lib_disc/function_spaces/approximation_space.h"            // for ApproximationSpace
#include "lib_disc/function_spaces/dof_position_util.h"              // for DoFPosition
#include "lib_disc/spatial_disc/ass_tuner.h"                         // for ConstraintType...
#include "lib_disc/spatial_disc/constraints/constraint_interface.h"  // for IConstraint
#include "lib_disc/spatial_disc/disc_util/conv_shape.h"              // for Convecti...
#include "lib_disc/spatial_disc/disc_util/fv1_geom.h"                // for FV1Geometry
#include "lib_disc/spatial_disc/disc_util/fvho_geom.h"               // for DimFVGeometry
#include "lib_disc/spatial_disc/disc_util/geom_provider.h"           // for GeomProvider
#include "lib_disc/spatial_disc/disc_util/hfv1_geom.h"               // for HFV1Geometry
#include "lib_grid/algorithms/subset_dim_util.h"                     // for DimensionOfSubset
#include "lib_grid/grid/grid_base_object_traits.h"                   // for geometry_traits
#include "lib_grid/multi_grid.h"                                     // for MultiGrid
#include "pcl/pcl_base.h"                                            // for NumProcs
#ifdef UG_PARALLEL
	#include "pcl/pcl_process_communicator.h"                        // for ProcessCommunicator
#endif

#include "vtk_export_ho.h"                                           // for vtk_export_ho
#include "choose_fvgeom.h"

namespace ug {

// forward declarations
class Octahedron;
class Pyramid;
template <int TDim> class VTKOutput;
template <int dim> class IConvectionShapes;

namespace nernst_planck {



template <typename TGridFunction>
FluxExporter<TGridFunction>::FluxExporter
(
	SmartPtr<TGridFunction> u,
	const std::string& cmp_name_spec,
	const std::string& cmp_name_pot

)
: m_u(u),
  m_diffConst(std::numeric_limits<number>::quiet_NaN()),
  m_convConst(std::numeric_limits<number>::quiet_NaN()),
  m_quadOrder(-1),
  m_bWriteFluxMatrix(false),
  m_hangingConstraintFlux(SPNULL),
  m_hangingConstraintVol(SPNULL),
  m_spConvShape(new ConvectionShapesNoUpwind<dim>)
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
void FluxExporter<TGridFunction>::set_hanging_constraint(SmartPtr<IDomainConstraint<dom_type, algebra_type> > constr)
{
	typedef OneSideP1Constraints<dom_type, algebra_type> AsymHC;
	AsymHC* hca = dynamic_cast<AsymHC*>(constr.get());
	if (hca)
	{
		m_hangingConstraintFlux = make_sp(new OneSideP1Constraints<dom_type, CPUBlockAlgebra<dim> >());
		m_hangingConstraintVol = make_sp(new OneSideP1Constraints<dom_type, CPUAlgebra>());
		m_hangingConstraintBoxFlux = make_sp(new OneSideP1Constraints<dom_type, CPUAlgebra>());
		return;
	}

	typedef SymP1Constraints<dom_type, algebra_type> SymHC;
	SymHC* hcs = dynamic_cast<SymHC*>(constr.get());
	if (hcs)
	{
		m_hangingConstraintFlux = make_sp(new SymP1Constraints<dom_type, CPUBlockAlgebra<dim> >());
		m_hangingConstraintVol = make_sp(new SymP1Constraints<dom_type, CPUAlgebra>());
		m_hangingConstraintBoxFlux = make_sp(new SymP1Constraints<dom_type, CPUAlgebra>());
		return;
	}

	UG_THROW("Only instances of SymP1Constraints or OneSideP1Constraints are allowed as hanging constraint,\n"
		"but given constraint is of neither type.");
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
				UG_COND_THROW(ssdim != dim , "Only subsets of full domain dim (" << dim << ") allowed, "
						"but subset " << m_sh->get_subset_name(si) << " has dim = " << ssdim << ".");
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
	UG_COND_THROW(m_diffConst != m_diffConst, "No diffusion constant specified.");
	UG_COND_THROW(m_convConst != m_convConst, "No convection constant (!= 0) specified.");

	// calculate flux as vector-valued grid function
	SmartPtr<GridFunction<dom_type, CPUBlockAlgebra<dim> > > flux;
	try {flux = calc_flux(scale_factor);}
	UG_CATCH_THROW("Error during flux calculation.")

	// export to vtk
	std::vector<std::string> flux_cmp_names(dim);
	if (dim >= 1) flux_cmp_names[0] = std::string("flux_x");
	if (dim >= 2) flux_cmp_names[1] = std::string("flux_y");
	if (dim >= 3) flux_cmp_names[2] = std::string("flux_z");
	vtkOutput->clear_selection();
	vtkOutput->select(flux_cmp_names, fluxName.c_str());

	try {vtk_export_ho<GridFunction<dom_type, CPUBlockAlgebra<dim> >, dim>
		(flux, flux_cmp_names, m_lfeid.order(), vtkOutput, filename.c_str(), step, time, m_sg);}
	UG_CATCH_THROW("Error during vtk export.")

	//size_t sgSz = m_sg.size();
	//for (size_t s = 0; s < sgSz; ++s)
	//	vtkOutput->print_subset(filename.c_str(), *flux, m_sg[s], step, time);

	// make sure that flux grid function is deleted!
	//vtkOutput->clear_data_selection();
}



template <typename TGridFunction>
void FluxExporter<TGridFunction>::write_box_fluxes
(
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
	UG_COND_THROW(m_diffConst != m_diffConst, "No diffusion constant specified.");
	UG_COND_THROW(m_convConst != m_convConst, "No convection constant (!= 0) specified.");

	// set flag for writing to flux map
	m_bWriteFluxMatrix = true;

	// resize matrix
	SmartPtr<ApproximationSpace<dom_type> > approxFlux =
		make_sp(new ApproximationSpace<dom_type>(m_u->approx_space()->domain(), AlgebraType(AlgebraType::CPU, 1)));
	std::vector<std::string> vol_cmp_names(1);
	vol_cmp_names[0] = "flux";
	try {approxFlux->add(vol_cmp_names, m_lfeid, m_vSubset);}
	UG_CATCH_THROW("Failed to add functions to approximation space on specified subsets.");
	approxFlux->init_top_surface();
	m_ddBoxFluxes = approxFlux->dd(GridLevel());
	m_fluxMatrix.resize_and_clear(m_ddBoxFluxes->num_indices(), m_ddBoxFluxes->num_indices());

	// calculate matrix entries
	calc_flux(scale_factor);

	// apply hanging constraint
	if (m_hangingConstraintBoxFlux.valid())
	{
		m_hangingConstraintBoxFlux->set_approximation_space(approxFlux);
		m_hangingConstraintBoxFlux->set_ass_tuner(make_sp(new AssemblingTuner<CPUAlgebra>())); // dummy
		GridFunction<dom_type, CPUAlgebra> dummy(approxFlux, false); // not needed in adjust_jacobian
		m_hangingConstraintBoxFlux->adjust_jacobian(m_fluxMatrix, dummy, m_ddBoxFluxes, CT_HANGING);
	}

#ifdef UG_PARALLEL
	// make consistent
	MatMakeConsistentOverlap0(m_fluxMatrix);
#endif

// save flux map contents to .csv file

	// first: associate a vertex with each DoF
	std::map<size_t, Vertex*> dof2Vrt;
	std::vector<DoFIndex> vInd;
	SmartPtr<MultiGrid> mg = m_u->approx_space()->domain()->grid();
	const size_t nSS = m_sg.size();
	for (size_t s = 0; s < nSS; ++s)
	{
		const int si = m_sg[s];
		it_type it = m_ddBoxFluxes->template begin<elem_type>(si);
		it_type itEnd = m_ddBoxFluxes->template end<elem_type>(si);
		for (; it != itEnd; ++it)
		{
			elem_type* elem = *it;
			const size_t nVrt = elem->num_vertices();
			for (size_t v = 0; v < nVrt; ++v)
			{
				m_ddBoxFluxes->inner_dof_indices(elem->vertex(v), 0, vInd, true);
				UG_COND_THROW(vInd.size() != 1, "Not exactly one DoF index for vertex, "
					"instead " << vInd.size() << ".");
				dof2Vrt[vInd[0][0]] = elem->vertex(v);
			}
		}
	}

	// get position accessor
	const typename dom_type::position_accessor_type& aaPos =
		m_u->approx_space()->domain()->position_accessor();

	// construct name
	AppendCounterToString(filename, "_t", step);
	filename.append(".csv");

	std::string header;

#ifdef UG_PARALLEL
	if (pcl::NumProcs() > 1)
	{
		// FIXME: write larger chunks
		// (possibly only one per proc and in a collective call)
		pcl::ProcessCommunicator pc;
		MPI_Status status;
		MPI_Comm m_mpiComm = pc.get_mpi_communicator();
		MPI_File fh;

		// open file
		if (MPI_File_open(m_mpiComm, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh))
			UG_THROW("Unable to open " << filename << ".");

		// calculate offsets for each proc
		const size_t precision = std::numeric_limits<number>::digits10;

		std::ostringstream testBuf;
		testBuf << std::scientific << std::setprecision(precision) << std::showpos << 1.0;
		const size_t doubleValueSize = testBuf.str().size();

		size_t nEntries = 0;
		const size_t nRow = m_fluxMatrix.num_rows();
		for (size_t r = 0; r < nRow; ++r)
		{
			typename matrix_type::row_iterator cit = m_fluxMatrix.begin_row(r);
			typename matrix_type::row_iterator citEnd = m_fluxMatrix.end_row(r);
			for (; cit != citEnd && cit.index() < r; ++cit) // only sub-diagonal values count
				if (cit.value()) // only non-zero values count
					++nEntries;
		}

		unsigned long mySize = nEntries * 2*dim * (doubleValueSize + 1);
		if (pcl::ProcRank() == 0)
		{
			// file head
			header = std::string("coordX");
			if (dim >= 2) header.append(std::string(",coordY"));
			if (dim >= 3) header.append(std::string(",coordZ"));
			header.append(std::string(",fluxX"));
			if (dim >= 2) header.append(std::string(",fluxY"));
			if (dim >= 3) header.append(std::string(",fluxZ"));
			header.append("\n");
			mySize += header.size();
		}

		unsigned long allSizesUpToMine;
		MPI_Scan(&mySize, &allSizesUpToMine, 1, MPI_LONG_LONG, MPI_SUM, m_mpiComm);
		unsigned long offset = allSizesUpToMine - mySize;

 		// write data at correct offset
		MPI_File_seek(fh, offset, MPI_SEEK_SET);

		if (pcl::ProcRank() == 0)
			MPI_File_write(fh, header.c_str(), header.size(), MPI_BYTE, &status);

		MathVector<dim> coords;
		MathVector<dim> fluxVec;
		for (size_t r = 0; r < nRow; ++r)
		{
			typename matrix_type::row_iterator cit = m_fluxMatrix.begin_row(r);
			typename matrix_type::row_iterator citEnd = m_fluxMatrix.end_row(r);
			for (; cit != citEnd && cit.index() < r; ++cit) // only sub-diagonal values count
			{
				if (cit.value()) // only sub-diagonal and non-zero values count
				{
					const size_t c = cit.index();

					std::ostringstream buf;
					buf << std::scientific << std::setprecision(precision) << std::showpos;

					UG_COND_THROW(dof2Vrt.find(r) == dof2Vrt.end(), "From vertex not found.");
					UG_COND_THROW(dof2Vrt.find(c) == dof2Vrt.end(), "To vertex not found.");
					const MathVector<dim>& coFrom = aaPos[dof2Vrt[r]];
					const MathVector<dim>& coTo = aaPos[dof2Vrt[c]];
					VecScaleAdd(coords, 0.5, coFrom, 0.5, coTo);
					VecScaleAdd(fluxVec, 1.0, coTo, -1.0, coFrom);
					VecNormalize(fluxVec, fluxVec);
					VecScale(fluxVec, fluxVec, cit.value());

					for (size_t d = 0; d < (size_t) dim; ++d)
						buf << coords[d] << ",";
					for (size_t d = 0; d < (size_t) dim-1; ++d)
						buf << fluxVec[d]*scale_factor << ",";
					buf << fluxVec[dim-1]*scale_factor << std::endl;

					MPI_File_write(fh, buf.str().c_str(), buf.str().size(), MPI_BYTE, &status);
				}
			}
		}

		// close file
		MPI_File_close(&fh);
	}
	else
#endif
	{
		std::ofstream outFile;
		outFile.open(filename.c_str(), std::ios_base::out);

		// write fluxes
		try
		{
			// write file header
			header = std::string("coordX");
			if (dim >= 2) header.append(std::string(",coordY"));
			if (dim >= 3) header.append(std::string(",coordZ"));
			header.append(std::string(",fluxX"));
			if (dim >= 2) header.append(std::string(",fluxY"));
			if (dim >= 3) header.append(std::string(",fluxZ"));
			outFile << header << std::endl;

			// write vectors
			const size_t precision = std::numeric_limits<number>::digits10;
			outFile << std::scientific << std::setprecision(precision) << std::showpos;
			MathVector<dim> coords;
			MathVector<dim> fluxVec;
			const size_t nRow = m_fluxMatrix.num_rows();
			for (size_t r = 0; r < nRow; ++r)
			{
				typename matrix_type::row_iterator cit = m_fluxMatrix.begin_row(r);
				typename matrix_type::row_iterator citEnd = m_fluxMatrix.end_row(r);
				for (; cit != citEnd && cit.index() < r; ++cit) // only sub-diagonal values count
				{
					if (cit.value()) // only sub-diagonal and non-zero values count
					{
						const size_t c = cit.index();

						std::ostringstream buf;
						buf << std::scientific << std::setprecision(precision) << std::showpos;

						const MathVector<dim>& coFrom = aaPos[dof2Vrt[r]];
						const MathVector<dim>& coTo = aaPos[dof2Vrt[c]];
						VecScaleAdd(coords, 0.5, coFrom, 0.5, coTo);
						VecScaleAdd(fluxVec, 1.0, coTo, -1.0, coFrom);
						VecNormalize(fluxVec, fluxVec);
						VecScale(fluxVec, fluxVec, cit.value());

						for (size_t d = 0; d < (size_t) dim; ++d)
							outFile << coords[d] << ",";
						for (size_t d = 0; d < (size_t) dim-1; ++d)
							outFile << fluxVec[d]*scale_factor << ",";
						outFile << fluxVec[dim-1]*scale_factor << std::endl;
					}
				}
			}
		}
		UG_CATCH_THROW("Output file" << filename << " could not be written to.");

		outFile.close();
	}

	// reset flag
	m_fluxMatrix.clear_and_free();
	m_bWriteFluxMatrix = false;
}



template <typename TGridFunction>
template <typename TFVGeom, typename Dummy>
FluxExporter<TGridFunction>::prep_elem_loop<TFVGeom, Dummy>::prep_elem_loop
(
	FluxExporter<TGridFunction>* flEx,
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
	FluxExporter<TGridFunction>* flEx,
	const ReferenceObjectID roid
)
{
	//	init upwind for element type
	if (!TFV1Geom<TElem, dim>::usesHangingNodes)
	{
		TFV1Geom<TElem, dim>& geo = GeomProvider<TFV1Geom<TElem, dim> >::get();
		if (!flEx->m_spConvShape->template set_geometry_type<TFV1Geom<TElem, dim> >(geo))
			UG_THROW("Cannot init upwind for element type.");
	}
}



template <typename TGridFunction>
template <typename TFVGeom, typename Dummy>
FluxExporter<TGridFunction>::prep_elem<TFVGeom, Dummy>::prep_elem
(
	FluxExporter<TGridFunction>* flEx,
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
	FluxExporter<TGridFunction>* flEx,
	GridObject* elem,
	const std::vector<MathVector<FluxExporter<TGridFunction>::dim> >& vCornerCoords
)
{
	TFV1Geom<TElem, dim>& geo = GeomProvider<TFV1Geom<TElem, dim> >::get();
	try {geo.update(elem, &vCornerCoords[0], flEx->m_sh.get());}
	UG_CATCH_THROW("Failed updating Finite Volume Geometry for elem.");

	if (TFV1Geom<TElem, dim>::usesHangingNodes)
	{
		if (flEx->m_spConvShape.valid())
			if (!flEx->m_spConvShape->template set_geometry_type<TFV1Geom<TElem, dim> >(geo))
				UG_THROW("Cannot init upwind for element type.");
	}

}



template <typename TGridFunction>
template<typename TFVGeom, typename Dummy>
FluxExporter<TGridFunction>::assemble_flux_elem<TFVGeom, Dummy>::assemble_flux_elem
(
	FluxExporter<TGridFunction>* flEx,
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

			// multiply by half the vector from box center to box center and add to flux
			size_t to = scvf.to();
			size_t from = scvf.from();
			VecScaleAdd(grad, 0.5, vDoFPos[to], -0.5, vDoFPos[from]);
			for (size_t i = 0; i < dim; ++i)
			{
				f(i, from) += scvf.weight(ip) * scalarFluxIP * grad[i];
				f(i, to)   += scvf.weight(ip) * scalarFluxIP * grad[i];
			}
		}
	}
}

template <typename TGridFunction>
template<typename TElem, int dim, template <class, int> class TFV1Geom, typename Dummy>
FluxExporter<TGridFunction>::assemble_flux_elem<TFV1Geom<TElem, dim>, Dummy>::assemble_flux_elem
(
	FluxExporter<TGridFunction>* flEx,
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

	// compute convection shapes
	size_t nIP = geo.num_scvf();
	MathVector<dim>* vVelAtIP = new MathVector<dim>[nIP];
	for (size_t i = 0; i < nIP; ++i)
	{
		VecSet(vVelAtIP[i], 0.0);
		for (size_t sh = 0; sh < geo.scvf(i).num_sh(); ++sh)
			VecScaleAppend(vVelAtIP[i], -flEx->m_convConst*u(_P_,sh), geo.scvf(i).global_grad(sh));
	}

	if (!flEx->m_spConvShape->update(&geo, vVelAtIP, NULL, false))
		UG_THROW("Cannot compute convection shapes.");
	IConvectionShapes<dim>& convShape = *flEx->m_spConvShape.get();

	//delete[] vDiffAtIP;
	delete[] vVelAtIP;

	// prepare local indices for box fluxes if they are written
	LocalIndices indBoxFlux;
	if (flEx->m_bWriteFluxMatrix)
		flEx->m_ddBoxFluxes->indices(elem, indBoxFlux, true);

	// loop SCVFs
	for (size_t s = 0; s < geo.num_scvf(); ++s)
	{
		const typename TFV1Geom<TElem, dim>::SCVF& scvf = geo.scvf(s);

		number scalarFluxIP;

		// diffusive term
		MathVector<dim> grad(0.0);
		for (size_t sh = 0; sh < scvf.num_sh(); ++sh)
			VecScaleAppend(grad, u(_C_,sh), scvf.global_grad(sh));
		scalarFluxIP = -flEx->m_diffConst * VecDot(grad, scvf.normal());

		// convective term
		number conv_flux = 0.0;
		for (size_t sh = 0; sh < convShape.num_sh(); ++sh)
			conv_flux += u(_C_, sh) * convShape(s, sh);
		scalarFluxIP += conv_flux;

		// multiply by vector between box centers and add to flux
		size_t from = scvf.from();
		const typename TFV1Geom<TElem, dim>::SCV& scvFrom = geo.scv(from);
		size_t to = scvf.to();
		const typename TFV1Geom<TElem, dim>::SCV& scvTo = geo.scv(to);

		VecScaleAdd(grad, 0.5, scvTo.global_ip(), -0.5, scvFrom.global_ip());
		for (size_t i = 0; i < (size_t) dim; ++i)
		{
			f(i, from) += scalarFluxIP * grad[i];
			f(i, to) += scalarFluxIP * grad[i];
		}

		// fill flux map if needed
		if (flEx->m_bWriteFluxMatrix)
		{
			size_t indFrom = indBoxFlux.index(_C_, from);
			size_t indTo = indBoxFlux.index(_C_, to);
			if (flEx->m_fluxMatrix.has_connection(indFrom, indTo))
				BlockRef(flEx->m_fluxMatrix(indFrom, indTo), 0, 0) += scalarFluxIP;
			else
				BlockRef(flEx->m_fluxMatrix(indFrom, indTo), 0, 0) = scalarFluxIP;
			if (flEx->m_fluxMatrix.has_connection(indTo, indFrom))
				BlockRef(flEx->m_fluxMatrix(indTo, indFrom), 0, 0) -= scalarFluxIP;
			else
				BlockRef(flEx->m_fluxMatrix(indTo, indFrom), 0, 0) = -scalarFluxIP;
		}
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



template <typename TGridFunction>
template <bool hanging, int order, typename TElem>
void FluxExporter<TGridFunction>::assemble
(
	SmartPtr<GridFunction<dom_type, CPUBlockAlgebra<dim> > > flux,
	SmartPtr<GridFunction<dom_type, CPUAlgebra> > vol,
	int si
)
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
	SmartPtr<GridFunction<dom_type, CPUBlockAlgebra<dim> > > flux,
	SmartPtr<GridFunction<dom_type, CPUAlgebra> > vol
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
				boost::mpl::for_each<ElemList>(AssembleWrapper<false, order>(this, flux, vol, si));
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
	SmartPtr<GridFunction<dom_type, CPUBlockAlgebra<dim> > > flux,
	SmartPtr<GridFunction<dom_type, CPUAlgebra> > vol,
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
SmartPtr<GridFunction<typename TGridFunction::domain_type, CPUBlockAlgebra<TGridFunction::domain_type::dim> > >
FluxExporter<TGridFunction>::calc_flux(number scale_factor)
{
	// set up approx space with flux function of the same order as solution grid function
	SmartPtr<ApproximationSpace<dom_type> > approxFlux =
		make_sp(new ApproximationSpace<dom_type>(m_u->approx_space()->domain(), AlgebraType(AlgebraType::CPU, dim)));

	std::vector<std::string> flux_cmp_names(dim);
	if (dim >= 1) flux_cmp_names[0] = std::string("flux_x");
	if (dim >= 2) flux_cmp_names[1] = std::string("flux_y");
	if (dim >= 3) flux_cmp_names[2] = std::string("flux_z");

	try {approxFlux->add(flux_cmp_names, m_lfeid, m_vSubset);}
	UG_CATCH_THROW("Failed to add functions to approximation space on specified subsets.");
	approxFlux->init_top_surface();

	// set up approx space for box volumes
	SmartPtr<ApproximationSpace<dom_type> > approxVol =
		make_sp(new ApproximationSpace<dom_type>(m_u->approx_space()->domain(), AlgebraType(AlgebraType::CPU, 1)));
	std::vector<std::string> vol_cmp_names(1);
	vol_cmp_names[0] = std::string("vol");
	try {approxVol->add(vol_cmp_names, m_lfeid, m_vSubset);}
	UG_CATCH_THROW("Failed to add functions to approximation space on specified subsets.");
	approxVol->init_top_surface();


	// create flux and volume grid functions from their respective approx spaces
	SmartPtr<GridFunction<dom_type, CPUBlockAlgebra<dim> > > flux
		= make_sp(new GridFunction<dom_type, CPUBlockAlgebra<dim> >(approxFlux));
	SmartPtr<GridFunction<dom_type, CPUAlgebra> > vol
		= make_sp(new GridFunction<dom_type, CPUAlgebra>(approxVol));
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
	if (m_hangingConstraintFlux.valid())
	{
		m_hangingConstraintFlux->set_approximation_space(approxFlux);
		m_hangingConstraintFlux->set_ass_tuner(make_sp(new AssemblingTuner<CPUBlockAlgebra<dim> >())); // dummy
		GridFunction<dom_type, CPUBlockAlgebra<dim> > dummy(approxFlux, false); // not needed in adjust_defect
		m_hangingConstraintFlux->adjust_defect(*flux, dummy, flux->dd(), CT_HANGING);
	}
	if (m_hangingConstraintVol.valid())
	{
		m_hangingConstraintVol->set_approximation_space(approxVol);
		m_hangingConstraintVol->set_ass_tuner(make_sp(new AssemblingTuner<CPUAlgebra>())); // dummy
		GridFunction<dom_type, CPUAlgebra> dummy(approxVol, false); // not needed in adjust_defect
		m_hangingConstraintVol->adjust_defect(*vol, dummy, vol->dd(), CT_HANGING);
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
#ifdef UG_CPU_5
	#ifdef UG_DIM_1
		template class FluxExporter<GridFunction<Domain1d, CPUBlockAlgebra<5> > >;
	#endif
	#ifdef UG_DIM_2
		template class FluxExporter<GridFunction<Domain2d, CPUBlockAlgebra<5> > >;
	#endif
	#ifdef UG_DIM_3
		template class FluxExporter<GridFunction<Domain3d, CPUBlockAlgebra<5> > >;
	#endif
#endif


} // namespace nernst_planck
} // namespace ug
