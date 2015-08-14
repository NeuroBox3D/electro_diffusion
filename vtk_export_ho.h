/*
 * vtk_export_ho.h
 *
 *  Created on: 04.05.2015
 *      Author: mbreit
 */

#ifndef UG__PLUGINS__EXPERIMENTAL__NERNST_PLANCK__VTK_EXPORT_HO_H_
#define UG__PLUGINS__EXPERIMENTAL__NERNST_PLANCK__VTK_EXPORT_HO_H_

#include "lib_disc/function_spaces/grid_function_global_user_data.h"			// GlobalGridFunctionNumberData
#include "lib_disc/function_spaces/dof_position_util.h"							// DoFPosition
#include "lib_grid/parallelization/parallel_refinement/parallel_refinement.h"	// ParallelGlobalRefiner_MultiGrid
#include "lib_disc/io/vtkoutput.h"												// VTKOutput

namespace ug {
namespace nernst_planck {


template <typename TDomain, class TElem>
inline void CopySelectedElements
(
	SmartPtr<TDomain> destDom,
	SmartPtr<TDomain> srcDom,
	Selector& sel,
	AVertex& aNewVrt
)
{
	typedef typename TDomain::subset_handler_type SH_type;

	MultiGrid& srcGrid = *srcDom->grid();
	MultiGrid& destGrid = *destDom->grid();

	ConstSmartPtr<SH_type> srcSH = srcDom->subset_handler();
	SmartPtr<SH_type> destSH = destDom->subset_handler();

	Grid::VertexAttachmentAccessor<AVertex> aaNewVrt(srcGrid, aNewVrt);

	CustomVertexGroup vrts;
	typedef typename Grid::traits<TElem>::iterator iter_t;

	for (iter_t eiter = sel.begin<TElem>();	eiter != sel.end<TElem>(); ++eiter)
	{
		TElem* e = *eiter;
		vrts.resize(e->num_vertices());
		for (size_t iv = 0; iv < e->num_vertices(); ++iv)
			vrts.set_vertex(iv, aaNewVrt[e->vertex(iv)]);

		TElem* ne;
		try {ne = *destGrid.create_by_cloning(e, vrts);}
		UG_CATCH_THROW("New element could not be created.");
		destSH->assign_subset(ne, srcSH->get_subset_index(e));
	}
}

template <typename TDomain>
inline void CopySelected
(
	SmartPtr<TDomain> destDom,
	SmartPtr<TDomain> srcDom,
	Selector& sel
)
{
	typedef typename TDomain::position_attachment_type APos_type;
	typedef typename TDomain::subset_handler_type SH_type;

	MultiGrid& srcGrid = *srcDom->grid();
	MultiGrid& destGrid = *destDom->grid();

	ConstSmartPtr<SH_type> srcSH = srcDom->subset_handler();
	SmartPtr<SH_type> destSH = destDom->subset_handler();

	APos_type& aPos = srcDom->position_attachment();
	UG_COND_THROW(!srcGrid.has_vertex_attachment(aPos), "Position attachment required.");
	Grid::VertexAttachmentAccessor<APos_type> aaPosSrc(srcGrid, aPos);
	if (!destGrid.has_vertex_attachment(aPos))
		destGrid.attach_to_vertices(aPos);
	Grid::VertexAttachmentAccessor<APos_type> aaPosDest(destGrid, aPos);

	AVertex aNewVrt;
	srcGrid.attach_to_vertices(aNewVrt);
	Grid::VertexAttachmentAccessor<AVertex> aaNewVrt(srcGrid, aNewVrt);

	for (int si = destSH->num_subsets(); si < srcSH->num_subsets(); ++si)
		destSH->subset_info(si) = srcSH->subset_info(si);

	SelectAssociatedGridObjects(sel);

	// create new vertices in destGrid
	for (VertexIterator viter = sel.begin<Vertex>(); viter != sel.end<Vertex>(); ++viter)
	{
		Vertex* v = *viter;
		Vertex* nv;
		try {nv = *destGrid.create_by_cloning(v);}
		UG_CATCH_THROW("New vertex could not be created.");
		aaNewVrt[v] = nv;
		aaPosDest[nv] = aaPosSrc[v];
		destSH->assign_subset(nv, srcSH->get_subset_index(v));
	}

	CopySelectedElements<TDomain,Edge>(destDom, srcDom, sel, aNewVrt);
	CopySelectedElements<TDomain,Face>(destDom, srcDom, sel, aNewVrt);
	CopySelectedElements<TDomain,Volume>(destDom, srcDom, sel, aNewVrt);

	srcGrid.detach_from_vertices(aNewVrt);
}


template <typename TElem, typename TGridFunction>
inline void interpolate_from_original_fct
(
	SmartPtr<TGridFunction> u_new,
	const GlobalGridFunctionNumberData<TGridFunction>& u_orig,
	size_t fct,
	const LFEID& lfeid
)
{
	typedef typename TGridFunction::domain_type dom_type;
	static const int dim = dom_type::dim;
	typedef typename DoFDistribution::traits<TElem>::const_iterator const_iter_type;

	const_iter_type elem_iter = u_new->template begin<TElem>();
	const_iter_type iterEnd = u_new->template end<TElem>();

	std::vector<DoFIndex> ind;
	for (; elem_iter != iterEnd; ++elem_iter)
	{
		u_new->inner_dof_indices(*elem_iter, fct, ind);

		// get dof positions
		std::vector<MathVector<dim> > globPos;
		InnerDoFPosition<dom_type>(globPos, *elem_iter, *u_new->domain(), lfeid);

		UG_ASSERT(globPos.size() == ind.size(), "#DoF mismatch");

		// write values in new grid function
		for (size_t dof = 0; dof < ind.size(); ++dof)
		{
			if (!u_orig.evaluate(DoFRef(*u_new, ind[dof]), globPos[dof]))
			{
				DoFRef(*u_new, ind[dof]) = std::numeric_limits<number>::quiet_NaN();
				//UG_THROW("Interpolation onto new grid did not succeed.\n"
				//		 "DoF with coords " << globPos[dof] << " is out of range.");
			}
		}
	}
}


/// export a solution from a high-order ansatz space to vtk file(s)
/**
 *  This function will create a temporary domain, copy all elements from the domain
 *  which the grid function u is defined on to the temporary domain and then refine
 *  the resulting grid until it has at least as many vertices as the original grid
 *  functions has unknowns (e.g. a grid for a function of order 2 would be refined
 *  once, a grid for a function of order 4 would be refined twice, and so on).
 *  After refinement, a temporary grid function of order 1 (Lagrange) is defined on
 *  the refined grid and its values interpolated from the original function u.
 *  The temporary first-order grid function is then exported to vtk using the usual
 *  mechanisms.
 *
 * @param u			original high-order grid function to be exported
 * @param vFct		vector of function names (contained in grid function) to be exported
 * @param order		order to be used
 * @param vtkOutput	VTKOutput object to use for export of linearized function
 * @param filename	file name to be used in export
 *
 * @todo The order parameter might be left out and determined automatically from the
 * 		 grid function.
 */
template <typename TGridFunction>
void vtk_export_ho
(
	SmartPtr<TGridFunction> u,
	const std::vector<std::string>& vFct,
	size_t order,
	SmartPtr<VTKOutput<TGridFunction::domain_type::dim> > vtkOutput,
	const char* filename,
	size_t step,
	number time
)
{
	typedef typename TGridFunction::domain_type dom_type;
	typedef typename dom_type::position_attachment_type position_attachment_type;
	typedef typename TGridFunction::approximation_space_type approx_space_type;
	typedef typename TGridFunction::element_type elem_type;

	SmartPtr<approx_space_type> srcApproxSpace = u->approx_space();
	SmartPtr<dom_type> srcDom = srcApproxSpace->domain();

	// select surface elements in old grid
	MultiGrid& srcGrid = *srcDom->grid();
	Selector srcSel(srcGrid);
	srcSel.select(u->template begin<elem_type>(), u->template end<elem_type>());

	// create new domain
	dom_type* dom_ptr;
	try	{dom_ptr = new dom_type();}
	UG_CATCH_THROW("Temporary domain could not be created.");
	SmartPtr<dom_type> destDom = make_sp(dom_ptr);

	// copy grid from old domain to new domain
	try	{CopySelected(destDom, srcDom, srcSel);}
	UG_CATCH_THROW("Temporary grid could not be created.");

	// refine
#ifdef UG_PARALLEL
	ParallelGlobalRefiner_MultiGrid refiner(*destDom->distributed_grid_manager());
#else
	MultiGrid& destGrid = *destDom->grid();
	GlobalMultiGridRefiner refiner(destGrid);
#endif
	size_t numRefs = ceil(log2(order));
	for (size_t iref = 0; iref < numRefs; ++iref)
	{
		try	{refiner.refine();}
		UG_CATCH_THROW("Refinement step " << iref << " could not be carried out.");
	}

	// retain function group for functions being exported
	FunctionGroup fg(srcApproxSpace->dof_distribution_info(), vFct);

	// create approx space and add functions
	approx_space_type* approx_ptr;
	try {approx_ptr = new approx_space_type(destDom);}
	UG_CATCH_THROW("Temporary approximation space could not be created.");
	SmartPtr<approx_space_type> destApproxSpace = make_sp(approx_ptr);

	for (size_t fct = 0; fct < fg.size(); ++fct)
	{
		if (fg.function_pattern()->is_def_everywhere(fg.unique_id(fct)))
			destApproxSpace->add(fg.name(fct), "Lagrange", 1);
		else
		{
			int num_subsets = fg.function_pattern()->num_subsets();
			std::string subsets;
			for (int si = 0; si < num_subsets; ++si)
			{
				if (fg.function_pattern()->is_def_in_subset(fg.unique_id(fct), si))
					subsets += std::string(",") + fg.function_pattern()->subset_name(si);
			}
			if (!subsets.empty())
				subsets.erase(0,1);	// delete leading ","

			destApproxSpace->add(fg.name(fct), "Lagrange", 1, subsets.c_str());
		}
	}
	destApproxSpace->init_top_surface();

	TGridFunction* gridFct_ptr;
	try {gridFct_ptr = new TGridFunction(destApproxSpace);}
	UG_CATCH_THROW("Temporary grid function could not be created.");
	SmartPtr<TGridFunction> u_new = make_sp(gridFct_ptr);

	// interpolate onto new grid
	for (size_t fct = 0; fct < fg.size(); ++fct)
	{
		const LFEID lfeid = u_new->dof_distribution()->lfeid(fg[fct]);

		GlobalGridFunctionNumberData<TGridFunction> ggfnd =
			GlobalGridFunctionNumberData<TGridFunction>(u, fg.name(fct));

		// iterate over DoFs in new function and evaluate
		// should be vertices only for Lagrange-1
		if (u_new->max_dofs(VERTEX))
			interpolate_from_original_fct<Vertex, TGridFunction>(u_new, ggfnd, fct, lfeid);
		if (u_new->max_dofs(EDGE))
			interpolate_from_original_fct<Edge, TGridFunction>(u_new, ggfnd, fct, lfeid);
		if (u_new->max_dofs(FACE))
			interpolate_from_original_fct<Face, TGridFunction>(u_new, ggfnd, fct, lfeid);
		if (u_new->max_dofs(VOLUME))
			interpolate_from_original_fct<Volume, TGridFunction>(u_new, ggfnd, fct, lfeid);
	}

#ifdef UG_PARALLEL
	// copy storage type and layouts
	u_new->set_storage_type(u->get_storage_mask());
	u_new->set_layouts(u->layouts());
#endif

	// return new grid function
	vtkOutput->print(filename, *u_new, step, time);
}


template <typename TGridFunction>
void vtk_export_ho
(
	SmartPtr<TGridFunction> u,
	const std::vector<std::string>& vFct,
	size_t order,
	SmartPtr<VTKOutput<TGridFunction::domain_type::dim> > vtkOutput,
	const char* filename
	)
{
	vtk_export_ho(u, vFct, order, vtkOutput, filename, 0, 0.0);
}


} // end namespace ug
} // end namespace nernst_planck


#endif // UG__PLUGINS__EXPERIMENTAL__NERNST_PLANCK__VTK_EXPORT_HO_H_
