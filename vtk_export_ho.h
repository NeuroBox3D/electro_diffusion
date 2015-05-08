/*
 * vtk_export_ho.cpp
 *
 *  Created on: 04.05.2015
 *      Author: mbreit
 */

#include "lib_disc/function_spaces/grid_function_global_user_data.h"			// GlobalGridFunctionNumberData
#include "lib_disc/function_spaces/dof_position_util.h"							// DoFPosition
#include "lib_grid/parallelization/parallel_refinement/parallel_refinement.h"	// ParallelGlobalRefiner_MultiGrid
#include "lib_disc/io/vtkoutput.h"												// VTKOutput

namespace ug {
namespace nernst_planck {


/** Make sure that aNewVrt is attached to srcMesh->grid() and contains a
 * pointer to a valid vertex in destMesh for each selected vertex in srcMesh.
 * Also make sure that all vertices belonging to a selected element have been
 * selected, too.*/
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

		TElem* ne = *destGrid.create_by_cloning(e, vrts);
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
		Vertex* nv = *destGrid.create_by_cloning(v);
		aaNewVrt[v] = nv;
		aaPosDest[nv] = aaPosSrc[v];
		destSH->assign_subset(nv, srcSH->get_subset_index(v));
	}

	CopySelectedElements<TDomain,Edge>(destDom, srcDom, sel, aNewVrt);
	CopySelectedElements<TDomain,Face>(destDom, srcDom, sel, aNewVrt);
	CopySelectedElements<TDomain,Volume>(destDom, srcDom, sel, aNewVrt);

	srcGrid.detach_from_vertices(aNewVrt);
}



template <typename TGridFunction>
//SmartPtr<TGridFunction> vtk_export_ho
void vtk_export_ho
(
	SmartPtr<TGridFunction> u,
	std::vector<std::string> vFct,
	size_t order,
	SmartPtr<VTKOutput<TGridFunction::domain_type::dim> > vtkOutput,
	const char* filename
)
{
	typedef typename TGridFunction::domain_type dom_type;
	typedef typename dom_type::position_attachment_type position_attachment_type;
	typedef typename TGridFunction::approximation_space_type approx_space_type;
	static const int dim = dom_type::dim;
	typedef typename TGridFunction::element_type elem_type;

	SmartPtr<approx_space_type> srcApproxSpace = u->approx_space();
	SmartPtr<dom_type> srcDom = srcApproxSpace->domain();

	// select surface elements in old grid
	MultiGrid& srcGrid = *srcDom->grid();
	Selector srcSel(srcGrid);
	srcSel.select(u->template begin<elem_type>(), u->template end<elem_type>());

	// create new domain
	SmartPtr<dom_type> destDom = make_sp(new dom_type());

	// copy grid from old domain to new domain
	CopySelected(destDom, srcDom, srcSel);

	// refine
#ifdef UG_PARALLEL
	ParallelGlobalRefiner_MultiGrid refiner(*destDom->distributed_grid_manager());
#else
	MultiGrid& destGrid = *destDom->grid();
	GlobalMultiGridRefiner refiner(destGrid);
#endif
	size_t numRefs = std::ceil(std::log2(order));
	for (size_t iref = 0; iref < numRefs; ++iref)
		refiner.refine();

	// retain function group for functions being exported
	FunctionGroup fg(srcApproxSpace->dof_distribution_info(), vFct);

	// create approx space and add functions
	SmartPtr<approx_space_type> destApproxSpace =
		make_sp(new approx_space_type(destDom));

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

	SmartPtr<TGridFunction> u_new = make_sp(new TGridFunction(destApproxSpace));

	// interpolate onto new grid
	for (size_t fct = 0; fct < fg.size(); ++fct)
	{
		UG_LOG("fct: " << fct << "\n");
		GlobalGridFunctionNumberData<TGridFunction> ggfnd =
			GlobalGridFunctionNumberData<TGridFunction>(u, fg.name(fct));

		// iterate over DoFs in new function and evaluate
		typedef typename DoFDistribution::dim_traits<dim>::grid_base_object elem_type;
		typedef typename DoFDistribution::traits<elem_type>::const_iterator const_iter_type;

		const_iter_type elem_iter = u_new->template begin<elem_type>();
		const_iter_type iterEnd = u_new->template end<elem_type>();

		// TODO: better not iterate over elements, but over vertices etc. like in Dirichlet bnd
		std::vector<DoFIndex> ind;
		for (; elem_iter != iterEnd; ++elem_iter)
		{
			u_new->dof_indices(*elem_iter, fct, ind);

			// get dof positions
			const LFEID lfeid = u_new->dof_distribution()->lfeid(fg[fct]);
			std::vector<MathVector<dim> > globPos;
			DoFPosition(globPos, *elem_iter, *destDom, lfeid);

			UG_ASSERT(globPos.size() == ind.size(), "#DoF mismatch");

			// write values in new grid function
			for (size_t dof = 0; dof < ind.size(); ++dof)
			{
				if (!ggfnd.evaluate(DoFRef(*u_new, ind[dof]), globPos[dof]))
				{
					DoFRef(*u_new, ind[dof]) = std::numeric_limits<number>::quiet_NaN();
					//UG_THROW("Interpolation onto new grid did not succeed.\n"
					//		 "DoF with coords " << globPos[dof] << " is out of range.");
				}
				if (fct == 4 && std::fabs(DoFRef(*u_new, ind[dof])) > 1.0)
					UG_LOG("phi " << globPos[dof] << " = " << DoFRef(*u_new, ind[dof]) << "\n");
			}
		}
	}

#ifdef UG_PARALLEL
	// copy storage type and layouts
	u_new->set_storage_type(u->get_storage_mask());
	u_new->set_layouts(u->layouts());
#endif

	// return new grid function
	vtkOutput->print(filename, *u_new);
}

} // end namespace ug
} // end namespace nernst_planck
