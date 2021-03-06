/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Markus Breit
 * Creation date: 2014-05-28
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

#include "bridge/util.h"                                    // for RegisterCommon, RegisterDimensionDependent
#include "bridge/util_domain_algebra_dependent.h"           // for RegisterDomainAlgebraDependent
#include "lib_algebra/operator/interface/preconditioner.h"  // for IPreconditioner
#include "lib_disc/function_spaces/grid_function.h"         // for GridFunction

#include "np_config.h"                                      // for #defines
#include "extension/domain1d_solution_adjuster.h"           // for Domain1dSolutionAdjuster
#include "util/edl_1d.h"                                    // for EDLSimulation
#include "extension/extension_refMarkAdjuster.h"            // for ExtensionRefMarkAdjuster
#include "util/flux_exporter.h"                             // for FluxExporter
#include "interface/interface1d_fv.h"                       // for IInterface1D, Interface1D, AdditiveInterface1D
#ifdef UG_PARALLEL
	#include "interface/intf_distro_adjuster.h"             // for PNPDistroManager, set_distro_adjuster
#endif
#include "interface/intf_refMarkAdjuster.h"                 // for InterfaceRefMarkAdjuster
#include "morpho_gen.h"                                     // for MorphoGen
#include "util/nernst_planck_util.h"                        // for adjust_geom_after_refinement, exportSolution
#include "util/neck_recorder.h"								// for current and voltage recordings across spine neck
#include "util/order.h"                                     // for reorder_dof_distros_lex, reorder_dofs
#include "pnp_smoother.h"                                   // for PNPSmoother, PNP_ILU
#include "pnp_upwind.h"                                     // for PNPUpwind
#include "util/refinement_error_estimator.h"                // for RefinementErrorEstimator


using namespace std;
using namespace ug::bridge;


namespace ug {
namespace nernst_planck {


/**
 *  \defgroup plugin_nernst_planck Plugin nernst_planck
 *  \ingroup plugins_experimental
 *  This is a plugin for Poisson-Nernst-Planck functionality.
 *  \{
 */

/**
 * Class exporting the functionality. All functionality that is to
 * be used in scripts or visualization must be registered here.
 */
struct Functionality
{

/**
 * Function called for the registration of Domain and Algebra dependent parts.
 * All Functions and Classes depending on both Domain and Algebra
 * are to be placed here when registering. The method is called for all
 * available Domain and Algebra types, based on the current build options.
 *
 * @param reg		registry
 * @param grp		group for sorting of functionality
 */
template <typename TDomain, typename TAlgebra>
static void DomainAlgebra(Registry& reg, string grp)
{
	string suffix = GetDomainAlgebraSuffix<TDomain,TAlgebra>();
	string tag = GetDomainAlgebraTag<TDomain,TAlgebra>();

	static const int dim = TDomain::dim;
	typedef GridFunction<TDomain, TAlgebra> TGridFunction;

	// write residuals to file function for param optimization
	{
		reg.add_function("write_residuals_to_file", &writeResidualsToFile<TGridFunction>, grp.c_str(),
			"", "solution#reference solution#outFileName",
			"writes residuals in every dof to file as required by Ivo's optimization routines "
			"and returns squared 2-norm of residual vector");
	}

	// interface (in the sense of programming) for the nD/1D interface (in the sense of manifold) class
	{
		typedef Interface1D<TDomain, TAlgebra> T;
		typedef IDomainConstraint<TDomain, TAlgebra> TBase;
		typedef IInterface1D TBase2;
		string name = string("Interface1D").append(suffix);
		reg.add_class_<T, TBase, TBase2>(name, grp)
			.add_method("update", static_cast<void (T::*) ()>(&T::update))
			.add_method("determine_subset_indices", &T::determine_subset_indices)
			.add_method("set_approx_space", &T::set_approx_space);
			//.add_method("check_values_at_interface", &T::check_values_at_interface, "", "solution grid function", "", "");
		reg.add_class_to_group(name, "Interface1D", tag);
	}

	// additive nD/1D interface
	{
		typedef AdditiveInterface1D<TDomain, TAlgebra> T;
		typedef Interface1D<TDomain, TAlgebra> TBase;
		string name = string("AdditiveInterface1D").append(suffix);
		reg.add_class_<T, TBase >(name, grp)
			.template add_constructor<void (*)(const char*, const char*, const char*, const char*, std::vector<number>)>
				("function(s)#constrained subset#high-dim interface node subset#"
				 "one-dim interface node subset")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "AdditiveInterface1D", tag);
	}

	// multiplicative nD/1D interface
	{
		typedef MultiplicativeInterface1D<TDomain, TAlgebra> T;
		typedef Interface1D<TDomain, TAlgebra> TBase;
		string name = string("MultiplicativeInterface1D").append(suffix);
		reg.add_class_<T, TBase >(name, grp)
			.template add_constructor<void (*)(const char*, const char*, const char*, const char*, std::vector<number>)>
				("function(s)#constrained subset#high-dim interface node subset#"
				 "one-dim interface node subset")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "MultiplicativeInterface1D", tag);
	}

	// Domain1dSolutionAdjuster
	{
		typedef Domain1dSolutionAdjuster<TDomain, TAlgebra> T;
		string name = string("Domain1dSolutionAdjuster").append(suffix);
		reg.add_class_<T>(name, grp)
			.template add_constructor<void (*)()> ()
			.add_method("add_constrained_subset", &T::add_constrained_subset, "", "subset name", "", "")
			.add_method("add_constrainer_subset", &T::add_constrainer_subset, "", "subset name", "", "")
			.add_method("set_sorting_direction", &T::set_sorting_direction, "", "direction vector", "", "")
			.add_method("adjust_solution", &T::adjust_solution, "", "solution grid function", "", "")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "Domain1dSolutionAdjuster", tag);
	}

	// flux field export
	{
		typedef FluxExporter<TGridFunction> T;
		string name = string("FluxExporter").append(suffix);
		reg.add_class_<T>(name, grp)
			.template add_constructor<void (*)(SmartPtr<TGridFunction>, std::string, std::string)> ("grid function, species name, potential name")
			.add_method("set_diff_const", &T::set_diff_const, "", "diffusion constant", "", "")
			.add_method("set_conv_const", &T::set_conv_const, "", "convection constant", "", "")
			.add_method("set_quad_order", &T::set_quad_order, "", "order", "", "")
			.add_method("set_hanging_constraint", &T::set_hanging_constraint, "", "constraint for hanging nodes", "", "")
			.add_method("set_upwind", &T::set_upwind, "", "upwind scheme", "")
			.add_method("set_subsets", static_cast<void (T::*) (const std::vector<std::string>&)>(&T::set_subsets), "", "subsets as vector of string", "", "")
			.add_method("set_subsets", static_cast<void (T::*) (const char* cSubsets)>(&T::set_subsets), "", "subsets as comma-separated c-string", "", "")
			.add_method("write_flux",
						static_cast<void (T::*)
						(
							SmartPtr<VTKOutput<dim> > vtkOutput,
							std::string filename,
							size_t step,
							number time,
							std::string fluxName,
							number scale_factor
						)>(&T::write_flux),
						"", "vtkOutput object#filename#step#time#flux name#scale factor", "", "")
			.add_method("write_flux",
						static_cast<void (T::*)
						(
							SmartPtr<VTKOutput<dim> > vtkOutput,
							std::string filename,
							size_t step,
							number time,
							std::string fluxName
						)>(&T::write_flux),
						"", "vtkOutput object#filename#step#time#flux name", "", "")
			.add_method("write_flux",
						static_cast<void (T::*)
						(
							SmartPtr<VTKOutput<dim> > vtkOutput,
							std::string filename,
							size_t step,
							number time,
							number scale_factor
						)>(&T::write_flux),
						"", "vtkOutput object#filename#step#time#scale factor", "", "")
			.add_method("write_flux",
						static_cast<void (T::*)
						(
							SmartPtr<VTKOutput<dim> > vtkOutput,
							std::string filename,
							size_t step,
							number time
						)>(&T::write_flux),
						"", "vtkOutput object#filename#step#time", "", "")
			.add_method("write_box_fluxes",
						static_cast<void (T::*)
						(
							std::string filename,
							size_t step,
							number time,
							std::string fluxName,
							number scale_factor
						)>(&T::write_box_fluxes),
						"", "filename#step#time#flux name#scale factor", "", "")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "FluxExporter", tag);
	}

	// PNP smoothers
	{
		typedef PNPSmoother<TDomain, TAlgebra, ILU> TILU;
		typedef IPreconditioner<TAlgebra> TBase;
		string name = string("PNP_ILU").append(suffix);
		reg.add_class_<TILU, TBase>(name, grp)
			.template add_constructor<void (*)(SmartPtr<ApproximationSpace<TDomain> >)> ("approx space")
			.add_method("set_method", &TILU::set_method, "", "blocking method: 0 for non-overlapping blocks", "", "")
			.add_method("add_charge_surface_pair", &TILU::add_charge_surface_pair, "", "charged surface subset name#volume subset name", "", "")
			.add_method("set_parallelization_strategy", &TILU::set_parallelization_strategy, "",
				"0 for unique matrix and defect; 1 for consistent matrix and additive defect", "", "")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "PNP_ILU", tag);

		typedef PNPSmoother<TDomain, TAlgebra, GaussSeidel> TGS;
		typedef IPreconditioner<TAlgebra> TBase;
		name = string("PNP_GS").append(suffix);
		reg.add_class_<TGS, TBase>(name, grp)
			.template add_constructor<void (*)(SmartPtr<ApproximationSpace<TDomain> >)> ("approx space")
			.add_method("set_method", &TGS::set_method, "", "blocking method: 0 for non-overlapping blocks", "", "")
			.add_method("add_charge_surface_pair", &TGS::add_charge_surface_pair, "", "charged surface subset name#volume subset name", "", "")
			.add_method("set_parallelization_strategy", &TGS::set_parallelization_strategy, "",
				"0 for unique matrix and defect; 1 for consistent matrix and additive defect", "", "")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "PNP_GS", tag);
	}

	// neck recorder
	{
		typedef NeckRecorder<TDomain, TAlgebra> T;
		std::string name = std::string("NeckRecorder");
		reg.add_class_<T>(name+suffix, grp)
			.template add_constructor<void (*)(SmartPtr<ApproximationSpace<TDomain> >)> ("approx space")
			.add_method("add_measurement_zone", &T::add_measurement_zone, "", "z coordinate # name",
				"add a measurement manifold at given z location", "")
			.add_method("set_cytosolic_subset", &T::set_cytosolic_subset, "", "cytosolic subset name",
				"set the subset that recordings are to be taken in", "")
			.add_method("set_upwind", &T::set_upwind, "", "upwind", "", "")
			.add_method("set_diffusion_constants", &T::set_diffusion_constants, "",
				"constants in the order K, Na, Cl, A", "set diffusion constants", "")
			.add_method("set_convection_constants", &T::set_convection_constants, "",
				"constants in the order K, Na, Cl, A", "set convection constants", "")
			.add_method("set_temperature", &T::set_temperature, "", "temperature", "set temperature", "")
			.add_method("set_record_individual_currents", &T::set_record_individual_currents, "", "", "", "")
			.add_method("record_current", &T::record_current, "", "output file name#time#solution#scale factor",
				"record current through measurement manifold to file (each measurement zone separately)", "")
			.add_method("record_potential", &T::record_potential, "", "output file name#time#solution",
				"record averaged potential to file (each measurement zone separately)", "")
			.add_method("record_concentrations", &T::record_concentrations, "", "output file name#time#solution",
				"record averaged concentrations to file (each measurement zone separately)", "")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name+suffix, name, tag);
	}

	// RefinementErrorEstimator
	{
		typedef RefinementErrorEstimator<TDomain, TAlgebra> T;
		string nameBase = string("RefinementErrorEstimator");
		string name = nameBase;	name.append(suffix);
		reg.add_class_<T>(name, grp)
			.add_constructor()
			.add_method("compute_elementwise_errors", &T::compute_elementwise_errors, "", "uFine # uCoarse # cmp", "", "")
			.add_method("mark_with_strategy", &T::mark_with_strategy, "", "refiner # strategy", "", "")
#ifdef UG_CPU_1
			.add_method("error_grid_function", &T::error_grid_function, "", "domain", "", "")
#endif
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, nameBase, tag);
	}
}

/**
 * Function called for the registration of Domain dependent parts.
 * All Functions and Classes depending on the Domain
 * are to be placed here when registering. The method is called for all
 * available Domain types, based on the current build options.
 *
 * @param reg		registry
 * @param grp		group for sorting of functionality
 */
template <typename TDomain>
static void Domain(Registry& reg, string grp)
{
	//const int dim = TDomain::dim;
	string suffix = GetDomainSuffix<TDomain>();
	string tag = GetDomainTag<TDomain>();

	/*
	* // general syntax:
	* reg.add_function( "NameAppearingInScript", PointerToFunction, "group", "returnInformation", "paramInformation1#paramInformation2#...", "tooltip", "help" )
	* reg.add_class_<ClassName>("ClassName", grp, "tooltip")
	* .add_method("method_name", &ClassName::method_name, "returnInformation", "paramInformation1#paramInformation2#...", "tooltip", "help")
	*
	* // parameter information for VRL
	* "name1 | style1 | options1 # name2 | style2 | options2 # ... "
	*/

	// adjust geometry after refinement
	reg.add_function("adjust_geom_after_refinement", &adjust_geom_after_refinement<TDomain>, grp.c_str(),
					 "", "approximation space#full-dim interface node subset name"
					 "#1d interface node subset name", "adjusts location of full-dim interface node after"
					 "global refinement of the interface");
	reg.add_function("reorder_dofs", &reorder_dofs<TDomain>, grp.c_str(),
					 "", "approximation space#constrained subsets", "Re-orders DoFs in such a way that the "
					 "constrained indices are last in the order.");
	reg.add_function("reorder_dof_distros_lex", &reorder_dof_distros_lex<TDomain>, grp.c_str(),
					 "", "approximation space", "Re-orders DoFs in all the approxSpace's DoFDistros "
					 "according to lexicographical ordering as proposed by Rose et al.");

#if 0
	// 1D PNP FV1
	{
		typedef PNP1D_FV1<TDomain> T;
		typedef IElemDisc<TDomain> TBase;
		string name = string("PNP1D_FV1").append(suffix);
		reg.add_class_<T, TBase >(name, grp)
			.template add_constructor<void (*)(SmartPtr<ApproximationSpace<TDomain> >, const char*, const char*)>("function(s)#subset(s)")
			.template add_constructor<void (*)(SmartPtr<ApproximationSpace<TDomain> >, std::vector<std::string>, std::vector<std::string>)>
				("function(s)#subset(s)")
			.add_method("set_ions", &T::set_ions, "", "vector of ion species names", "", "")
			.add_method("set_valencies", &T::set_valencies, "", "vector of valencies", "", "")
			.add_method("set_reversal_potential", &T::set_reversal_potentials, "", "vector of reversal potentials", "", "")
			.add_method("set_specific_conductances", &T::set_specific_conductances, "", "vector of specific conductances", "", "")
			.add_method("set_specific_capacities", &T::set_specific_capacities, "", "vector of specific capacities", "", "")
			.add_method("set_diffusion_constants", &T::set_diffusion_constants, "", "vector of ion diffusion constants", "", "")
			.add_method("set_permettivities", &T::set_permettivities, "", "permettivity value dendrite (eps_dend*eps_0)#permettivity value membrane (eps_mem*eps_0)", "", "")
			.add_method("set_membrane_thickness", &T::set_membrane_thickness, "", "membrane thickness", "", "")
			.add_method("set_dendritic_radius", &T::set_dendritic_radius, "", "radius", "", "")
			.add_method("set_rtf", &T::set_rtf, "", "universal gas constant#temperature#Faraday constant", "Set natural constants in dimensions of choice.", "")
			.add_method("set_represented_dimension", &T::set_represented_dimension, "", "represented dimension (either 2 or 3)", "")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "PNP1D_FV1", tag);
	}
	// 1D PNP FV
	{
		typedef PNP1D_FV<TDomain> T;
		typedef IElemDisc<TDomain> TBase;
		string name = string("PNP1D_FV").append(suffix);
		reg.add_class_<T, TBase >(name, grp)
			.template add_constructor<void (*)(SmartPtr<ApproximationSpace<TDomain> >, const char*, const char*)>("function(s)#subset(s)")
			.template add_constructor<void (*)(SmartPtr<ApproximationSpace<TDomain> >, std::vector<std::string>, std::vector<std::string>)>
				("function(s)#subset(s)")
			.add_method("set_ions", &T::set_ions, "", "vector of ion species names", "", "")
			.add_method("set_valencies", &T::set_valencies, "", "vector of valencies", "", "")
			.add_method("set_reversal_potential", &T::set_reversal_potentials, "", "vector of reversal potentials", "", "")
			.add_method("set_specific_conductances", &T::set_specific_conductances, "", "vector of specific conductances", "", "")
			.add_method("set_specific_capacities", &T::set_specific_capacities, "", "vector of specific capacities", "", "")
			.add_method("set_diffusion_constants", &T::set_diffusion_constants, "", "vector of ion diffusion constants", "", "")
			.add_method("set_permettivities", &T::set_permettivities, "", "permettivity value dendrite (eps_dend*eps_0)#permettivity value membrane (eps_mem*eps_0)", "", "")
			.add_method("set_membrane_thickness", &T::set_membrane_thickness, "", "membrane thickness", "", "")
			.add_method("set_dendritic_radius", &T::set_dendritic_radius, "", "radius", "", "")
			.add_method("set_rtf", &T::set_rtf, "", "universal gas constant#temperature#Faraday constant", "Set natural constants in dimensions of choice.", "")
			.add_method("set_represented_dimension", &T::set_represented_dimension, "", "represented dimension (either 2 or 3)", "")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "PNP1D_FV", tag);
	}
#endif

	// PNPUpwind
	{
		typedef PNPUpwind<TDomain::dim> T;
		typedef IConvectionShapes<TDomain::dim> TBase;
		string nameBase = string("PNPUpwind");
		string name = nameBase;	name.append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.add_method("set_alpha", &T::set_alpha, "", "exponential weighing factor", "default is 1e2", "")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, nameBase, tag);
	}

#ifdef UG_PARALLEL
	// PNPDistroAdjuster
	{
		typedef PNPDistroManager<TDomain> T;
		string nameBase = string("PNPDistroManager");
		string name = nameBase;	name.append(suffix);
#ifdef NPParmetis
		typedef parmetis::AnisotropyUnificator<TDomain, typename grid_dim_traits<TDomain::dim>::grid_base_object> TBase;
		reg.add_class_<T, TBase>(name, grp)
#else
		reg.add_class_<T>(name, grp)
#endif
			.template add_constructor<void (*)(SmartPtr<ApproximationSpace<TDomain> >)>("approx space")
			.add_method("add_interface", &T::add_interface, "", "interface", "", "")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, nameBase, tag);
	}

	reg.add_function("set_distro_adjuster", &set_distro_adjuster<TDomain>, grp.c_str(), "", "", "");
#endif


	// extension refinement mark adjuster
	{
		typedef ExtensionRefMarkAdjuster<TDomain> T;
		string nameBase = string("ExtensionRefMarkAdjuster");
		string name = nameBase;	name.append(suffix);
		reg.add_class_<T>(name, grp)
		    .template add_constructor<void (*)(SmartPtr<TDomain>, std::vector<number>, const std::string)>
				("subset handler#direction (as vector)#useless subset")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, nameBase, tag);
	}

	reg.add_function("add_extension_ref_mark_adjuster", &add_extension_ref_mark_adjuster<TDomain>, grp.c_str(), "", "", "");
}

/**
 * Function called for the registration of Dimension dependent parts.
 * All Functions and Classes depending on the Dimension
 * are to be placed here when registering. The method is called for all
 * available Dimension types, based on the current build options.
 *
 * @param reg		registry
 * @param grp		group for sorting of functionality
 */
template <int dim>
static void Dimension(Registry& reg, string grp)
{
	string suffix = GetDimensionSuffix<dim>();
	string tag = GetDimensionTag<dim>();
}

/**
 * Function called for the registration of Algebra dependent parts.
 * All Functions and Classes depending on Algebra
 * are to be placed here when registering. The method is called for all
 * available Algebra types, based on the current build options.
 *
 * @param reg		registry
 * @param grp		group for sorting of functionality
 */
template <typename TAlgebra>
static void Algebra(Registry& reg, string grp)
{
	string suffix = GetAlgebraSuffix<TAlgebra>();
	string tag = GetAlgebraTag<TAlgebra>();
}

/**
 * Function called for the registration of Domain and Algebra independent parts.
 * All Functions and Classes not depending on Domain and Algebra
 * are to be placed here when registering.
 *
 * @param reg		registry
 * @param grp		group for sorting of functionality
 */
static void Common(Registry& reg, string grp)
{
	// EDL simulation
	{
		typedef EDLSimulation T;
		reg.add_class_<T>("EDLSimulation", grp)
			.add_constructor()
			.add_constructor<void (*)(number)>("membrane surface charge")
			.add_method("set_constants", &T::set_constants)
			.add_method("set_bnd_cond", &T::set_bnd_cond)
			.add_method("set_geom_specs", &T::set_geom_specs)
			.add_method("compute_solution", &T::compute_solution)
			.add_method("set_verbosity_level", &T::set_verbosity_level)
			.add_method("set_output_file", &T::set_output_file)
			.add_method("set_min_reduction", &T::set_min_reduction)
			.add_method("set_min_defect", &T::set_min_defect)
			.add_method("set_max_iter", &T::set_max_iter)
			.add_method("calc_EDL_ions", &T::calc_EDL_ions)
			.set_construct_as_smart_pointer(true);
	}

	// interface base class
	{
		typedef IInterface1D T;
		string name = string("IInterface1D");
		reg.add_class_<T>(name, grp);
	}

	// interface refinement mark adjuster
	{
		typedef InterfaceRefMarkAdjuster T;
		reg.add_class_<T>("InterfaceRefMarkAdjuster", grp)
		    .add_constructor()
		    .add_method("add_interfaces", &T::add_interfaces, "", "vector of interfaces", "")
			.add_method("set_subset_handler", &T::set_subset_handler, "", "subset handler", "")
			.set_construct_as_smart_pointer(true);
	}

	reg.add_function("add_interface_ref_mark_adjuster", &add_interface_ref_mark_adjuster, grp.c_str(), "", "", "");

	// morphology generator
	{
		typedef MorphoGen T;
		string name = string("MorphoGen");
		reg.add_class_<T>(name, grp)
			.add_constructor()
			.add_method("set_num_neck_filaments", &T::set_num_neck_filaments, "", "", "")
			.add_method("set_num_filaments", &T::set_num_filaments, "", "", "")
			.add_method("set_fil_anisotropic", &T::set_fil_anisotropic, "", "", "")
			.add_method("set_with_refinement_projector", &T::set_with_refinement_projector, "", "", "")
			.add_method("set_seed", &T::set_seed, "", "", "")
			.add_method("set_randomized", &T::set_randomized, "", "", "")
			.add_method("set_membrane_envelope_radius", &T::set_membrane_envelope_radius, "", "", "")
			.add_method("set_filament_envelope_radius", &T::set_filament_envelope_radius, "", "", "")
			.add_method("set_resolution", &T::set_resolution, "", "", "")
			.add_method("create_dendrite", &T::create_dendrite, "", "", "")
			.add_method("create_dendrite_2d", &T::create_dendrite_2d, "", "", "")
			.add_method("create_dendrite_1d", &T::create_dendrite_1d, "", "", "")
			.set_construct_as_smart_pointer(true);
	}

}


}; // end Functionality

// end group plugin_nernst_planck
/// \}

} // end namespace nernst_planck


/**
 * This function is called when the plugin is loaded.
 */
extern "C" void
InitUGPlugin_nernst_planck(Registry* reg, string grp)
{
	grp.append("/nernst_planck");
	typedef nernst_planck::Functionality Functionality;

	// we only register with algebra types CPU1 and CPU5,
	// as we explicitly instantiate templates only with these in cpp files
	typedef boost::mpl::list
	<
		#ifdef UG_CPU_1
		CPUAlgebra,
		#endif
		#ifdef UG_CPU_5
		CPUBlockAlgebra<5>,
		#endif
		#ifdef UG_CPU_6
		CPUBlockAlgebra<6>,
		#endif
		end_boost_list
	> MyCompileAlgebraList;

	try
	{
		RegisterCommon<Functionality>(*reg,grp);
		RegisterDimensionDependent<Functionality>(*reg,grp);
		RegisterDomainDependent<Functionality>(*reg,grp);
		//RegisterAlgebraDependent<Functionality, MyCompileAlgebraList>(*reg,grp);
		RegisterDomainAlgebraDependent<Functionality, CompileDomainList, MyCompileAlgebraList>(*reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

} // namespace ug



