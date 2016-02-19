/*
 * nernst_planck_plugin.cpp
 *
 *  Created on: 28.05.2014
 *      Author: mbreit
 */

#include "bridge/util.h"
#include "bridge/util_domain_algebra_dependent.h"

#include "nernst_planck_util.h"
#include "copy_neighbor_value_constraint.h"
#include "electric_circuit.h"
#include "interface1d_fv.h"
#include "pnp1d_fv1.h"
#include "pnp1d_fv.h"
#include "vtk_export_ho.h"


using namespace std;
using namespace ug::bridge;

namespace ug{
namespace nernst_planck{

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

	typedef typename TAlgebra::vector_type vector_type;
	typedef GridFunction<TDomain, TAlgebra> TGridFunction;

	// write residuals to file function for param optimization
	reg.add_function("write_residuals_to_file", &writeResidualsToFile<TGridFunction>, grp.c_str(),
						 "", "solution#reference solution#outFileName",
						 "writes residuals in every dof to file as required by Ivo's optimization routines "
						 "and returns squared 2-norm of residual vector");

	// export solution
	reg.add_function("export_solution", &exportSolution<GridFunction<TDomain, TAlgebra> >, grp.c_str(),
					 "", "solution#time#subsetNames#functionNames#outFileName",
					 "outputs solutions to file");

	// import solution
	reg.add_function("import_solution", &importSolution<GridFunction<TDomain, TAlgebra> >, grp.c_str(),
					 "", "solution#subset names#function name#input file name",
					 "writes values for the given function and on the given subsets "
					 "from the given file to the given solution vector "
					 "(using the value of the nearest neighbor for each vertex)");

	// interface (in the sense of programming) for the nD/1D interface (in the sense of manifold) class
	{
		typedef IInterface1D<TDomain, TAlgebra> T;
		typedef IDomainConstraint<TDomain, TAlgebra> TBase;
		string name = string("IInterface1D").append(suffix);
		reg.add_class_<T, TBase >(name, grp)
			.add_method("update", &T::update);
			//.add_method("check_values_at_interface", &T::check_values_at_interface, "", "solution grid function", "", "");
		reg.add_class_to_group(name, "IInterface1D", tag);
	}

	// additive nD/1D interface
	{
		typedef AdditiveInterface1D<TDomain, TAlgebra> T;
		typedef IInterface1D<TDomain, TAlgebra> TBase;
		string name = string("AdditiveInterface1D").append(suffix);
		reg.add_class_<T, TBase >(name, grp)
			.template add_constructor<void (*)(const char*, const char*, const char*, const char*)>
				("function(s)#constrained subset#high-dim interface node subset#"
				 "one-dim interface node subset")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "AdditiveInterface1D", tag);
	}

	// multiplicative nD/1D interface
	{
		typedef MultiplicativeInterface1D<TDomain, TAlgebra> T;
		typedef IInterface1D<TDomain, TAlgebra> TBase;
		string name = string("MultiplicativeInterface1D").append(suffix);
		reg.add_class_<T, TBase >(name, grp)
			.template add_constructor<void (*)(const char*, const char*, const char*, const char*)>
				("function(s)#constrained subset#high-dim interface node subset#"
				 "one-dim interface node subset")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "MultiplicativeInterface1D", tag);
	}

	// InterfaceMapper
	{
		typedef typename IInterface1D<TDomain, TAlgebra>::Interface1DMapper T;
		typedef ILocalToGlobalMapper<TAlgebra> TBase;
		string name = string("Interface1DMapper").append(suffix);
		reg.add_class_<T, TBase >(name, grp)
			.template add_constructor<void (*)(SmartPtr<IAssemble<TAlgebra> >)>
				("SmartPtr to domainDisc")
			.add_method("add_interface", &T::add_interface, "", "SmartPtr to class of IInterface1D", "", "")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "Interface1DMapper", tag);

	}

	// constraint for useless DoFs
	{
		typedef CopyNeighborValueConstraint<TDomain, TAlgebra> T;
		typedef IDomainConstraint<TDomain, TAlgebra> TBase;
		string name = string("CopyNeighborValueConstraint").append(suffix);
		reg.add_class_<T, TBase >(name, grp)
			.template add_constructor<void (*)(const char*, const char*)>
				("function(s)#constrained subset")
				.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "CopyNeighborValueConstraint", tag);
	}

	// vtk export for higher order grid functions
	{
		reg.add_function("vtk_export_ho",
			static_cast<void (*) (SmartPtr<TGridFunction>, const std::vector<std::string>&,
						size_t, SmartPtr<VTKOutput<TGridFunction::domain_type::dim> >, const char*)>
				(&vtk_export_ho<GridFunction<TDomain, TAlgebra> >),
			grp.c_str(), "new grid function", "input grid function#functions to be exported#order",
			"creates a grid function of order 1 containing interpolated values from high-order input grid function on a refined grid");
		reg.add_function("vtk_export_ho",
			static_cast<void (*) (SmartPtr<TGridFunction>, const std::vector<std::string>&, size_t,
						SmartPtr<VTKOutput<TGridFunction::domain_type::dim> >, const char*, size_t, number)>
				(&vtk_export_ho<GridFunction<TDomain, TAlgebra> >),
			grp.c_str(), "new grid function", "input grid function#functions to be exported#order",
			"creates a grid function of order 1 containing interpolated values from high-order input grid function on a refined grid");
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
					 "", "approximation space#inner subset name#full-dim interface node subset name"
					 "#1d interface node subset name", "adjusts location of full-dim interface node after"
					 "global refinement of the interface");

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
			.add_method("set_represented_dimension", &T::set_represented_dimension, "", "represented dimension (either 2 or 3)", "")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "PNP1D_FV", tag);
	}

	//reg.add_function("test_positions", &TestPositions<TDomain>, grp.c_str(), "", "", "");
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
	typedef ElectricCircuit T;
	reg.add_class_<T>("ElectricCircuit", grp)
		.add_constructor()
		.add_method("add_capacitor", &T::add_capacitor)
		.add_method("add_resistor", &T::add_resistor)
		.add_method("add_voltage_source", &T::add_voltage_source)
		.add_method("add_current_source", &T::add_current_source)
		.add_method("add_initial_solution", &T::add_initial_solution)
		.add_method("init", &T::init)
		.add_method("update_rhs", &T::update_rhs)
		.add_method("write_to_file", &T::write_to_file)
		.add_method("close_file", &T::close_file)
		.add_method("solve_stationary", &T::solve_stationary)
		.add_method("solve_euler", &T::solve_euler)
		.add_method("solve_trapezoid", &T::solve_trapezoid)
		.add_method("get_solution", &T::get_solution)
		.add_method("get_rhs", &T::get_rhs)
		.set_construct_as_smart_pointer(true);
}

}; // end Functionality

// end group plugin_nernst_planck
/// \}

} // end namespace calciumDynamics


/**
 * This function is called when the plugin is loaded.
 */
extern "C" void
InitUGPlugin_nernst_planck(Registry* reg, string grp)
{
	grp.append("/nernst_planck");
	typedef nernst_planck::Functionality Functionality;

	try{
		RegisterCommon<Functionality>(*reg,grp);
		//RegisterDimensionDependent<Functionality>(*reg,grp);
		RegisterDomainDependent<Functionality>(*reg,grp);
		//RegisterAlgebraDependent<Functionality>(*reg,grp);
		RegisterDomainAlgebraDependent<Functionality>(*reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

}// namespace ug



