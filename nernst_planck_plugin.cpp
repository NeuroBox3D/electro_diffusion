/*
 * nernst_planck_plugin.cpp
 *
 *  Created on: 28.05.2014
 *      Author: mbreit
 */

#include "bridge/util.h"
#include "bridge/util_domain_algebra_dependent.h"

#include "nernst_planck_util.h"
#include "PNP_1D.h"
#include "interface1d_fv1.h"
//#include "constrained_ilu.h"
#include "electric_circuit.h"

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

	// write residuals to file function for param optimization
	reg.add_function("write_residuals_to_file", &writeResidualsToFile<GridFunction<TDomain, TAlgebra> >, grp.c_str(),
						 "", "solution#reference solution#outFileName",
						 "writes residuals in every dof to file as required by Ivo's optimization routines "
						 "and returns squared 2-norm of residual vector");

	// interface (in the sense of programming) for the nD/1D interface (in the sense of manifold) class
	{
		typedef IInterface1DFV1<TDomain, TAlgebra> T;
		typedef IDomainConstraint<TDomain, TAlgebra> TBase;
		string name = string("IInterface1DFV1").append(suffix);
		reg.add_class_<T, TBase >(name, grp);
			//.add_method("check_values_at_interface", &T::check_values_at_interface, "", "solution grid function", "", "");
		reg.add_class_to_group(name, "IInterface1DFV1", tag);
	}

	// additive nD/1D interface
	{
		typedef AdditiveInterface1DFV1<TDomain, TAlgebra> T;
		typedef IInterface1DFV1<TDomain, TAlgebra> TBase;
		string name = string("AdditiveInterface1DFV1").append(suffix);
		reg.add_class_<T, TBase >(name, grp)
			.template add_constructor<void (*)(const char*, const char*, const char*, const char*, const char*)>
					  ("function(s)#high-dim constrained subset#one-dim constrained subset#"
					   "high-dim interface node subset#one-dim extension subset")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "AdditiveInterface1DFV1", tag);
	}

	// multiplicative nD/1D interface
	{
		typedef MultiplicativeInterface1DFV1<TDomain, TAlgebra> T;
		typedef IInterface1DFV1<TDomain, TAlgebra> TBase;
		string name = string("MultiplicativeInterface1DFV1").append(suffix);
		reg.add_class_<T, TBase >(name, grp)
			.template add_constructor<void (*)(const char*, const char*, const char*, const char*, const char*)>
					  ("function(s)#high-dim constrained subset#one-dim constrained subset#"
					   "high-dim interface node subset#one-dim extension subset")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "MultiplicativeInterface1DFV1", tag);
	}

	/*
	//	ILUC preconditioner (no longer needed)
	{
		typedef ILUC<TDomain, TAlgebra> T;
		typedef IPreconditioner<TAlgebra> TBase;
		string name = string("ILUC").append(suffix);
		reg.add_class_<T,TBase>(name, grp, "Incomplete LU Decomposition respecting constraints")
		.template add_constructor<void (*)(ConstSmartPtr<ApproximationSpace<TDomain> >)>("approximation space")
			.add_method("set_beta", &T::set_beta, "", "beta")
			.add_method("set_sort", &T::set_sort, "", "bSort", "if bSort=true, use a cuthill-mckey sorting to reduce fill-in. default false")
			.add_method("add_constraint", &T::add_constraint, "", "constraint")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ILUC", tag);
	}
	*/
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

	// 1D PNP
	{
		typedef PNP_1D<TDomain> T;
		typedef IElemDisc<TDomain> TBase;
		string name = string("PNP_1D").append(suffix);
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
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "PNP_1D", tag);
	}
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
	grp.append("/CalciumDynamics");
	typedef nernst_planck::Functionality Functionality;

	try{
		RegisterCommon<Functionality>(*reg,grp);
		RegisterDimensionDependent<Functionality>(*reg,grp);
		RegisterDomainDependent<Functionality>(*reg,grp);
		//RegisterAlgebraDependent<Functionality>(*reg,grp);
		RegisterDomainAlgebraDependent<Functionality>(*reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

}// namespace ug



