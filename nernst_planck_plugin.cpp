/*
 * nernst_planck_plugin.cpp
 *
 *  Created on: 28.05.2014
 *      Author: mbreit
 */

#include "bridge/util.h"
#include "bridge/util_domain_algebra_dependent.h"

#include "interface1d_fv1.h"

using namespace std;
using namespace ug::bridge;

namespace ug{
namespace nernst_planck{

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
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TDomain, typename TAlgebra>
static void DomainAlgebra(Registry& reg, string grp)
{
	string suffix = GetDomainAlgebraSuffix<TDomain,TAlgebra>();
	string tag = GetDomainAlgebraTag<TDomain,TAlgebra>();

	// interface (in the sense of programming) for the nD/1D interface (in the sense of manifold) class
	{
		typedef IInterface1DFV1<TDomain, TAlgebra> T;
		typedef IDomainConstraint<TDomain, TAlgebra> TBase;
		string name = string("IInterface1DFV1").append(suffix);
		reg.add_class_<T, TBase >(name, grp);
		reg.add_class_to_group(name, "IInterface1DFV1", tag);
	}

	// additive nD71D interface
	{
		typedef AdditiveInterface1DFV1<TDomain, TAlgebra> T;
		typedef IInterface1DFV1<TDomain, TAlgebra> TBase;
		string name = string("AdditiveInterface1DFV1").append(suffix);
		reg.add_class_<T, TBase >(name, grp)
			.template add_constructor<void (*)(const char*, const char*, const char*)>
					  ("function(s)#high-dim constrained subset#one-dim constrained subset")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "AdditiveInterface1DFV1", tag);
	}

	// additive nD71D interface
	{
		typedef MultiplicativeInterface1DFV1<TDomain, TAlgebra> T;
		typedef IInterface1DFV1<TDomain, TAlgebra> TBase;
		string name = string("MultiplicativeInterface1DFV1").append(suffix);
		reg.add_class_<T, TBase >(name, grp)
			.template add_constructor<void (*)(const char*, const char*, const char*)>
					  ("function(s)#high-dim constrained subset#one-dim constrained subset")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "MultiplicativeInterface1DFV1", tag);
	}
}

/**
 * Function called for the registration of Domain dependent parts.
 * All Functions and Classes depending on the Domain
 * are to be placed here when registering. The method is called for all
 * available Domain types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
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
}

/**
 * Function called for the registration of Dimension dependent parts.
 * All Functions and Classes depending on the Dimension
 * are to be placed here when registering. The method is called for all
 * available Dimension types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
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
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
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
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
static void Common(Registry& reg, string grp)
{

}

}; // end Functionality
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



