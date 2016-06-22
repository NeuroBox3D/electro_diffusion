################################################################################
# NERNST PLANCK
#
# The electro-diffusion plugin.
################################################################################

cmake_minimum_required(VERSION 2.6)

project(NERNST_PLANCK)


# include the definitions and dependencies for ug-plugins.
include(${UG_ROOT_CMAKE_PATH}/ug_plugin_includes.cmake)


# sources
set(SOURCES nernst_planck_plugin.cpp
            pnp1d_fv1.cpp
            pnp1d_fv.cpp
            electric_circuit.cpp
            edl_1d.cpp
            order.cpp
            charge_marking.cpp
            intf_refMarkAdjuster.cpp
   )

if (PARALLEL)
	set(SOURCES ${SOURCES} intf_distro_adjuster.cpp) 
endif(PARALLEL)


if(buildEmbeddedPlugins)
	# add the sources to ug4's sources
	EXPORTSOURCES(${CMAKE_CURRENT_SOURCE_DIR} ${SOURCES})
else(buildEmbeddedPlugins)
	# create a shared library from the sources and link it against ug4 and external libraries
	add_library(nernst_planck SHARED ${SOURCES})
	target_link_libraries(nernst_planck ug4)
endif(buildEmbeddedPlugins)
