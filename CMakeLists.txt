################################################################################
# NERNST PLANCK
#
# The electro-diffusion plugin.
################################################################################

cmake_minimum_required(VERSION 2.6)

project(NERNST_PLANCK)


# include the definitions and dependencies for ug-plugins.
include("../../../cmake/ug_plugin_includes.cmake")


# sources
set(SOURCES nernst_planck_plugin.cpp
            interface_1d_fv1.cpp)

message(STATUS "-------------------information regarding calciumDynamics-------------------")


# store libs to be linked against here
set(linkLibraries)

if(buildEmbeddedPlugins)
	# add the sources to ug4's sources
	EXPORTSOURCES(plugins/experimental/nernst_planck ${SOURCES})
	#exportDependencies(${ANN_LIBRARY})
else(buildEmbeddedPlugins)
	# create a shared library from the sources and link it against ug4 and external libraries
	add_library(nernst_planck SHARED ${SOURCES})
	target_link_libraries(nernst_planck ug4 ${linkLibraries})
endif(buildEmbeddedPlugins)

message(STATUS "---------------------------------------------------------------------------") 
