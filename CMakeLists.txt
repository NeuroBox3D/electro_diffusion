##
 # Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 #
 # Author: Markus Breit
 # Creation date: 2013
 #
 # This file is part of NeuroBox, which is based on UG4.
 #
 # NeuroBox and UG4 are free software: You can redistribute it and/or modify it
 # under the terms of the GNU Lesser General Public License version 3
 # (as published by the Free Software Foundation) with the following additional
 # attribution requirements (according to LGPL/GPL v3 §7):
 #
 # (1) The following notice must be displayed in the appropriate legal notices
 # of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 #
 # (2) The following notice must be displayed at a prominent place in the
 # terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 #
 # (3) The following bibliography is recommended for citation and must be
 # preserved in all covered files:
 # "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 #   parallel geometric multigrid solver on hierarchically distributed grids.
 #   Computing and visualization in science 16, 4 (2013), 151-164"
 # "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 #   flexible software system for simulating PDE based models on high performance
 #   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 # "Stepniewski, M., Breit, M., Hoffer, M. and Queisser, G.
 #   NeuroBox: computational mathematics in multiscale neuroscience.
 #   Computing and visualization in science (2019).
 # "Breit, M. et al. Anatomically detailed and large-scale simulations studying
 #   synapse loss and synchrony using NeuroBox. Front. Neuroanat. 10 (2016), 8"
 #
 # This program is distributed in the hope that it will be useful,
 # but WITHOUT ANY WARRANTY; without even the implied warranty of
 # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 # GNU Lesser General Public License for more details.
 ##


################################################################################
# NERNST PLANCK
#
# The electro-diffusion plugin.
################################################################################

cmake_minimum_required(VERSION 2.6)

project(NERNST_PLANCK)
include(CheckCXXSymbolExists)

# include the definitions and dependencies for ug-plugins.
include(${UG_ROOT_CMAKE_PATH}/ug_plugin_includes.cmake)



message(STATUS "")
message(STATUS "Additional options (-DNP) for the nernst_planck plugin")


# compile with tetgen 1.5 support
option(NPTetgen "Tetgen")
message(STATUS "      Tetgen:   ${NPTetgen} (options are: ON, OFF)")


# sources
set(SOURCES extension/domain1d_solution_adjuster.cpp
            extension/extension_refMarkAdjuster.cpp
            interface/interface1d_fv.cpp
            interface/intf_refMarkAdjuster.cpp
            util/edl_1d.cpp
            util/flux_exporter.cpp
            util/neck_recorder.cpp
			util/nernst_planck_util.cpp
            util/order.cpp
            util/refinement_error_estimator.cpp
            morpho_gen.cpp
            nernst_planck_plugin.cpp
            pnp_smoother.cpp
            pnp_upwind.cpp
   )


if (PARALLEL)
	set(SOURCES ${SOURCES} interface/intf_distro_adjuster.cpp)
endif (PARALLEL)

if (Parmetis)
	set (NPParmetis ON)
endif (Parmetis)


if (MembranePotentialMapping)
	set(NPWithMPM 1)
endif (MembranePotentialMapping)


# error if tetgen libs/includes are not found
if (NPTetgen)
	UGInstall("tetgen" TETGEN15_PATH)

	unset(TETGEN15_LIBS CACHE)
	find_library(TETGEN15_LIBS NAMES tet PATHS ${TETGEN15_PATH})
	if (NOT TETGEN15_LIBS)
		message(FATAL_ERROR "ERROR: TetGen-1.5 library not found in ${TETGEN15_PATH}.")
	endif (NOT TETGEN15_LIBS)
	
	find_path(TETGEN15_INCLUDE_PATH tetgen.h PATHS ${TETGEN15_PATH})
	if (NOT TETGEN15_INCLUDE_PATH)
		message(FATAL_ERROR "ERROR: TetGen-1.5 include file tetgen.h not found in ${TETGEN15_PATH}.")
	endif (NOT TETGEN15_INCLUDE_PATH)
		
	check_cxx_symbol_exists(fsymself ${TETGEN15_INCLUDE_PATH}/tetgen.h FSMYSELF_EXISTS)
	if (NOT FSMYSELF_EXISTS)
		message(FATAL_ERROR "ERROR: Though TetGen was found, it does not seem to be the required version 1.5.")
		unset(FSMYSELF_EXISTS CACHE)
	endif()
endif (NPTetgen)


# This will copy the configurable file config_cmake.h to
# the final configuration header np_config.h,
# setting the #defines as given by the current CMake variables.
configure_file(config_cmake.h ${CMAKE_BINARY_DIR}/np_config.h)


if(buildEmbeddedPlugins)
	# add the sources to ug4's sources
	ExportSources(${CMAKE_CURRENT_SOURCE_DIR} ${SOURCES})
	
	if (NPTetgen)		
		# add include path for tetgen.h
		ExportIncludes(${TETGEN15_INCLUDE_PATH})
		
		# add library to link against
		ExportDependencies(${TETGEN15_LIBS})
	endif (NPTetgen)
	
else(buildEmbeddedPlugins)
	# create a shared library from the sources and link it against ug4 and external libraries
	add_library(nernst_planck SHARED ${SOURCES})

	if (NPTetgen)
		# add include path for tetgen.h
		include_directories(${TETGEN15_INCLUDE_PATH})
	
		# add tetgen library to link against
		set(additionalLibs ${additionalLibs} ${TETGEN15_LIBS})
	endif (NPTetgen)
	
	# make plugin link against MPM
	if (MembranePotentialMapping)
		set(additionalLibs MembranePotentialMapping ${additionalLibs})
	endif (MembranePotentialMapping)
	
	# add Parmetis to linker (if built alongside)
	if (Parmetis)
		set (additionalLibs ${additionalLibs} Parmetis)
	endif (Parmetis)
	
	target_link_libraries(nernst_planck ug4 ${additionalLibs})
endif(buildEmbeddedPlugins)
