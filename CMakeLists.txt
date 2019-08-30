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
set(SOURCES nernst_planck_plugin.cpp
			nernst_planck_util.cpp
            #pnp1d_fv1.cpp
            #pnp1d_fv.cpp
            electric_circuit.cpp
            edl_1d.cpp
            interface1d_fv.cpp
            order.cpp
            intf_refMarkAdjuster.cpp
            extension_refMarkAdjuster.cpp
            domain1d_solution_adjuster.cpp
            flux_exporter.cpp
            pnp_upwind.cpp
            morpho_gen.cpp
            pnp_smoother.cpp
            neck_recorder.cpp
            refinement_error_estimator.cpp
   )
   
if (APPLE)
	set(SOURCES ${SOURCES} mem_info.cpp)
endif (APPLE)

if (PARALLEL)
	set(SOURCES ${SOURCES} intf_distro_adjuster.cpp redistribution_util.cpp)
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
		# add definition for tetgen usage to morpho_gen.cpp
		# Unfortunately, this still results in the whole target (i.e., all of ug)
		# to be recompiled. This is a known issue when creating MakeFiles with CMake.
		#ExportSingleFileCompileFlags(${CMAKE_CURRENT_SOURCE_DIR}/morpho_gen.cpp "TETGEN_15_ENABLED;TETLIBRARY")
		
		# add include path for tetgen.h
		ExportIncludes(${TETGEN15_INCLUDE_PATH})
		
		# add library to link against
		ExportDependencies(${TETGEN15_LIBS})
	endif (NPTetgen)
	
else(buildEmbeddedPlugins)
	# create a shared library from the sources and link it against ug4 and external libraries
	add_library(nernst_planck SHARED ${SOURCES})

	if (NPTetgen)
		# add definition for tetgen usage to morpho_gen.cpp
		#set_property(SOURCE morpho_gen.cpp APPEND PROPERTY COMPILE_DEFINITIONS "TETGEN_15_ENABLED;TETLIBRARY")
			
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
