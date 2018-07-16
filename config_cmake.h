/*
 * tetgen_config_cmake.h
 *
 *  Created on: 22.09.2016
 *      Author: mbreit */

#ifndef UG__PLUGINS__NERNST_PLANCK__TETGEN_CONFIG_H
#define UG__PLUGINS__NERNST_PLANCK__TETGEN_CONFIG_H

#cmakedefine NPTetgen
#ifdef NPTetgen
	#ifndef TETGEN_15_ENABLED
		#define TETGEN_15_ENABLED
	#endif
	#ifndef TETLIBRARY
		#define TETLIBRARY
	#endif
#endif

#cmakedefine NPParmetis



#endif // UG__PLUGINS__NERNST_PLANCK__TETGEN_CONFIG_H
