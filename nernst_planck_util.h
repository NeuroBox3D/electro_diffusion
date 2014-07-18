/*
 * nernst_planck_util.h
 *
 *  Created on: 17.07.2014
 *      Author: mbreit
 */

#ifndef _NERNST_PLANCK__NERNST_PLANCK_UTIL_H_
#define _NERNST_PLANCK__NERNST_PLANCK_UTIL_H_


#include "common/common.h"
#include "lib_disc/function_spaces/approximation_space.h"
#include "lib_disc/function_spaces/grid_function.h"

#include <iostream>
#include <fstream>


namespace ug{
namespace nernst_planck{


template <typename TGridFunction>
number writeResidualsToFile(SmartPtr<TGridFunction> sol1, SmartPtr<TGridFunction> sol2, const char* fileName);



} // namspace calciumDynamics
} // namespace ug

#include "nernst_planck_util_impl.h"

#endif // _NERNST_PLANCK__NERNST_PLANCK_UTIL_H_
