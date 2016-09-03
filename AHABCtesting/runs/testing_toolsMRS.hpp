/* MCT - Markov Chains on Trees.


   Copyright (C) 2011, 2012 Jennifer Harlow

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor,
   Boston, MA 02110-1301, USA.  */

/*! \file
\brief Tools for doing various routines using MRS types

*/

#ifndef _INC_TESTINGTOOLS_MRS_H
#define _INC_TESTINGTOOLS_MRS_H


#include "ahabc_values.hpp"
//#include "parameter_prior_set.hpp"
#include "summary_statistic_set.hpp"

#include <mscplusplus/population_structure.hpp>

#include "piecewise_constant_function.hpp"

#include <boost/shared_ptr.hpp>

#include <vector>
#include <string>



// make some values
boost::shared_ptr < const mct::AHABCValues  > makeValuesPtr(
		size_t n_reps,
		size_t nloci,
		size_t nsites,
		boost::shared_ptr< hudson_ms::PopulationStructure > pop,
		size_t n_reps_per_particle,	
		const subpavings::PiecewiseConstantFunction& pcfPosterior,
		const std::vector < double >& scale,
		const std::vector < double >& shift,
		const std::vector < size_t >& reqStats,
		int seed,
		const std::string& desStatsFilename);

	
#endif
