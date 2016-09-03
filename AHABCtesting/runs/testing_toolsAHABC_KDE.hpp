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
\brief Tools for doing various AHABC routines using KDE

*/

#ifndef _INC_TESTINGTOOLS_AHABC_KDE_H
#define _INC_TESTINGTOOLS_AHABC_KDE_H


#include "ahabc_values.hpp"
#include "parameter_prior_set.hpp"
#include "summary_statistic_set.hpp"

#include <mscplusplus/population_structure.hpp>

#include "piecewise_constant_function.hpp"

#include <boost/shared_ptr.hpp>

#include <vector>
#include <utility>
#include <string>



void doAHABC_KDEHeterozygosity(
				const std::string& prefix,
				const std::string& spec,
				size_t n_obs,
				size_t n_reps,
				size_t nloci,
				size_t nsites,
				const std::vector< std::pair<double, double> >& priorLims,
				boost::shared_ptr< hudson_ms::PopulationStructure > pop,
				size_t n_reps_per_particle,
				const std::vector < double >& trueParams,
				int seed,
				cxsc::real maxMeasure,
				bool reducePriorLims);

void doAHABC_KDEHeterozygosity(
				const std::string& prefix,
				const std::string& spec,
				size_t n_obs,
				size_t n_reps,
				size_t nloci,
				size_t nsites,
				const std::vector< std::pair<double, double> >& priorLims,
				boost::shared_ptr< hudson_ms::PopulationStructure > pop,
				size_t n_reps_per_particle,
				const std::vector < double >& trueParams,
				const std::vector < std::vector < double > >& otherParams,
				int seed,
				cxsc::real maxMeasure,
				bool reducePriorLims);
				
void doAHABC_KDEHeterozygosityAndSegStats(
				const std::string& prefix,
				const std::string& spec,
				size_t n_obs,
				size_t n_reps,
				size_t nloci,
				size_t nsites,
				const std::vector< std::pair<double, double> >& priorLims,
				boost::shared_ptr< hudson_ms::PopulationStructure > pop,
				size_t n_reps_per_particle,
				const std::vector < double >& trueParams,
				int seed,
				cxsc::real maxMeasure,
				bool reducePriorLims);

void doAHABC_KDEHeterozygosityAndSegStats(
				const std::string& prefix,
				const std::string& spec,
				size_t n_obs,
				size_t n_reps,
				size_t nloci,
				size_t nsites,
				const std::vector< std::pair<double, double> >& priorLims,
				boost::shared_ptr< hudson_ms::PopulationStructure > pop,
				size_t n_reps_per_particle,
				const std::vector < double >& trueParams,
				const std::vector < std::vector < double > >& otherParams,
				int seed,
				cxsc::real maxMeasure,
				bool reducePriorLims);


void doAHABC_KDE(
				const std::string& prefix,
				const std::string& spec,
				size_t n_obs,
				size_t n_reps,
				size_t nloci,
				size_t nsites,
				const std::vector< std::pair<double, double> >& priorLims,
				boost::shared_ptr< hudson_ms::PopulationStructure > pop,
				size_t n_reps_per_particle,
				const std::vector < double >& trueParams,
				const std::vector < std::vector < double > >& otherParams,
				int seed,
				cxsc::real maxMeasure,
				const std::vector < size_t >& reqStats,
				bool reducePriorLims);
	
#endif
