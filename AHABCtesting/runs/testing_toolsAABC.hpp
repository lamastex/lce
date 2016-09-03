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
\brief Tools for doing various AABC routines

*/

#ifndef _INC_TESTINGTOOLS_AABC_H
#define _INC_TESTINGTOOLS_AABC_H




#include "config.h" // mct

#include "ahabc_values.hpp"

#include "summary_statistic_set.hpp"

#include "pop_mc_abc_simulator_adaptive.hpp"
#include "parameter_particle_validator.hpp"
#include "multi_loci_polytable_set.hpp"
#include "multi_loci_polytable.hpp"

#include "prng.hpp"


#include <mscplusplus/population_structure.hpp>

#include <boost/shared_ptr.hpp>

#include <boost/function.hpp>

#include <string>
#include <vector>
#include <stack>
#include <utility>

boost::shared_ptr < aabc::ParameterParticleValidator > 
			makeValidatorHeterozygosity(boost::shared_ptr < mct::PRNGen > rr_sampler,
			const boost::shared_ptr < const mct::MultiLociPolyTable > refPtr,
			size_t n_reps_per_particle);

boost::shared_ptr < aabc::ParameterParticleValidator > 
			makeValidatorHeterozygosityAndSegSites(
			boost::shared_ptr < mct::PRNGen > rr_sampler,
			const boost::shared_ptr < const mct::MultiLociPolyTable > refPtr,
			size_t n_reps_per_particle);

	
boost::shared_ptr < aabc::PopMCabcSimulatorAdaptive > doSimulationMultLociPolytable(
				boost::shared_ptr < mct::PRNGen > r_loops,
				size_t n_reps_per_particle,	
				const std::vector< std::pair<double, double> >& priorLims,
				boost::shared_ptr < const mct::MultiLociPolyTable > refPtr,
				size_t iterations,
				double baseEpsilonPercentile,
				std::stack < double > epsilonMults,
				boost::function < 
					boost::shared_ptr < aabc::ParameterParticleValidator > 
					(boost::shared_ptr < mct::PRNGen > rr_sampler,
					const boost::shared_ptr < const mct::MultiLociPolyTable > refPtr,
					size_t n_reps_per_particle) 
								> validatorMaker,
				const std::string& logfilename); 

				
	
void outputAABCLogFileString(const std::string& s, 
								const std::string& line, 
								bool append = true);


#endif
