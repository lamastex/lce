/* MCT - Markov Chains on Trees.


   Copyright (C) 2012 Jennifer Harlow

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
\brief Tools for AABC runs
*/

#include "testing_toolsAABC.hpp"

#include "testing_tools.hpp"

//#include "test_routines.hpp"


//#include "piecewise_constant_function.hpp"
//#include "spsnode_measure_obj.hpp"

#include "config.h" // mct

//#include "ahabc_values.hpp"

//#include "descriptive_stats.hpp"
//#include "summary_statistic_set.hpp"
//#include "summary_statistic.hpp"

#include "parameter_particle_validator_multi_loci_polytable_heterozygosity.hpp"
#include "parameter_particle_validator_multi_loci_polytable.hpp"
#include "multi_loci_simdata_sampler.hpp"
#include "multi_loci_polytable_sampler.hpp"

#include "parameter_particle_set.hpp"

#include "parameter_prior_uniform.hpp"
#include "parameter_prior_set.hpp"
#include "parameter_prior.hpp"

#include "prng_gsl.hpp"

#include "utilities.hpp"


#include <iostream>
#include <sstream>

using namespace std;
using namespace mct;
using namespace aabc;
using namespace subpavings;

// really, these are 'internal' headers, which do not need to be publically declared



boost::shared_ptr < PopMCabcSimulatorAdaptive > doSimulationMultLociPolytable(
				boost::shared_ptr < PRNGen > r_loops,
				size_t n_reps_per_particle,	
				const std::vector< std::pair<double, double> >& priorLims,
				boost::shared_ptr < const mct::MultiLociPolyTable > refPtr,
				size_t iterations,
				double baseEpsilonPercentile,
				std::stack < double > epsilonMults, // pass by value
				boost::function < 
					boost::shared_ptr < aabc::ParameterParticleValidator > 
					(boost::shared_ptr < mct::PRNGen > rr_sampler,
					const boost::shared_ptr < const mct::MultiLociPolyTable > refPtr,
					size_t n_reps_per_particle) 
								> validatorMaker,
				const std::string& logfilename) 
{
	boost::shared_ptr < PRNGen > rr_sampler = r_loops->spawnAnother();
	boost::shared_ptr < PRNGen > rr_simulator = r_loops->spawnAnother();

	size_t n_params = priorLims.size(); // number of parameters we are using
	
	ParameterPriorSet priors(n_params); // prior set with n_params element
	/* set that element to a uniform prior on -10, 10
	 * and give the prior its random number generator */
	for (size_t j = 0; j < n_params; ++j) {
		boost::shared_ptr < PRNGen > rr_prior = r_loops->spawnAnother();
		priors.at(j) = boost::shared_ptr < ParameterPrior > 
					( new ParameterPriorUniform(priorLims[j].first,
												priorLims[j].second,
												rr_prior));
	}
	
	boost::shared_ptr < ParameterParticleValidator > validator = 
			validatorMaker(rr_sampler,
			refPtr,
			n_reps_per_particle);

	// Finally create the simulator itself
	boost::shared_ptr < PopMCabcSimulatorAdaptive > simulatorPtr( new 
			PopMCabcSimulatorAdaptive(priors, 
							validator,
							baseEpsilonPercentile,
							rr_simulator));
	
	//iterate
	while (!epsilonMults.empty()) {
		
		simulatorPtr->iterate(epsilonMults.top(), iterations, logfilename);

		// pop the stack			
		epsilonMults.pop();
		
	}
	
	return simulatorPtr;
}

boost::shared_ptr < ParameterParticleValidator > 
			makeValidatorHeterozygosity(boost::shared_ptr < PRNGen > rr_sampler,
			const boost::shared_ptr < const mct::MultiLociPolyTable > refPtr,
			size_t n_reps_per_particle)
{
	// make a simdata (infinite sites model) sampler with its prng
	boost::shared_ptr < MultiLociPolyTableSampler > 
			sampler_ptr ( new MultiLociSimDataSampler(rr_sampler) );

	/* Setting up the validator.  This is the thing that says how
	 * close a given particle is to what wanted
	 * 
	 * This one only looks at heterozygosity*/
	boost::shared_ptr < ParameterParticleValidator > 
		validator( new ParameterParticleValidatorMultiLociPolyTableHeterozygosity(
				refPtr, sampler_ptr, n_reps_per_particle) );
	
	return validator;
}

boost::shared_ptr < ParameterParticleValidator > 
			makeValidatorHeterozygosityAndSegSites(
			boost::shared_ptr < PRNGen > rr_sampler,
			const boost::shared_ptr < const mct::MultiLociPolyTable > refPtr,
			size_t n_reps_per_particle)
{
	// make a simdata (infinite sites model) sampler with its prng
	boost::shared_ptr < MultiLociPolyTableSampler > 
			sampler_ptr ( new MultiLociSimDataSampler(rr_sampler) );

	/* Setting up the validator.  This is the thing that says how
	 * close a given particle is to what wanted
	 * 
	 * This one only looks at heterozygosity and seg sites*/
	boost::shared_ptr < ParameterParticleValidator > 
		validator( new ParameterParticleValidatorMultiLociPolyTable(
				refPtr, sampler_ptr, n_reps_per_particle) );
	
	return validator;
}


// Method to add a line to a file
// Output goes to file named according to argument s
// does not add any new lines
void outputAABCLogFileString(const std::string& s, 
								const std::string& line, 
								bool append)
{
	std::ofstream os;
	if (append) os.open(s.c_str(), ios::app);         // append
	else os.open(s.c_str()); // don't append
	if (os.is_open()) {
		os << "\n" << line << flush;
		os.close();
	}
	else {
		std::cerr << "Error: could not open file named "
			<< s << std::endl << std::endl;
	}
}
