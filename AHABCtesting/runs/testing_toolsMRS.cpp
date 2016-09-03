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

#include "testing_toolsMRS.hpp"

#include "sptypes.hpp"

#include "testing_toolsSimdata.hpp"
#include "testing_tools.hpp"
#include "log_memory.hpp"

//#include "test_routines.hpp"

//#include "config.h" // mct

#include "parameter_prior_uniform.hpp"
#include "parameter_prior_set.hpp"
#include "parameter_prior.hpp"

#include "multi_loci_simdata_sampler.hpp"
#include "multi_loci_polytable_sampler.hpp"

#include "parameter_particle_set.hpp"
#include "descriptive_stats.hpp"
#include "prng.hpp"
#include "prng_gsl.hpp"

#include "utilities.hpp"




#include <real.hpp>


#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <cassert>
#include <iterator>

#include <boost/lambda/bind.hpp>
#include <boost/lambda/lambda.hpp>

#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_rng.h>


using namespace std;
using namespace mct;
using namespace aabc;
using namespace subpavings;


// make some values, specifying required stats
boost::shared_ptr < const AHABCValues  > makeValuesPtr(
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
		const std::string& logFilename)
{
	
	boost::shared_ptr < mct::PRNGen > rr ( new mct::PRNGenGSL(seed) );
	
	//vector of nloci copies of nsites
	std::vector < size_t > nloci_nsites(nloci, nsites);
		
	// prng for sampler
	boost::shared_ptr < PRNGen > rr_sampler = rr->spawnAnother();
	
	// make a simdata (infinite sites model) sampler with its prng
	boost::shared_ptr < MultiLociPolyTableSampler > 
					sampler_ptr ( new MultiLociSimDataSampler(rr_sampler) );
	
			
	cout << "about to do particles and sampling" << endl;
	
	/* memory problems if our samples are too large: 100,000 uses < 10%
	 * on mathxeon3 which is very reasonable, so limit to this and sample
	 * in sections.  Note that 6,000,000 reps with Simdata2d
	 * seems to use ~ 5.7g = 18% on 
	 * mathxeon even using this approach, ie when we have the histogram */
	
	size_t toDo = n_reps;
	size_t maxToDo = 100000;
	
	// a set for the parameter particles
	boost::shared_ptr< ParameterParticleSet > 
			pset_ptr ( new ParameterParticleSet(n_reps) );
	
	boost::shared_ptr < mct::SummaryStatisticSet > 
		sset_ptr ( new SummaryStatisticSet(n_reps) );
		
	size_t nParams = pcfPosterior.getDimensions();
	
	
	RVecData paramContainer;
	paramContainer.reserve(n_reps);
	
	/* simulate values in transformed world */
	pcfPosterior.simulateData(paramContainer,
					n_reps, seed);
	
	RVecData::iterator rit = paramContainer.begin();
	RVecData::iterator erit = paramContainer.begin();
	
			
	
	while (toDo > 0) {
		
		
		size_t doThisLoop = (toDo <= maxToDo ? toDo : maxToDo);
				
		cout << "\tsampling " << doThisLoop << " from posterior" << endl;
	
		advance(erit, doThisLoop);
		boost::shared_ptr< ParameterParticleSet > 
			pset_ptr_tmp(new ParameterParticleSet(doThisLoop));
		
		/* take rvectors in paramContainer and convert to parameter
		 * particles in pset_ptr_tmp */
		for ( ; rit < erit; ++rit) {
			rvector rv = *rit;
			vector < double > values(nParams, 0.0);
			// untransform again as we go
			for (size_t i = 0; i < nParams; ++i) {
				values[i] = _double(rv[i+1]) * scale[i] + shift[i];
			}
			boost::shared_ptr< ParameterParticle > ptmp(
				new ParameterParticle(values) );
			pset_ptr_tmp->add(ptmp);
		}
		
		// get a set of values from the sampler, using the particles
		boost::shared_ptr < MultiLociPolyTableSet > sampledValues_ptr = 
				sampler_ptr->sample(n_reps_per_particle, 
					pop,
					nloci_nsites,
					pset_ptr_tmp);
		
		cout << "\tgetting summary stats set" << endl;
	
		// and get a set of summary statistics from the values
		boost::shared_ptr < mct::SummaryStatisticSet > sset_ptr_tmp =
						sampledValues_ptr->getSummaryStatisticSet();
	
		// add the req stats from this lot into the set
		cout << "now getting required stats" << endl;
		if (!reqStats.empty()) {
			sset_ptr->add(sset_ptr_tmp->makeSelectedSummaryStatisticSet(reqStats));
		}
		else {
			sset_ptr->add(sset_ptr_tmp);
		}
	
		// make sure we add these parameters into to the total
		pset_ptr->add(pset_ptr_tmp);
		
		assert (doThisLoop <= toDo);
		toDo -= doThisLoop;
		cout << "\ttoDo is now " << toDo << endl;
	}
	

	
	assert(pset_ptr->size() == n_reps);
	assert(sset_ptr->size() == n_reps);
	
	{
		//summary stats from particles, just so we can do descriptive stats
		boost::shared_ptr < mct::SummaryStatisticSet > ss_tmp =
					pset_ptr->makeSummaryStatisticSet();
		
		{
			std::cout << "\nTry to get descriptive stats from particles for posterior predictive" << std::endl;
			outputAHABCLogFileString(logFilename, "\nparameters for posterio predictive");
			
			for (size_t i = 0; i < ss_tmp->at(0)->size() ; ++i) {
					
				cout << "Descriptive stats for parameters index " << i << endl;
				DescriptiveStats ds(ss_tmp, i);
				std::cout << (ds.toString()) << std::endl;
				outputAHABCLogFileString(logFilename, ds.toString());
			}
		}
	}
					
	{
		std::cout << "\nTry to get descriptive stats from summary stats for posterior predictive" << std::endl;
		outputAHABCLogFileString(logFilename, "\nsummary stats for posterior predictive");
		
		for (size_t i = 0; i < sset_ptr->at(0)->size() ; ++i) {
				
			cout << "Descriptive stats for summary stats index " << i << endl;
			DescriptiveStats ds(sset_ptr, i);
			std::cout << (ds.toString()) << std::endl;
			outputAHABCLogFileString(logFilename, ds.toString());
		}
	}
	
	// make an ahabcvalues object from the particles and summary stats
	
	boost::shared_ptr < const AHABCValues  > ahabcValuesPtr
		= boost::shared_ptr < const AHABCValues  >(
			new AHABCValues(pset_ptr, sset_ptr));

	return ahabcValuesPtr;
}


