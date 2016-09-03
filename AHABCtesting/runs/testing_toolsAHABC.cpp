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
\brief Tools for doing various AHABC routines

*/

#include "testing_toolsAHABC.hpp"
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

#include "adaptivehistogram.hpp"



#include <real.hpp>


#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <cassert>


//#define DEBUG_GRIDDING
//#define DEBUG_GRIDDING_OUTPUT

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
		const std::vector< std::pair<double, double> >& priorLims,
		const std::vector < size_t >& reqStats,
		int seed,
		const std::string& logFilename)
{
	
	boost::shared_ptr < mct::PRNGen > rr ( new mct::PRNGenGSL(seed) );
	
	size_t n_params = priorLims.size(); // number of parameters we are using
	
	ParameterPriorSet priors(n_params); // prior set with n_params element
	/* set that element to a uniform prior
	 * and give the prior its random number generator */
	for (size_t j = 0; j < n_params; ++j) {
		boost::shared_ptr < PRNGen > rr_prior = rr->spawnAnother();
		priors.at(j) = boost::shared_ptr < ParameterPrior > 
					( new ParameterPriorUniform(priorLims[j].first,
												priorLims[j].second,
												rr_prior));
	}
	
	//summary of priors to logfile
	outputAHABCLogFileString(logFilename, priors.toString());
	
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
	
	while (toDo > 0) {
		
		size_t doThisLoop = (toDo <= maxToDo ? toDo : maxToDo);
				
		cout << "\tsampling " << doThisLoop << endl;
	
		// get a set of parameter particles from the priors
		boost::shared_ptr< ParameterParticleSet > 
		pset_ptr_tmp = priors.generateParameterParticleSet(doThisLoop);
		
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
			std::cout << "\nTry to get descriptive stats from particles" << std::endl;
			outputAHABCLogFileString(logFilename, "\nparameters");
			
			for (size_t i = 0; i < ss_tmp->at(0)->size() ; ++i) {
					
				cout << "Descriptive stats for parameters index " << i << endl;
				DescriptiveStats ds(ss_tmp, i);
				std::cout << (ds.toString()) << std::endl;
				outputAHABCLogFileString(logFilename, ds.toString());
			}
		}
	}
					
	{
		std::cout << "\nTry to get descriptive stats from summary stats" << std::endl;
		outputAHABCLogFileString(logFilename, "\nsummary stats");
		
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





std::vector< std::pair<double, double> > getPriorLims(
		size_t n_obs,
		size_t nloci,
		size_t nsites,
		boost::shared_ptr< hudson_ms::PopulationStructure > pop,
		size_t n_reps_test,	
		const std::vector < double >& trueParams,
		const std::vector< std::pair<double, double> >& priorLims,
		const std::vector < size_t >& reqStats,
		int seed,
		size_t gridSegments,
		double percentClosest)
{
	
	int trueseed = seed + 1;
	
	boost::shared_ptr < mct::SummaryStatisticSet > refset_ptr
					= makeAHABCRefSet(n_obs, nloci, nsites,
								pop, trueParams, reqStats,
								trueseed);
	#ifdef DEBUG_GRIDDING
		cout << "size of refset is " << refset_ptr->size() << endl;
	#endif
	
	return getPriorLims(
		refset_ptr,
		nloci,
		nsites,
		pop,
		n_reps_test,	
		priorLims,
		reqStats,
		seed,
		gridSegments,
		percentClosest);

}

std::vector< std::pair<double, double> > getPriorLims(
		boost::shared_ptr < mct::SummaryStatisticSet > refset_ptr,
		size_t nloci,
		size_t nsites,
		boost::shared_ptr< hudson_ms::PopulationStructure > pop,
		size_t n_reps_test,	
		const std::vector< std::pair<double, double> >& priorLims,
		const std::vector < size_t >& reqStats,
		int seed,
		size_t gridSegments,
		double percentClosest)
{
	
	std::vector < std::vector < double > > refValues;
	size_t nrs = refset_ptr->size();
	
	for(size_t i = 0; i < nrs; ++i) {
		refValues.push_back(refset_ptr->at(i)->getValues());
	} 
	
	boost::shared_ptr < mct::PRNGen > rr ( new mct::PRNGenGSL(seed) );
	
	size_t n_params = priorLims.size(); // number of parameters we are using
	
	std::vector < double >tmp;
	tmp.reserve(gridSegments+1);
	std::vector < std::vector < double > > 
				priorGridSegments(n_params, tmp);
	std::vector < double > priorGridSegmentSizes(n_params, 0.0);
	
	#ifdef DEBUG_GRIDDING
		cout << "About to do the gridding" << endl;
	#endif
	for (size_t i = 0; i < n_params; ++i) {
		
		// get side space and divvie it up
		double min = priorLims[i].first;
		double max = priorLims[i].second;
		double seg = (max - min)/gridSegments;
		priorGridSegmentSizes[i] = seg;
		priorGridSegments[i].push_back(min);
		for (size_t j = 1; j < gridSegments ; ++j) {
			priorGridSegments[i].push_back(min + j*seg);
		}
		priorGridSegments[i].push_back(max);
		
	}
	
	#ifdef DEBUG_GRIDDING
		cout << "About to build the points" << endl;
	#endif
	std::vector < std::vector < double > >	gridPoints;
	gridPoints.reserve(
		static_cast<size_t>(std::pow(gridSegments+1,n_params)));
	std::vector < double > gridPoint;
	gridPoint.reserve(n_params);
	buildGridPoint(gridPoints,
		gridPoint,
		priorGridSegments,
		0);
		
	#ifdef DEBUG_GRIDDING
		cout << "Grid points size = " << gridPoints.size() << endl;
			
		{
			cout << "Priors grid is" << endl;
			for (size_t i = 0; i < n_params; ++i) {
				ostream_iterator<double> out_it (cout,"\t");
				copy ( priorGridSegments[i].begin(), priorGridSegments[i].end(), out_it );
				cout <<  endl;
			}
		}
		{
			cout << "Grid points are is" << endl;
			for (size_t i = 0; i < gridPoints.size(); ++i) {
				ostream_iterator<double> out_it (cout," ");
				copy ( gridPoints[i].begin(), gridPoints[i].end(), out_it );
				cout <<  endl;
			}
		}
		
	#endif
	
	
	//vector of nloci copies of nsites
	std::vector < size_t > nloci_nsites(nloci, nsites);
		
	// prng for sampler
	boost::shared_ptr < PRNGen > rr_sampler = rr->spawnAnother();
	
	// make a simdata (infinite sites model) sampler with its prng
	boost::shared_ptr < MultiLociPolyTableSampler > 
					sampler_ptr ( new MultiLociSimDataSampler(rr_sampler) );
	
	size_t npts = gridPoints.size();
	
	/* Overall scores for each point*/ 
	std::vector < int > overallPtScores(npts, 0);
	
	//for each ob 
	for(size_t i = 0; i < nrs; ++i) {
		//ob is refValues[i]
		
		//point scores for this ob
		std::vector < double> distances;
		distances.reserve(n_reps_test*npts);
		std::vector < std::vector < double> > 
		distanceByPt(npts, std::vector < double>(n_reps_test, 0.0));
		
		//for each particle
		for (size_t j = 0; j < npts; ++j) {
		
			boost::shared_ptr < ParameterParticle > 
					pPtr( new ParameterParticle(gridPoints[j]));
			
			// get a set of values from the sampler, using the particles
			boost::shared_ptr < MultiLociPolyTableSet > sampledValues_ptr = 
					sampler_ptr->sample(n_reps_test, 
						pop,
						nloci_nsites,
						pPtr);
			
			// and get a set of summary statistics from the values
			boost::shared_ptr < mct::SummaryStatisticSet > sset_ptr;
			boost::shared_ptr < mct::SummaryStatisticSet > sset_ptr_tmp =
							sampledValues_ptr->getSummaryStatisticSet();
		
			// add the req stats from this lot into the set
			if (!reqStats.empty()) {
				sset_ptr = sset_ptr_tmp->makeSelectedSummaryStatisticSet(reqStats);
			}
			else {
				sset_ptr= sset_ptr_tmp;
			}
			
			// for each stat in the set
			for (size_t k = 0; k < n_reps_test; k++) {
				double dis = getDistance(refValues[i],	sset_ptr->at(k));
				distances.push_back(dis);
				distanceByPt[j][k] = dis;
			}
		}	
		
		/* then we would find the percentClosest point cutoff*/
		double cutoff = getPercentile(distances, percentClosest);
		#ifdef DEBUG_GRIDDING
			cout << "cutoff = " << cutoff << endl;
		#endif
		
		/* and go through pt by pt in distancesByPt
		 * then put a 1 in overallPtScores for every pt where
		 * at least one of the nss distancesByPt[i] values is <= cutoff */
		
		for (size_t j = 0; j < npts; ++j) {
			if (!overallPtScores[j]) {
				for (size_t k = 0; k < n_reps_test; ++k) {
					if (!(distanceByPt[j][k] > cutoff)) {
						overallPtScores[j] = 1;
						break;
					}
				}
			}
		
		}

	} // end loop through obs
	
	#ifdef DEBUG_GRIDDING
		cout << "Grid points and grid point scores are is" << endl;
			for (size_t i = 0; i < gridPoints.size(); ++i) {
				ostream_iterator<double> out_it (cout," ");
				copy ( gridPoints[i].begin(), gridPoints[i].end(), out_it );
				cout << "\tscore\t"<< overallPtScores[i] << endl;
			}
		
	#endif
	
	/*set up pair min and max for each param and set initial values 
	for each to the elements of the first point in gridPoints that has
	score 1.*/
	
	/* now go through overallPtScores and adjust min and max for 
	 * any point with score 1 */
	
	
	std::vector< std::pair<double, double> > 
			newPriorLims(n_params, std::pair<double, double>(0.0,0.0));
	
	size_t ptsIn = 0;
	for (size_t j = 0; j < npts; ++j) {
		
		if (overallPtScores[j]) {
			for (size_t k = 0; k < n_params; ++k) {
				
				if (ptsIn) {
					if (gridPoints[j][k] < newPriorLims[k].first)  
						newPriorLims[k].first = gridPoints[j][k];
					else if (gridPoints[j][k] > newPriorLims[k].second)
						newPriorLims[k].second = gridPoints[j][k];
					
				}
				else { // first one
					newPriorLims[k].first = gridPoints[j][k];
					newPriorLims[k].second = gridPoints[j][k];
				
				}
			}
			++ptsIn;
		}
	}
	
	if (ptsIn) {
		/* and now go out one segment in each direction */
		for (size_t i = 0; i < n_params; ++i) {
			if (!(newPriorLims[i].first - priorGridSegmentSizes[i]
				< priorLims[i].first)) {
					newPriorLims[i].first -= priorGridSegmentSizes[i];
			}
			else newPriorLims[i].first = priorLims[i].first;
			
			if (!(newPriorLims[i].second + priorGridSegmentSizes[i]
				> priorLims[i].second)) {
					newPriorLims[i].second += priorGridSegmentSizes[i];
			}
			else newPriorLims[i].second = priorLims[i].second;
			
		}
	}
	/* If 0 of the gridpoints in priors is close enough, 
	 * just return the old priors */
	else newPriorLims = priorLims;
	
	#ifdef DEBUG_GRIDDING_OUTPUT
		cout << "\nNew prior limits are " << endl;
		for (size_t i = 0; i < n_params; ++i) {
			cout << newPriorLims[i].first << "\t" << newPriorLims[i].second
				<< "\t(" << priorLims[i].first << "\t"
					<< priorLims[i].second << ")" << endl;
		}
	#endif
	
	return newPriorLims;
	
}


