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
\brief Basic testing of AHABC using NormalSD sampler, 
1 parameter, 1 summary stat.


*/

#include "automcmc_factory.hpp"
#include "testing_tools.hpp"

#include "test_routines.hpp"

#include "piecewise_constant_function.hpp"
#include "spsnode_measure_obj.hpp"

#include "config.h" // mct

#include "ahabc_values.hpp"

#include "descriptive_stats.hpp"
#include "normalsd_test_sampler.hpp"
#include "normalsd_test.hpp"

#include "generic_data_object_set.hpp"
#include "parameter_particle_set.hpp"
#include "summary_statistic_set.hpp"

#include "parameter_prior_uniform.hpp"
#include "parameter_prior.hpp"
#include "parameter_prior_set.hpp"


#include "prng.hpp"
#include "prng_gsl.hpp"

#include "utilities.hpp"

#include <mscplusplus/population_structure.hpp>

//#include <boost/lambda/bind.hpp>
//#include <boost/lambda/lambda.hpp>


#include <boost/shared_ptr.hpp>

#include <real.hpp>


#include <iostream>
#include <vector>
#include <iterator>
#include <cassert>

using namespace std;
using namespace mct;
using namespace aabc;
using namespace subpavings;


boost::shared_ptr < mct::SummaryStatisticSet > makeRefSet(
		size_t n_obs,
		double trueSD,
		int seed);

boost::shared_ptr < const AHABCValues  > makeValuesPtr(
		size_t n_reps,
		int seed);
		

void carverMCMCGR(
				boost::shared_ptr < const AHABCValues  > ahabcValuesPtr,
				boost::shared_ptr < mct::SummaryStatisticSet > refset_ptr,
				const std::vector < double >& trueParams,
				const std::string& prefix,
				const std::string& spec);

void carverStandardisedMCMCGR(
				boost::shared_ptr < const AHABCValues  > ahabcValuesPtr,
				boost::shared_ptr < mct::SummaryStatisticSet > refset_ptr,
				const std::vector < double >& trueParams,
				const std::string& prefix,
				const std::string& spec);
				

int main()
{
	int myPid = getpid();
	std::cout << "basic ahabc test for normal sd test " << std::endl;
	std::cout << "This process id is " << myPid << std::endl;
	
	size_t n_reps = 10000; // number of joint observations to generate
	
	size_t n_obs = 50; // number of 'observed values' to have
	
	int seed = 2345;
	
	//true sd
	double trueSD = 1.0;
	
	std::vector < double > trueParams (1, trueSD);
	
	cout << "\n\nn_reps = " << n_reps << ", n_obs = " << n_obs << endl;
	
	cout << "true sd " << trueSD << endl;
	
	string spec;
	{
		std::ostringstream stm;
		stm << "nr_" << n_reps;
		stm.precision(mct_utilities::fitPrecision(trueSD));
		stm <<   fixed;
		stm << "_sd_" << trueSD;
				
		spec = stm.str();
	}
	
	boost::shared_ptr < const AHABCValues  > ahabcValuesPtr 
					= makeValuesPtr(n_reps, seed);

	assert (ahabcValuesPtr->size() == n_reps);
		

	boost::shared_ptr < mct::SummaryStatisticSet > refset_ptr
					= makeRefSet(n_obs, trueSD, seed);

	string refstatsFilename;
	{
		std::ostringstream stm;
		stm << "../figures/NormalSD/normalsdRefStats_" << spec << ".txt";
		refstatsFilename = stm.str();
	}
	refset_ptr->outputToFile(refstatsFilename, false); // overwrite, not append
	
	string prefix("../figures/NormalSD/normalsdGR");
	
	carverMCMCGR(ahabcValuesPtr, refset_ptr, 
					trueParams, prefix, spec);
	//carverMCMC();
	//carverStartPoint(ahabcValuesPtr, refset_ptr, prefix);
				
	//basicMCMC();
	
	
	return 0;
	
}


	

// make some values
boost::shared_ptr < const AHABCValues  > makeValuesPtr(
		size_t n_reps,
		int seed)
{
	
	boost::shared_ptr < mct::PRNGen > rr ( new mct::PRNGenGSL(seed) );
	// one prng each for each prior and sampler
	boost::shared_ptr < PRNGen > rr_prior_sd = rr->spawnAnother();
	boost::shared_ptr < PRNGen > rr_sampler = rr->spawnAnother();
	
	// number of repetitions in each comparison sample generated for each parameter particle
	size_t n_reps_per_particle = 1; 
	
	size_t n_params = 1; // number of parameters we are using
	
	ParameterPriorSet priors(n_params); // prior set with n_params element
	/* set that element to a uniform prior on -10, 10
	 * and give the prior its random number generator */
	priors.at(0) = boost::shared_ptr < ParameterPrior > 
				( new ParameterPriorUniform(0.001,3,
											rr_prior_sd));
	
	// make asampler with its prng
	boost::shared_ptr < NormalSDTestSampler > 
					sampler_ptr ( new NormalSDTestSampler(rr_sampler) );
	
	// get a set of parameter particles from the priors
	boost::shared_ptr< ParameterParticleSet > 
			pset_ptr = priors.generateParameterParticleSet(n_reps);
			
	cout << "Size of particles is " << pset_ptr->size() << endl;

	{
		//summary stats from particles, just so we can do descriptive stats
		boost::shared_ptr < mct::SummaryStatisticSet > ss_tmp =
					pset_ptr->makeSummaryStatisticSet();
		
		{
			std::cout << "\nTry to get descriptive stats from particles" << std::endl;
			
			for (size_t i = 0; i < ss_tmp->at(0)->size() ; ++i) {
					
				cout << "Descriptive stats for parameters index " << i << endl;
				DescriptiveStats ds(ss_tmp, i);
				std::cout << (ds.toString()) << std::endl;
			}
		}
	}
	
	// get a set of values from the sampler, using the particles
	boost::shared_ptr < GenericDataObjectSet > sampledValues_ptr = 
			sampler_ptr->sample(n_reps_per_particle, 
				pset_ptr);
	
	// and get a set of summary statistics from the values
	boost::shared_ptr < mct::SummaryStatisticSet > sset_ptr =
					sampledValues_ptr->getSummaryStatisticSet();
					
	cout << "Size of summary stats is " << sset_ptr->size() << endl;
	
	assert(!sset_ptr->empty());
	
	{
		std::cout << "\nTry to get descriptive stats from summary stats" << std::endl;
		
		for (size_t i = 0; i < sset_ptr->at(0)->size() ; ++i) {
				
			cout << "Descriptive stats for summary stats index " << i << endl;
			DescriptiveStats ds(sset_ptr, i);
			std::cout << (ds.toString()) << std::endl;
		}
	}
	
	// make an ahabcvalues object from the particles and summary stats
	
	boost::shared_ptr < const AHABCValues  > ahabcValuesPtr
		= boost::shared_ptr < const AHABCValues  >(
			new AHABCValues(pset_ptr, sset_ptr));

	return ahabcValuesPtr;
}

// make some reference statistics
boost::shared_ptr < mct::SummaryStatisticSet > makeRefSet(
		size_t n_obs,
		double trueSD,
		int seed)
{
	
	int trueseed = seed+1;
	boost::shared_ptr < mct::PRNGen > rr_sampler ( new mct::PRNGenGSL(trueseed) );
	
	// make a  sampler with its prng
	boost::shared_ptr < NormalSDTestSampler > 
					sampler_ptr ( new NormalSDTestSampler(rr_sampler) );
	
	double realvalues[] = {trueSD}; 
	
	boost::shared_ptr < aabc::ParameterParticle > 
			p( new aabc::ParameterParticle(realvalues, 1) );

	// get a set of values from the sampler, using the particles
	boost::shared_ptr < GenericDataObjectSet > refValues_ptr = 
			sampler_ptr->sample(n_obs, p);
	
	// and get a set of reference summary statistics from the values
	boost::shared_ptr < mct::SummaryStatisticSet > refset_ptr =
					refValues_ptr->getSummaryStatisticSet();

	cout << "Observed summary statistics are:" << endl;
	cout << refset_ptr->toString() << endl;

	
	return refset_ptr;
}

// test Carver with GR convergence
void carverMCMCGR(
				boost::shared_ptr < const AHABCValues  > ahabcValuesPtr,
				boost::shared_ptr < mct::SummaryStatisticSet > refset_ptr,
				const std::vector < double >& trueParams,
				const std::string& prefix,
				const std::string& spec)
{
	try {
		
		cout << "Getting the values " << endl;
		//fill a container with the parameter particle and summary statistic values
		std::vector < std::vector < double > > container;
		ahabcValuesPtr->fillDataContainer(container);
		
		std::string containerFilename = prefix + "Data_" + spec + ".txt"; 
		mct_utilities::outputVectorVectorToFile(container,
						containerFilename, false); // overwrite
		
		std::vector < std::vector < double > > slicePtsVec = 
						getSlicePoints(refset_ptr);
		{
			cout << "\nref stats are = " << endl;
			for (std::vector< std::vector<double> >::iterator it = slicePtsVec.begin();
							it < slicePtsVec.end();
							++it) {
				ostream_iterator<double> out_it (cout,"\t");
				copy ( it->begin(), it->end(), out_it );
				cout << endl;
			}
			cout << endl;
		}
		
		std::vector < int > sliceDims = ahabcValuesPtr->getSummaryStatisticsDims();
		
		size_t maxLeavesCarving = 1000; // leaves in carver queue
	
		size_t critSEB = 5;  // max points in a node to stop posterior queue
			
		size_t maxLeavesSEB = 2000000; // maximum leaves we'll let SEB queue go to
									// even if critSEB not met
			
		size_t minPoints = 1;
		size_t keepBest = 2; // keep two best
			
		std::string postFileName = prefix + "LogPosteriorsForStarts.txt";
		std::string checkPostFileNameBase = prefix + "CheckLogPosteriors";
		
		string scalarsFileName = prefix;
		{
			ostringstream oss;
			oss << "_" << spec << ".log";
			scalarsFileName += oss.str();
			
		}
		
		unsigned int maxLoops = 5000000;
		unsigned int thinout = 500;
		int samplesNeeded = 100;
		cxsc::real tolerance(0.1);
		
		cxsc::real padding = 0.000005;
		ivector box = getBox(container, slicePtsVec, sliceDims, padding);
		
		ahabc::AutoMCMCFactory::AUTOTYPE autotype = ahabc::AutoMCMCFactory::ALLTHREE;
	
		int rep = 1;
		
			
		boost::shared_ptr < subpavings::PiecewiseConstantFunction >
			mcmcAvPtr = testCarverMCMCGR(
											box,
											container,
											maxLeavesCarving,
											critSEB,	
											maxLeavesSEB,	
											minPoints,
											keepBest,
											postFileName,
											checkPostFileNameBase,
											scalarsFileName,
											maxLoops,
											thinout,
											samplesNeeded,
											tolerance,
											rep,
											autotype);
		
		
		int prec = 5;
		bool confirm = true;
		size_t posteriorCaptureInterval = 10;
		
		
		processAverage(
			*mcmcAvPtr,
			ahabcValuesPtr,
			prefix,
			spec,
			sliceDims,
			slicePtsVec,
			trueParams,
			prec,
			confirm,
			posteriorCaptureInterval);
		
		
		
	}
	catch (std::exception& e) {
	
		std::cout << "std::exception::\n" + std::string(e.what()) << std::endl;
		throw;
	}
	catch (...) {
	
		std::cout << "Unknown error\n" << std::endl;
		throw;
	}
}



#if(0)

// test Carver with GR convergence, standardised data
void carverStandardisedMCMCGR(
				boost::shared_ptr < const AHABCValues  > ahabcValuesPtr,
				boost::shared_ptr < mct::SummaryStatisticSet > refset_ptr,
				const std::vector < double >& trueParams,
				const std::string& prefix,
				const std::string& spec)
{
	try {
		
		cout << "Getting the standardised values " << endl;
		//fill a container with the parameter particle and summary statistic values
		std::vector < std::vector < double > > container;
		ahabcValuesPtr->fillStandardisedDataContainer(container);
		
		std::string containerFilename = prefix + "Data" + spec + ".txt"; 
		mct_utilities::outputVectorVectorToFile(container,
						containerFilename, false); // overwrite
		
		std::vector < std::vector < double > > slicePtsVec = 
		getStandardisedSlicePoints(refset_ptr, ahabcValuesPtr);
		{
			cout << "\nstandardised ref stats are = " << endl;
			for (std::vector< std::vector<double> >::iterator it = slicePtsVec.begin();
							it < slicePtsVec.end();
							++it) {
				ostream_iterator<double> out_it (cout,"\t");
				copy ( it->begin(), it->end(), out_it );
				cout << endl;
			}
			cout << endl;
		}
		
		std::vector < int > sliceDims = ahabcValuesPtr->getSummaryStatisticsDims();
		
		size_t maxLeavesCarving = 1000; // leaves in carver queue
	
		size_t critSEB = 5;  // max points in a node to stop posterior queue
			
		size_t maxLeavesSEB = 2000000; // maximum leaves we'll let SEB queue go to
									// even if critSEB not met
			
		size_t minPoints = 1;
		size_t keepBest = 2; // keep two best
			
		std::string postFileName = prefix + "LogPosteriorsForStarts.txt";
		std::string checkPostFileNameBase = prefix + "CheckLogPosteriors";
		
		string scalarsFileName = prefix;
		{
			ostringstream oss;
			oss << "_" << spec << ".log";
			scalarsFileName += oss.str();
			
		}
		
		unsigned int maxLoops = 5000000;
		unsigned int thinout = 500;
		int samplesNeeded = 100;
		cxsc::real tolerance(0.1);
		
		cxsc::real padding = 0.000005;
		ivector box = getBox(container, slicePtsVec, sliceDims, padding);
		
		ahabc::AutoMCMCFactory::AUTOTYPE autotype = ahabc::AutoMCMCFactory::ALLTHREE;
	
		int rep = 1;
			
		boost::shared_ptr < subpavings::PiecewiseConstantFunction >
			mcmcAvPtr = testCarverMCMCGR(
										box,
										container,
										maxLeavesCarving,
										critSEB,	
										maxLeavesSEB,	
										minPoints,
										keepBest,
										postFileName,
										checkPostFileNameBase,
										scalarsFileName,
										maxLoops,
										thinout,
										samplesNeeded,
										tolerance,
										rep,
										autotype);
		
		int prec = 5;
		bool confirm = true;
		
		processStandardisedAverage(
			*mcmcAvPtr,
			ahabcValuesPtr,
			prefix,
			spec,
			sliceDims,
			refset_ptr,
			slicePtsVec,
			trueParams,
			prec,
			confirm);
		
		
		
	}
	catch (std::exception& e) {
	
		std::cout << "std::exception::\n" + std::string(e.what()) << std::endl;
		throw;
	}
	catch (...) {
	
		std::cout << "Unknown error\n" << std::endl;
		throw;
	}
}
#endif



