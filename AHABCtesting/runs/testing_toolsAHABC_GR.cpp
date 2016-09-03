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
\brief Tools for doing various AHABC routines using MCMC histogram average

*/

#include "testing_toolsAHABC_GR.hpp"
#include "testing_toolsAHABC.hpp"
#include "testing_toolsGR.hpp"
#include "testing_toolsMCMCGR.hpp"
#include "testing_toolsSimdata.hpp"
#include "testing_toolsMRS.hpp"
#include "testing_tools.hpp"
#include "log_memory.hpp"



#include "config.h" // mct

//#include "parameter_prior_uniform.hpp"

//#include "parameter_prior.hpp"

//#include "multi_loci_simdata_sampler.hpp"
//#include "multi_loci_polytable_sampler.hpp"

//#include "parameter_particle_set.hpp"
//#include "descriptive_stats.hpp"
//#include "prng.hpp"
//#include "prng_gsl.hpp"

#include "utilities.hpp"

//#include "adaptivehistogram.hpp"

#include "toolz.hpp" // for realVolume

#include <real.hpp>


#include <iostream>
//#include <fstream>
//#include <sstream>
//#include <ctime>
//#include <cmath>
//#include <numeric>
//#include <algorithm>
#include <cassert>

//#include <boost/lambda/bind.hpp>
//#include <boost/lambda/lambda.hpp>

//#include <gsl/gsl_statistics_double.h>

//#define MINVOLMAIN
#define MINDEPTHMAIN 27
#define MINVOLPRED
#define MINDEPTHPRED 20


using namespace std;
using namespace mct;
using namespace aabc;
using namespace subpavings;


void doAHABCHeterozygosity(
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
				int simSeed,
				ahabc::AutoMCMCFactory::AUTOTYPE autotype,
				bool reducePriorLims)
{
	//dummy
	std::vector < std::vector < double > > otherParams;
	doAHABCHeterozygosity(
				prefix,
				spec,
				n_obs,
				n_reps,
				nloci,
				nsites,
				priorLims,
				pop,
				n_reps_per_particle,
				trueParams,
				otherParams,
				seed,
				simSeed,
				autotype,
				reducePriorLims);
}

void doAHABCHeterozygosity(
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
				int simSeed,
				ahabc::AutoMCMCFactory::AUTOTYPE autotype,
				bool reducePriorLims)
{
	
	// only want to use heterozygosity stat
	std::vector < size_t > reqStats(1, 2);
	
	doAHABC(
				prefix,
				spec,
				n_obs,
				n_reps,
				nloci,
				nsites,
				priorLims,
				pop,
				n_reps_per_particle,
				trueParams,
				otherParams,
				seed,
				simSeed,
				reqStats,
				autotype,
				reducePriorLims);
	
}
void doAHABCHeterozygosityAndSegStats(
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
				int simSeed,
				ahabc::AutoMCMCFactory::AUTOTYPE autotype,
				bool reducePriorLims)
{
	//dummy
	std::vector < std::vector < double > > otherParams;
	doAHABCHeterozygosityAndSegStats(
				prefix,
				spec,
				n_obs,
				n_reps,
				nloci,
				nsites,
				priorLims,
				pop,
				n_reps_per_particle,
				trueParams,
				otherParams,
				seed,
				simSeed,
				autotype,
				reducePriorLims);
}

void doAHABCHeterozygosityAndSegStats(
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
				int simSeed,
				ahabc::AutoMCMCFactory::AUTOTYPE autotype,
				bool reducePriorLims)
{
	
	// dummy
	std::vector < size_t > reqStats;
	
	doAHABC(
				prefix,
				spec,
				n_obs,
				n_reps,
				nloci,
				nsites,
				priorLims,
				pop,
				n_reps_per_particle,
				trueParams,
				otherParams,
				seed,
				simSeed,
				reqStats,
				autotype,
				reducePriorLims);
	
}

void doAHABC(
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
				int simSeed,
				const std::vector < size_t >& reqStats,
				ahabc::AutoMCMCFactory::AUTOTYPE autotype,
				bool reducePriorLims)
{
	string logFilename;
	{
		std::ostringstream stm;
		stm << prefix << "AHABCLog_" << spec << ".txt";
		logFilename = stm.str();
	}
	outputAHABCLogFileStart(logFilename);

	string memLogFilename;
	{
		std::ostringstream stm;
		stm << prefix << "MemLog_" << spec << ".log";
		memLogFilename = stm.str();
	}
	
	logUsage(memLogFilename, "Start doAHABC", false); // new mem log file
	
	clock_t startTime = clock();
		
	/*stuff I know I'll need */
	boost::shared_ptr < subpavings::PiecewiseConstantFunction >
			mcmcAvPtr;
	std::vector < double > sds;
	std::vector < double > means;
	int np;
	int ns;
	std::vector < std::vector < double > > slicePtsVec;
	std::vector < int > sliceDims;
	
	int trueseed = seed + 1;
	
	{
		std::ostringstream stm;
		stm << "\n\nSeed for 'true' values (ref stats)\t" << trueseed;
		outputAHABCLogFileString(logFilename, stm.str());
	}
	
	boost::shared_ptr < mct::SummaryStatisticSet > refset_ptr
					= makeAHABCRefSet(n_obs, nloci, nsites,
								pop, trueParams, reqStats,
								trueseed);
	
	string refstatsFilename;
	{
		std::ostringstream stm;
		stm << prefix << "RefStats_" << spec << ".txt";
		refstatsFilename = stm.str();
	}
	refset_ptr->outputToFile(refstatsFilename, false); // overwrite, not append
	
	
		
	cout << "\nFirst 10 ref stats are = " << endl;
	for (size_t i = 0; (i < 10) && (i < refset_ptr->size()); ++i) {
		cout << refset_ptr->at(i)->toString() << endl;
		
	}
	
	logUsage(memLogFilename, "Made refstats");
	
	clock_t endTime = clock();
	{
		double timing = ((static_cast<double>(endTime - startTime)) / CLOCKS_PER_SEC);
		std::ostringstream stm;
		stm << "\nTiming to make ref stats\t" << timing;
		outputAHABCLogFileString(logFilename, stm.str());
	}
	
	std::vector< std::pair<double, double> > newPriorLims = priorLims;
	if (reducePriorLims) {
		
		
	
		double percentClosest = 0.2;
		size_t n_reps_test = 100;
		size_t gridSegments = 10;
		cout << "\nReducing prior limits, n_reps_test = " << n_reps_test
			<< ", gridSegments = " << gridSegments 
			<< ", percentClosest = " << percentClosest << endl;
		
		newPriorLims = getPriorLims(refset_ptr,
									nloci,
									nsites,
									pop,
									n_reps_test,	
									priorLims,
									reqStats,
									seed,
									gridSegments,
									percentClosest);
		
		cout << "\nNew prior limits are " << endl;
		for (size_t i = 0; i < priorLims.size(); ++i) {
			cout << newPriorLims[i].first << "\t" << newPriorLims[i].second
				<< "\t(" << priorLims[i].first << "\t"
					<< priorLims[i].second << ")" << endl;
		}
		endTime = clock();
		{
			double timing = ((static_cast<double>(endTime - startTime)) / CLOCKS_PER_SEC);
			std::ostringstream stm;
			stm << "\n\nTiming to reduce prior lims\t" << timing;
			outputAHABCLogFileString(logFilename, stm.str());
		}
		startTime = clock();
	
	}
		
	/* put a lot of the work in it's own scope so that anything
	 * we do not need goes out of scope when we've done the heavy lifting... */
	{
		//fill a container with the parameter particle and summary statistic values
		boost::shared_ptr < std::vector < std::vector < double > > >
			containerPtr (new std::vector < std::vector < double > >() );
		
		{
			
			{
				std::ostringstream stm;
				stm << "\n\nSeed for simulations\t" << simSeed;
				outputAHABCLogFileString(logFilename, stm.str());
			}
			// pass in the newPriorLims
			boost::shared_ptr < const AHABCValues  > ahabcValuesPtr 
						= makeValuesPtr(n_reps, nloci, nsites, 
										pop, n_reps_per_particle,
										newPriorLims, reqStats,
										simSeed,
										logFilename);

			endTime = clock();
			{
				double timing = ((static_cast<double>(endTime - startTime)) / CLOCKS_PER_SEC);
				std::ostringstream stm;
				stm << "\nTiming to make joint samples\t" << timing;
				outputAHABCLogFileString(logFilename, stm.str());
			}
			startTime = clock();
		
			assert (ahabcValuesPtr->size() == n_reps);
			
			logUsage(memLogFilename, "Made AHABC values");
				
			means = ahabcValuesPtr->getMeans();
			
			sds = ahabcValuesPtr->getSDs();
		
			np = ahabcValuesPtr->getNumberParameters();
			ns = ahabcValuesPtr->getNumberSummaryStatistics();
			
			sliceDims = ahabcValuesPtr->getSummaryStatisticsDims();
			
			/* get and output the unstandardised data */
			std::vector < std::vector < double > > containerTmp;
			ahabcValuesPtr->fillDataContainer(containerTmp);
			std::string containerFilenameTmp = prefix + "UnstandardisedData_" + spec + ".txt"; 
			mct_utilities::outputVectorVectorToFile(containerTmp,
						containerFilenameTmp, false); // overwrite
			
			ahabcValuesPtr->fillStandardisedDataContainer(*containerPtr);
			
			endTime = clock();
			{
				double timing = ((static_cast<double>(endTime - startTime)) / CLOCKS_PER_SEC);
				std::ostringstream stm;
				stm << "\n\nTiming to standardise data\t" << timing;
				outputAHABCLogFileString(logFilename, stm.str());
			}
			startTime = clock();
			
			logUsage(memLogFilename, "Filled container");
	
		}
		cout << "ahabc values should go out of scope" << endl;
		logUsage(memLogFilename, "xxAHABC values");
			
		std::string containerFilename = prefix + "StandardisedData_" + spec + ".txt"; 
		mct_utilities::outputVectorVectorToFile(*containerPtr,
					containerFilename, false); // overwrite
	
		slicePtsVec = getStandardisedSlicePoints(refset_ptr, np, means, sds);
	
		cxsc::real padding = 0.000005;
		ivector box = getBox(*containerPtr, slicePtsVec, sliceDims, padding);
		
		cout << "\nGetting the joint distn\n" << endl;
		#ifdef MINVOLMAIN
			cxsc::real boxVol = realVolume(box);
			cxsc::real minVol = 0.0;
			if (ns > 1) minVol = boxVol/cxsc::power(2.0, MINDEPTHMAIN);
			
			mcmcAvPtr = getMCMCGRaverage(
						*containerPtr,
						box,
						minVol,
						np,
						ns,
						prefix,
						spec,
						autotype,
						logFilename,
						memLogFilename);
		#else
			mcmcAvPtr = getMCMCGRaverage(
						*containerPtr,
						box,
						np,
						ns,
						prefix,
						spec,
						autotype,
						logFilename,
						memLogFilename);
		#endif
		endTime = clock();
		{
			double timing = ((static_cast<double>(endTime - startTime)) / CLOCKS_PER_SEC);
			std::ostringstream stm;
			stm << "\n\nTiming to make joint rmrp\t" << timing;
			outputAHABCLogFileString(logFilename, stm.str());
		}
		startTime = clock();
		
		logUsage(memLogFilename, "xxhistograms");
	
	}
	logUsage(memLogFilename, "xxcontainer");
		
	string meanAndSdFilename = prefix + "_MeanAndSd_" + spec + ".txt";

	printAndOutputMeansAndSDs(means,sds, meanAndSdFilename);
			
	{
		std::vector< std::vector<double> >::iterator eit = slicePtsVec.begin();
		advance(eit, 10);
		
		cout << "\nFirst 10 standardised ref stats are = " << endl;
		for (std::vector< std::vector<double> >::iterator it = slicePtsVec.begin();
						(it < eit) && (it < slicePtsVec.end());
						++it) {
			ostream_iterator<double> out_it (cout,"\t");
			copy ( it->begin(), it->end(), out_it );
			cout << endl;
		}
		cout << endl;
	}
	
	int prec = 5;
	bool confirm = true;
	size_t posteriorCaptureInterval = 10;
	
	std::vector < 
		boost::shared_ptr < PiecewiseConstantFunction > 
	> posteriorPtrsVec = processStandardisedAverage(
		*mcmcAvPtr,
		sds,
		means,
		np,
		ns,
		n_reps,
		prefix,
		spec,
		sliceDims,
		refset_ptr,
		slicePtsVec,
		trueParams,
		otherParams,
		prec,
		confirm,
		posteriorCaptureInterval);
	
	
	
	logUsage(memLogFilename, "At end of main ahabc routine");
	
	endTime = clock();
	
	double timing = ((static_cast<double>(endTime - startTime)) / CLOCKS_PER_SEC);
	
	{
		std::ostringstream stm;
		stm << "\n\nTiming to do processing\t" << timing;
		outputAHABCLogFileString(logFilename, stm.str());
	}
	
	cout << "\n\n\nEnd of main AHABC:  time = " << timing << "\n" <<endl;
	
	assert(!posteriorPtrsVec.empty());
	cout << "\nDo posterior predictive" <<endl;
	
	startTime = clock();
	
	/*make the posterior predictive by redoing the whole business  of 
	 * making the joint, and then marginalise over the parameters
	 * to get marginal on summary stats */
	
	/* I am just going to do it for the first and last posterior here, ie
	 * the one using posterior capture interval of the ref stats
	 * and the one using all of all the refstats.*/
	
	std::vector < size_t > posteriorObsVec;
	int ip = 1;
	while (ip*posteriorCaptureInterval < n_obs) {
		posteriorObsVec.push_back(ip*posteriorCaptureInterval);
		++ip;
	}
	posteriorObsVec.push_back(n_obs);
	size_t nPosts = posteriorObsVec.size();
	assert(nPosts == posteriorPtrsVec.size());
	std::vector < size_t > posteriorIndexVec(1, 0);
	if (nPosts > 1) posteriorIndexVec.push_back(nPosts - 1);
	for (size_t i = 0; i < posteriorIndexVec.size(); ++i) {
		
		size_t index = posteriorIndexVec[i];
		
		boost::shared_ptr < PiecewiseConstantFunction > pcfPtr 
				= posteriorPtrsVec.at(index);
		size_t nrefstats = posteriorObsVec.at(index);
		
		{
			ostringstream oss;
			oss << "Doing posterior predictive for nrefstats = " 
				<< nrefstats;
			logUsage(memLogFilename, oss.str());
			cout << "\n\n" << oss.str() << "\n" << endl;
		}		
		
		boost::shared_ptr < subpavings::PiecewiseConstantFunction >
				mcmcAvPostPtr;
		std::vector < double > sdsPost;
		std::vector < double > meansPost;
	
		std::vector < std::vector < double > > slicePtsPostVec;
		
		// only want some of the refset 
		boost::shared_ptr < mct::SummaryStatisticSet > refset_ptr_post
			(new mct::SummaryStatisticSet(nrefstats) );
		for (size_t i = 0; i < nrefstats; ++i ) {
			refset_ptr_post->add(refset_ptr->at(i));
			
		}
		//output them
		string refstatsFilename;
		{
			std::ostringstream stm;
			stm << prefix << "_PredictiveRefStats_nrefstats" 
					<< nrefstats << spec << ".txt";
			refstatsFilename = stm.str();
		}
		refset_ptr_post->outputToFile(refstatsFilename, false); // overwrite, not append
			
			
		{
			//fill a container with the parameter particle and summary statistic values
			boost::shared_ptr < std::vector < std::vector < double > > >
				containerPostPtr (new std::vector < std::vector < double > >() );
			
			
			{
				// pass in the posterior
				boost::shared_ptr < const AHABCValues  > ahabcValuesPtr 
							= makeValuesPtr(n_reps, nloci, nsites, 
											pop, n_reps_per_particle,
											*pcfPtr,
											sds, // orignal sds
											means, // and mean
											reqStats,
											seed,
											logFilename);

				assert (ahabcValuesPtr->size() == n_reps);
				
				logUsage(memLogFilename, "Made AHABC values for posterior");
					
				meansPost = ahabcValuesPtr->getMeans();
				
				sdsPost = ahabcValuesPtr->getSDs();
			
				ahabcValuesPtr->fillStandardisedDataContainer(*containerPostPtr);
				
				logUsage(memLogFilename, "Filled container for posterior");
		
			}
			logUsage(memLogFilename, "xxAHABC values for posterior");
				
			slicePtsPostVec = getStandardisedSlicePoints(refset_ptr_post, np, meansPost, sdsPost);
		
			cxsc::real padding = 0.000005;
			ivector box = getBox(*containerPostPtr, slicePtsPostVec, sliceDims, padding);
			
			string prefixPosterior;
			{
				ostringstream oss;
				oss << prefix << "Predictive" << nrefstats;
				prefixPosterior = oss.str();
			}
			
			
			cout << "\nGetting the posterior joint distn\n" << endl;
			#ifdef MINVOLPRED
				cxsc::real boxVol = realVolume(box);
				cxsc::real minVol = 0.0;
				if (ns > 1) minVol = boxVol/cxsc::power(2.0, MINDEPTHPRED);
				mcmcAvPostPtr = getMCMCGRaverage(
							*containerPostPtr,
							box,
							minVol,
							np,
							ns,
							prefixPosterior,
							spec,
							autotype,
							logFilename,
							memLogFilename);
			#else
				mcmcAvPostPtr = getMCMCGRaverage(
							*containerPostPtr,
							box,
							np,
							ns,
							prefixPosterior,
							spec,
							autotype,
							logFilename,
							memLogFilename);
			#endif
			
			logUsage(memLogFilename, "xxhistograms");
		
		}
		logUsage(memLogFilename, "xxcontainer");
			
		string meanAndSdPostFilename;
		{
			ostringstream oss;
			oss << prefix << "_Posterior" << nrefstats <<
				"JointMeanAndSd_" + spec + ".txt";
			meanAndSdPostFilename = oss.str();
		}

		printAndOutputMeansAndSDs(meansPost,sdsPost, meanAndSdPostFilename);
		
		processStandardisedPredictivePosteriorAverage(
			*mcmcAvPostPtr,
			sdsPost,
			meansPost,
			np,
			ns,
			prefix,
			spec,
			refset_ptr_post,
			slicePtsPostVec,
			prec,
			confirm);
		
		
		
	}	
	
	logUsage(memLogFilename, "At end of ahabc routine");
	
	clock_t finalEndTime = clock();
	
	double finalTiming = ((static_cast<double>(finalEndTime - startTime)) / CLOCKS_PER_SEC);
	
	{
		std::ostringstream stm;
		stm << "\n\nTiming\t" << finalTiming;
		outputAHABCLogFileString(logFilename, stm.str());
	}
	
	cout << "\n\n\nEnd of AHABC:  time = " << finalTiming << "\n\n\n" <<endl;
	
	
}

void doFixedParamsHeterozygosity(
				const std::string& prefix,
				const std::string& spec,
				size_t n_obs,
				size_t n_reps,
				size_t nloci,
				size_t nsites,
				boost::shared_ptr< hudson_ms::PopulationStructure > pop,
				const std::vector < double >& trueParams,
				const std::vector < double >& testParams,
				int seed,
				int simSeed,
				ahabc::AutoMCMCFactory::AUTOTYPE autotype,
				int maxSplits)
{
	
	// only want to use heterozygosity stat
	std::vector < size_t > reqStats(1, 2);
	
	
	doFixedParams(
				prefix,
				spec,
				n_obs,
				n_reps,
				nloci,
				nsites,
				pop,
				trueParams,
				testParams,
				seed,
				simSeed,
				reqStats,
				autotype,
				maxSplits);
	
}


void doFixedParamsHeterozygosityAndSegStats(
				const std::string& prefix,
				const std::string& spec,
				size_t n_obs,
				size_t n_reps,
				size_t nloci,
				size_t nsites,
				boost::shared_ptr< hudson_ms::PopulationStructure > pop,
				const std::vector < double >& trueParams,
				const std::vector < double >& testParams,
				int seed,
				int simSeed,
				ahabc::AutoMCMCFactory::AUTOTYPE autotype,
				int maxSplits)
{
	
	// dummy
	std::vector < size_t > reqStats;
	
	doFixedParams(
				prefix,
				spec,
				n_obs,
				n_reps,
				nloci,
				nsites,
				pop,
				trueParams,
				testParams,
				seed,
				simSeed,
				reqStats,
				autotype,
				maxSplits);
	
}


void doFixedParams(
				const std::string& prefix,
				const std::string& spec,
				size_t n_obs,
				size_t n_reps,
				size_t nloci,
				size_t nsites,
				boost::shared_ptr< hudson_ms::PopulationStructure > pop,
				const std::vector < double >& trueParams,
				const std::vector < double >& testParams,
				int seed,
				int simSeed,
				const std::vector < size_t >& reqStats,
				ahabc::AutoMCMCFactory::AUTOTYPE autotype,
				int maxSplits)
{
	string logFilename;
	{
		std::ostringstream stm;
		stm << prefix << "FixedLog_" << spec << ".txt";
		logFilename = stm.str();
	}
	outputAHABCLogFileStart(logFilename);

	string memLogFilename;
	{
		std::ostringstream stm;
		stm << prefix << "MemLog_" << spec << ".log";
		memLogFilename = stm.str();
	}
	
	
	
	logUsage(memLogFilename, "Start doFixedParams", false); // new mem log file
	
	clock_t startTime = clock();
		
	/*stuff I know I'll need */
	boost::shared_ptr < subpavings::PiecewiseConstantFunction >
			mcmcAvPtr;
	std::vector < double > sds;
	std::vector < double > means;
	int ns = 0;
	
	
	int trueseed = seed + 1;
	
	{
		std::ostringstream stm;
		stm << "\n\nSeed for 'true' values (ref stats)\t" << trueseed;
		outputAHABCLogFileString(logFilename, stm.str());
	}
	
	boost::shared_ptr < mct::SummaryStatisticSet > refset_ptr
					= makeAHABCRefSet(n_obs, nloci, nsites,
								pop, trueParams, reqStats,
								trueseed);
	
	string refstatsFilename;
	{
		std::ostringstream stm;
		stm << prefix << "RefStats_" << spec << ".txt";
		refstatsFilename = stm.str();
	}
	refset_ptr->outputToFile(refstatsFilename, false); // overwrite, not append
	
	
		
	cout << "\nFirst 10 ref stats are = " << endl;
	for (size_t i = 0; (i < 10) && (i < refset_ptr->size()); ++i) {
		cout << refset_ptr->at(i)->toString() << endl;
		
	}
	
	logUsage(memLogFilename, "Made refstats");
	
	{
		std::ostringstream stm;
		stm << "\n\nSeed for simulations\t" << simSeed;
		outputAHABCLogFileString(logFilename, stm.str());
	}
	
	boost::shared_ptr < mct::SummaryStatisticSet > ssset_ptr
					= makeAHABCRefSet(n_reps, nloci, nsites,
								pop, testParams, reqStats,
								simSeed);
	
	string statsFilename;
	{
		std::ostringstream stm;
		stm << prefix << "SumStats_" << spec << ".txt";
		statsFilename = stm.str();
	}
	ssset_ptr->outputToFile(statsFilename, false); // overwrite, not append
	
	
	logUsage(memLogFilename, "Made stats");
	
	clock_t endTime = clock();
	{
		double timing = ((static_cast<double>(endTime - startTime)) / CLOCKS_PER_SEC);
		std::ostringstream stm;
		stm << "\nTiming to make stats\t" << timing;
		outputAHABCLogFileString(logFilename, stm.str());
	}
	
	std::vector < std::vector < double > > slicePtsVec;
		
	/* put a lot of the work in it's own scope so that anything
	 * we do not need goes out of scope when we've done the heavy lifting... */
	{
		//fill a container with the parameter particle and summary statistic values
		boost::shared_ptr < std::vector < std::vector < double > > >
			containerPtr 
					= ssset_ptr->getAllValuesStandardisedSumStatLayout(
						means,
						sds);

		ns = containerPtr->at(0).size();
		int np = 0;
		slicePtsVec = getStandardisedSlicePoints(refset_ptr, np, means, sds);

		std::string containerFilename = prefix + "StandardisedData_" + spec + ".txt"; 
		mct_utilities::outputVectorVectorToFile(*containerPtr,
					containerFilename, false); // overwrite
	
		cxsc::real padding = 0.000005;
		ivector box = subpavings::makeBox(*containerPtr, padding);
		
		cxsc::real minVol = 0.0;
		if (maxSplits > 0) {
			cxsc::real boxVol = realVolume(box);
			
			minVol = boxVol/cxsc::power(2.0, maxSplits);
			minVol = pred(minVol);
			
		}
		mcmcAvPtr = getMCMCGRaverage(
				*containerPtr,
				box,
				minVol,
				ns,
				prefix,
				spec,
				autotype,
				logFilename,
				memLogFilename);
		
		
		endTime = clock();
		{
			double timing = ((static_cast<double>(endTime - startTime)) / CLOCKS_PER_SEC);
			std::ostringstream stm;
			stm << "\n\nTiming to make summary stat density estimate rmrp\t" << timing;
			outputAHABCLogFileString(logFilename, stm.str());
		}
		startTime = clock();
		
		logUsage(memLogFilename, "xxhistograms");
	
	}
	logUsage(memLogFilename, "xxcontainer");
		
	string meanAndSdFilename = prefix + "_MeanAndSd_" + spec + ".txt";

	printAndOutputMeansAndSDs(means,sds, meanAndSdFilename);
			
	int prec = 5;
	bool confirm = true;
	
	processStandardisedSummaryStatsAverage(	*mcmcAvPtr,
											n_reps,
											sds,
											means,
											ns,
											prefix,
											spec,
											slicePtsVec,
											prec,
											confirm);
	
	
	
	
	
	logUsage(memLogFilename, "At end of ahabc routine");
	
	clock_t finalEndTime = clock();
	
	double finalTiming = ((static_cast<double>(finalEndTime - startTime)) / CLOCKS_PER_SEC);
	
	{
		std::ostringstream stm;
		stm << "\n\nTiming\t" << finalTiming;
		outputAHABCLogFileString(logFilename, stm.str());
	}
	
	cout << "\n\n\nEnd of AHABC:  time = " << finalTiming << "\n\n\n" <<endl;
	
	
}
