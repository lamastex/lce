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
\brief Try to get some sort of rmrpkde for data using Beaumont etc toy model in 2D

*/


//#include "testing_toolsAHABC.hpp"
//#include "testing_toolsGR.hpp"
#include "testing_toolsKDE.hpp"

//#include "testing_toolsMRS.hpp"
#include "testing_tools.hpp"
#include "log_memory.hpp"


#include "automcmc_factory.hpp"


#include "piecewise_constant_function.hpp"
#include "spsnode_measure_obj.hpp"
#include "sptools.hpp"

#include "config.h" // mct

#include "ahabc_values.hpp"

#include "mixture_test_sampler.hpp"
#include "mixture_test.hpp"

#include "generic_data_object_set.hpp"
#include "parameter_particle_set.hpp"
#include "summary_statistic_set.hpp"

#include "parameter_prior_uniform.hpp"
#include "parameter_prior.hpp"
#include "parameter_prior_set.hpp"

#include "descriptive_stats.hpp"

#include "prng.hpp"
#include "prng_gsl.hpp"


#include "utilities.hpp"

#include "toolz.hpp" // for realVolume




#include <boost/shared_ptr.hpp>


#include <real.hpp>

#include <iostream>
#include <vector>
#include <iterator>
#include <cassert>
#include <utility>
#include <string>

using namespace std;
using namespace mct;
using namespace aabc;
using namespace subpavings;


#define REDUCE_PRIOR_LIMS

// make some reference statistics
boost::shared_ptr < mct::SummaryStatisticSet > makeAHABCRefSet(
		size_t n_obs,
		const std::vector < double >& trueParams,
		int trueseed);

std::vector< std::pair<double, double> > getPriorLimsMixTest(
		boost::shared_ptr < mct::SummaryStatisticSet > refset_ptr,
		size_t n_reps_test,	
		const std::vector< std::pair<double, double> >& priorLims,
		int seed,
		size_t gridSegments,
		double percentClosest);

// make some values
boost::shared_ptr < const AHABCValues  > makeValuesPtr(
		size_t n_reps,
		size_t n_reps_per_particle,	
		std::vector< std::pair<double, double> >& priorLims,
		int seed,
		const std::string& logFilename);

boost::shared_ptr < const AHABCValues  > makeValuesPtr(
		size_t n_reps,
		size_t n_reps_per_particle,	
		const subpavings::PiecewiseConstantFunction& pcfPosterior,
		int seed,
		const std::string& logFilename);

void doAHABCMixtureKDE(
				const std::string& prefix,
				const std::string& spec,
				size_t n_obs,
				size_t n_reps,
				const std::vector< std::pair<double, double> >& priorLims,
				size_t n_reps_per_particle,
				const std::vector < double >& trueParams,
				const std::vector < std::vector < double > >& otherParams,
				int seed,
				int simSeed,
				cxsc::real maxMeasure,
				bool reducePriorLims);
				

int main()
{
	int myPid = getpid();
	std::cout << "basic ahabc test for mixture testwith small prior range" << std::endl;
	std::cout << "This process id is " << myPid << std::endl;
	
	std::string baseOutputDir("../figures/");
	
	std::string thisDir("MixtureTest2DKDETest");
	#ifdef REDUCE_PRIOR_LIMS
		thisDir += "Rdcd";
	#endif
	
	std::string descriptionString("Mixture2DKDE");
	
	// make the output dir and get back path
	std::string path = mct_utilities::makeDir(baseOutputDir, thisDir);
	
	string prefix = path + "/" + descriptionString;
	
	
	size_t n_reps = 2000; // number of joint observations to generate
	
	size_t n_obs = 50; // number of 'observed values' to have
	
	cxsc::real maxMeasure(0.00001);
	
	int seed = 2345;
	
	#if(1)
		int simSeed = seed;
	#else
		int simSeed  = seed + 10;
	#endif
	
	//true mu1 and mu2
	double mu1 = 0.0;
	
	std::vector < double > trueParams (1, mu1);
	
	cout << "\n\nn_reps = " << n_reps << ", n_obs = " << n_obs << endl;
	
	cout << "true mu " << mu1 << endl;
	
	
	string spec;
	{
		std::ostringstream stm;
		stm << "_no_" << n_obs << "_nr_" << n_reps;
		stm.precision(mct_utilities::fitPrecision(mu1));
		stm <<   fixed;
		stm << "_mu_" << mu1;
		
		spec = stm.str();
	}
	
	
	
	#ifdef REDUCE_PRIOR_LIMS
		std::vector< std::pair<double, double> > priorLims;
				
			priorLims.push_back( pair<double, double>
					//for mu1
					(-10.0, 10.0));
			
		
		bool reducePriorLims = true;
	#else
		std::vector< std::pair<double, double> > priorLims;
		priorLims.push_back( pair<double, double>
					//for mu1
					(-10.0, 10.0));
			
		
		bool reducePriorLims = false;
	#endif
	
	// a container for some other 'interesting' sets of parameters
	std::vector < std::vector < double > > otherParams;
	{
		
			double ptl = priorLims.front().first;
			double ptt = priorLims.front().second;
			double t1 = 0.25;
			if (!(t1 > ptl)) t1 = ptl + (mu1 - ptl)*0.8;
			double t2 = 0.5;
			if (!(t2 < ptt)) t2 = mu1 + (ptt - mu1)*0.2;
			double t3 = 1.0;
			if (!(t3 > ptl)) t3 = ptl + (mu1 - ptl)*0.9;
			double t4 = 1.5;
			if (!(t4 < ptt)) t4 = mu1 + (ptt - mu1)*0.1;
			
			
		{
			std::vector < double > otherP(1, t1);
			otherParams.push_back(otherP);
		}
		
		{
			std::vector < double > otherP(1, t2);
			otherParams.push_back(otherP);
		}
		
		{
			std::vector < double > otherP(1, t3);
			otherParams.push_back(otherP);
		}
		
		{
			std::vector < double > otherP(1, t4);
			otherParams.push_back(otherP);
		}
		
		
	}

	
	
	// number of repetitions in each comparison sample generated for each parameter particle
	size_t n_reps_per_particle = 1; 
	
	doAHABCMixtureKDE(
				prefix,
				spec,
				n_obs,
				n_reps,
				priorLims,
				n_reps_per_particle,
				trueParams,
				otherParams,
				seed,
				simSeed,
				maxMeasure,
				reducePriorLims);

	
	return 0;
	
}







void doAHABCMixtureKDE(
				const std::string& prefix,
				const std::string& spec,
				size_t n_obs,
				size_t n_reps,
				const std::vector< std::pair<double, double> >& priorLims,
				size_t n_reps_per_particle,
				const std::vector < double >& trueParams,
				const std::vector < std::vector < double > >& otherParams,
				int seed,
				int simSeed,
				cxsc::real maxMeasure,
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
	
	logUsage(memLogFilename, "Start doAHABC_KDE", false); // new mem log file
	
	clock_t startTime = clock();
		
	/*stuff I know I'll need */
	boost::shared_ptr < subpavings::PiecewiseConstantFunction >
			kdeEstPtr;
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
					= makeAHABCRefSet(n_obs, 
								trueParams,
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
		
		newPriorLims = getPriorLimsMixTest(refset_ptr,
									n_reps_test,	
									priorLims,
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
						= makeValuesPtr(n_reps, 
										n_reps_per_particle,
										newPriorLims,
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
				
			np = ahabcValuesPtr->getNumberParameters();
			ns = ahabcValuesPtr->getNumberSummaryStatistics();
			
			sliceDims = ahabcValuesPtr->getSummaryStatisticsDims();
			
			/* get and output the unstandardised data */
			ahabcValuesPtr->fillDataContainer(*containerPtr);
			std::string containerFilenameTmp = prefix + "UnstandardisedData_" + spec + ".txt"; 
			mct_utilities::outputVectorVectorToFile(*containerPtr,
						containerFilenameTmp, false); // overwrite
			
			
			endTime = clock();
			{
				double timing = ((static_cast<double>(endTime - startTime)) / CLOCKS_PER_SEC);
				std::ostringstream stm;
				stm << "\n\nTiming to get data\t" << timing;
				outputAHABCLogFileString(logFilename, stm.str());
			}
			startTime = clock();
			
			logUsage(memLogFilename, "Filled container");
	
		}
		cout << "ahabc values should go out of scope" << endl;
		logUsage(memLogFilename, "xxAHABC values");
			
		
		
		
		slicePtsVec = getSlicePoints(refset_ptr);
	
		cxsc::real padding = 0.000005;
		ivector box = getBox(*containerPtr, slicePtsVec, sliceDims, padding);
		
		cout << "\nGetting the joint distn\n" << endl;
		
		kdeEstPtr = getKDEapprox(
				*containerPtr,
				box,
				prefix,
				spec,
				maxMeasure,
				logFilename);
		
		logUsage(memLogFilename, "xxKDE");
		
		endTime = clock();
		{
			double timing = ((static_cast<double>(endTime - startTime)) / CLOCKS_PER_SEC);
			std::ostringstream stm;
			stm << "\n\nTiming to make joint rmrp\t" << timing;
			outputAHABCLogFileString(logFilename, stm.str());
		}
		startTime = clock();
		
		
	}
	logUsage(memLogFilename, "xxcontainer");
		
	int prec = 5;
	bool confirm = true;
	size_t posteriorCaptureInterval = 10;
	
	std::vector < 
		boost::shared_ptr < PiecewiseConstantFunction > 
	> posteriorPtrsVec = processAverage(
		*kdeEstPtr,
		np,
		ns,
		n_reps,
		prefix,
		spec,
		sliceDims,
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
				kdeEstPostPtr;
		
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
							= makeValuesPtr(n_reps, 
											n_reps_per_particle,
											*pcfPtr,
											seed,
											logFilename);

				assert (ahabcValuesPtr->size() == n_reps);
				
				logUsage(memLogFilename, "Made AHABC values for posterior");
					
				/* get and output the unstandardised data */
				ahabcValuesPtr->fillDataContainer(*containerPostPtr);
				
				std::string containerFilenameTmp = prefix + "PredictiveUnstandardisedData_" + spec + ".txt"; 
				mct_utilities::outputVectorVectorToFile(*containerPostPtr,
							containerFilenameTmp, false); // overwrite
			
				
				logUsage(memLogFilename, "Filled container for posterior");
		
			}
			logUsage(memLogFilename, "xxAHABC values for posterior");
				
			slicePtsPostVec = getSlicePoints(refset_ptr_post);
		
			cxsc::real padding = 0.000005;
			ivector box = getBox(*containerPostPtr, slicePtsPostVec, sliceDims, padding);
			
			string prefixPosterior;
			{
				ostringstream oss;
				oss << prefix << "Predictive" << nrefstats;
				prefixPosterior = oss.str();
			}
			
			
			cout << "\nGetting the posterior joint distn\n" << endl;
			
			kdeEstPostPtr = getKDEapprox(
				*containerPostPtr,
				box,
				prefixPosterior,
				spec,
				maxMeasure,
				logFilename);
					
			
			logUsage(memLogFilename, "xxkde");
		
		}
		logUsage(memLogFilename, "xxcontainer");
			
		processPredictivePosteriorAverage(
			*kdeEstPostPtr,
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




// make some reference statistics
boost::shared_ptr < mct::SummaryStatisticSet > makeAHABCRefSet(
		size_t n_obs,
		const std::vector < double >& trueParams,
		int trueseed)
{
	
	boost::shared_ptr < mct::PRNGen > rr_sampler ( new mct::PRNGenGSL(trueseed) );
	
	
	// make a MixtureTest sampler with its prng
	boost::shared_ptr < MixtureTestSampler > 
					sampler_ptr ( new MixtureTestSampler(rr_sampler) );
	
	boost::shared_ptr < aabc::ParameterParticle > 
			p( new aabc::ParameterParticle(trueParams) );

	cout << "sampling 'true' values " << ( p->toString() )<< endl;
	

	// get a set of values from the sampler, using the particles
	boost::shared_ptr < GenericDataObjectSet > sampledValues_ptr = 
			sampler_ptr->sample(n_obs, p);
	
	// and get a set of summary statistics from the values
	boost::shared_ptr < mct::SummaryStatisticSet > refValues_ptr =
					sampledValues_ptr->getSummaryStatisticSet();

	return refValues_ptr;
	
}


std::vector< std::pair<double, double> > getPriorLimsMixTest(
		boost::shared_ptr < mct::SummaryStatisticSet > refset_ptr,
		size_t n_reps_test,	
		const std::vector< std::pair<double, double> >& priorLims,
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
	
	
	// prng for sampler
	boost::shared_ptr < PRNGen > rr_sampler = rr->spawnAnother();
	
	// make a simdata (infinite sites model) sampler with its prng
	boost::shared_ptr < MixtureTestSampler > 
					sampler_ptr ( new MixtureTestSampler(rr_sampler) );
	
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
			boost::shared_ptr < GenericDataObjectSet > sampledValues_ptr = 
					sampler_ptr->sample(n_reps_test, 
						pPtr);
			
			// and get a set of summary statistics from the values
			boost::shared_ptr < mct::SummaryStatisticSet > sset_ptr =
							sampledValues_ptr->getSummaryStatisticSet();
		
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

boost::shared_ptr < const AHABCValues  > makeValuesPtr(
		size_t n_reps,
		size_t n_reps_per_particle,	
		std::vector< std::pair<double, double> >& priorLims,
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
	
	// prng for sampler
	boost::shared_ptr < PRNGen > rr_sampler = rr->spawnAnother();
	
	// make a simdata (infinite sites model) sampler with its prng
	boost::shared_ptr < MixtureTestSampler > 
					sampler_ptr ( new MixtureTestSampler(rr_sampler) );
	
	
	// get a set of parameter particles from the priors
	boost::shared_ptr< ParameterParticleSet > 
			pset_ptr = priors.generateParameterParticleSet(n_reps);
			
	cout << "Size of particles is " << pset_ptr->size() << endl;

	// get a set of values from the sampler, using the particles
	boost::shared_ptr < GenericDataObjectSet > sampledValues_ptr = 
			sampler_ptr->sample(n_reps_per_particle, pset_ptr);
	
	// and get a set of summary statistics from the values
	boost::shared_ptr < mct::SummaryStatisticSet > sset_ptr =
					sampledValues_ptr->getSummaryStatisticSet();

	cout << "Size of summary stats is " << sset_ptr->size() << endl;

	// make an ahabcvalues object from the particles and summary stats
	boost::shared_ptr < const AHABCValues  > ahabcValuesPtr
		= boost::shared_ptr < const AHABCValues  >(new AHABCValues(pset_ptr, sset_ptr));
	return ahabcValuesPtr;
}


// make some values, specifying required stats
boost::shared_ptr < const AHABCValues  > makeValuesPtr(
		size_t n_reps,
		size_t n_reps_per_particle,	
		const subpavings::PiecewiseConstantFunction& pcfPosterior,
		int seed,
		const std::string& logFilename)
{
	
	boost::shared_ptr < mct::PRNGen > rr ( new mct::PRNGenGSL(seed) );
	
	// prng for sampler
	boost::shared_ptr < PRNGen > rr_sampler = rr->spawnAnother();
	
	boost::shared_ptr < MixtureTestSampler > 
					sampler_ptr ( new MixtureTestSampler(rr_sampler) );
	
			
	cout << "about to do particles and sampling" << endl;
	
	
	// a set for the parameter particles
	boost::shared_ptr< ParameterParticleSet > 
			pset_ptr ( new ParameterParticleSet(n_reps) );
	
	size_t nParams = pcfPosterior.getDimensions();
	
	RVecData paramContainer;
	paramContainer.reserve(n_reps);
	
	/* simulate values in transformed world */
	pcfPosterior.simulateData(paramContainer,
					n_reps, seed);
	
	
	/* take rvectors in paramContainer and convert to parameter
	 * particles in pset_ptr */
	for (RVecData::iterator rit = paramContainer.begin();
		rit < paramContainer.end(); ++rit) {
		rvector rv = *rit;
		vector < double > values(nParams, 0.0);
		// untransform again as we go
		for (size_t i = 0; i < nParams; ++i) {
			values[i] = _double(rv[i+1]);
		}
		boost::shared_ptr< ParameterParticle > ptmp(
			new ParameterParticle(values) );
		pset_ptr->add(ptmp);
	}
	
		
	// get a set of values from the sampler, using the particles
	boost::shared_ptr < GenericDataObjectSet > sampledValues_ptr = 
		sampler_ptr->sample(n_reps_per_particle, pset_ptr);

	// and get a set of summary statistics from the values
	boost::shared_ptr < mct::SummaryStatisticSet > sset_ptr =
				sampledValues_ptr->getSummaryStatisticSet();

	
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


