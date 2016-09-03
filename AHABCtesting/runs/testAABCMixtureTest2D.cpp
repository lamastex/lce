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
\brief Do a number of AABC runs using MixtureTest2D.
* 
Aim is to compare 'product likelihood' of separate results with 
posterior for all obs together.
* 
*/


//#include "testing_toolsAABC.hpp"

#include "testing_toolsGR.hpp"

#include "testing_tools.hpp"

#include "test_routines.hpp"





#include "automcmc_factory.hpp"

#include "pop_mc_abc_simulator_adaptive.hpp"
#include "parameter_particle_validator_generic_data_object.hpp"
#include "parameter_particle_validator.hpp"

#include "piecewise_constant_function.hpp"
#include "spsnode_measure_obj.hpp"
#include "sptools.hpp"

#include "config.h" // mct

#include "ahabc_values.hpp"

#include "mixture_test_2d_sampler.hpp"
#include "mixture_test_2d.hpp"

#include "generic_data_object_set.hpp"
#include "data_object.hpp"

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



#include <boost/shared_ptr.hpp>


#include <iostream>
#include <vector>
#include <stack>
#include <utility>

using namespace std;
using namespace mct;
using namespace aabc;
using namespace subpavings;

void runAABC_Mixture();

// make some reference statistics
boost::shared_ptr < GenericDataObjectSet > makeRefSet(
		size_t n_obs,
		const std::vector < double >& trueParams,
		int trueseed);

boost::shared_ptr < PopMCabcSimulatorAdaptive > doSimulationMixture(
				boost::shared_ptr < PRNGen > r_loops,
				size_t n_reps_per_particle,	
				const std::vector< std::pair<double, double> >& priorLims,
				boost::shared_ptr < const DataObject > refPtr,
				size_t iterations,
				double baseEpsilonPercentile,
				std::stack < double > epsilonMults, // pass by value
				const std::string& logfilename);

void doAABCMixture(const string& prefix,
				const string& spec,
				boost::shared_ptr < PRNGen > r,
				size_t iterations,
				size_t n_obs,
				const std::vector< std::pair<double, double> >& priorLims,
				boost::shared_ptr < mct::SummaryStatisticSet > refSetSumStatsPtr,
				boost::shared_ptr < GenericDataObjectSet > refObjSetPtr,
				size_t n_reps_per_particle,
				double baseEpsilonPercentile,
				const std::stack < double >& epsilonMults,
				ahabc::AutoMCMCFactory::AUTOTYPE autotype,
				bool doIndividual);
				
// test Carver with GR convergence
void carverMCMCGR(
				const ivector& box,
				boost::shared_ptr < std::vector < std::vector < double > > >
							containerPtr,
				int setIndex,
				const std::string& prefix,
				const std::string& spec,
				ahabc::AutoMCMCFactory::AUTOTYPE autotype,
				const std::string& logfilename);

void doJointForSlice(
				const AHABCValues& ahabcValues,
				boost::shared_ptr < mct::SummaryStatisticSet > refset_ptr,
				const std::string& prefix,
				const std::string& spec,
				ahabc::AutoMCMCFactory::AUTOTYPE autotype);
			
int main()
{
	try {
		runAABC_Mixture();
	
	
		return 0;
	}
	catch (std::exception& e) {
		std::cout << "Exception " << e.what() << std::endl; 
		throw;
	}
	catch (...) {
		std::cout << "Unknown error " << std::endl; 
		throw;
	}

	
}

	

void runAABC_Mixture()
{
	
		
	int myPid = getpid();
	std::cout << "MixtureTest2D aabc" << std::endl;
	std::cout << "This process id is " << myPid << std::endl;
	
	std::string baseOutputDir("../figures/");
	
	std::string thisDir("AABCMixture4DNewAvOldPriorSupport50Obs");
	
	std::string descriptionString("AABCMixture4DGR");
	
	int seed = 2345;
	boost::shared_ptr < PRNGen > r ( new PRNGenGSL(seed) );
	
	/* How many iterations we want in each generation - ie how many members of the 
	 * population we are aiming for in each generation */
	size_t iterations = 10000;
	
	bool doIndividual = false; // whether to do individual posteriors
	size_t n_obs = 50; // number of 'observed values' to have
	
	// make the output dir and get back path
	if (!doIndividual) thisDir += "Conly";
	std::string path = mct_utilities::makeDir(baseOutputDir, thisDir);
	string prefix = path + "/" + descriptionString;
	
			
	//true mu1 and mu2
	double mu1 = 0.0;
	double mu2 = 0.0;
	
	std::vector < double > trueParams (1, mu1);
	trueParams.push_back(mu2);
	
	cout << "\n\nn_reps = " << iterations << ", n_obs = " << n_obs << endl;
	
	cout << "true mu1 " << mu1 << ", mu2 " << mu2 << endl;
	
	
	std::vector< std::pair<double, double> > priorLims;
	priorLims.push_back( pair<double, double>
					//for mu1
					(-10.0, 10.0));
			// prior limits for mu2
			priorLims.push_back( pair<double, double>
				(-10.0, 10.0));
	
	
	
	
	int trueseed = seed + 1;
	
	boost::shared_ptr < GenericDataObjectSet > refObjSetPtr
							= makeRefSet(
										n_obs,
										trueParams,
										trueseed);
	
	//get the summary statistics for all the observed values
	boost::shared_ptr < SummaryStatisticSet > refSetSumStatsPtr = 
						refObjSetPtr->getSummaryStatisticSet();
	
	// stuff for the aabc simulation
	// number of repetitions in each comparison sample generated for each parameter particle
	size_t n_reps_per_particle = 1; 
	
	// the percentile distance on which to base calculations of epsilon
	double baseEpsilonPercentile = 0.1;
	
	/* In each generation, we get tighter and tighter in our 
	 * idea of 'close enough'.  These multipliers are applied 
	 * progressively, ie 1.0 * the value we get from the base 
	 * epsilon criteria for the first loop, then 0.75 * that, then
	 * 0.75 * that, then 0.5 * that */
	double mydoubles[] = {0.5, 0.75, 0.75, 1.0}; 
	std::stack < double > epsilonMults (std::deque <double>(
			mydoubles, mydoubles + sizeof(mydoubles) / sizeof(double) ));
	
	ahabc::AutoMCMCFactory::AUTOTYPE autotype = ahabc::AutoMCMCFactory::ALLTHREE;
	
	string spec;
	{
		std::ostringstream stm;
		#if(0)
			stm << ahabc::AutoMCMCFactory::getDiagnosticName(autotype);
		#endif
		stm << "_no_" << n_obs;
		stm.precision(mct_utilities::fitPrecision(mu1));
		stm <<   fixed;
		stm << "_mu1_" << mu1;
		stm.precision(mct_utilities::fitPrecision(mu2));
		stm << "_mu2_" << mu2;
		stm << "_ni_" << iterations;
		
		spec = stm.str();
	}
	
	doAABCMixture(prefix,
				spec,
				r,
				iterations,
				n_obs,
				priorLims,
				refSetSumStatsPtr,
				refObjSetPtr,
				n_reps_per_particle,
				baseEpsilonPercentile,
				epsilonMults,
				autotype,
				doIndividual);
	
}



// make some reference statistics
boost::shared_ptr < GenericDataObjectSet > makeRefSet(
		size_t n_obs,
		const std::vector < double >& trueParams,
		int trueseed)
{
	
	boost::shared_ptr < mct::PRNGen > rr_sampler ( new mct::PRNGenGSL(trueseed) );
	
	
	// make a MixtureTest sampler with its prng
	boost::shared_ptr < MixtureTest2DSampler > 
					sampler_ptr ( new MixtureTest2DSampler(rr_sampler) );
	
	boost::shared_ptr < aabc::ParameterParticle > 
			p( new aabc::ParameterParticle(trueParams) );

	cout << "sampling 'true' values " << ( p->toString() )<< endl;
	

	// get a set of values from the sampler, using the particles
	boost::shared_ptr < GenericDataObjectSet > sampledValues_ptr = 
			sampler_ptr->sample(n_obs, p);
	
	
	return sampledValues_ptr;
	
}


void doAABCMixture(const string& prefix,
				const string& spec,
				boost::shared_ptr < PRNGen > r,
				size_t iterations,
				size_t n_obs,
				const std::vector< std::pair<double, double> >& priorLims,
				boost::shared_ptr < mct::SummaryStatisticSet > refSetSumStatsPtr,
				boost::shared_ptr < GenericDataObjectSet > refObjSetPtr,
				size_t n_reps_per_particle,
				double baseEpsilonPercentile,
				const std::stack < double >& epsilonMults,
				ahabc::AutoMCMCFactory::AUTOTYPE autotype,
				bool doIndividual)
{
	
	
	/* use this prng to spawn others */
	boost::shared_ptr < PRNGen > r_loops = r->spawnAnother();

	
	// output the stats values for observed values to file
	{
		bool statsAppend = false;
		ostringstream oss;
		oss << prefix << "_Stats_" << spec << "_obs" << ".txt";
		refSetSumStatsPtr->outputToFile(oss.str(), statsAppend);
	}
	
	std::string refSetSumStatsString = refSetSumStatsPtr->plainString();
	cout << "\nThe 'observed' values are\n" << refSetSumStatsString << endl;
	
	boost::shared_ptr < GenericDataObjectSet > allRefObjsPtr
		(new GenericDataObjectSet() );
	
	// a container for all the final parameter particle sets
	std::vector < boost::shared_ptr < ParameterParticleSet > > paramsSetPtrsVec;
	
	
	
	// stuff for recording filenames for params and stats for each obs
	bool appendParamsAndStats = false;
	string paramsAndStatsFilenames;
	{
		ostringstream oss;
		oss << prefix << "_ParamsAndStats_" << spec << ".txt";
		paramsAndStatsFilenames = oss.str();
	}	
		
	std::vector <std::string > logfilenames;
	
	for (size_t i = 0; i < n_obs; ++i) {
		
		cout << "\nStarting obs " << (i+1) << endl;
				
		boost::shared_ptr < const DataObject > 
			refPtr = refObjSetPtr->at(i);
		
		allRefObjsPtr->add(refPtr);
		
		if(doIndividual) { // if we are doing the individual stuff
			cout << "\nStarting aabc simulations for obs " << (i+1) << endl;
		
			string logfilename;
			{
				ostringstream oss;
				oss << prefix << "_SimulationLog_" << spec << "_" << (i+1) << ".txt";
				logfilename = oss.str();		
			}
			logfilenames.push_back(logfilename);
		
			cout << "simulating" << endl;
			boost::shared_ptr < PopMCabcSimulatorAdaptive > simulatorPtr 
							 = doSimulationMixture(
								r_loops,
								n_reps_per_particle,	
								priorLims,
								refPtr,
								iterations,
								baseEpsilonPercentile,
								epsilonMults,
								logfilename);
								
			// get the parameters
			boost::shared_ptr < ParameterParticleSet > paramsPtr
						= simulatorPtr->getCurrentParticles();
					
			cout << "simulation ended: paramsPtr.size() = " << paramsPtr->size() << endl;
			paramsSetPtrsVec.push_back(paramsPtr);
			
			/*get the summary stats corresponding to the particles */
			boost::shared_ptr < SummaryStatisticSet > statsPtr
					= simulatorPtr->getCurrentSummaryStatistics();
				
			//make an ahabc out of them
			AHABCValues ahabc(paramsPtr, statsPtr);
			// use this to get the values and output
			std::vector < std::vector < double > > container;
			ahabc.fillDataContainer(container);
			{
				bool ppAppend = false;
				ostringstream oss;
				oss << prefix << "_ParamsAndStats_" << spec << "_" << (i+1) << ".txt";
				string fname = oss.str();
				mct_utilities::outputVectorVectorToFile(container, fname, 
														ppAppend);
				outputFileLine(paramsAndStatsFilenames, stripPath(fname),
														appendParamsAndStats);
				if (!appendParamsAndStats) appendParamsAndStats = true;
			}
		}
	} // end of loop through separate obs
	
	// paramsSetPtrsVec has the paramPtr for each observation IF ANY
	
	assert((!doIndividual && paramsSetPtrsVec.empty()) || 
		(doIndividual && (paramsSetPtrsVec.size() == n_obs)));
	
	string logfilenameOverall;
	{
		ostringstream oss;
		oss << prefix << "_SimulationLog_" << spec << "_" << 0 << ".txt";
		logfilenameOverall = oss.str();		
	}
	
	cout << "\nDo simulation for all ref objs together" << endl;
	// do the simulation for all the ref stats together
	boost::shared_ptr < PopMCabcSimulatorAdaptive > simulatorPtr
							 = doSimulationMixture(
									r_loops,
									n_reps_per_particle*n_obs, // need to take account of n_obs in average	
									priorLims,
									allRefObjsPtr,
									iterations,
									baseEpsilonPercentile,
									epsilonMults,
									logfilenameOverall);
	
	// a 'set' (just one) of summary statistics for the combined ref data
	boost::shared_ptr < mct::SummaryStatisticSet > 
		allRefObjsRefSet_ptr ( new SummaryStatisticSet() );
	
	allRefObjsRefSet_ptr->add(allRefObjsPtr->getSummaryStatistic());
		
	/* allRefObjsRefSet_ptr should now contain the summary
	 *  stats for the combined data - we will use it later */

	boost::shared_ptr < ParameterParticleSet > allRefStatsParamsPtr
				= simulatorPtr->getCurrentParticles();		
	/*get the summary stats corresponding to the particles*/
	boost::shared_ptr < SummaryStatisticSet >
		allRefStatsStatsPtr = simulatorPtr->getCurrentSummaryStatistics();
	
	/* an ahabc ptr with params and ref stats using all the pseudo-observed values
	 * as one reference object*/
	AHABCValues ahabcAllRefStats(allRefStatsParamsPtr,	allRefStatsStatsPtr);
	
	{		
		// use this to get the values and output
		std::vector < std::vector < double > > container;
		ahabcAllRefStats.fillDataContainer(container);
		
		{
			bool ppAppend = false;
			ostringstream oss;
			oss << prefix << "_ParamsAndStats_" << spec << "_" << 0 << ".txt";
			string fname = oss.str();
			mct_utilities::outputVectorVectorToFile(container, fname, 
												ppAppend);
			outputFileLine(paramsAndStatsFilenames, stripPath(fname), 
												appendParamsAndStats);
			if (!appendParamsAndStats) appendParamsAndStats = true;
		}
			
	}
	
	/* and turn these param particles into a summary stats set form */
	boost::shared_ptr < mct::SummaryStatisticSet > allRefStatsParamsStatsSet
				= allRefStatsParamsPtr->makeSummaryStatisticSet();
	
	//get a container of the params values for all ref stats together
	// for overall means and sds
	std::vector < double > allRefStatsParamsSds;
	std::vector < double > allRefStatsParamsMeans;
	boost::shared_ptr < std::vector < std::vector < double > > >
		allRefStatsParamsContainerPtr =
		allRefStatsParamsStatsSet->getAllValuesSumStatLayout();
	
	{ // output the data to file
		bool ppAppend = false;
		ostringstream oss;
		oss << prefix << "_ParamsData_" << spec << "_" << 0 << ".txt";
		string fname = oss.str();
		mct_utilities::outputVectorVectorToFile(
				*allRefStatsParamsContainerPtr, fname, ppAppend);
	}
		
	
	std::vector < boost::shared_ptr < std::vector < std::vector < double > > > >
		paramsStatsSetPtrsVec;
	
	ivector box; 
	// we are aiming to make the box and collect some other useful things...
	{
		std::vector < std::vector < double > > 
						boxContainer = *allRefStatsParamsContainerPtr;

		/* for each separate set of values, do the same thing
		 * and scoop all that data into one container */
		for (size_t i = 0; i <paramsSetPtrsVec.size(); ++i) {
		
			boost::shared_ptr < mct::SummaryStatisticSet > paramsStatsSet
					= paramsSetPtrsVec[i]->makeSummaryStatisticSet();
					
			boost::shared_ptr < std::vector < std::vector < double > > >
				paramsContainerPtr =
					paramsStatsSet->getAllValuesSumStatLayout();
			
			paramsStatsSetPtrsVec.push_back(paramsContainerPtr);
			
			{ // output the standardised data to file
				bool ppAppend = false;
				ostringstream oss;
				oss << prefix << "_ParamsData_" << spec << "_" << (i+1) << ".txt";
				string fname = oss.str();
				mct_utilities::outputVectorVectorToFile(*paramsContainerPtr, fname, 
													ppAppend);
			}
			
			
			
			boxContainer.insert(boxContainer.end(), 
								paramsContainerPtr->begin(),
								paramsContainerPtr->end() ) ;
		}
		
		//make the box using all the params values
		cxsc::real padding = 0.000005;
		box = subpavings::makeBox(boxContainer, padding);
	}
	
	//do one histogram for each separate run	
	for (size_t i = 0; i < paramsStatsSetPtrsVec.size(); ++i) {
	
		cout << "\nStarting histograms for set " << (i+1) << endl;
		carverMCMCGR(
				box,
				paramsStatsSetPtrsVec[i],
				i+1,
				prefix,
				spec,
				autotype,
				logfilenames[i]);
	}
	
	//now do a histogram for the results for all ref stats together
	{
		cout << "\nStarting histograms for all ref objs together" << endl;
		carverMCMCGR(
				box,
				allRefStatsParamsContainerPtr,
				0, // 'index' for the combined lot 
				prefix,
				spec,
				autotype,
				logfilenameOverall);
	}
	
	/* now do a histogram of the joint using params and summary stats
	 * for the combined data, in order to take a slice of this
	 * across the combined ref stats */
	doJointForSlice(
				ahabcAllRefStats,
				allRefObjsRefSet_ptr,
				prefix,
				spec,
				autotype);
	
	
	cout << "\nEnd of all simulations and histograms" << endl;
}

boost::shared_ptr < PopMCabcSimulatorAdaptive > doSimulationMixture(
				boost::shared_ptr < PRNGen > r_loops,
				size_t n_reps_per_particle,	
				const std::vector< std::pair<double, double> >& priorLims,
				boost::shared_ptr < const DataObject > refPtr,
				size_t iterations,
				double baseEpsilonPercentile,
				std::stack < double > epsilonMults, // pass by value
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
	
	boost::shared_ptr < GenericDataObjectSampler > sampler ( new MixtureTest2DSampler(rr_sampler) );
	
	boost::shared_ptr < ParameterParticleValidator > 
			validator( 
			new ParameterParticleValidatorGenericDataObject(refPtr,
													sampler,
													n_reps_per_particle) 
					);
					
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

// test Carver with GR convergence
void carverMCMCGR(
				const ivector& box,
				boost::shared_ptr < std::vector < std::vector < double > > >
							containerPtr,
				int setIndex,
				const std::string& prefix,
				const std::string& spec,
				ahabc::AutoMCMCFactory::AUTOTYPE autotype,
				const std::string& logfilename)
{
	try {
		
		assert(!containerPtr->empty());
		
		size_t n = containerPtr->size();
		size_t d = containerPtr->at(0).size();
		
		
		/* some guesses for leaves in carver queue */
		double mult = 2.0;
		if ( d > 1) mult = std::pow(d, 2.0);
		size_t maxLeavesCarving = static_cast<size_t>(
					mult *
					std::pow(static_cast<double>(n), 0.5) );
		if (maxLeavesCarving > n) maxLeavesCarving = n;
		
		
		/* some guesses for max points in a node to stop posterior queue */
		size_t critSEB = static_cast<size_t>(std::log(static_cast<double>(n)));
		
			
		/* some guesses for maximum leaves we'll let SEB queue go to */
		size_t maxLeavesSEB = n / critSEB; // integer division
		
			
		size_t minPoints = 1;
		size_t keep = 3; // keep three histograms to start chains with
			
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
		double percentSpread = 0.95;
		
		cxsc::real minVol(0.0);
		
		int chooseStarts = 5; // number of starts to pick from (min 2)
		if (static_cast<int>(d) > chooseStarts) chooseStarts = d;
		
		// log values used for carver to file
		{
			std::ostringstream stm;
			stm << "\n\nStarts chosen using\n\tmaxLeavesCarving\t" << maxLeavesCarving;
			stm << "\n\tcritSEB\t" << critSEB;
			stm << "\n\tmaxLeavesSEB\t" << maxLeavesSEB;
			stm << "\n\tpercentSpread\t" << percentSpread;
			stm << "\n\tmaxLoops\t" << maxLoops;
			stm << "\tthinout\t" << thinout;
			stm << "\tsamplesNeeded\t" << samplesNeeded;
			outputAHABCLogFileString(logfilename, stm.str());
		}
		
		clock_t startTime = clock();
		
		boost::shared_ptr < subpavings::PiecewiseConstantFunction >
		  mcmcAvPtr = testCarverMCMCGR(
			box,
			*containerPtr,
			maxLeavesCarving,
			critSEB,	
			maxLeavesSEB,
			percentSpread,	
			minPoints,
			minVol,
			keep,
			postFileName,
			checkPostFileNameBase,
			scalarsFileName,
			maxLoops,
			thinout,
			samplesNeeded,
			tolerance,
			setIndex,
			autotype,
			chooseStarts);
		
		int prec = 5;
		bool confirm = true;
		
		// append false if setIndex == 1, true if setIndex != 1
		bool setAppend = (setIndex != 1);
		
		processAverage(
			*mcmcAvPtr,
			setAppend,
			setIndex,
			prefix,
			spec,
			prec,
			confirm);
		
		clock_t endTime = clock();
		double timing = ((static_cast<double>(endTime - startTime)) / CLOCKS_PER_SEC);
	
		{
			std::ostringstream stm;
			stm << "\n\nhistogram timing\t" << timing;
			outputAHABCLogFileString(logfilename, stm.str());
		}
		
	}
	catch (std::exception& e) {
	
		std::cout << "In carverStandardisedMCMCGR - std::exception::\n" + std::string(e.what()) << std::endl;
		throw;
	}
	catch (...) {
	
		std::cout << "In carverStandardisedMCMCGR - Unknown error\n" << std::endl;
		throw;
	}
}

void doJointForSlice(
				const AHABCValues& ahabcValues,
				boost::shared_ptr < mct::SummaryStatisticSet > refset_ptr,
				const std::string& prefix,
				const std::string& spec,
				ahabc::AutoMCMCFactory::AUTOTYPE autotype)
{
	std::string prefixJoint = prefix + "Joint";
	
	string logFilename;
	{
		std::ostringstream stm;
		stm << prefixJoint << "AHABCLog_" << spec << ".txt";
		logFilename = stm.str();
	}
	outputAHABCLogFileStart(logFilename);

	string memLogFilename;
	{
		std::ostringstream stm;
		stm << prefixJoint << "MemLog_" << spec << ".log";
		memLogFilename = stm.str();
	}
	
	clock_t startTime = clock();
		
	/*stuff I know I'll need */
	boost::shared_ptr < subpavings::PiecewiseConstantFunction >
			mcmcAvPtr;
	int np;
	int ns;
	std::vector < std::vector < double > > slicePtsVec;
	std::vector < int > sliceDims;
	
	string refstatsFilename;
	{
		std::ostringstream stm;
		stm << prefixJoint << "RefStats_" << spec << ".txt";
		refstatsFilename = stm.str();
	}
	refset_ptr->outputToFile(refstatsFilename, false); // overwrite, not append
	
	
		
	cout << "\nref stats for joint are "
		<< refset_ptr->at(0)->toString() << endl;
		
		
		
	/* put a lot of the work in it's own scope so that anything
	 * we do not need goes out of scope when we've done the heavy lifting... */
	{
		//fill a container with the parameter particle and summary statistic values
		boost::shared_ptr < std::vector < std::vector < double > > >
			containerPtr (new std::vector < std::vector < double > >() );
		
		
		np = ahabcValues.getNumberParameters();
		ns = ahabcValues.getNumberSummaryStatistics();
		
		sliceDims = ahabcValues.getSummaryStatisticsDims();
		
		ahabcValues.fillDataContainer(*containerPtr);
		
		std::string containerFilename = prefixJoint + "Data_" + spec + ".txt"; 
		mct_utilities::outputVectorVectorToFile(*containerPtr,
					containerFilename, false); // overwrite
	
		slicePtsVec = getSlicePoints(refset_ptr);
	
		cxsc::real padding = 0.000005;
		ivector box = getBox(*containerPtr, slicePtsVec, sliceDims, padding);
		
		mcmcAvPtr = getMCMCGRaverage(
					*containerPtr,
					box,
					np,
					ns,
					prefixJoint,
					spec,
					autotype,
					logFilename,
					memLogFilename);
		
		
	
	}
		
			
	int prec = 5;
	bool confirm = true;
	size_t posteriorCaptureInterval = 0;
	std::vector < double > trueParams; // dummy
	std::vector < std::vector < double > > otherParams; // dummy
	size_t n_reps = ahabcValues.size();
	
	std::vector < 
		boost::shared_ptr < PiecewiseConstantFunction > 
	> posteriorPtrsVec = processAverage(
		*mcmcAvPtr,
		np,
		ns,
		n_reps,
		prefixJoint,
		spec,
		sliceDims,
		slicePtsVec,
		trueParams,
		otherParams,
		prec,
		confirm,
		posteriorCaptureInterval);
	
	
	
	clock_t endTime = clock();
	
	double timing = ((static_cast<double>(endTime - startTime)) / CLOCKS_PER_SEC);
	
	{
		std::ostringstream stm;
		stm << "\n\nTiming\t" << timing;
		outputAHABCLogFileString(logFilename, stm.str());
	}
	
	cout << "\n\n\nEnd of joint AHABC:  time = " << timing << "\n" <<endl;
	
	
}


