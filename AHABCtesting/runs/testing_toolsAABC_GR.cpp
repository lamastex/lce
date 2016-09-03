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
\brief Tools for AABC runs with GR histograms
*/

#include "testing_toolsAABC_GR.hpp"

#include "testing_toolsAABC.hpp"

#include "testing_toolsGR.hpp"

#include "testing_tools.hpp"

#include "test_routines.hpp"



#include "testing_toolsKDE.hpp"

#include "piecewise_constant_function.hpp"
#include "spsnode_measure_obj.hpp"

#include "config.h" // mct

#include "ahabc_values.hpp"

#include "descriptive_stats.hpp"
#include "summary_statistic_set.hpp"
#include "summary_statistic.hpp"

//#include "parameter_particle_validator_multi_loci_polytable_heterozygosity.hpp"
//#include "parameter_particle_validator_multi_loci_polytable.hpp"
//#include "multi_loci_simdata_sampler.hpp"
//#include "multi_loci_polytable_sampler.hpp"

#include "parameter_particle_set.hpp"

#include "parameter_prior_uniform.hpp"
#include "parameter_prior_set.hpp"
#include "parameter_prior.hpp"

#include "prng_gsl.hpp"

#include "utilities.hpp"


#include <iostream>
#include <sstream>
#include <ctime>

using namespace std;
using namespace mct;
using namespace aabc;
using namespace subpavings;


void doJointForSlice(
				const AHABCValues& ahabcValuesPtr,
				boost::shared_ptr < mct::SummaryStatisticSet > refset_ptr,
				const std::string& prefix,
				const std::string& spec,
				ahabc::AutoMCMCFactory::AUTOTYPE autotype);
				
				
void carverStandardisedMCMCGR(
		const cxsc::ivector& box,
		boost::shared_ptr < std::vector < std::vector < double > > >
							containerPtr,
		const std::vector < double >& sds,
		const std::vector < double >& means,
		int setIndex,
		const std::string& prefix,
		const std::string& spec,
		ahabc::AutoMCMCFactory::AUTOTYPE autotype,
		const std::string& logfilename);


void doAABCSimsMultLociPolytableHeterozygosity(const string& prefix,
				const string& spec,
				boost::shared_ptr < PRNGen > r,
				size_t iterations,
				size_t n_obs,
				size_t nloci,
				size_t nsites,
				const std::vector< std::pair<double, double> >& priorLims,
				boost::shared_ptr< hudson_ms::PopulationStructure > pop,
				boost::shared_ptr < mct::SummaryStatisticSet > refSetSumStatsPtr,
				boost::shared_ptr < MultiLociPolyTableSet > refObjSetPtr,
				size_t n_reps_per_particle,
				double baseEpsilonPercentile,
				const std::stack < double >& epsilonMults,
				ahabc::AutoMCMCFactory::AUTOTYPE autotype,
				bool doIndividual)
{
	boost::function < 
					boost::shared_ptr < aabc::ParameterParticleValidator > 
					(boost::shared_ptr < mct::PRNGen > rr_sampler,
					const boost::shared_ptr < const mct::MultiLociPolyTable > refPtr,
					size_t n_reps_per_particle) 
								> validatorMaker = &makeValidatorHeterozygosity;
	
	// only want to use heterozygosity stat
	std::vector < size_t > reqStats(1, 2);
		
	doAABCSimsMultLociPolytable(prefix,
				spec,
				r,
				iterations,
				n_obs,
				nloci,
				nsites,
				priorLims,
				pop,
				refSetSumStatsPtr,
				refObjSetPtr,
				n_reps_per_particle,
				baseEpsilonPercentile,
				epsilonMults,
				validatorMaker,
				reqStats,
				autotype,
				doIndividual);
}

void doAABCSimsMultLociPolytableHeterozygosityAndSegSites(const string& prefix,
				const string& spec,
				boost::shared_ptr < mct::PRNGen > r,
				size_t iterations,
				size_t n_obs,
				size_t nloci,
				size_t nsites,
				const std::vector< std::pair<double, double> >& priorLims,
				boost::shared_ptr< hudson_ms::PopulationStructure > pop,
				boost::shared_ptr < mct::SummaryStatisticSet > refSetSumStatsPtr,
				boost::shared_ptr < mct::MultiLociPolyTableSet > refObjSetPtr,
				size_t n_reps_per_particle,
				double baseEpsilonPercentile,
				const std::stack < double >& epsilonMults,
				ahabc::AutoMCMCFactory::AUTOTYPE autotype,
				bool doIndividual)
{
	
	// want to use both heterozygosity and seg sites stats
	std::vector < size_t > reqStats;
	
	boost::function < 
					boost::shared_ptr < aabc::ParameterParticleValidator > 
					(boost::shared_ptr < mct::PRNGen > rr_sampler,
					const boost::shared_ptr < const mct::MultiLociPolyTable > refPtr,
					size_t n_reps_per_particle) 
								> validatorMaker = &makeValidatorHeterozygosityAndSegSites;
	
	doAABCSimsMultLociPolytable(prefix,
				spec,
				r,
				iterations,
				n_obs,
				nloci,
				nsites,
				priorLims,
				pop,
				refSetSumStatsPtr,
				refObjSetPtr,
				n_reps_per_particle,
				baseEpsilonPercentile,
				epsilonMults,
				validatorMaker,
				reqStats,
				autotype,
				doIndividual);
}

void doAABCSimsMultLociPolytable(const string& prefix,
				const string& spec,
				boost::shared_ptr < PRNGen > r,
				size_t iterations,
				size_t n_obs,
				size_t nloci,
				size_t nsites,
				const std::vector< std::pair<double, double> >& priorLims,
				boost::shared_ptr< hudson_ms::PopulationStructure > pop,
				boost::shared_ptr < mct::SummaryStatisticSet > refSetSumStatsPtr,
				boost::shared_ptr < MultiLociPolyTableSet > refObjSetPtr,
				size_t n_reps_per_particle,
				double baseEpsilonPercentile,
				const std::stack < double >& epsilonMults,
				boost::function < 
					boost::shared_ptr < aabc::ParameterParticleValidator > 
					(boost::shared_ptr < mct::PRNGen > rr_sampler,
					const boost::shared_ptr < const mct::MultiLociPolyTable > refPtr,
					size_t n_reps_per_particle) 
								> validatorMaker,
				const std::vector < size_t >& reqStats,
				ahabc::AutoMCMCFactory::AUTOTYPE autotype,
				bool doIndividual)
{
	
	//vector of nloci copies of nsites
	std::vector < size_t > nloci_nsites(nloci, nsites);
	
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
	
	
	boost::shared_ptr < mct::MultiLociPolyTable > allRefObjsPtr
		(new mct::MultiLociPolyTable(pop) );
	
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
		
		cout << "\nStarting loci " << (i+1) << endl;
				
		boost::shared_ptr < const mct::MultiLociPolyTable > 
			refPtr = refObjSetPtr->at(i);
		
		// add this into the large-loci table we are building up	
		for (size_t j = 0; j < nloci; ++j) {
			allRefObjsPtr->add(refPtr->at(j), nsites);
		}
		
		if(doIndividual) { // if we are doing the individual stuff
			cout << "\nStarting aabc simulations for loci " << (i+1) << endl;
		
			string logfilename;
			{
				ostringstream oss;
				oss << prefix << "_SimulationLog_" << spec << "_" << (i+1) << ".txt";
				logfilename = oss.str();		
			}
			logfilenames.push_back(logfilename);
		
			cout << "simulating" << endl;
			boost::shared_ptr < PopMCabcSimulatorAdaptive > simulatorPtr 
							 = doSimulationMultLociPolytable(
								r_loops,
								n_reps_per_particle,	
								priorLims,
								refPtr,
								iterations,
								baseEpsilonPercentile,
								epsilonMults,
								validatorMaker,
								logfilename);
			
			/*The simulation gives stats as specified by the type of validator
			 * from the validator maker */
		
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
							 = doSimulationMultLociPolytable(
									r_loops,
									n_reps_per_particle,	
									priorLims,
									allRefObjsPtr,
									iterations,
									baseEpsilonPercentile,
									epsilonMults,
									validatorMaker,
									logfilenameOverall);
	
	/*The simulation gives stats as specified by the type of validator
		 * from the validator maker */
		
		
	// a 'set' (just one) of summary statistics for the combined ref data
	boost::shared_ptr < mct::SummaryStatisticSet > 
		allRefObjsRefSet_ptr ( new SummaryStatisticSet() );
	
	{
		boost::shared_ptr < mct::SummaryStatisticSet > 
		sset_ptr_tmp ( new SummaryStatisticSet() );
		
		// and get a set of summary statistics from the values
		boost::shared_ptr < const mct::SummaryStatistic > ss_ptr_tmp =
						allRefObjsPtr->getSummaryStatistic();

		sset_ptr_tmp->add(ss_ptr_tmp);
		// add the req stats from this into the set we actually want
		if (!reqStats.empty()) {
			allRefObjsRefSet_ptr->add(sset_ptr_tmp->makeSelectedSummaryStatisticSet(reqStats));
		}
		else {
			allRefObjsRefSet_ptr->add(sset_ptr_tmp);
		}
	}
	/* allRefObjsRefSet_ptr should now contain the summary
	 *  stats for the combined data - we will use it later */

	boost::shared_ptr < ParameterParticleSet > allRefStatsParamsPtr
				= simulatorPtr->getCurrentParticles();		
	/*get the summary stats corresponding to the particles*/
	boost::shared_ptr < SummaryStatisticSet >
		allRefStatsStatsPtr = simulatorPtr->getCurrentSummaryStatistics();
	
	/* an ahabc ptr with params and ref stats using all the pseudo-observed values
	 * as one reference object, ie assumin linked (non-independent) loci */
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
		allRefStatsParamsStatsSet->getAllValuesStandardisedSumStatLayout(
				allRefStatsParamsMeans, allRefStatsParamsSds);
	
	{ // output the standardised data to file
		bool ppAppend = false;
		ostringstream oss;
		oss << prefix << "_StandardisedParamsData_" << spec << "_" << 0 << ".txt";
		string fname = oss.str();
		mct_utilities::outputVectorVectorToFile(
				*allRefStatsParamsContainerPtr, fname, ppAppend);
	}
	//print and output means
	{ 
		ostringstream oss;
		oss << prefix << "_MeansAndSds_" << spec << "_" << 0 << ".txt";
		printAndOutputMeansAndSDs(allRefStatsParamsMeans, 
								allRefStatsParamsSds, oss.str());
	}
	
	std::vector < std::vector < double > > sdsVec;
	std::vector < std::vector < double > > meansVec;

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
					
			// for overall means and sds
			std::vector < double > sds;
			std::vector < double > means;
		
			boost::shared_ptr < std::vector < std::vector < double > > >
				paramsContainerPtr =
					paramsStatsSet->getAllValuesStandardisedSumStatLayout(
							means, sds);
			sdsVec.push_back(sds);
			meansVec.push_back(means);
			paramsStatsSetPtrsVec.push_back(paramsContainerPtr);
			
			{ // output the standardised data to file
				bool ppAppend = false;
				ostringstream oss;
				oss << prefix << "_StandardisedParamsData_" << spec << "_" << (i+1) << ".txt";
				string fname = oss.str();
				mct_utilities::outputVectorVectorToFile(*paramsContainerPtr, fname, 
													ppAppend);
			}
			
			//print and output means
			{ 
				ostringstream oss;
				oss << prefix << "_MeansAndSds_" << spec << "_" << (i+1) << ".txt";
				printAndOutputMeansAndSDs(means, sds, oss.str());
			}
			
			boxContainer.insert(boxContainer.end(), 
								paramsContainerPtr->begin(),
								paramsContainerPtr->end() ) ;
		}
		
		//make the box using all the standardised params values
		cxsc::real padding = 0.000005;
		box = subpavings::makeBox(boxContainer, padding);
	}
	
	//do one histogram for each separate run	
	for (size_t i = 0; i < paramsStatsSetPtrsVec.size(); ++i) {
	
		cout << "\nStarting histograms for set " << (i+1) << endl;
		carverStandardisedMCMCGR(
				box,
				paramsStatsSetPtrsVec[i],
				sdsVec[i],
				meansVec[i],
				i+1,
				prefix,
				spec,
				autotype,
				logfilenames[i]);
	}
	
	//now do a histogram for the results for all ref stats together
	{
		cout << "\nStarting histograms for all ref objs together" << endl;
		carverStandardisedMCMCGR(
				box,
				allRefStatsParamsContainerPtr,
				allRefStatsParamsSds,
				allRefStatsParamsMeans,
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




// test Carver with GR convergence
void carverStandardisedMCMCGR(
				const ivector& box,
				boost::shared_ptr < std::vector < std::vector < double > > >
							containerPtr,
				const std::vector < double >& sds,
				const std::vector < double >& means,
				int setIndex,
				const std::string& prefix,
				const std::string& spec,
				ahabc::AutoMCMCFactory::AUTOTYPE autotype,
				const std::string& logfilename)
{
	try {
		
		assert(!containerPtr->empty());
		
		size_t n = containerPtr->size();
		size_t d = means.size();
		
		
		/* some guesses for leaves in carver queue */
		double mult = 2.0;
		if ( d > 1) mult = std::pow(d, 2.0);
		size_t maxLeavesCarving = static_cast<size_t>(
					mult *
					std::pow(static_cast<double>(n), 0.5) );
		if (maxLeavesCarving > n) maxLeavesCarving = n;
		
		#if(0)
			size_t maxLeavesCarving = 1000; 
		#endif
		
		/* some guesses for max points in a node to stop posterior queue */
		size_t critSEB = static_cast<size_t>(std::log(static_cast<double>(n)));
		
		#if(0)
			size_t critSEB = 5;  // max points in a node to stop posterior queue
		#endif
			
		/* some guesses for maximum leaves we'll let SEB queue go to */
		size_t maxLeavesSEB = n / critSEB; // integer division
		
		#if(0)
			size_t maxLeavesSEB = 2000000; // maximum leaves we'll let SEB queue go to
									// even if critSEB not met
		#endif
			
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
			outputAABCLogFileString(logfilename, stm.str());
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
		
		processStandardisedAverage(
			*mcmcAvPtr,
			setAppend,
			setIndex,
			sds,
			means,
			prefix,
			spec,
			prec,
			confirm);
		
		clock_t endTime = clock();
		double timing = ((static_cast<double>(endTime - startTime)) / CLOCKS_PER_SEC);
	
		{
			std::ostringstream stm;
			stm << "\n\nhistogram timing\t" << timing;
			outputAABCLogFileString(logfilename, stm.str());
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
	std::vector < double > sds;
	std::vector < double > means;
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
		
		
		means = ahabcValues.getMeans();
		
		sds = ahabcValues.getSDs();
	
		np = ahabcValues.getNumberParameters();
		ns = ahabcValues.getNumberSummaryStatistics();
		
		sliceDims = ahabcValues.getSummaryStatisticsDims();
		
		ahabcValues.fillStandardisedDataContainer(*containerPtr);
		
		std::string containerFilename = prefixJoint + "Data_" + spec + ".txt"; 
		mct_utilities::outputVectorVectorToFile(*containerPtr,
					containerFilename, false); // overwrite
	
		slicePtsVec = getStandardisedSlicePoints(refset_ptr, np, means, sds);
	
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
		
	string meanAndSdFilename = prefixJoint + "_MeanAndSd_" + spec + ".txt";

	printAndOutputMeansAndSDs(means,sds, meanAndSdFilename);
			
	int prec = 5;
	bool confirm = true;
	size_t posteriorCaptureInterval = 0;
	std::vector < double > trueParams; // dummy
	std::vector < std::vector < double > > otherParams; // dummy
	size_t n_reps = ahabcValues.size();
	
	std::vector < 
		boost::shared_ptr < PiecewiseConstantFunction > 
	> posteriorPtrsVec = processStandardisedAverage(
		*mcmcAvPtr,
		sds,
		means,
		np,
		ns,
		n_reps,
		prefixJoint,
		spec,
		sliceDims,
		refset_ptr,
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

