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
\brief Tools for AABC runs with KDE
*/

#include "testing_toolsAABC_KDE.hpp"

#include "testing_toolsAABC.hpp"

#include "testing_tools.hpp"

#include "testing_toolsKDE.hpp"


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

using namespace std;
using namespace mct;
using namespace aabc;
using namespace subpavings;



void doKDE_AABCSimsMultLociPolytableHeterozygosity(const string& prefix,
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
				cxsc::real maxMeasure)
{
	boost::function < 
					boost::shared_ptr < aabc::ParameterParticleValidator > 
					(boost::shared_ptr < mct::PRNGen > rr_sampler,
					const boost::shared_ptr < const mct::MultiLociPolyTable > refPtr,
					size_t n_reps_per_particle) 
								> validatorMaker = &makeValidatorHeterozygosity;
	
	doKDE_AABCSimsMultLociPolytable(prefix,
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
				maxMeasure);
}

void doKDE_AABCSimsMultLociPolytableHeterozygosityAndSegSites(const string& prefix,
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
				cxsc::real maxMeasure)
{
	boost::function < 
					boost::shared_ptr < aabc::ParameterParticleValidator > 
					(boost::shared_ptr < mct::PRNGen > rr_sampler,
					const boost::shared_ptr < const mct::MultiLociPolyTable > refPtr,
					size_t n_reps_per_particle) 
								> validatorMaker = &makeValidatorHeterozygosityAndSegSites;
	
	doKDE_AABCSimsMultLociPolytable(prefix,
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
				maxMeasure);
}



/* do sims and just make a kde */
void doKDE_AABCSimsMultLociPolytable(const string& prefix,
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
				cxsc::real maxMeasure)
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
		
		cout << "\nStarting aabc simulations for set " << (i+1) << endl;
				
		boost::shared_ptr < const mct::MultiLociPolyTable > 
			refPtr = refObjSetPtr->at(i);
		
		// add this into the large-loci table we are building up	
		for (size_t j = 0; j < nloci; ++j) {
			allRefObjsPtr->add(refPtr->at(j), nsites);
		}
		
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
	} // end of loop through separate obs
	
	// paramsSetPtrsVec has the paramPtr for each observation
	
	assert(paramsSetPtrsVec.size() == n_obs);
	
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
	
	boost::shared_ptr < ParameterParticleSet > allRefStatsParamsPtr
				= simulatorPtr->getCurrentParticles();		
	/*get the summary stats corresponding to the particles*/
	boost::shared_ptr < SummaryStatisticSet >
		allRefStatsStatsPtr = simulatorPtr->getCurrentSummaryStatistics();
	
	{//make an ahabc out of them
	
		AHABCValues ahabc(allRefStatsParamsPtr,	allRefStatsStatsPtr);
		// use this to get the values and output
		std::vector < std::vector < double > > container;
		ahabc.fillDataContainer(container);
		
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
		for (size_t i = 0; i < n_obs; ++i) {
		
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
	
	//do one kde for each separate run	
	for (size_t i = 0; i < n_obs; ++i) {
	
		cout << "\nStarting kdes for set " << (i+1) << endl;
		
		doKDE(
				*(paramsStatsSetPtrsVec[i]),
				box,
				prefix,
				spec,
				maxMeasure,
				i+1,
				logfilenames[i]);
	}
	
	//now do a kde for the results for all ref stats together
	{
		cout << "\nStarting kde for all ref objs together" << endl;
		doKDE(
				*allRefStatsParamsContainerPtr,
				box,
				prefix,
				spec,
				maxMeasure,
				0, // 'index' for the combined lot 
				logfilenameOverall);
	}
	cout << "\nEnd of all simulations and kdes" << endl;
}


