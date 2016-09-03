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
\brief Tools for doing various AHABC routines using KDE

*/

#include "testing_toolsAHABC_KDE.hpp"
#include "testing_toolsKDE.hpp"
#include "testing_toolsAHABC.hpp"
#include "testing_toolsSimdata.hpp"
#include "testing_tools.hpp"
#include "log_memory.hpp"

#include "test_routines.hpp"

#include "config.h" // mct

#include "parameter_prior_uniform.hpp"

#include "parameter_prior.hpp"

#include "multi_loci_simdata_sampler.hpp"
#include "multi_loci_polytable_sampler.hpp"

#include "parameter_particle_set.hpp"
#include "descriptive_stats.hpp"
#include "prng.hpp"
#include "prng_gsl.hpp"

#include "utilities.hpp"

//#include "adaptivehistogram.hpp"



#include <real.hpp>


#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <cassert>

#include <boost/lambda/bind.hpp>
#include <boost/lambda/lambda.hpp>

#include <gsl/gsl_statistics_double.h>

//#define DEBUG_GRIDDING
//#define DEBUG_GRIDDING_OUTPUT

using namespace std;
using namespace mct;
using namespace aabc;
using namespace subpavings;



void doAHABC_KDEHeterozygosity(
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
				cxsc::real maxMeasure,
				bool reducePriorLims)
{
	//dummy
	std::vector < std::vector < double > > otherParams;
	doAHABC_KDEHeterozygosity(
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
				maxMeasure,
				reducePriorLims);
}

void doAHABC_KDEHeterozygosity(
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
				cxsc::real maxMeasure,
				bool reducePriorLims)
{
	
	// only want to use heterozygosity stat
	std::vector < size_t > reqStats(1, 2);
	
	doAHABC_KDE(
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
				maxMeasure,
				reqStats,
				reducePriorLims);
	
}
void doAHABC_KDEHeterozygosityAndSegStats(
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
				cxsc::real maxMeasure,
				bool reducePriorLims)
{
	//dummy
	std::vector < std::vector < double > > otherParams;
	doAHABC_KDEHeterozygosityAndSegStats(
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
				maxMeasure,
				reducePriorLims);
}

void doAHABC_KDEHeterozygosityAndSegStats(
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
				cxsc::real maxMeasure,
				bool reducePriorLims)
{
	
	// dummy
	std::vector < size_t > reqStats;
	
	doAHABC_KDE(
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
				maxMeasure,
				reqStats,
				reducePriorLims);
	
}


void doAHABC_KDE(
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
				cxsc::real maxMeasure,
				const std::vector < size_t >& reqStats,
				bool reducePriorLims)
{
	string logFilename;
	{
		std::ostringstream stm;
		stm << prefix << "AHABC_Log_" << spec << ".txt";
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
	std::vector < double > sds;
	std::vector < double > means;
	int np;
	int ns;
	std::vector < std::vector < double > > slicePtsVec;
	std::vector < int > sliceDims;
	
	int trueseed = seed + 1;
	
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
	
	
	std::vector< std::vector<double> >::iterator eit = slicePtsVec.begin();
		advance(eit, 10);
		
	cout << "\nFirst 10 ref stats are = " << endl;
	for (size_t i = 0; (i < 10) && (i < refset_ptr->size()); ++i) {
		cout << refset_ptr->at(i)->toString() << endl;
		
	}
	
	#if(0)
	
		for (std::vector< std::vector<double> >::iterator it = slicePtsVec.begin();
						(it < eit) && (it < slicePtsVec.end());
						++it) {
			ostream_iterator<double> out_it (cout,"\t");
			copy ( it->begin(), it->end(), out_it );
			cout << endl;
		}
		cout << endl;

		
		if (refset_ptr->size() <= 10) {
			std::string refSetSumStatsString = refset_ptr->plainString();
			cout << "\nThe 'observed' values are\n" << refSetSumStatsString << endl;
		}
	#endif
	logUsage(memLogFilename, "Made refstats");
	
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
	
	}
		
		
	/* put a lot of the work in it's own scope so that anything
	 * we do not need goes out of scope when we've done the heavy lifting... */
	{
		//fill a container with the parameter particle and summary statistic values
		boost::shared_ptr < std::vector < std::vector < double > > >
			containerPtr (new std::vector < std::vector < double > >() );
		
		{
			// pass in the newPriorLims
			boost::shared_ptr < const AHABCValues  > ahabcValuesPtr 
						= makeValuesPtr(n_reps, nloci, nsites, 
										pop, n_reps_per_particle,
										newPriorLims, reqStats,
										seed,
										logFilename);

			assert (ahabcValuesPtr->size() == n_reps);
			
			logUsage(memLogFilename, "Made AHABC values");
				
			means = ahabcValuesPtr->getMeans();
			
			sds = ahabcValuesPtr->getSDs();
		
			np = ahabcValuesPtr->getNumberParameters();
			ns = ahabcValuesPtr->getNumberSummaryStatistics();
			
			sliceDims = ahabcValuesPtr->getSummaryStatisticsDims();
			
			ahabcValuesPtr->fillStandardisedDataContainer(*containerPtr);
			
			logUsage(memLogFilename, "Filled container");
	
		}
		cout << "ahabc values should go out of scope" << endl;
		logUsage(memLogFilename, "xxAHABC values");
			
		std::string containerFilename = prefix + "Data_" + spec + ".txt"; 
		mct_utilities::outputVectorVectorToFile(*containerPtr,
					containerFilename, false); // overwrite
	
		slicePtsVec = getStandardisedSlicePoints(refset_ptr, np, means, sds);
	
		cxsc::real padding = 0.000005;
		ivector box = getBox(*containerPtr, slicePtsVec, sliceDims, padding);
		
		kdeEstPtr = getKDEapprox(
				*containerPtr,
				box,
				prefix,
				spec,
				maxMeasure,
				logFilename);
		
		logUsage(memLogFilename, "xxKDE");
	
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
	
	processStandardisedAverage(
		*kdeEstPtr,
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
	
	logUsage(memLogFilename, "At end");
	
	clock_t endTime = clock();
	
	double timing = ((static_cast<double>(endTime - startTime)) / CLOCKS_PER_SEC);
	
	{
		std::ostringstream stm;
		stm << "\n\nTiming\t" << timing;
		outputAHABCLogFileString(logFilename, stm.str());
	}
	
	cout << "\n\n\nEnd of AHABC:  time = " << timing << "\n\n\n" <<endl;
	
	
}

