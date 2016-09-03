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
\brief Tools for doing various mcmc and ahabc routines

*/

#include "testing_toolsMCMCGR.hpp"
#include "testing_toolsGR.hpp"
#include "testing_tools.hpp" // for logfilestring
#include "test_routines.hpp"

//#include "MCMCGRAutoNew.hpp"


#include "cxsc.hpp"

#include "config.h" // mct


#include <iostream>
#include <fstream>
#include <iterator>
#include <iomanip>
#include <cassert>


using namespace std;
using namespace subpavings;

// get a histogram
boost::shared_ptr < subpavings::PiecewiseConstantFunction > getMCMCGRaverage(
				const std::vector < std::vector < double > >& container,
				const cxsc::ivector& box,
				int np,
				int ns,
				const std::string& prefix,
				const std::string& spec,
				ahabc::AutoMCMCFactory::AUTOTYPE autotype,
				const std::string& logFilename,
				const std::string& memLogFilename)
{
	real minVol(0.0);
	return getMCMCGRaverage(
				container,
				box,
				minVol,
				np,
				ns,
				prefix,
				spec,
				autotype,
				logFilename,
				memLogFilename);
}





// get a histogram
boost::shared_ptr < subpavings::PiecewiseConstantFunction > getMCMCGRaverage(
				const std::vector < std::vector < double > >& container,
				const cxsc::ivector& box,
				real minVol,
				int ns,
				const std::string& prefix,
				const std::string& spec,
				ahabc::AutoMCMCFactory::AUTOTYPE autotype,
				const std::string& logFilename,
				const std::string& memLogFilename)
{
	int np = 0;
	
	return getMCMCGRaverage(
				container,
				box,
				minVol,
				np,
				ns,
				prefix,
				spec,
				autotype,
				logFilename,
				memLogFilename);
}

// get a histogram
boost::shared_ptr < subpavings::PiecewiseConstantFunction > getMCMCGRaverage(
				const std::vector < std::vector < double > >& container,
				const cxsc::ivector& box,
				real minVol,
				int np,
				int ns,
				const std::string& prefix,
				const std::string& spec,
				ahabc::AutoMCMCFactory::AUTOTYPE autotype,
				const std::string& logFilename,
				const std::string& memLogFilename)
{
	try {
		
		size_t d = np + ns;
		size_t n = container.size();
		
		/* some guesses for leaves in carver queue */
		size_t maxLeavesCarving = static_cast<size_t>(
					std::pow(d, 2.0)*
					std::pow(static_cast<double>(n), 0.5) );
		if (maxLeavesCarving > n) maxLeavesCarving = n;
		
		
		/* some guesses for max points in a node to stop posterior queue */
		size_t critSEB = static_cast<size_t>(std::log(static_cast<double>(n)));
		
		/* some guesses for maximum leaves we'll let SEB queue go to */
		size_t maxLeavesSEB = n / critSEB; // integer division
		
		size_t keep = 3; // keep three histograms to start chains with
			
		std::string postFileName = prefix + "LogPosteriorsForStarts.txt";
		std::string checkPostFileNameBase = prefix + "CheckLogPosteriors";
		
		string scalarsFileName = prefix;
		{
			ostringstream oss;
			oss << "_" << spec << ".log";
			scalarsFileName += oss.str();
			
		}
		
		size_t minPoints = 1;
		
		unsigned int maxLoops = 5000000;
		if (d > 2) maxLoops *= d;
		unsigned int thinout = 1000;
		int samplesNeeded = 100;
		cxsc::real tolerance(0.1);
		
		size_t samplingInterval = 100;
		size_t maxStatesForCalcs = 1000000;
		
		if (d > 3) {
			tolerance = 0.2;
			samplingInterval = 1000;
			
			#if(0)
			if (n > 50000) {
				maxStatesForCalcs = 500000;
			}
			#endif
			
		}
		
		
		
		int rep = 1; // repeats in the MCMC
		int chooseStarts = 5; // number of starts to pick from (min 2)
		if (d > static_cast<size_t>(chooseStarts)) chooseStarts = d;
				
		double percentSpread = 0.95;
		
		if (d > 3) percentSpread = 0.98;
		if (d > 3 && n > 50000) percentSpread = 0.99;
		
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
			stm << "\ntolerance\t" << tolerance;
			stm << "\tsamplingInterval\t" << samplingInterval;
			stm << "\tmaxStatesForCalcs\t" << maxStatesForCalcs;
			stm << "\nminVol\t" << minVol;
			outputAHABCLogFileString(logFilename, stm.str());
			
			cout << stm.str() << endl;	
		}
		
		
		
		boost::shared_ptr < subpavings::PiecewiseConstantFunction >
			mcmcAvPtr = testCarverMCMCGR(
										box,
										container,
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
										rep,
										autotype,
										chooseStarts,
										memLogFilename);
		
		return mcmcAvPtr;
		
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
