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



#include "config.h" // mct

#include "testing_toolsMCMCGR.hpp"

#include "piecewise_constant_function.hpp"
#include "adaptivehistogram.hpp"

#include "carver_seb.hpp"

#include <vector>
#include <string>
#include <iostream>
#include <iterator>
#include <cassert>
#include <stdexcept>

using namespace std;
using namespace subpavings;


PiecewiseConstantFunction testCarverMCMC(
			const std::vector < std::vector < double > >& container,
			size_t maxLeavesCarving,
			size_t critSEB,	
			size_t maxLeavesSEB,	
			size_t minPoints,
			const std::string& postFileName,
			const std::string& checkPostFileNameBase,
			const std::string& bestoneFilename,
			const std::string& sliceFilename,
			unsigned int loops,
			unsigned int burnin,
			unsigned int thinout,
			LOGGING_LEVEL logging,				
			const std::string& mcmcAvFilename,
			const std::string& mcmcSliceFilename,
			const std::vector < int >& sliceDims,
			const std::vector < double >& slicePts)
{
	AdaptiveHistogram adh;
	bool successfulInsertion = 
			adh.insertRvectorsFromVectorOfVecDbls(container);
	
	if (!successfulInsertion) 
					throw std::runtime_error("Failed to insert data");
					
	cout << "histogram has " << adh.getRootCounter() << " values" << endl;

	// container for hist pointers
	std::vector< subpavings::AdaptiveHistogram* > hists;
	
	
	SPSNodeMeasureVolMassMinus compCarving(adh.getRootCounter());
	
	AdaptiveHistogram::PrioritySplitQueueEvaluator 
			evaluatorCarving( compCarving, maxLeavesCarving);
	
	SPSNodeMeasureCount compSEB;
	
	AdaptiveHistogram::PrioritySplitQueueEvaluator 
			evaluatorSEB( compSEB, critSEB, maxLeavesSEB);
	
	LogCatalanPrior logPrior; // prior
	int chooseStarts = 10; // try to find 10 launch points
	size_t keepBest = 1; // keep just one best
	bool stopOnMaxPosterior = true; // gives a little latitude
	int precPQ = 5;
	
	CarverSEB::findStartingPoints(adh, hists,
					evaluatorCarving,
					evaluatorSEB,
					logPrior, minPoints,
					chooseStarts, keepBest,
					stopOnMaxPosterior,
					postFileName, checkPostFileNameBase,	
					precPQ);
	
	
	int prec = 10;
	bool confirm = true;
		
	assert(!hists.empty());
	
	// at least one there, and we only wanted one, so the single one is the best
	
	// make it into a pcf
	PiecewiseConstantFunction bestOne(**hists.begin());
	

	// get the best as a pcf
	bestOne.outputToTxtTabs(bestoneFilename, prec, confirm);
	
	{
		PiecewiseConstantFunction slice = bestOne.makeSlice(
											sliceDims, slicePts);

		slice.outputToTxtTabs(sliceFilename, prec, confirm);
	}
	
		
	PiecewiseConstantFunction mcmcAv = makeMCMCAv(**hists.begin(), 
							loops,burnin, thinout,
							 minPoints, logging);
	{	
		// get the average as a pcf
		mcmcAv.outputToTxtTabs(mcmcAvFilename, prec, confirm);
		

		PiecewiseConstantFunction slice = mcmcAv.makeSlice(
											sliceDims, slicePts);

		// output slice
		slice.outputToTxtTabs(mcmcSliceFilename, prec, confirm);
	}	
	try {
		vector<AdaptiveHistogram*>::iterator ait;
		for (ait = hists.begin(); ait < hists.end(); ++ait) {
			if (NULL != *ait) delete (*ait);
			*ait = NULL;
		}
	}
	catch(...) {} //catch and swallow
	
	return mcmcAv;
	

}
