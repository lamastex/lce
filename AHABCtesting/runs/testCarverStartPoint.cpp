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



std::vector< subpavings::AdaptiveHistogram* >& testCarverStartPoint(
			std::vector< subpavings::AdaptiveHistogram* >& hists,
			const std::vector < std::vector < double > >& container,
			size_t maxLeavesCarving,
			size_t critSEB,	
			size_t maxLeavesSEB,	
			size_t minPoints,
			size_t keepBest,
			int bestReps,
			const std::string& postFileName,
			const std::string& checkPostFileNameBase,
			const std::string& bestoneFilenameBase,
			const std::string& sliceFilenameBase,
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
	
	
	SPSNodeMeasureVolMassMinus compCarving(adh.getRootCounter());
	
	AdaptiveHistogram::PrioritySplitQueueEvaluator 
			evaluatorCarving( compCarving, maxLeavesCarving);
	
	SPSNodeMeasureCount compSEB;
	
	AdaptiveHistogram::PrioritySplitQueueEvaluator 
			evaluatorSEB( compSEB, critSEB, maxLeavesSEB);
	
	LogCatalanPrior logPrior; // prior
	int chooseStarts = 10; // try to find 10 launch points
	bool stopOnMaxPosterior = false; // don't stop when we find the first best
	int precPQ = 5;
	
	unsigned long int seed = 1234;
	
	for (int rep = 0; rep < bestReps; ++rep) {
			
		unsigned long int seedRep = seed + rep;
			
		cout << "\n\nBest rep loop " << rep+1 << "with seed " << seedRep << endl;
		
		//copy the base histogram for this rep
		AdaptiveHistogram adhRep(adh);
		
		// best one from each rep should get kept
		CarverSEB::findStartingPoints(adhRep, hists,
						evaluatorCarving,
						evaluatorSEB,
						logPrior, minPoints,
						chooseStarts, keepBest,
						stopOnMaxPosterior,
						postFileName, checkPostFileNameBase,	
						precPQ,
						seedRep);
	}
	
	
	int prec = 10;
	bool confirm = true;
		
	assert(!hists.empty());
	
	// at least one there
	
	
	size_t nHists = hists.size();
	for (size_t h = 1; h <= nHists; ++h) {
		
		// make it into a pcf
		PiecewiseConstantFunction bestOne(**hists.begin());
		
		string bestoneFilename = bestoneFilenameBase + ".txt";
		if (nHists > 1) { 
			ostringstream oss;
			oss << bestoneFilenameBase << "_" << h << ".txt";
			bestoneFilename = oss.str();
		}
		// output
		bestOne.outputToTxtTabs(bestoneFilename, prec, confirm);
		
		PiecewiseConstantFunction slice = bestOne.makeSlice(
										sliceDims, slicePts);

		
		string sliceFilename = sliceFilenameBase + ".txt";
		if (nHists > 1) { 
			ostringstream oss;
			oss << sliceFilenameBase << "_" << h << ".txt";
			sliceFilename = oss.str();
		}

		slice.outputToTxtTabs(sliceFilename, prec, confirm);
		
	}
	
	return hists;
	
}

std::vector< subpavings::AdaptiveHistogram* >& testCarverStartPoint(
			const ivector& box,
			std::vector< subpavings::AdaptiveHistogram* >& hists,
			const std::vector < std::vector < double > >& container,
			size_t maxLeavesCarving,
			size_t critSEB,	
			size_t maxLeavesSEB,	
			size_t minPoints,
			size_t keepBest,
			int bestReps,
			const std::string& postFileName,
			const std::string& checkPostFileNameBase)
				
{
	AdaptiveHistogram adh(box);
	
	bool successfulInsertion = 
			adh.insertRvectorsFromVectorOfVecDbls(container);
	
	if (!successfulInsertion) 
					throw std::runtime_error("Failed to insert data");
	
	cout << "histogram has " << adh.getRootCounter() << " values" << endl;

		
	// container for hist pointers
	
	
	SPSNodeMeasureVolMassMinus compCarving(adh.getRootCounter());
	
	AdaptiveHistogram::PrioritySplitQueueEvaluator 
			evaluatorCarving( compCarving, maxLeavesCarving);
	
	SPSNodeMeasureCount compSEB;
	
	AdaptiveHistogram::PrioritySplitQueueEvaluator 
			evaluatorSEB( compSEB, critSEB, maxLeavesSEB);
	
	LogCatalanPrior logPrior; // prior
	int chooseStarts = 10; // try to find 10 launch points
	bool stopOnMaxPosterior = false; // don't stop when we find the first best
	int precPQ = 5;
	
	unsigned long int seed = 1234;
	
	for (int rep = 0; rep < bestReps; ++rep) {
			
		unsigned long int seedRep = seed + rep;
			
		cout << "\n\nBest rep loop " << rep+1 << "with seed " << seedRep << endl;
		
		//copy the base histogram for this rep
		AdaptiveHistogram adhRep(adh);
		
		// best one from each rep should get kept
		CarverSEB::findStartingPoints(adhRep, hists,
						evaluatorCarving,
						evaluatorSEB,
						logPrior, minPoints,
						chooseStarts, keepBest,
						stopOnMaxPosterior,
						postFileName, checkPostFileNameBase,	
						precPQ,
						seedRep);
	}
	
	
		
	return hists;
	
}

