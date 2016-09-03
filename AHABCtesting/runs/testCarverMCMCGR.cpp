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
\brief Try to get some sort of histogram for data using Beaumont etc toy model

*/

#include "config.h" // mct

#include "test_routines.hpp"
#include "testing_tools.hpp"
#include "testing_toolsMCMCGR.hpp"
#include "log_memory.hpp"
#include "MCMCGRAutoLeaves.hpp"

#include "piecewise_constant_function.hpp"
#include "adaptivehistogram.hpp"

#include "carver_seb.hpp"

#include <vector>
#include <string>
#include <iostream>
#include <iterator>
#include <cassert>
#include <stdexcept>

#include <boost/shared_ptr.hpp>

using namespace std;
using namespace subpavings;

boost::shared_ptr < subpavings::PiecewiseConstantFunction > testCarverMCMCGR(
			const ivector& box,
			const std::vector < std::vector < double > >& container,
			size_t maxLeavesCarving,
			size_t critSEB,	
			size_t maxLeavesSEB,	
			size_t minPoints,
			size_t keep,
			const std::string& postFileName,
			const std::string& checkPostFileNameBase,
			const std::string& scalarsFileName,
			unsigned int maxLoops,
			unsigned int thinout,
			int samplesNeeded,
			real tolerance,
			int rep,
			ahabc::AutoMCMCFactory::AUTOTYPE autotype,
			const std::string& memLogFilename)
{
	double percentSpread = 0.95;
	int chooseStarts = 10; // try to find 10 launch points
	real minVol(0.0);
	
	return testCarverMCMCGR(
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
	
}

boost::shared_ptr < subpavings::PiecewiseConstantFunction > testCarverMCMCGR(
			const ivector& box,
			const std::vector < std::vector < double > >& container,
			size_t maxLeavesCarving,
			size_t critSEB,	
			size_t maxLeavesSEB,	
			size_t minPoints,
			size_t keep,
			const std::string& postFileName,
			const std::string& checkPostFileNameBase,
			const std::string& scalarsFileName,
			unsigned int maxLoops,
			unsigned int thinout,
			int samplesNeeded,
			real tolerance,
			int rep,
			ahabc::AutoMCMCFactory::AUTOTYPE autotype,
			int chooseStarts,
			const std::string& memLogFilename)
{
	real minVol(0.0);
	double percentSpread = 0.95;
	
	return testCarverMCMCGR(
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
}

boost::shared_ptr < subpavings::PiecewiseConstantFunction > testCarverMCMCGR(
			const ivector& box,
			const std::vector < std::vector < double > >& container,
			size_t maxLeavesCarving,
			size_t critSEB,	
			size_t maxLeavesSEB,	
			double percentSpread,	
			size_t minPoints,
			size_t keep,
			const std::string& postFileName,
			const std::string& checkPostFileNameBase,
			const std::string& scalarsFileName,
			unsigned int maxLoops,
			unsigned int thinout,
			int samplesNeeded,
			real tolerance,
			int rep,
			ahabc::AutoMCMCFactory::AUTOTYPE autotype,
			int chooseStarts,
			const std::string& memLogFilename)
{
	real minVol(0.0);
	
	return testCarverMCMCGR(
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
}


boost::shared_ptr < subpavings::PiecewiseConstantFunction > testCarverMCMCGR(
			const ivector& box,
			const std::vector < std::vector < double > >& container,
			size_t maxLeavesCarving,
			size_t critSEB,	
			size_t maxLeavesSEB,	
			double percentSpread,
			size_t minPoints,
			real minVol,
			size_t keep,
			const std::string& postFileName,
			const std::string& checkPostFileNameBase,
			const std::string& scalarsFileName,
			unsigned int maxLoops,
			unsigned int thinout,
			int samplesNeeded,
			real tolerance,
			int rep,
			ahabc::AutoMCMCFactory::AUTOTYPE autotype,
			int chooseStarts,
			const std::string& memLogFilename)
{
	if (chooseStarts < 0) 
					throw std::invalid_argument("chooseStarts < 0");
	
	
	cout << "\nFind " << keep << " overdispersed starting points with maxLeavesCarving = " << maxLeavesCarving
	<< ", critSEB = " << critSEB
	<< ", maxLeavesSEB = " << maxLeavesSEB
	<< ", percentSpread = " << percentSpread << endl;

	boost::shared_ptr < AdaptiveHistogram  > adhPtr(new AdaptiveHistogram(box));
	
	bool successfulInsertion = 
			adhPtr->insertRvectorsFromVectorOfVecDbls(container);
	
	if (!successfulInsertion) 
					throw std::runtime_error("Failed to insert data");
	
	cout << "histogram has " << adhPtr->getRootCounter() << " values"
		" and variance-covariance matrix is" << endl;

	subpavings::RealVec varcovar = adhPtr->getDataCollectionVarCovar();
	printVarCovar(varcovar, adhPtr->getDimensions());
	
	cout << endl;
	
	// container for hist pointers
	std::vector< subpavings::AdaptiveHistogram* > hists;
	
	SPSNodeMeasureVolMassMinus compCarving(adhPtr->getRootCounter());
	
	AdaptiveHistogram::PrioritySplitQueueEvaluator 
			evaluatorCarving( compCarving, maxLeavesCarving);
	
	SPSNodeMeasureCount compSEB;
	
	AdaptiveHistogram::PrioritySplitQueueEvaluator 
			evaluatorSEB( compSEB, critSEB, maxLeavesSEB);
	
	LogCatalanPrior logPrior; // prior
		
	bool stopOnMaxPosterior = true;
	int precPQ = 5;

	double minVolDbl = _double(minVol);

	/* find start points overdispersed */
	CarverSEB::findStartingPointsOverdispersed(
					*adhPtr,
					hists,
					evaluatorCarving,
					evaluatorSEB,
					logPrior,
					minPoints,
					minVolDbl,
					chooseStarts,
					keep,
					percentSpread,
					stopOnMaxPosterior, 
					postFileName,
					checkPostFileNameBase,	
					precPQ);
	
	assert(!hists.empty());
	
	if (!memLogFilename.empty()) logUsage(memLogFilename, "Made starts");
	
	// at least one there
	
	boost::shared_ptr < PiecewiseConstantFunction > mcmcAvPtr
										= doAutoMCMC(maxLoops, 
											samplesNeeded,
											thinout, 
											tolerance,
											minPoints,
											minVol,
											hists,
											scalarsFileName,
											autotype,
											rep);
	
	
	
	if (!memLogFilename.empty()) logUsage(memLogFilename, "Made average");
	
	try {
		vector<AdaptiveHistogram*>::iterator ait;
		for (ait = hists.begin(); ait < hists.end(); ++ait) {
			if (NULL != *ait) delete (*ait);
			*ait = NULL;
		}
	}
	catch(...) {} //catch and swallow
	
	if (!memLogFilename.empty()) logUsage(memLogFilename, "xxStarts");
	
	
	return mcmcAvPtr;

}

