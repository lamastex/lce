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



#ifndef _INC_TEST_ROUTINES_H
#define _INC_TEST_ROUTINES_H

	#include "automcmc_factory.hpp"

	#include "piecewise_constant_function.hpp"
	#include "sptypes.hpp"
	
	#include <vector>
	#include <string>
	
	#include <boost/shared_ptr.hpp>


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
				const std::vector < double >& slicePts);
	
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
			const std::string& checkPostFileNameBase);
				
	subpavings::PiecewiseConstantFunction testCarverMCMC(
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
				subpavings::LOGGING_LEVEL logging,				
				const std::string& mcmcAvFilename,
				const std::string& mcmcSliceFilename,
				const std::vector < int >& sliceDims,
				const std::vector < double >& slicePts);
				
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
			const std::string& memLogFilename = "");
	
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
			const std::string& memLogFilename = "");
	
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
			const std::string& memLogFilename = "");
	
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
			const std::string& memLogFilename = "");
			
	subpavings::PiecewiseConstantFunction testBasicMCMC(
			const std::vector < std::vector < double > >& container,
			size_t minPoints,
			unsigned int loops,
			unsigned int burnin,
			unsigned int thinout,
			subpavings::LOGGING_LEVEL logging,				
			const std::string& mcmcAvFilename,
			const std::string& mcmcSliceFilename,
			const std::vector < int >& sliceDims,
			const std::vector < double >& slicePts);

#endif
