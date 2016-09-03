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

#ifndef _INC_TESTINGTOOLS_MCMCGR_H
#define _INC_TESTINGTOOLS_MCMCGR_H


#include "adaptivehistogram.hpp"

#include "piecewise_constant_function.hpp"

#include "automcmc_factory.hpp"


#include <boost/shared_ptr.hpp>


#include <string>
#include <vector>

	const subpavings::PiecewiseConstantFunction makeMCMCAv(
			subpavings::AdaptiveHistogram& adh, 
			unsigned int loops, unsigned int burnin, unsigned int thinout,
			int minPoints, subpavings::LOGGING_LEVEL logging);
			
	

	boost::shared_ptr < subpavings::PiecewiseConstantFunction > doAutoMCMC(
					int maxLoops, int samplesNeeded,
					int thinout, 
					real tolerance,
					size_t minPoints,
					std::vector< subpavings::AdaptiveHistogram* >& hists,
					const string& scalarsFileName,
					int rep = 1);
	
	boost::shared_ptr < subpavings::PiecewiseConstantFunction >
			doAutoMCMC(int maxLoops, int samplesNeeded,
				int thinout, 
				real tolerance,
				size_t minPoints,
				real minVol,
				std::vector< subpavings::AdaptiveHistogram* >& hists,
				const string& scalarsFileName,
				ahabc::AutoMCMCFactory::AUTOTYPE autotype,
				int rep = 1);
	
#endif
