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

#include <vector>
#include <string>
#include <iostream>
#include <iterator>
#include <cassert>
#include <stdexcept>

using namespace std;
using namespace subpavings;


PiecewiseConstantFunction testBasicMCMC(
			const std::vector < std::vector < double > >& container,
			size_t minPoints,
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

	
	int prec = 10;
	bool confirm = true;
		

	PiecewiseConstantFunction mcmcAv = makeMCMCAv(adh, 
						loops,burnin, thinout,
						 minPoints, logging);
	

	// get the average as a pcf
	mcmcAv.outputToTxtTabs(mcmcAvFilename, prec, confirm);
	
	PiecewiseConstantFunction slice = mcmcAv.makeSlice(
										sliceDims, slicePts);

	// output slice
	slice.outputToTxtTabs(mcmcSliceFilename, prec, confirm);
	
	return mcmcAv;

}

