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


#include "MCMCGRAutoNew.hpp"


#include "cxsc.hpp"

#include "config.h" // mct


#include <iostream>
#include <fstream>
#include <iterator>
#include <iomanip>
#include <cassert>


using namespace std;
using namespace subpavings;


const PiecewiseConstantFunction makeMCMCAv(AdaptiveHistogram& adh, 
		unsigned int loops, unsigned int burnin, unsigned int thinout,
		int minPoints, LOGGING_LEVEL logging)
{
	
			
	// set up proposal distribution object
	UniformProposal proposal;
	// set up prior distribution object
	LogCatalanPrior logPrior;
		
	cout << endl << endl;
    cout << "\n\nStart mcmc\n" << endl;

			
	/* make a generator for this mcmc run */
	int seed = 1234;
	gsl_rng * rgsl = gsl_rng_alloc (gsl_rng_mt19937);
		
	gsl_rng_set(rgsl, seed);

	clock_t start, end;
	start = clock();

	PiecewiseConstantFunction mcmcAv = adh.MCMCNEW(
											loops, burnin, thinout,
											proposal, logPrior,
											minPoints, logging,
											rgsl);
	
	end = clock();

	cout << "Computing time : "
	   << ((static_cast<double>(end - start)) / CLOCKS_PER_SEC) << " s."<< endl;
	
	cout << "Finished MCMC get average" << endl;
	
	try {
		if (rgsl != NULL) gsl_rng_free (rgsl);
		rgsl = NULL;
	}
	catch(...) {}// catch and swallow
	
	return mcmcAv;

}


boost::shared_ptr < subpavings::PiecewiseConstantFunction >
			doAutoMCMC(int maxLoops, int samplesNeeded,
				int thinout, 
				real tolerance,
				size_t minPoints,
				std::vector< subpavings::AdaptiveHistogram* >& hists,
				const string& scalarsFileName,
				int rep)
{
	ahabc::AutoMCMCFactory::AUTOTYPE autotype = ahabc::AutoMCMCFactory::LEAVES;
	
	real minVol(0.0);
	
	return doAutoMCMC(maxLoops, samplesNeeded,
				thinout, 
				tolerance,
				minPoints,
				minVol,
				hists,
				scalarsFileName,
				autotype,
				rep);

	
}


boost::shared_ptr < subpavings::PiecewiseConstantFunction >
			doAutoMCMC(int maxLoops, int samplesNeeded,
				int thinout, 
				real tolerance,
				size_t minPoints,
				real minVol,
				std::vector< subpavings::AdaptiveHistogram* >& hists,
				const string& scalarsFileName,
				ahabc::AutoMCMCFactory::AUTOTYPE autotype,
				int rep)
{
	int chains = hists.size();		

	cout << "\ndoMCMCGRAuto with " << chains << " chains" << endl;
	cout << "maxLoops = " << maxLoops
			<< ", thinout\t" << thinout
			<< ", samplesNeeded\t" << samplesNeeded << endl;
	
	
	std::vector< MCMCGRAuto::Diagnostic * > diagPtrs;
	
	// fill in diagPtrs
	ahabc::AutoMCMCFactory::makeAutoMCMCDiags(autotype,
					diagPtrs,
					tolerance);
		
	{
		std::cout << "Diagnostics are:";
		// bolt together a list of the diagnostics we are using
		for (size_t d = 0; d < diagPtrs.size(); ++d) {
			string sName = diagPtrs[d]->getScalarsName();
			double tol = _double(diagPtrs[d]->getTolerance());
			cout << " " << sName << "(" << tol << ")";
		}
		cout << "\n" << endl;
	}
	
	unsigned long int seedMCMC = 1234+rep;
	
	
	MCMCGRAutoNew autoMCMC(diagPtrs, seedMCMC);
			
	
	int flagThreshold = static_cast<int>(diagPtrs.size());
	
	
	//---------------=Gelman-Rubin method------------------------=//
	 
	 // Gelman and Carlin [1996] recommend ten or more chains if the target
	 // distribution is unimodal (Martinez and Martinez, 2000, p. 453)
	 // the starting points are chosen to be widely dispersed

	clock_t startMCMC = clock();
	
	boost::shared_ptr < subpavings::PiecewiseConstantFunction > 
		samplesAvgPtr( new PiecewiseConstantFunction(
			autoMCMC.doMCMCGRAuto(hists, maxLoops, 
								samplesNeeded, thinout, 
								minPoints,
								minVol,
								flagThreshold,
								rep,
								scalarsFileName) ) );
	
	// the ptrs are newed -arrghh - so need to delete them
	// this is NOT GOOD design
	for (size_t i = 0; i < diagPtrs.size(); ++i) delete diagPtrs[i];
	
	// stop recording time here
	clock_t endMCMC = clock();	
	double timingMCMC = ((static_cast<double>(endMCMC - startMCMC)) / CLOCKS_PER_SEC);
	cout << "Computing time for MCMC: " << timingMCMC << " s."<< endl;

	cxsc::real avIntegral = samplesAvgPtr->getTotalIntegral();
	
	/* Check integral is close enough to 1. 
	This is Knuth's essentiallyEqual test using epsilon 10*cxsc::Epsilon
	It's not always the best 
	way to compare anything represented in the machine as a float though so not in toolz.  See
	http://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
	*/
	if(! ( cxsc::abs(avIntegral - cxsc::real(1.0)) 
		<= ( cxsc::Epsilon * 10 * (cxsc::abs(avIntegral) > cxsc::abs(cxsc::real(1.0)) ?
			cxsc::abs(cxsc::real(1.0)) : cxsc::abs(avIntegral)) ) ) ) {
		
		cerr << "integral of samples average = " << avIntegral << endl;
		cerr << "cxsc::Epsilon = " << cxsc::Epsilon << endl;
		cerr << "integral - 1 = " << (avIntegral - cxsc::real(1.0)) << endl;
		cerr << "cxsc::Epsilon * 10 * cxsc::abs(cxsc::real(1.0)) = " << (cxsc::Epsilon * 10 * cxsc::abs(cxsc::real(1.0))) << endl;
		cerr << "cxsc::Epsilon * 10 * cxsc::abs(avIntegral) = " << (cxsc::Epsilon * 10 * cxsc::abs(avIntegral)) << endl;
		throw std::runtime_error("samples average integral");
		
	}
	
	
	return samplesAvgPtr;
	
	
} // end of doAutoMCMC

