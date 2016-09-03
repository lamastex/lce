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

#include "testing_tools.hpp"

#include "cxsc.hpp"

#include "config.h" // mct

#include "utilities.hpp"

#include <iostream>
#include <fstream>
#include <iterator>
#include <iomanip>
#include <cassert>
#include <numeric>
#include <algorithm>


#include <boost/lambda/bind.hpp>
#include <boost/lambda/lambda.hpp>

#include <gsl/gsl_statistics_double.h>


using namespace std;
using namespace mct;
using namespace subpavings;


// Method to put opening line into a log file
void outputAHABCLogFileStart(const std::string& s)
{
	// Make a string with filename and timestamp to start log file
	time_t rawtime;
	struct tm * timeinfo;
	time ( &rawtime );
	timeinfo = localtime ( &rawtime );
	std::ofstream os(s.c_str());         // replace data
	if (os.is_open()) {
		os << "File " << s << " created " <<  asctime (timeinfo) << std::endl;
		os.close();
	}
	else {
		std::cerr << "Error: could not open file named "
			<< s << std::endl << std::endl;
	}
}


// Method to add a line to a file
// Output goes to file named according to argument s
// does not add any new lines
void outputAHABCLogFileString(const std::string& s, 
								const std::string& line, 
								bool append)
{
	std::ofstream os;
	if (append) os.open(s.c_str(), ios::app);         // append
	else os.open(s.c_str()); // don't append
	if (os.is_open()) {
		os << "\n" << line << flush;
		os.close();
	}
	else {
		std::cerr << "Error: could not open file named "
			<< s << std::endl << std::endl;
	}
}


void printVarCovar(const subpavings::RealVec& varcovar, int dim)
{
	
	// loop through the elements and print as matrix
	for (int i = 0; i < dim; i++) {
		for (int j = 0; j < dim; j++) {
			std::cout << "  " << varcovar[(i*dim)+j];
		}
		std::cout << std::endl;
	}
}

void printHouseholderMatrix(
		boost::shared_ptr < const AHABCValues  > ahabcValuesPtr)
{
	std::vector < std::vector < double > > hm
						= ahabcValuesPtr->getHouseholderMatrix();
	
	cout << "\nHouseholder matrix is = " << endl;
	for (std::vector< std::vector<double> >::iterator it = hm.begin();
					it < hm.end();
					++it) {
		ostream_iterator<double> out_it (cout,"\t");
		copy ( it->begin(), it->end(), out_it );
		cout << endl;
	}
	cout << endl;
	
}

void printAndOutputMeansAndSDs(const std::vector < double >& means,
				const std::vector < double >& sds,
				const std::string& filename)

{
	
	cout << "\nMeans are " << endl;
	{
		ostream_iterator<double> out_it (cout,"\t");
		copy ( means.begin(), means.end(), out_it );
		cout << endl;
	}
	cout << endl;
	
	cout << "\nSDs are " << endl;
	{
		ostream_iterator<double> out_it (cout,"\t");
		copy ( sds.begin(), sds.end(), out_it );
		cout << endl;
	}
	cout << endl;
	
	size_t ms = means.size();
	size_t ss = sds.size();
	
	if (ms != ss) throw std::invalid_argument(
		"printAndOutputMeansAndSDs(...)");
	
	std::ofstream os;
	
	os.open(filename.c_str());         // overwrite 
	
	if (os.is_open()) {
	
		for (size_t i = 0; i < ms; ++i) {
			os << sds[i] << "\t" << means[i] << "\n";
		}
		os.close();
	}
	else {
		std::cerr << "Error: could not open file named "
		<< filename << std::endl << std::endl;
	}
	
}

void printMeansAndSDs(boost::shared_ptr < const AHABCValues  > ahabcValuesPtr)
{
	std::vector < double > means = ahabcValuesPtr->getMeans();
	std::vector < double > sds = ahabcValuesPtr->getSDs();
	
	cout << "\nMeans are " << endl;
	{
		ostream_iterator<double> out_it (cout,"\t");
		copy ( means.begin(), means.end(), out_it );
		cout << endl;
	}
	cout << endl;
	
	cout << "\nSDs are " << endl;
	{
		ostream_iterator<double> out_it (cout,"\t");
		copy ( sds.begin(), sds.end(), out_it );
		cout << endl;
	}
	cout << endl;
	
	
}

void outputStats(const std::string& s,
						const std::vector < double >& stats)
{
	ofstream os(s.c_str());         // replace data
	if (os.is_open()) {
		{
			ostream_iterator<double> out_it (os,"\n");
			copy ( stats.begin(), stats.end(), out_it );
		}
		os.close();
	}
	else {
		std::cerr << "Error: could not open file named "
			<< s << std::endl << std::endl;
	}
	
}

void printMeansAndSDsAndHouseholderMatrixStandardised(
		boost::shared_ptr < const AHABCValues  > ahabcValuesPtr)
{
	
	printMeansAndSDs(ahabcValuesPtr);
	
	printHouseholderMatrix(ahabcValuesPtr);
	
	
	
}

std::vector < std::vector < double > > getSlicePoints(
		boost::shared_ptr < mct::SummaryStatisticSet > refset_ptr)

{
	size_t n = refset_ptr->size();
	
	std::vector < std::vector < double > > slicePts;
	
	for (size_t i = 0; i < n; ++i) {
		
		slicePts.push_back( refset_ptr->at(i)->getValues() );
	}
	return slicePts;
		
}



/* standardise */
std::vector < double > standardise(
							const std::vector < double >& invec,
							const std::vector < double >& means,
							const std::vector < double >& sds)
{
	size_t d = means.size();
	
	if (invec.size() != d) throw std::invalid_argument("standardised(...)");
	
	std::vector < double > standardised = invec;
	for (size_t i = 0; i < d; ++i) {
		standardised[i] = (invec[i] - means[i]) / (sds[i]);
	}
	
	
	return standardised;
}

std::vector < double > getStandardisedSlicePoint(
		boost::shared_ptr < const mct::SummaryStatistic  > refStat,
		size_t np,
		const std::vector < double >& means,
		const std::vector < double >& sds)
{
	
	std::vector < double > slicePtsDummy(np, 0.0);
	
	std::vector < double > refValues = refStat->getValues();
	
	slicePtsDummy.insert(slicePtsDummy.end(), 
				refValues.begin(), refValues.end());
	
	std::vector < double > allSlicePts 
			= standardise(slicePtsDummy, means, sds);
			
	return std::vector < double > (allSlicePts.begin() + np, allSlicePts.end());
}

std::vector < double > getStandardisedParamSlicePoint(
		const std::vector < double >& paramPts,
		size_t np,
		size_t ns,
		const std::vector < double >& means,
		const std::vector < double >& sds)
{
	std::vector < double > slicePtsDummy(paramPts);
	
	std::vector < double > tmp(ns, 0.0);
	
	slicePtsDummy.insert(slicePtsDummy.end(), 
				tmp.begin(), tmp.end());
	
	std::vector < double > allSlicePts 
			= standardise(slicePtsDummy, means, sds);
			
	return std::vector < double > (allSlicePts.begin(), allSlicePts.begin() + np);
}




std::vector < std::vector < double > > getStandardisedSlicePoints(
		boost::shared_ptr < mct::SummaryStatisticSet > refset_ptr,
		size_t np,
		const std::vector < double >& means,
		const std::vector < double >& sds)

{
	size_t n = refset_ptr->size();
	
	std::vector < std::vector < double > > slicePts;
	
	for (size_t i = 0; i < n; ++i) {
		
		std::vector < double > slicePtsDummy;
		if (np) { // stuff np elements with 0s
			slicePtsDummy = std::vector < double >(np, 0.0);
		}
	
		std::vector < double > refValues = refset_ptr->at(i)->getValues();
		
		slicePtsDummy.insert(slicePtsDummy.end(), 
					refValues.begin(), refValues.end());
		
		std::vector < double > allSlicePts 
				= standardise(slicePtsDummy, means, sds);
				
		slicePts.push_back(
			std::vector < double > (allSlicePts.begin() + np, allSlicePts.end()) );
	}
	return slicePts;
		
}

std::vector < double > getStandardisedSlicePoint(
		boost::shared_ptr < const mct::SummaryStatistic  > refStat,
		boost::shared_ptr < const mct::AHABCValues  > ahabcValuesPtr)
{
	
	size_t np = ahabcValuesPtr->getNumberParameters();
	
	std::vector < double > slicePtsDummy(np, 0.0);
	
	std::vector < double > refValues = refStat->getValues();
	
	slicePtsDummy.insert(slicePtsDummy.end(), 
				refValues.begin(), refValues.end());
	
	std::vector < double > allSlicePts 
			= ahabcValuesPtr->standardise(slicePtsDummy);
			
	return std::vector < double > (allSlicePts.begin() + np, allSlicePts.end());
}

std::vector < double > getStandardisedParamSlicePoint(
		const std::vector < double >& paramPts,
		boost::shared_ptr < const mct::AHABCValues  > ahabcValuesPtr)
{
	size_t np = ahabcValuesPtr->getNumberParameters();
	assert( np == paramPts.size() );
	
	size_t ns = ahabcValuesPtr->getNumberSummaryStatistics();
	
	std::vector < double > slicePtsDummy(paramPts);
	
	std::vector < double > tmp(ns, 0.0);
	
	slicePtsDummy.insert(slicePtsDummy.end(), 
				tmp.begin(), tmp.end());
	
	std::vector < double > allSlicePts 
			= ahabcValuesPtr->standardise(slicePtsDummy);
			
	return std::vector < double > (allSlicePts.begin(), allSlicePts.begin() + np);
}

std::vector < std::vector < double > > getStandardisedSlicePoints(
		boost::shared_ptr < mct::SummaryStatisticSet > refset_ptr,
		boost::shared_ptr < const mct::AHABCValues  > ahabcValuesPtr)

{
	size_t np = ahabcValuesPtr->getNumberParameters();
	size_t n = refset_ptr->size();
	
	std::vector < std::vector < double > > slicePts;
	
	for (size_t i = 0; i < n; ++i) {
		
		std::vector < double > slicePtsDummy(np, 0.0);
	
		std::vector < double > refValues = refset_ptr->at(i)->getValues();
		
		slicePtsDummy.insert(slicePtsDummy.end(), 
					refValues.begin(), refValues.end());
		
		std::vector < double > allSlicePts 
				= ahabcValuesPtr->standardise(slicePtsDummy);
		slicePts.push_back(
			std::vector < double > (allSlicePts.begin() + np, allSlicePts.end()) );
	}
	return slicePts;
		
}

std::vector < double > getStandardisedTransformedSlicePoint(
		boost::shared_ptr < const mct::SummaryStatistic  > refStat,
		boost::shared_ptr < const mct::AHABCValues  > ahabcValuesPtr)
{
	
	size_t np = ahabcValuesPtr->getNumberParameters();
	
	std::vector < double > slicePtsDummy(np, 0.0);
	
	std::vector < double > refValues = refStat->getValues();
	
	slicePtsDummy.insert(slicePtsDummy.end(), 
				refValues.begin(), refValues.end());
	
	std::vector < double > allSlicePts 
			= ahabcValuesPtr->standardiseAndHouseholderTransform
										(slicePtsDummy);
	return std::vector < double > (allSlicePts.begin() + np, allSlicePts.end());
}


std::vector < std::vector < double > > getStandardisedTransformedSlicePoints(
		boost::shared_ptr < mct::SummaryStatisticSet > refset_ptr,
		boost::shared_ptr < const mct::AHABCValues  > ahabcValuesPtr)

{
	size_t np = ahabcValuesPtr->getNumberParameters();
	size_t n = refset_ptr->size();
	
	std::vector < std::vector < double > > slicePts;
	
	for (size_t i = 0; i < n; ++i) {
		
		std::vector < double > slicePtsDummy(np, 0.0);
	
		std::vector < double > refValues = refset_ptr->at(i)->getValues();
		
		slicePtsDummy.insert(slicePtsDummy.end(), 
					refValues.begin(), refValues.end());
		
		std::vector < double > allSlicePts 
				= ahabcValuesPtr->standardiseAndHouseholderTransform
											(slicePtsDummy);
		slicePts.push_back(
			std::vector < double > (allSlicePts.begin() + np, allSlicePts.end()) );
	}
	return slicePts;
		
}

cxsc::ivector getBox(const std::vector < std::vector < double > >& container,
			const std::vector < std::vector < double > >& slicePtsVec,
			const std::vector < int >& sliceDims,
			cxsc::real padding)
{
	
	
	cxsc::ivector box = subpavings::makeBox(container, padding);
	
	std::cout << "The box to fit just the data is" << std::endl;
    for (int i = 1; i <=VecLen(box); i++) {
        std::cout << box[i] << "  ";    
    }
    std::cout << std::endl;
	
	if (!slicePtsVec.empty()) {
	
		std::vector < double > maxs = slicePtsVec.front();
		std::vector < double > mins = maxs;
		
		for (size_t i = 0; i < sliceDims.size(); ++i) {
			int d = sliceDims[i];
			for (size_t j = 0; j < slicePtsVec.size(); ++j) {
				double v = slicePtsVec[j][i];
				if (v > maxs[i]) maxs[i] = v;
				if (v < mins[i]) mins[i] = v;
			}
			
			if (Inf(box[d]) > mins[i] - padding) {
				SetInf(box[d], mins[i] - padding);
			}
			if (Sup(box[d]) < maxs[i] + padding) {
				SetSup(box[d], maxs[i] + padding);
			}
		}
		
		std::cout << "The box to fit data and reference stats is" << std::endl;
    for (int i = 1; i <=VecLen(box); i++) {
        std::cout << box[i] << "  ";    
    }
    std::cout << std::endl;
	}
	return box;
}



void processStandardisedTransformedAverage(
			const PiecewiseConstantFunction&  pcf,
			boost::shared_ptr < const mct::AHABCValues  > ahabcValuesPtr,
			const std::string& pcfFilename,
			const std::string& pcfSliceFilename,
			const std::string& pcfUntransformedFilename,
			const std::string& pcfSliceUntransformedFilename,
			const std::vector < int >& sliceDims,
			const std::vector < double >& slicePts,
			int prec,
			bool confirm)
{
	// output
	pcf.outputToTxtTabs(pcfFilename, prec, confirm);
	
	//not bothering about sliceCapture here
	PiecewiseConstantFunction slice = pcf.makeSlice(
										sliceDims, slicePts);

	// output slice
	slice.outputToTxtTabs(pcfSliceFilename, prec, confirm);

	
	//untransform the average
	
	std::vector < double > scale = ahabcValuesPtr->getSDs();
	std::vector < double > shift = ahabcValuesPtr->getMeans();
	
	printMeansAndSDs(ahabcValuesPtr);
					
	std::vector < std::vector < double > > hm
			= ahabcValuesPtr->getHouseholderMatrixStandardised();
	
	PolygonPiecewiseConstantFunction
			ppcfAv (pcf, hm, scale, shift);
	
	ppcfAv.outputToTxtTabs(pcfUntransformedFilename, prec, confirm);

	//untransform the slice
	size_t np = ahabcValuesPtr->getNumberParameters();
	scale.resize(np);
	shift.resize(np);
	hm.resize(np);
	for (std::vector < std::vector < double > >::iterator it = hm.begin();
				it < hm.end(); ++it) {
		it->resize(np);			
					
	}
	PolygonPiecewiseConstantFunction
			ppcfSlice (slice, hm, scale, shift);
	
	ppcfSlice.outputToTxtTabs(pcfSliceUntransformedFilename, prec, confirm);
	
}

void processAverage(
			const PiecewiseConstantFunction&  pcf,
			boost::shared_ptr < const mct::AHABCValues  > ahabcValuesPtr,
			const std::string& prefix,
			const std::string& spec,
			const std::vector < int >& sliceDims,
			const std::vector < std::vector < double > >& slicePtsVec,
			const std::vector < double >& trueParams,
			int prec,
			bool confirm,
			size_t posteriorCaptureInterval)
{
	string pcfFilename = prefix + "_Av_" + spec + ".txt";
				
	string pcfPosteriorFilenameBase = prefix + "_Posterior_" + spec;
	
	string pcfParamMargFilename = prefix + "_Prior_" + spec + ".txt";
				
	string pcfMargFilenameBase = prefix + "_" + spec;
	
	string pcfSliceFilenameBase = prefix + "_" + spec;
	
	// output pcf
	pcf.outputToTxtTabs(pcfFilename, prec, confirm);
	
	// number of ref stats
	size_t nrefstats = slicePtsVec.size();
	
	int np = ahabcValuesPtr->getNumberParameters();
	int ns = ahabcValuesPtr->getNumberSummaryStatistics();
		
	// marginalisations
	{
		if (np > 1 || ns > 1) {
			
			string filenameBit("_Marg");
			string pcfMargFilenames = prefix + filenameBit + "s.txt";
			bool append = false;
			
			// for each param-ns pair do a marginal
			for (int p = 1; p <= np; ++p) {
				
				for (int s = 1; s <= ns; ++s) {
					std::vector<int> reqDims(1, p);
					reqDims.push_back(np+s);
					
					//for output filename
					std::vector < int > margDims(1,p);
					margDims.push_back(s);
					
					PiecewiseConstantFunction pcfMarg 
										= pcf.makeMarginal(reqDims);
										
					outputMarginal(	pcfMarg,
									pcfMargFilenameBase,
									pcfMargFilenames,
									filenameBit,
									margDims,
									prec,
									confirm,
									append);
			
					if (!append) append = true;
					
					/* Also want to make a slice of that marginal
					 * at the appropriate point of the ref stats
					 * - all our marginals are 2-d and have the ref
					 * stat on the second dim, so the slice dim
					 * in every case is 2 */
					std::vector<int> thisSliceDim(1, 2); 
					
					for (size_t i = 0; i < nrefstats; ++i) {
						vector < real > thisSlicePt(1, slicePtsVec[i][s-1]);
						PiecewiseConstantFunction tmp
							= pcfMarg.makeSlice(thisSliceDim, thisSlicePt);
					}
				
				}
				
			}
			
		}
		
		/* do marginals just on ss */
		if (ns > 0) {
			string filenameBit("_SummaryStatMarg");
			
			outputSummaryStatMarginals(pcf,
					ns,
					np,
					pcfMargFilenameBase,
					prefix,
					filenameBit,
					prec,
					confirm);
		}
	} // end marginals
	
	
	// Get a 'smeared' version of the pcf to use for slices and posteriors
	PiecewiseConstantFunction pcfSmearedAverage 
			= pcf.makeSmearZeroValues(1/(1.0*ahabcValuesPtr->size()));
			
	/* recover the prior as a histogram, as the marginal on the parameters */
	PiecewiseConstantFunction pcfPrior;
	{
		std::vector<int> reqDims(np, 0.0);
		for (int p = 1; p <= np; ++p) reqDims[p-1]=p;
		pcfPrior = pcfSmearedAverage.makeMarginal(reqDims);
		pcfPrior.normalise();// should not be needed ...
		// output pcfPrior
		pcfPrior.outputToTxtTabs(pcfParamMargFilename, prec, confirm);
	}
	
	PiecewiseConstantFunction pcfPosterior;
	
	//capture slice boxes only for 2-d data
	string sliceBoxesFilenameBase("");
	string sliceBoxesFilenames("");
	if (np == 1 && ns == 1) {
		sliceBoxesFilenameBase = prefix + "_SliceBoxes";
		sliceBoxesFilenames = prefix + "_SliceBoxes.txt";
	}
	bool appendSliceBoxes = false;
	
	string pcfSliceFilenames = prefix + "_Slices.txt";
	bool append = false;
	
	string pcfPosteriorFilenames = prefix + "_Posteriors.txt";
	bool appendLik = false;
			
	/* For each ref stat, now represented as a slice point
	 * in slicePtsVec, make  slice and output it*/
	for (size_t i = 0; i < nrefstats; ++i) {
		
		{
			cout << "This standardised slice point is\t";
			ostream_iterator<double> out_it(std::cout, "\t");
			copy ( slicePtsVec[i].begin(), slicePtsVec[i].end(), out_it );
		}
		std::cout << std::endl;	
		
		std::string sliceBoxesFilename("");
		if (!sliceBoxesFilenameBase.empty()) sliceBoxesFilename 
					= setSliceBoxesName(sliceBoxesFilenameBase,
								sliceBoxesFilenames,
								slicePtsVec[i],
								i,
								appendSliceBoxes);
		if (!appendSliceBoxes) appendSliceBoxes = true;
						
			
		PiecewiseConstantFunction slice 
					= pcfSmearedAverage.makeSlice(sliceDims, slicePtsVec[i],
										sliceBoxesFilename);

		outputSlice(slice,
					pcfSliceFilenameBase,
					pcfSliceFilenames,
					i,
					slicePtsVec[i],
					prec,
					confirm,
					append);
		if (!append) append = true;
				
		if (i == 0) pcfPosterior = slice;
		else {
			pcfPosterior *= slice;
			pcfPosterior/=pcfPrior;
		}
		
		if (posteriorCaptureInterval 
					&& ((i+1)%posteriorCaptureInterval == 0)
					&& (i < nrefstats-1)) {
			
			PiecewiseConstantFunction 
					pcfPosteriorNorm = pcfPosterior.makeNormalised();
						
			ostringstream oss;
			oss << pcfPosteriorFilenameBase << "_" << (i+1) << ".txt";
			pcfPosteriorNorm.outputToTxtTabs(oss.str(), prec, confirm);
			// output filename to file
			std::vector <int> tmp(1, i+1);
			outputFileLineAndInts(pcfPosteriorFilenames, 
						stripPath(oss.str()), tmp, appendLik);
			if (!appendLik) appendLik = true;		
		}
		
	}
	
	//normalise posterior
	pcfPosterior.normalise();
	
	// output posterior
	{
		ostringstream oss;
		oss << pcfPosteriorFilenameBase << "_" << nrefstats << ".txt";
		pcfPosterior.outputToTxtTabs(oss.str(), prec, confirm);
		std::vector <int> tmp(1, nrefstats);
		outputFileLineAndInts(pcfPosteriorFilenames, 
							stripPath(oss.str()), tmp, appendLik);
	}
	
	/* For true parameter, make slice and output marginals of it*/
	if (!trueParams.empty()) {
		
		std::vector < int > paramSliceDims;
		for (int i = 1; i <= np; ++i) paramSliceDims.push_back(i);
		
		std::vector < double > standardisedParamSlicePt 
			= getStandardisedParamSlicePoint(trueParams, ahabcValuesPtr);
			
		{
			cout << "The standardised true parameter values are\t";
			ostream_iterator<double> out_it(std::cout, "\t");
			copy ( standardisedParamSlicePt.begin(),
					standardisedParamSlicePt.end(), out_it );
		}
		
		PiecewiseConstantFunction slice 
					= pcfSmearedAverage.makeSlice(paramSliceDims, 
												standardisedParamSlicePt);

		
		/* do marginals just on ss */
		if (ns > 0) {
			string filenameBit("_TrueParamSliceSummaryStatMarg");
			
			int npNow = 0;
			outputSummaryStatMarginals(slice,
					ns,
					npNow,
					pcfMargFilenameBase,
					prefix,
					filenameBit,
					prec,
					confirm);
			
		}
		slice.normalise();

		/* do marginals just on ss */
		if (ns > 0) {
			string filenameBit("_TrueParamNormalisedSliceSummaryStatMarg");
			
			int npNow = 0;
			outputSummaryStatMarginals(slice,
					ns,
					npNow,
					pcfMargFilenameBase,
					prefix,
					filenameBit,
					prec,
					confirm);
			
		}

	}

}

//process average for mixture test
std::vector < boost::shared_ptr < PiecewiseConstantFunction > > 
		processAverage(
			const PiecewiseConstantFunction&  pcf,
			int np,
			int ns,
			size_t n_reps,
			const std::string& prefix,
			const std::string& spec,
			const std::vector < int >& sliceDims,
			const std::vector < std::vector < double > >& slicePtsVec,
			const std::vector < double >& trueParams,
			const std::vector < std::vector < double > >& otherParams,
			int prec,
			bool confirm,
			size_t posteriorCaptureInterval)
{
	std::vector < real > covs(1, 1.0);
	covs.push_back(0.95);
	covs.push_back(0.80);
	std::vector < std::string > covStrs;
	for (size_t ci = 0; ci < covs.size(); ++ci) {
				
		ostringstream oss;
		oss << "_cov" << _double(covs[ci]);
		covStrs.push_back(oss.str());
	}
	
	
	
	string pcfFilename = prefix + "_Av_" + spec + ".txt";
				
	string pcfFilenameBase = prefix + "_Av_" + spec;
		
	string pcfParamMargFilename = prefix + "_Prior_" + spec + ".txt";
		
	string pcfPosteriorFilenameBase = prefix + "_Posterior_" + spec;
			
	string pcfPosteriorCoverageFilenames = prefix + "PosteriorCoverages.txt";
	
	string pcfMargFilenameBase = prefix + "_" + spec;
	
	string pcfSliceFilenameBase = prefix + "_" + spec;
	
	
	// output pcf
	pcf.outputToTxtTabs(pcfFilename, prec, confirm);
	
	//coverages
	
	std::vector < double > scaleDummy;
	std::vector < double > shiftDummy;
		
	{
		string pcfCoverageFilenames = prefix + "Coverages.txt";
		bool append = false;
	
		outputCoverages(pcf,
					scaleDummy,
					shiftDummy,
					covs,
					covStrs,
					pcfFilenameBase,
					"", //dummy
					pcfCoverageFilenames,
					"", // dummy
					append,
					prec,
					confirm);
	}
	
	// number of ref stats
	size_t nrefstats = slicePtsVec.size();
	
	// marginalisations
	{
		if (np > 1 || ns > 1) {
			
			string filenameBit("_Marg");
			string pcfMargFilenames = prefix + filenameBit + "s.txt";
			bool append = false;
			
			std::string sliceBoxesFilenameBaseBase = prefix + filenameBit + "SliceBoxes";
			std::string sliceBoxesFilenames = prefix + filenameBit + "SliceBoxes.txt";
			bool appendSliceBoxes = false;
			
			// for each param-ns pair do a marginal
			for (int p = 1; p <= np; ++p) {
				
				for (int s = 1; s <= ns; ++s) {
					
					string sliceBoxesFilenameBase("");
					{
						std::ostringstream oss;
							
						oss << sliceBoxesFilenameBaseBase << "_" << p << "_" << s << "_";
						sliceBoxesFilenameBase = oss.str();
					}
					
					std::vector<int> reqDims(1, p);
					reqDims.push_back(np+s);
					
					//for output filename
					std::vector < int > margDims(1,p);
					margDims.push_back(s);
					
					PiecewiseConstantFunction pcfMarg 
										= pcf.makeMarginal(reqDims);
										
					outputMarginal(	pcfMarg,
									pcfMargFilenameBase,
									pcfMargFilenames,
									filenameBit,
									margDims,
									prec,
									confirm,
									append);
					
					/* Also want to make a slice of that marginal
					 * at the appropriate point of the ref stats
					 * - all our marginals are 2-d and have the ref
					 * stat on the second dim, so the slice dim
					 * in every case is 2 */
					std::vector<int> thisSliceDim(1, 2); 
					
					for (size_t i = 0; i < nrefstats; ++i) {
						
						vector < double > thisSlicePt(1, slicePtsVec[i][s-1]);
						
						std::string sliceBoxesFilename 
							= setSliceBoxesName(sliceBoxesFilenameBase,
								sliceBoxesFilenames,
								thisSlicePt,
								i,
								appendSliceBoxes);
						
						if (!appendSliceBoxes) appendSliceBoxes = true;
												
						PiecewiseConstantFunction tmp
							= pcfMarg.makeSlice(thisSliceDim, thisSlicePt,
										sliceBoxesFilename);
					}
					
					if (!append) append = true;
					
				}
				
			}
			
		}
		
		/* do marginals just on ss */
		if (ns > 0) {
			string filenameBit("_SummaryStatMarg");
			
			outputSummaryStatMarginals(pcf,
					ns,
					np,
					pcfMargFilenameBase,
					"", // dummy
					prefix,
					filenameBit,
					scaleDummy,
					shiftDummy,
					prec,
					confirm);
			
			
		}
	} // end marginals
	
	
	// Get a 'smeared' version of the pcf to use for slices and posteriors
	
	PiecewiseConstantFunction pcfSmearedAverage 
			= pcf.makeSmearZeroValues(1/(1.0*n_reps));
	
	/* recover the prior as a histogram, as the marginal on the parameters */
	PiecewiseConstantFunction pcfPrior;
	{
		std::vector<int> reqDims(np, 0.0);
		for (int p = 1; p <= np; ++p) reqDims[p-1]=p;
		pcfPrior = pcfSmearedAverage.makeMarginal(reqDims);
		pcfPrior.normalise();// should not be needed ...
		// output pcfPrior
		pcfPrior.outputToTxtTabs(pcfParamMargFilename, prec, confirm);
	}
	
	PiecewiseConstantFunction pcfPosterior;
	std::vector < boost::shared_ptr < PiecewiseConstantFunction > > posteriorPtrsVec;
	
	
	//capture slice boxes only for 2-d data
	string sliceBoxesFilenameBase("");
	string sliceBoxesFilenames("");
	if (np == 1 && ns == 1) {
		sliceBoxesFilenameBase = prefix + "_SliceBoxes";
		sliceBoxesFilenames = prefix + "_SlicesBoxes.txt";
	}
	bool appendSliceBoxes = false;
	
	string pcfSliceFilenames = prefix + "_Slices.txt";
	bool append = false;
	
	string pcfPosteriorFilenames = prefix + "_Posteriors.txt";
	bool appendLik = false;
	
	
	/* For each ref stat, now represented as a slice point
	 * in slicePtsVec, make  slice and output it*/
	
	cout << "starting slices: nrefstats = "  << nrefstats << endl;
	
	for (size_t i = 0; i < nrefstats; ++i) {
		
		{
			cout << "This standardised slice point is\t";
			ostream_iterator<double> out_it(std::cout, "\t");
			copy ( slicePtsVec[i].begin(), slicePtsVec[i].end(), out_it );
			std::cout << std::endl;
		}
		
		
		std::string sliceBoxesFilename("");
		if (!sliceBoxesFilenameBase.empty()) sliceBoxesFilename 
					= setSliceBoxesName(sliceBoxesFilenameBase,
								sliceBoxesFilenames,
								slicePtsVec[i],
								i,
								appendSliceBoxes);
		
		if (!appendSliceBoxes) appendSliceBoxes = true; 
		
		PiecewiseConstantFunction slice 
					= pcfSmearedAverage.makeSlice(sliceDims, slicePtsVec[i],
					sliceBoxesFilename);

		outputSlice(slice,
					pcfSliceFilenameBase,
					pcfSliceFilenames,
					i,
					slicePtsVec[i],
					prec,
					confirm,
					append);
		
		if (i == 0) pcfPosterior = slice;
		else {
			pcfPosterior *= slice;
			pcfPosterior/=pcfPrior;
		}
		
		/* capture the posterior every so many multiplications,
		 * ie normalise it and output as histogram and polygon
		 * and also get coverage region and that standardised and as polygon */
		if (posteriorCaptureInterval 
					&& ((i+1)%posteriorCaptureInterval == 0)
					&& (i < nrefstats-1)) { // not the final one
			PiecewiseConstantFunction 
					pcfPosteriorNorm = pcfPosterior.makeNormalised();
					
			boost::shared_ptr < PiecewiseConstantFunction > 
					postPtr( new PiecewiseConstantFunction(pcfPosteriorNorm));
			posteriorPtrsVec.push_back(postPtr);
					
			outputPosteriorsAndCoverage(
					pcfPosteriorNorm,
					pcfPosteriorFilenameBase,
					"",
					pcfPosteriorFilenames,
					"",
					pcfPosteriorCoverageFilenames,
					"",
					scaleDummy,
					shiftDummy,
					covs,
					covStrs,
					i+1,
					appendLik,
					prec,
					confirm);
					
					if (!appendLik) appendLik = true;
				
		}
		
			
		if (!append) append = true;
		
	}
	
	//normalise the final posterior
	pcfPosterior.normalise();
	
	{
		boost::shared_ptr < PiecewiseConstantFunction > 
					postPtr( new PiecewiseConstantFunction(pcfPosterior));
		posteriorPtrsVec.push_back(postPtr);
	}		
	
	outputPosteriorsAndCoverage(
			pcfPosterior,
			pcfPosteriorFilenameBase,
			"",
			pcfPosteriorFilenames,
			"",
			pcfPosteriorCoverageFilenames,
			"",
			scaleDummy,
			shiftDummy,
			covs,
			covStrs,
			nrefstats,
			appendLik,
			prec,
			confirm);
			
	cout << "starting true parameter slices" << endl;
	
	std::vector < int > paramSliceDims;
	for (int i = 1; i <= np; ++i) paramSliceDims.push_back(i);
	
	/* container for observed likelihoods given parameter values - 
	 * one outer vector for each observed summary stat */
	std::vector < std::vector < double > > obsLikelihoods(nrefstats);
	/* have to convert the slice points into rvectors */
	std::vector < rvector > rvectorSlicePtsVec(nrefstats, rvector(ns));
	if (!trueParams.empty() || !otherParams.empty() ) {
		for (size_t i = 0; i < nrefstats; ++i) {
				
			std::vector < double > slPt = slicePtsVec[i];
			for (int k = 1; k <= ns; ++k) rvectorSlicePtsVec[i][k] = slPt[k-1];
			
		}
	}
	
	/* For true parameter, make slice 
	 * and get 'likelihoods' for observed stats
	 * and output marginals of it*/
	if (!trueParams.empty()) {
		
		
		std::vector < double > paramSlicePt  = trueParams;
			
		PiecewiseConstantFunction slice 
					= pcfSmearedAverage.makeSlice(paramSliceDims, 
												paramSlicePt);

		cout << "About to do marginals on ss" << endl;

		/* do marginals just on ss */
		if (ns > 0) {
			string filenameBit("_TrueParamSliceSummaryStatMarg");
			
			int npNow = 0;
					
			outputSummaryStatMarginals(slice,
					ns,
					npNow,
					pcfMargFilenameBase,
					"",
					prefix,
					filenameBit,
					scaleDummy,
					shiftDummy,
					prec,
					confirm);
			
		}
		
		cout << "Finished marginals on ss: now normalise and get pointwise values for each ob " << endl;
		slice.normalise();
		/* get the pointwise values for each obs on normalised slice */
		for (size_t i = 0; i < nrefstats; ++i) {
			
			obsLikelihoods[i].push_back(
				_double(slice.pointwiseExtension(rvectorSlicePtsVec[i])) );
			
		}
		
		cout << "Do marginals on normalised slice" << endl;
		
		/* do marginals just on ss */
		if (ns > 0) {
			string filenameBit("_TrueParamNormalisedSliceSummaryStatMarg");
			
			int npNow = 0;
					
			outputSummaryStatMarginals(slice,
					ns,
					npNow,
					pcfMargFilenameBase,
					"",
					prefix,
					filenameBit,
					scaleDummy,
					shiftDummy,
					prec,
					confirm);
			
		}

	}
	
	cout << "Do slices on other parameter sets of interest" << endl;
		
	/* For other parameter sets of interest, make slice 
	 * and get 'likelihoods' for observed stats (also do marginals etc)*/
	for (size_t j = 0; j < otherParams.size(); ++j) {
		
		std::vector < double > paramSlicePt = otherParams[j];
		
		PiecewiseConstantFunction slice 
					= pcfSmearedAverage.makeSlice(paramSliceDims, 
												paramSlicePt);
		slice.normalise();
		
		cout << "Do marginals on normalised slice" << endl;
		
		/* do marginals just on ss */
		if (ns > 0) {
			
			ostringstream oss;
			oss << "_OtherParam";
			for (size_t k = 0; k < otherParams[j].size(); ++k) 
				oss << "_" << otherParams[j][k];
			oss << "_NormalisedSliceSummaryStatMarg";
			string filenameBit = oss.str();
			
			int npNow = 0;
					
			outputSummaryStatMarginals(slice,
					ns,
					npNow,
					pcfMargFilenameBase,
					"",
					prefix,
					filenameBit,
					scaleDummy,
					shiftDummy,
					otherParams[j],
					prec,
					confirm);
			
		}

		
		for (size_t i = 0; i < nrefstats; ++i) {
			
			obsLikelihoods[i].push_back(
				_double(slice.pointwiseExtension(rvectorSlicePtsVec[i])) );
			
		}

	}
	
	if (!trueParams.empty() || !otherParams.empty() ) {
		
		string eFilename = prefix + "_PtWiseExtns_" + spec + ".txt";
	
		outputPointwiseExtensions(eFilename,
			trueParams,
			otherParams,
			obsLikelihoods);
	}
	
	return posteriorPtrsVec;
}




std::vector < boost::shared_ptr < PiecewiseConstantFunction > > 
	processStandardisedAverage (
			const PiecewiseConstantFunction&  pcf,
			const std::vector < double >& scale,
			const std::vector < double >& shift,
			int np,
			int ns,
			size_t n_reps,
			const std::string& prefix,
			const std::string& spec,
			const std::vector < int >& sliceDims,
			boost::shared_ptr < mct::SummaryStatisticSet > refset_ptr,
			const std::vector < std::vector < double > >& slicePtsVec,
			const std::vector < double >& trueParams,
			int prec,
			bool confirm,
			size_t posteriorCaptureInterval)
{
	std::vector < std::vector < double > > otherParams;
	
	return processStandardisedAverage(
			pcf,
			scale,
			shift,
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
}

std::vector < boost::shared_ptr < PiecewiseConstantFunction > > 
		processStandardisedAverage(
			const PiecewiseConstantFunction&  pcf,
			const std::vector < double >& scale,
			const std::vector < double >& shift,
			int np,
			int ns,
			size_t n_reps,
			//boost::shared_ptr < const mct::AHABCValues  > ahabcValuesPtr,
			const std::string& prefix,
			const std::string& spec,
			const std::vector < int >& sliceDims,
			boost::shared_ptr < mct::SummaryStatisticSet > refset_ptr,
			const std::vector < std::vector < double > >& slicePtsVec,
			const std::vector < double >& trueParams,
			const std::vector < std::vector < double > >& otherParams,
			int prec,
			bool confirm,
			size_t posteriorCaptureInterval)
{
	std::vector < real > covs(1, 1.0);
	covs.push_back(0.95);
	covs.push_back(0.80);
	std::vector < std::string > covStrs;
	for (size_t ci = 0; ci < covs.size(); ++ci) {
				
		ostringstream oss;
		oss << "_cov" << _double(covs[ci]);
		covStrs.push_back(oss.str());
	}
	
	#if(0)
		string strCov;
		{
			ostringstream oss;
			oss << "_" << _double(cov) << "_";
			strCov = oss.str();
		}
	#endif
	
	string pcfFilename = prefix + "_Av_" + spec + ".txt";
				
	string pcfResizedFilename = prefix + "_AvPolygon_" + spec + ".txt";
	
	string pcfFilenameBase = prefix + "_Av_" + spec;
				
	string pcfResizedFilenameBase = prefix + "_AvPolygon_" + spec;
	
	string pcfParamMargFilename = prefix + "_Prior_" + spec + ".txt";
		
	string pcfPosteriorFilenameBase = prefix + "_Posterior_" + spec;
			
	string pcfPosteriorResizedFilenameBase = prefix + "_PosteriorPolygon_" + spec;
	
	string pcfPosteriorCoverageFilenames = prefix + "PosteriorCoverages.txt";
	
	string pcfPosteriorResizedCoverageFilenames = prefix + "PosteriorCoveragePolygons.txt";
	
	string pcfMargFilenameBase = prefix + "_" + spec;
	
	string pcfMargResizedFilenameBase = prefix + "_Polygon" + spec;
	
	string pcfSliceFilenameBase = prefix + "_" + spec;
	
	string pcfSliceResizedFilenameBase = prefix + "_Polygon" + spec;
	
	string meanAndSdFilename = prefix + "_MeanAndSd_" + spec + ".txt";
	
	
	// output pcf
	pcf.outputToTxtTabs(pcfFilename, prec, confirm);
	
	printAndOutputMeansAndSDs(shift,scale, meanAndSdFilename);
	
	//unstandardised the average
	
	PolygonPiecewiseConstantFunction
			ppcfAv (pcf, scale, shift);
	
	ppcfAv.outputToTxtTabs(pcfResizedFilename, prec, confirm);

	//coverages
	
	{
		string pcfCoverageFilenames = prefix + "Coverages.txt";
		string pcfResizedCoverageFilenames = prefix + "CoveragePolygons.txt";
		bool append = false;
	
		outputCoverages(pcf,
					scale,
					shift,
					covs,
					covStrs,
					pcfFilenameBase,
					pcfResizedFilenameBase,
					pcfCoverageFilenames,
					pcfResizedCoverageFilenames,
					append,
					prec,
					confirm);
	}
	
	// number of ref stats
	size_t nrefstats = slicePtsVec.size();
	
	// marginalisations
	{
		if (np > 1 || ns > 1) {
			
			string filenameBit("_Marg");
			string pcfMargFilenames = prefix + filenameBit + "s.txt";
			string pcfMargResizedFilenames = prefix + filenameBit + "Polygons.txt";
			bool append = false;
			
			std::string sliceBoxesFilenameBaseBase = prefix + filenameBit + "SliceBoxes";
			std::string sliceBoxesFilenames = prefix + filenameBit + "SliceBoxes.txt";
			bool appendSliceBoxes = false;
			
			// for each param-ns pair do a marginal
			for (int p = 1; p <= np; ++p) {
				
				for (int s = 1; s <= ns; ++s) {
					
					string sliceBoxesFilenameBase("");
					{
						std::ostringstream oss;
							
						oss << sliceBoxesFilenameBaseBase << "_" << p << "_" << s << "_";
						sliceBoxesFilenameBase = oss.str();
					}
					
					std::vector<int> reqDims(1, p);
					reqDims.push_back(np+s);
					
					//for output filename
					std::vector < int > margDims(1,p);
					margDims.push_back(s);
					
					std::vector<double> thisScale(1, scale[p-1]);
					thisScale.push_back(scale[np+s-1]);
					
					std::vector<double> thisShift(1, shift[p-1]);
					thisShift.push_back(shift[np+s-1]);
					
					PiecewiseConstantFunction pcfMarg 
										= pcf.makeMarginal(reqDims);
										
					outputMarginal(	pcfMarg,
									pcfMargFilenameBase,
									pcfMargFilenames,
									filenameBit,
									margDims,
									prec,
									confirm,
									append);
					
					/* Also want to make a slice of that marginal
					 * at the appropriate point of the ref stats
					 * - all our marginals are 2-d and have the ref
					 * stat on the second dim, so the slice dim
					 * in every case is 2 */
					std::vector<int> thisSliceDim(1, 2); 
					
					for (size_t i = 0; i < nrefstats; ++i) {
						
						vector < double > thisSlicePt(1, slicePtsVec[i][s-1]);
						
						std::string sliceBoxesFilename 
							= setSliceBoxesName(sliceBoxesFilenameBase,
								sliceBoxesFilenames,
								thisSlicePt,
								i,
								appendSliceBoxes);
						
						if (!appendSliceBoxes) appendSliceBoxes = true;
												
						PiecewiseConstantFunction tmp
							= pcfMarg.makeSlice(thisSliceDim, thisSlicePt,
										sliceBoxesFilename);
					}
					
					
					PolygonPiecewiseConstantFunction
						ppcfMarg (pcfMarg, thisScale, thisShift);
					
					outputPolygonMarginal(	ppcfMarg,
									pcfMargResizedFilenameBase,
									pcfMargResizedFilenames,
									filenameBit,
									margDims,
									prec,
									confirm,
									append);
					
					if (!append) append = true;
					
				}
				
			}
			
		}
		
		/* do marginals just on ss */
		if (ns > 0) {
			string filenameBit("_SummaryStatMarg");
			
			outputSummaryStatMarginals(pcf,
					ns,
					np,
					pcfMargFilenameBase,
					pcfMargResizedFilenameBase,
					prefix,
					filenameBit,
					scale,
					shift,
					prec,
					confirm);
			
			
		}
	} // end marginals
	
	
	// Get a 'smeared' version of the pcf to use for slices and posteriors
	
	PiecewiseConstantFunction pcfSmearedAverage 
			= pcf.makeSmearZeroValues(1/(1.0*n_reps));
	
	/* recover the prior as a histogram, as the marginal on the parameters */
	PiecewiseConstantFunction pcfPrior;
	{
		std::vector<int> reqDims(np, 0.0);
		for (int p = 1; p <= np; ++p) reqDims[p-1]=p;
		pcfPrior = pcfSmearedAverage.makeMarginal(reqDims);
		pcfPrior.normalise();// should not be needed ...
		// output pcfPrior
		pcfPrior.outputToTxtTabs(pcfParamMargFilename, prec, confirm);
	}
	
	PiecewiseConstantFunction pcfPosterior;
	std::vector < boost::shared_ptr < PiecewiseConstantFunction > > posteriorPtrsVec;
	
	//for untransforming the slices
	std::vector < double > scaleParams = scale;
	scaleParams.resize(np);
	std::vector < double > shiftParams = shift;
	shiftParams.resize(np);
	
	//capture slice boxes only for 2-d data
	string sliceBoxesFilenameBase("");
	string sliceBoxesFilenames("");
	if (np == 1 && ns == 1) {
		sliceBoxesFilenameBase = prefix + "_SliceBoxes";
		sliceBoxesFilenames = prefix + "_SlicesBoxes.txt";
	}
	bool appendSliceBoxes = false;
	
	string pcfSliceFilenames = prefix + "_Slices.txt";
	string pcfSliceResizedFilenames = prefix + "_SlicePolygons.txt";
	bool append = false;
	
	string pcfPosteriorFilenames = prefix + "_Posteriors.txt";
	string pcfPosteriorResizedFilenames = prefix + "_PosteriorPolygons.txt";
	bool appendLik = false;
	
	
	/* For each ref stat, now represented as a slice point
	 * in slicePtsVec, make  slice and output it*/
	
	cout << "starting slices: nrefstats = "  << nrefstats << endl;
	
	for (size_t i = 0; i < nrefstats; ++i) {
		
		std::string sliceBoxesFilename("");
		if (!sliceBoxesFilenameBase.empty()) sliceBoxesFilename 
					= setSliceBoxesName(sliceBoxesFilenameBase,
								sliceBoxesFilenames,
								slicePtsVec[i],
								i,
								appendSliceBoxes);
		
		if (!appendSliceBoxes) appendSliceBoxes = true; 
		
		PiecewiseConstantFunction slice 
					= pcfSmearedAverage.makeSlice(sliceDims, slicePtsVec[i],
					sliceBoxesFilename);

		outputSlice(slice,
					pcfSliceFilenameBase,
					pcfSliceFilenames,
					i,
					slicePtsVec[i],
					prec,
					confirm,
					append);
		
		if (i == 0) pcfPosterior = slice;
		else {
			pcfPosterior *= slice;
			pcfPosterior/=pcfPrior;
		}
		
		/* capture the posterior every so many multiplications,
		 * ie normalise it and output as histogram and polygon
		 * and also get coverage region and that standardised and as polygon */
		if (posteriorCaptureInterval 
					&& ((i+1)%posteriorCaptureInterval == 0)
					&& (i < nrefstats-1)) { // not the final one
			
			{
				ostringstream oss;
				oss << prefix << "_UnnormalisedPosteriorTmp" << i << ".txt";
				pcfPosterior.outputToTxtTabs(oss.str(), prec, confirm);
			}
			
			PiecewiseConstantFunction 
					pcfPosteriorNorm = pcfPosterior.makeNormalised();
					
			boost::shared_ptr < PiecewiseConstantFunction > 
					postPtr( new PiecewiseConstantFunction(pcfPosteriorNorm));
			posteriorPtrsVec.push_back(postPtr);
					
			outputPosteriorsAndCoverage(
					pcfPosteriorNorm,
					pcfPosteriorFilenameBase,
					pcfPosteriorResizedFilenameBase,
					pcfPosteriorFilenames,
					pcfPosteriorResizedFilenames,
					pcfPosteriorCoverageFilenames,
					pcfPosteriorResizedCoverageFilenames,
					scaleParams,
					shiftParams,
					covs,
					covStrs,
					i+1,
					appendLik,
					prec,
					confirm);
					
					if (!appendLik) appendLik = true;
				
		}
		
		/* finally make the polygon version of the slice and output that */
		PolygonPiecewiseConstantFunction
			ppcfSlice (slice, scaleParams, shiftParams);
	
		std::vector < double > slicePts = refset_ptr->at(i)->getValues();
	
		outputPolygonSlice(ppcfSlice,
								pcfSliceResizedFilenameBase,
								pcfSliceResizedFilenames,
								i,
								slicePts,
								prec,
								confirm,
								append);
			
		if (!append) append = true;
		
	}
	
	//normalise the final posterior
	pcfPosterior.normalise();
	
	{
		boost::shared_ptr < PiecewiseConstantFunction > 
					postPtr( new PiecewiseConstantFunction(pcfPosterior));
		posteriorPtrsVec.push_back(postPtr);
	}		
	
	outputPosteriorsAndCoverage(
			pcfPosterior,
			pcfPosteriorFilenameBase,
			pcfPosteriorResizedFilenameBase,
			pcfPosteriorFilenames,
			pcfPosteriorResizedFilenames,
			pcfPosteriorCoverageFilenames,
			pcfPosteriorResizedCoverageFilenames,
			scaleParams,
			shiftParams,
			covs,
			covStrs,
			nrefstats,
			appendLik,
			prec,
			confirm);
			
	cout << "starting true parameter slices" << endl;
	
	std::vector < int > paramSliceDims;
	for (int i = 1; i <= np; ++i) paramSliceDims.push_back(i);
	
	/* container for observed likelihoods given parameter values - 
	 * one outer vector for each observed summary stat */
	std::vector < std::vector < double > > obsLikelihoods(nrefstats);
	/* have to convert the slice points into rvectors */
	std::vector < rvector > rvectorSlicePtsVec(nrefstats, rvector(ns));
	if (!trueParams.empty() || !otherParams.empty() ) {
		for (size_t i = 0; i < nrefstats; ++i) {
				
			std::vector < double > slPt = slicePtsVec[i];
			for (int k = 1; k <= ns; ++k) rvectorSlicePtsVec[i][k] = slPt[k-1];
			
		}
	}
	
	/* For true parameter, make slice 
	 * and get 'likelihoods' for observed stats
	 * and output marginals of it*/
	if (!trueParams.empty()) {
		
		
		#if(0)
			std::vector < double > standardisedParamSlicePt 
				= getStandardisedParamSlicePoint(trueParams, ahabcValuesPtr);
		#endif
		std::vector < double > standardisedParamSlicePt 
				= getStandardisedParamSlicePoint(trueParams,
						np,
						ns,
						shift,
						scale);
			
		{
			cout << "The standardised true parameter values are\t";
			ostream_iterator<double> out_it(std::cout, "\t");
			copy ( standardisedParamSlicePt.begin(),
					standardisedParamSlicePt.end(), out_it );
			cout << endl;
		}
		
		PiecewiseConstantFunction slice 
					= pcfSmearedAverage.makeSlice(paramSliceDims, 
												standardisedParamSlicePt);

		cout << "About to do marginals on ss" << endl;

		/* do marginals just on ss */
		if (ns > 0) {
			string filenameBit("_TrueParamSliceSummaryStatMarg");
			
			int npNow = 0;
			const std::vector < double> scaleSS(scale.begin()+np,scale.end());
			const std::vector < double> shiftSS(shift.begin()+np,shift.end());
					
			outputSummaryStatMarginals(slice,
					ns,
					npNow,
					pcfMargFilenameBase,
					pcfMargResizedFilenameBase,
					prefix,
					filenameBit,
					scaleSS,
					shiftSS,
					prec,
					confirm);
			
		}
		
		cout << "Finished marginals on ss: now normalise and get pointwise values for each ob " << endl;
		slice.normalise();
		/* get the pointwise values for each obs on normalised slice */
		for (size_t i = 0; i < nrefstats; ++i) {
			
			obsLikelihoods[i].push_back(
				_double(slice.pointwiseExtension(rvectorSlicePtsVec[i])) );
			
		}
		
		cout << "Do marginals on normalised slice" << endl;
		
		/* do marginals just on ss */
		if (ns > 0) {
			string filenameBit("_TrueParamNormalisedSliceSummaryStatMarg");
			
			int npNow = 0;
			const std::vector < double> scaleSS(scale.begin()+np,scale.end());
			const std::vector < double> shiftSS(shift.begin()+np,shift.end());
					
			outputSummaryStatMarginals(slice,
					ns,
					npNow,
					pcfMargFilenameBase,
					pcfMargResizedFilenameBase,
					prefix,
					filenameBit,
					scaleSS,
					shiftSS,
					prec,
					confirm);
			
		}

	}
	
	cout << "Do slices on other parameter sets of interest" << endl;
		
	/* For other parameter sets of interest, make slice 
	 * and get 'likelihoods' for observed stats (also do marginals etc)*/
	for (size_t j = 0; j < otherParams.size(); ++j) {
		
		std::vector < double > standardisedParamSlicePt 
				= getStandardisedParamSlicePoint(otherParams[j],
						np,
						ns,
						shift,
						scale);
		
		PiecewiseConstantFunction slice 
					= pcfSmearedAverage.makeSlice(paramSliceDims, 
												standardisedParamSlicePt);
		slice.normalise();
		
				cout << "Do marginals on normalised slice" << endl;
		
		/* do marginals just on ss */
		if (ns > 0) {
			
			ostringstream oss;
			oss << "_OtherParam";
			for (size_t k = 0; k < otherParams[j].size(); ++k) 
				oss << "_" << otherParams[j][k];
			oss << "_NormalisedSliceSummaryStatMarg";
			string filenameBit = oss.str();
			
			int npNow = 0;
			const std::vector < double> scaleSS(scale.begin()+np,scale.end());
			const std::vector < double> shiftSS(shift.begin()+np,shift.end());
					
			outputSummaryStatMarginals(slice,
					ns,
					npNow,
					pcfMargFilenameBase,
					pcfMargResizedFilenameBase,
					prefix,
					filenameBit,
					scaleSS,
					shiftSS,
					otherParams[j],
					prec,
					confirm);
			
		}

		
		for (size_t i = 0; i < nrefstats; ++i) {
			
			obsLikelihoods[i].push_back(
				_double(slice.pointwiseExtension(rvectorSlicePtsVec[i])) );
			
		}

	}
	
	if (!trueParams.empty() || !otherParams.empty() ) {
		
		string eFilename = prefix + "_PtWiseExtns_" + spec + ".txt";
	
		outputPointwiseExtensions(eFilename,
			trueParams,
			otherParams,
			obsLikelihoods);
	}
	
	return posteriorPtrsVec;
}


void processStandardisedPredictivePosteriorAverage(
			const PiecewiseConstantFunction&  pcf,
			const std::vector < double >& scale,
			const std::vector < double >& shift,
			int np,
			int ns,
			const std::string& prefix,
			const std::string& spec,
			boost::shared_ptr < mct::SummaryStatisticSet > refset_ptr,
			const std::vector < std::vector < double > >& slicePtsVec,
			int prec,
			bool confirm)
{
	// number of ref stats
	size_t nrefstats = slicePtsVec.size();
	string strNrefstats;
	{
		ostringstream oss;
		oss << "_nrefstats" << nrefstats << "_";
		strNrefstats = oss.str();
	}
	
	string pcfPredictiveFilenameBase = prefix + "_Predictive" 
						+ strNrefstats + spec;
	string pcfPredictiveFilename = pcfPredictiveFilenameBase + ".txt";
				
	string pcfPredictiveResizedFilenameBase = prefix + "_PredictivePolygon" 
						+ strNrefstats + spec;
	string pcfPredictiveResizedFilename 
				= pcfPredictiveResizedFilenameBase + ".txt";
	
	string pcfPredictiveCoverageFilenames = prefix + 
			"_PredictiveCoverages" + strNrefstats + ".txt";
	
	string pcfPredictiveResizedCoverageFilenames = prefix + 
			"_PredictivePolygonCoverages" + strNrefstats + ".txt";
			
	string pcfPredictiveMargFilenameBase = prefix + "_Predictive" + spec;
	
	string pcfPredictiveMargResizedFilenameBase = prefix + "_PredictivePolygon" + spec;
	
	string ptCoveragesFilename = prefix + 
			"_PredictiveRefStatPointCoverages" + strNrefstats + ".txt";
	
	std::vector<double> predictiveScale(ns);
	std::vector<double> predictiveShift(ns);
	
	
	/* marginalise on parameters */
	PiecewiseConstantFunction posteriorPredictive;
	{
		std::vector < int > margDims(ns);
		for (int s1 = 1; s1 <= ns; ++s1) {
			margDims[s1-1] = np+s1;
			predictiveScale[s1-1] = scale[np + s1 -1];
			predictiveShift[s1-1] = shift[np + s1 -1];
		}
		posteriorPredictive = pcf.makeMarginal(margDims);
	}
	
	// output posterior predictive
	posteriorPredictive.outputToTxtTabs( 
						pcfPredictiveFilename, prec, confirm);
	
	//unstandardised the posterior predictive and output
	{
	
		PolygonPiecewiseConstantFunction
				ppcfPredictive (posteriorPredictive, 
					predictiveScale, predictiveShift);
		
		ppcfPredictive.outputToTxtTabs( 
						pcfPredictiveResizedFilename, prec, confirm);
						
	}
	
	// do coverages
	std::vector < real > covs(1, 1.0);
	covs.push_back(0.95);
	covs.push_back(0.80);
	std::vector < std::string > covStrs;
	for (size_t ci = 0; ci < covs.size(); ++ci) {
				
		ostringstream oss;
		oss << "_cov" << _double(covs[ci]);
		covStrs.push_back(oss.str());
	}
	{
		bool append = false;
		
		outputCoverages(
				posteriorPredictive,
				predictiveScale,
				predictiveShift,
				covs,
				covStrs,
				pcfPredictiveFilenameBase,
				pcfPredictiveResizedFilenameBase,
				pcfPredictiveCoverageFilenames,
				pcfPredictiveResizedCoverageFilenames,
				nrefstats,
				append,
				prec,
				confirm);
	}
	
	// do marginals on each pair of ss
	if (ns > 1) {	
		string filenameBit("_SummaryStatMarg");
		
		outputSummaryStatMarginals(pcf,
				ns,
				np,
				pcfPredictiveMargFilenameBase,
				pcfPredictiveMargResizedFilenameBase,
				prefix,
				filenameBit,
				scale,
				shift,
				prec,
				confirm);
		
			
	} // end marginals on predictive	
	
	
	cout << "Do pt coverages" << endl;
	
	/* container for points and coverage values  */
	std::vector < std::vector < double > > ptsAndCoverageValues(nrefstats);
	for (size_t i = 0; i < nrefstats; ++i) {
	
		std::vector < double > slPt = slicePtsVec[i];
		
		std::vector<double> values = refset_ptr->at(i)->getValues();
		
		rvector rv(ns);
		for (int k = 0; k < ns; ++k) rv[k+1] = slPt[k];
		
		values.push_back(_double(posteriorPredictive.findCoverage(rv)));
		ptsAndCoverageValues.push_back(values);
		
		
	}
	// output them
	bool append = false;
	mct_utilities::outputVectorVectorToFile(ptsAndCoverageValues,
						ptCoveragesFilename, append);
	
}


void processPredictivePosteriorAverage(
			const PiecewiseConstantFunction&  pcf,
			int np,
			int ns,
			const std::string& prefix,
			const std::string& spec,
			boost::shared_ptr < mct::SummaryStatisticSet > refset_ptr,
			const std::vector < std::vector < double > >& slicePtsVec,
			int prec,
			bool confirm)
{
	// number of ref stats
	size_t nrefstats = slicePtsVec.size();
	string strNrefstats;
	{
		ostringstream oss;
		oss << "_nrefstats" << nrefstats << "_";
		strNrefstats = oss.str();
	}
	
	string pcfPredictiveFilenameBase = prefix + "_Predictive" 
						+ strNrefstats + spec;
	string pcfPredictiveFilename = pcfPredictiveFilenameBase + ".txt";
				
	string pcfPredictiveCoverageFilenames = prefix + 
			"_PredictiveCoverages" + strNrefstats + ".txt";
	
	string pcfPredictiveMargFilenameBase = prefix + "_Predictive" + spec;
	
	string ptCoveragesFilename = prefix + 
			"_PredictiveRefStatPointCoverages" + strNrefstats + ".txt";
	
	std::vector<double> predictiveScaleDummy;
	std::vector<double> predictiveShiftDummy;
	
	
	/* marginalise on parameters */
	PiecewiseConstantFunction posteriorPredictive;
	{
		std::vector < int > margDims(ns);
		for (int s1 = 1; s1 <= ns; ++s1) {
			margDims[s1-1] = np+s1;
			
		}
		posteriorPredictive = pcf.makeMarginal(margDims);
	}
	
	// output posterior predictive
	posteriorPredictive.outputToTxtTabs( 
						pcfPredictiveFilename, prec, confirm);
	
	// do coverages
	std::vector < real > covs(1, 1.0);
	covs.push_back(0.95);
	covs.push_back(0.80);
	std::vector < std::string > covStrs;
	for (size_t ci = 0; ci < covs.size(); ++ci) {
				
		ostringstream oss;
		oss << "_cov" << _double(covs[ci]);
		covStrs.push_back(oss.str());
	}
	{
		bool append = false;
		
		outputCoverages(
				posteriorPredictive,
				predictiveScaleDummy,
				predictiveShiftDummy,
				covs,
				covStrs,
				pcfPredictiveFilenameBase,
				"",
				pcfPredictiveCoverageFilenames,
				"",
				nrefstats,
				append,
				prec,
				confirm);
	}
	
	// do marginals on each pair of ss
	if (ns > 1) {	
		string filenameBit("_SummaryStatMarg");
		
		outputSummaryStatMarginals(pcf,
				ns,
				np,
				pcfPredictiveMargFilenameBase,
				"",
				prefix,
				filenameBit,
				predictiveScaleDummy,
				predictiveShiftDummy,
				prec,
				confirm);
		
			
	} // end marginals on predictive	
	
	
	cout << "Do pt coverages" << endl;
	
	/* container for points and coverage values  */
	std::vector < std::vector < double > > ptsAndCoverageValues(nrefstats);
	for (size_t i = 0; i < nrefstats; ++i) {
	
		std::vector < double > slPt = slicePtsVec[i];
		
		std::vector<double> values = refset_ptr->at(i)->getValues();
		
		rvector rv(ns);
		for (int k = 0; k < ns; ++k) rv[k+1] = slPt[k];
		
		values.push_back(_double(posteriorPredictive.findCoverage(rv)));
		ptsAndCoverageValues.push_back(values);
		
		
	}
	// output them
	bool append = false;
	mct_utilities::outputVectorVectorToFile(ptsAndCoverageValues,
						ptCoveragesFilename, append);
	
}



// for processing results of AABC
void processStandardisedAverage(
			const PiecewiseConstantFunction&  pcf,
			bool setAppend,
			int setIndex,
			const std::vector < double >& sds,
			const std::vector < double >& means,
			const std::string& prefix,
			const std::string& spec,
			int prec,
			bool confirm)
{
	
	cout << "\n\nProcessing for set " << setIndex << endl;
	
	std::vector < real > covs(1, 1.0);
	covs.push_back(0.95);
	covs.push_back(0.80);
	std::vector < std::string > covStrs;
	for (size_t ci = 0; ci < covs.size(); ++ci) {
				
		ostringstream oss;
		oss << "_cov" << _double(covs[ci]);
		covStrs.push_back(oss.str());
	}
	string setIndexStr;
	{
		ostringstream oss;
		oss << "_" << setIndex;
		setIndexStr = oss.str();
	}
	
	string pcfFilename = prefix + "_Av_" + spec + setIndexStr + ".txt";
		
	string pcfResizedFilename = prefix + "_AvPolygon_" + spec + setIndexStr + ".txt";
	
	string pcfFilenameBase = prefix + "_Av_" + spec + setIndexStr;
		
	string pcfResizedFilenameBase = prefix + "_AvPolygon_" + spec + setIndexStr;
			
	string pcfMargFilenameBase = prefix + "_" + spec + setIndexStr;
	
	string pcfMargResizedFilenameBase = prefix + "_" + spec + "_Polygon" + setIndexStr;
		
	string pcfFilenames = prefix + "_" + spec + "_Avs.txt";
	string pcfResizedFilenames = prefix + "_" + spec + "_AvPolygons.txt";
	
	string pcfCoverageFilenames = prefix + "_Coverages.txt";
	string pcfResizedCoverageFilenames = prefix + "_CoveragePolygons.txt";
	
	// output pcf
	pcf.outputToTxtTabs(pcfFilename, prec, confirm);
	std::vector <int> tmp(1, setIndex);
	outputFileLineAndInts(pcfFilenames, stripPath(pcfFilename), 
				tmp, setAppend);
	
	//unstandardised the average
	
	PolygonPiecewiseConstantFunction
			ppcfAv (pcf, sds, means);
	
	ppcfAv.outputToTxtTabs(pcfResizedFilename, prec, confirm);
	outputFileLineAndInts(pcfResizedFilenames, stripPath(pcfResizedFilename),
				tmp, setAppend);

	
	outputCoverages(
			pcf,
			sds,
			means,
			covs,
			covStrs,
			pcfFilenameBase,
			pcfResizedFilenameBase,
			pcfCoverageFilenames,
			pcfResizedCoverageFilenames,
			setIndex,
			setAppend,
			prec,
			confirm);
		
	// marginalisations
	
	int np = sds.size();
	
	if (np > 1) {
		string filenameBit("_ParamMarg");
		
		outputAABCParamMarginals(
					pcf,
					np,
					pcfMargFilenameBase,
					pcfMargResizedFilenameBase,
					prefix,
					filenameBit,
					setAppend,
					setIndex,
					sds,
					means,
					covs,
					covStrs,
					prec,
					confirm);
		
		
	} // end marginals
	
}


// for processing results of AABC, non-standardised
void processAverage(
			const PiecewiseConstantFunction&  pcf,
			bool setAppend,
			int setIndex,
			const std::string& prefix,
			const std::string& spec,
			int prec,
			bool confirm)
{
	
	cout << "\n\nProcessing for set " << setIndex << endl;
	
	std::vector < real > covs(1, 1.0);
	covs.push_back(0.95);
	covs.push_back(0.80);
	std::vector < std::string > covStrs;
	for (size_t ci = 0; ci < covs.size(); ++ci) {
				
		ostringstream oss;
		oss << "_cov" << _double(covs[ci]);
		covStrs.push_back(oss.str());
	}
	string setIndexStr;
	{
		ostringstream oss;
		oss << "_" << setIndex;
		setIndexStr = oss.str();
	}
	
	string pcfFilename = prefix + "_Av_" + spec + setIndexStr + ".txt";
		
	string pcfFilenameBase = prefix + "_Av_" + spec + setIndexStr;
		
	string pcfMargFilenameBase = prefix + "_" + spec + setIndexStr;
	
	string pcfFilenames = prefix + "_" + spec + "_Avs.txt";
	string pcfCoverageFilenames = prefix + "_Coverages.txt";
	
	// output pcf
	pcf.outputToTxtTabs(pcfFilename, prec, confirm);
	std::vector <int> tmp(1, setIndex);
	outputFileLineAndInts(pcfFilenames, stripPath(pcfFilename), 
				tmp, setAppend);
	
	std::vector< double > meansDummy;
	std::vector< double > sdsDummy;
	
	outputCoverages(
			pcf,
			sdsDummy,
			meansDummy,
			covs,
			covStrs,
			pcfFilenameBase,
			"",
			pcfCoverageFilenames,
			"",
			setIndex,
			setAppend,
			prec,
			confirm);
		
	// marginalisations
	
	int np = pcf.getDimensions();
	
	if (np > 1) {
		string filenameBit("_ParamMarg");
		
		outputAABCParamMarginals(
					pcf,
					np,
					pcfMargFilenameBase,
					"",
					prefix,
					filenameBit,
					setAppend,
					setIndex,
					sdsDummy,
					meansDummy,
					covs,
					covStrs,
					prec,
					confirm);
		
		
	} // end marginals
	
}


/* process standardised average for summary stats only */
void processStandardisedSummaryStatsAverage(
			const PiecewiseConstantFunction&  pcf,
			size_t n_reps,
			const std::vector < double >& scale,
			const std::vector < double >& shift,
			int ns,
			const std::string& prefix,
			const std::string& spec,
			const std::vector < std::vector < double > >& slicePtsVec,
			int prec,
			bool confirm)
{
	std::vector < real > covs(1, 1.0);
	covs.push_back(0.95);
	covs.push_back(0.80);
	std::vector < std::string > covStrs;
	for (size_t ci = 0; ci < covs.size(); ++ci) {
				
		ostringstream oss;
		oss << "_cov" << _double(covs[ci]);
		covStrs.push_back(oss.str());
	}
	
	string pcfFilename = prefix + "_Av_" + spec + ".txt";
				
	string pcfResizedFilename = prefix + "_AvPolygon_" + spec + ".txt";
	
	string pcfFilenameBase = prefix + "_Av_" + spec;
				
	string pcfResizedFilenameBase = prefix + "_AvPolygon_" + spec;
	
	string pcfMargFilenameBase = prefix + "_" + spec;
	
	string pcfMargResizedFilenameBase = prefix + "_Polygon" + spec;
	
	string meanAndSdFilename = prefix + "_MeanAndSd_" + spec + ".txt";
	
	// output pcf
	pcf.outputToTxtTabs(pcfFilename, prec, confirm);
	
	printAndOutputMeansAndSDs(shift,scale, meanAndSdFilename);
	
	//unstandardised the average
	
	PolygonPiecewiseConstantFunction
			ppcfAv (pcf, scale, shift);
	
	ppcfAv.outputToTxtTabs(pcfResizedFilename, prec, confirm);

	//coverages
	
	{
		string pcfCoverageFilenames = prefix + "Coverages.txt";
		string pcfResizedCoverageFilenames = prefix + "CoveragePolygons.txt";
		bool append = false;
	
		outputCoverages(pcf,
					scale,
					shift,
					covs,
					covStrs,
					pcfFilenameBase,
					pcfResizedFilenameBase,
					pcfCoverageFilenames,
					pcfResizedCoverageFilenames,
					append,
					prec,
					confirm);
	}
	
		
	/* do marginals on ss */
	if (ns > 0) {
		string filenameBit("_SummaryStatMarg");
		
		outputSummaryStatMarginals(pcf,
				ns,
				pcfMargFilenameBase,
				pcfMargResizedFilenameBase,
				prefix,
				filenameBit,
				scale,
				shift,
				prec,
				confirm);
		
		
	}
	
	/* Do 'likelihoods' for each observed stat */
	
	PiecewiseConstantFunction pcfSmearedAverage 
			= pcf.makeSmearZeroValues(1/(1.0*n_reps));
	
	// number of ref stats
	size_t nrefstats = slicePtsVec.size();
	
	/* container for observed likelihoods given parameter values - 
	 * one outer vector for each observed summary stat */
	std::vector < double > obsLikelihoods;
	obsLikelihoods.reserve(nrefstats);
	
	for (size_t i = 0; i < nrefstats; ++i) {
		/* have to convert the slice points into rvectors */
		rvector rv(ns);	
		std::vector < double > slPt = slicePtsVec[i];
		for (int k = 1; k <= ns; ++k) rv[k] = slPt[k-1];
		
		#if(0)
			cout << "standardised points is " << rv << endl;
		#endif
		
		double lik = _double(pcfSmearedAverage.pointwiseExtension(rv));
		#if(0)
			cout << "lik is " << lik << "\n" << endl;
		#endif
		
		obsLikelihoods.push_back(lik) ;
		
	}
	string eFilename = prefix + "_PtWiseExtns_" + spec + ".txt";
	
	outputPointwiseExtensions(eFilename, obsLikelihoods);
	
}


void outputSummaryStatMarginals(const subpavings::PiecewiseConstantFunction& pcf,
					int ns,
					const std::string& pcfMargFilenameBase,
					const std::string& pcfMargResizedFilenameBase,
					const std::string& prefix,
					const std::string& filenameBit,
					const std::vector < double>& scale,
					const std::vector < double>& shift,
					int prec,
					bool confirm)
{
	int np =0;
	outputSummaryStatMarginals(pcf,
					ns,
					np,
					pcfMargFilenameBase,
					pcfMargResizedFilenameBase,
					prefix,
					filenameBit,
					scale,
					shift,
					prec,
					confirm);
}

void outputCoverages(
			const PiecewiseConstantFunction&  pcf,
			const std::vector < double >& scale,
			const std::vector < double >& shift,
			const std::vector < cxsc::real >& covs,
			const std::vector < std::string >& covStrs,
			const std::string& pcfFilenameBase,
			const std::string& pcfResizedFilenameBase,
			const std::string& pcfCoverageFilenames,
			const std::string& pcfResizedCoverageFilenames,
			bool append,
			int prec,
			bool confirm)
{
	std::vector < int > margDims;
	int setIndex = -1;
	outputCoverages(
			pcf,
			scale,
			shift,
			covs,
			covStrs,
			pcfFilenameBase,
			pcfResizedFilenameBase,
			pcfCoverageFilenames,
			pcfResizedCoverageFilenames,
			margDims,
			setIndex,
			append,
			prec,
			confirm);
			
}

void outputCoverages(
			const PiecewiseConstantFunction&  pcf,
			const std::vector < double >& scale,
			const std::vector < double >& shift,
			const std::vector < cxsc::real >& covs,
			const std::vector < std::string >& covStrs,
			const std::string& pcfFilenameBase,
			const std::string& pcfResizedFilenameBase,
			const std::string& pcfCoverageFilenames,
			const std::string& pcfResizedCoverageFilenames,
			const std::vector < int >& margDims,
			bool append,
			int prec,
			bool confirm)
{
	int setIndex = -1;
	outputCoverages(
			pcf,
			scale,
			shift,
			covs,
			covStrs,
			pcfFilenameBase,
			pcfResizedFilenameBase,
			pcfCoverageFilenames,
			pcfResizedCoverageFilenames,
			margDims,
			setIndex,
			append,
			prec,
			confirm);
			
}

void outputCoverages(
			const PiecewiseConstantFunction&  pcf,
			const std::vector < double >& scale,
			const std::vector < double >& shift,
			const std::vector < cxsc::real >& covs,
			const std::vector < std::string >& covStrs,
			const std::string& pcfFilenameBase,
			const std::string& pcfResizedFilenameBase,
			const std::string& pcfCoverageFilenames,
			const std::string& pcfResizedCoverageFilenames,
			int setIndex,
			bool append,
			int prec,
			bool confirm)
{
	std::vector < int > margDims;
	outputCoverages(
			pcf,
			scale,
			shift,
			covs,
			covStrs,
			pcfFilenameBase,
			pcfResizedFilenameBase,
			pcfCoverageFilenames,
			pcfResizedCoverageFilenames,
			margDims,
			setIndex,
			append,
			prec,
			confirm);
			
}

void outputCoverages(
			const PiecewiseConstantFunction&  pcf,
			const std::vector < double >& scale,
			const std::vector < double >& shift,
			const std::vector < cxsc::real >& covs,
			const std::vector < std::string >& covStrs,
			const std::string& pcfFilenameBase,
			const std::string& pcfResizedFilenameBase,
			const std::string& pcfCoverageFilenames,
			const std::string& pcfResizedCoverageFilenames,
			const std::vector < int >& margDims,
			int setIndex,
			bool append,
			int prec,
			bool confirm)
{
	
	for (size_t ci = 0; ci < covs.size(); ++ci) {
					
		real cov = covs[ci];
		
		string filename = pcfFilenameBase
				+ covStrs[ci] + ".txt";
		string filenameResized = pcfResizedFilenameBase
				+ covStrs[ci] + ".txt";
		
		pcf.outputCoverageRegion(filename, cov);
		std::vector < double > tmpdbl(1, _double(cov));
		
		if (!scale.empty()) {
			
			PolygonPiecewiseConstantFunction
				ppcf (pcf, scale, shift, cov);
			
			ppcf.outputToTxtTabs(filenameResized, prec, confirm);
		}
			
		if (setIndex < 0 && margDims.empty()) {
			outputFileLineAndDbls(pcfCoverageFilenames, 
						stripPath(filename), 
						tmpdbl, append);
			if (!scale.empty()) {
				outputFileLineAndDbls(pcfResizedCoverageFilenames, 
							stripPath(filenameResized), 
							tmpdbl, append);
			}
		}
		else {
			std::vector < int > tmp(margDims);
			if (setIndex >= 0) tmp.push_back(setIndex);
			outputFileLineAndDblsAndInts(pcfCoverageFilenames, 
						stripPath(filename), 
						tmpdbl, tmp, append);
			if (!scale.empty()) {
				outputFileLineAndDblsAndInts(pcfResizedCoverageFilenames, 
						stripPath(filenameResized), 
						tmpdbl, tmp, append);
			}
		}
		if (!append) append = true;
		
	}
}	

// coverage from AABC Marginal
void outputAABCMarginalCoverage(const subpavings::PiecewiseConstantFunction&  pcf,
			const std::vector < double >& scale,
			const std::vector < double >& shift,
			const std::vector < cxsc::real >& covs,
			const std::vector < std::string >& covStrs,
			std::string pcfFilenameBase,
			std::string pcfResizedFilenameBase,
			const std::string& pcfCoverageFilenames,
			const std::string& pcfResizedCoverageFilenames,
			const std::vector < int >& margDims,
			int setIndex,
			bool append,
			int prec,
			bool confirm)
{
	{
		ostringstream oss;
		oss << pcfFilenameBase;
		for (size_t i = 0; i < margDims.size(); ++i) {
			oss << "_" << margDims[i];
		} 
		pcfFilenameBase = oss.str();
	}
	{
		ostringstream oss;
		oss << pcfResizedFilenameBase;
		for (size_t i = 0; i < margDims.size(); ++i) {
			oss << "_" << margDims[i];
		} 
		pcfResizedFilenameBase = oss.str();
	}
		
	outputCoverages(pcf,
					scale,
					shift,
					covs,
					covStrs,
					pcfFilenameBase,
					pcfResizedFilenameBase,
					pcfCoverageFilenames,
					pcfResizedCoverageFilenames,
					margDims,
					setIndex,
					append,
					prec,
					confirm);
	
}

// marginals from AABC
void outputAABCParamMarginals(const subpavings::PiecewiseConstantFunction& pcf,
					int np,
					const std::string& pcfMargFilenameBase,
					const std::string& pcfMargResizedFilenameBase,
					const std::string& prefix,
					const std::string& filenameBit,
					bool append,
					int setIndex,
					const std::vector < double>& scale,
					const std::vector < double>& shift,
					const std::vector < cxsc::real >& covs,
					const std::vector < std::string >& covStrs,
					int prec,
					bool confirm)
{
	string pcfMargFilenames1D = prefix + filenameBit + "s1D.txt";
	string pcfMargResizedFilenames1D
					= prefix + filenameBit + "Polygons1D.txt";
	string pcfMargCoverageFilenames1D 
					= prefix  + filenameBit + "Coverages1D.txt";
	string pcfMargResizedCoverageFilenames1D 
					= prefix  + filenameBit + "CoveragePolygons1D.txt";
	bool append1D = append;
	string pcfMargFilenames2D = prefix + filenameBit + "s2D.txt";
	string pcfMargResizedFilenames2D
					= prefix + filenameBit + "Polygons2D.txt";
	string pcfMargCoverageFilenames2D 
					= prefix + filenameBit + "Coverages2D.txt";
	string pcfMargResizedCoverageFilenames2D 
					= prefix + filenameBit + "CoveragePolygons2D.txt";
	bool append2D = append;
				
	
	for (int p1 = 1; p1 <= np; ++p1) {
		
		std::vector<int> reqDims1D(1, p1);
		std::vector < int > margDims1D(1,p1);
		
		std::vector<double> thisScale1D;
		std::vector<double> thisShift1D;
		
		if (!scale.empty()) {
			thisScale1D = std::vector<double>(1, scale[p1-1]);
			thisShift1D = std::vector<double>(1, shift[p1-1]);
		}
	
		PiecewiseConstantFunction pcfMarg1D
								= pcf.makeMarginal(reqDims1D);
								
		outputMarginal(	pcfMarg1D,
						pcfMargFilenameBase,
						pcfMargFilenames1D,
						filenameBit,
						margDims1D,
						setIndex,
						prec,
						confirm,
						append1D);
		
		if (!scale.empty()) {
			PolygonPiecewiseConstantFunction
				ppcfMarg1D (pcfMarg1D, thisScale1D, thisShift1D);
			
			outputPolygonMarginal(	ppcfMarg1D,
							pcfMargResizedFilenameBase,
							pcfMargResizedFilenames1D,
							filenameBit,
							margDims1D,
							setIndex,
							prec,
							confirm,
							append1D);
							
		}
			
		string pcfMargCoverageFilenameBase = pcfMargFilenameBase 
				+ filenameBit;
		string pcfMargResizedCoverageFilenameBase = pcfMargResizedFilenameBase 
				+ filenameBit;
		
		outputAABCMarginalCoverage(pcfMarg1D,
						thisScale1D,
						thisShift1D,
						covs,
						covStrs,
						pcfMargCoverageFilenameBase,
						pcfMargResizedCoverageFilenameBase,
						pcfMargCoverageFilenames1D,
						pcfMargResizedCoverageFilenames1D,
						margDims1D,
						setIndex,
						append1D,
						prec,
						confirm);
			
		if (!append1D) append1D = true;
			
		
		for (int p2 = p1+1; p2 < np; ++p2) {
			std::vector<int> reqDims2D = reqDims1D;
			reqDims2D.push_back(p2);
			
			//for output filename
			std::vector < int > margDims2D =  margDims1D;
			margDims2D.push_back(p2);

			std::vector<double> thisScale2D = thisScale1D;
			std::vector<double> thisShift2D = thisShift1D;
			
			if (!scale.empty()) {
				
				thisScale2D.push_back(scale[p2-1]);
			
				thisShift2D.push_back(shift[p2-1]);
			}
			
			PiecewiseConstantFunction pcfMarg2D
								= pcf.makeMarginal(reqDims2D);
								
			outputMarginal(	pcfMarg2D,
							pcfMargFilenameBase,
							pcfMargFilenames2D,
							filenameBit,
							margDims2D,
							setIndex,
							prec,
							confirm,
							append2D);
			
			if (!scale.empty()) {
				PolygonPiecewiseConstantFunction
					ppcfMarg2D (pcfMarg2D, thisScale2D, thisShift2D);
				
				outputPolygonMarginal(	ppcfMarg2D,
								pcfMargResizedFilenameBase,
								pcfMargResizedFilenames2D,
								filenameBit,
								margDims2D,
								setIndex,
								prec,
								confirm,
								append2D);
			}
			
			outputAABCMarginalCoverage(pcfMarg2D,
						thisScale2D,
						thisShift2D,
						covs,
						covStrs,
						pcfMargCoverageFilenameBase,
						pcfMargResizedCoverageFilenameBase,
						pcfMargCoverageFilenames2D,
						pcfMargResizedCoverageFilenames2D,
						margDims2D,
						setIndex,
						append2D,
						prec,
						confirm);
			
			if (!append2D) append2D = true;
			
		}
		
	}
	
}

void outputSummaryStatMarginals(const subpavings::PiecewiseConstantFunction& pcf,
					int ns,
					int np,
					const std::string& pcfMargFilenameBase,
					const std::string& pcfMargResizedFilenameBase,
					const std::string& prefix,
					const std::string& filenameBit,
					const std::vector < double>& scale,
					const std::vector < double>& shift,
					int prec,
					bool confirm)
{
	
	std::vector < double> details;
					
	outputSummaryStatMarginals(pcf,
					ns,
					np,
					pcfMargFilenameBase,
					pcfMargResizedFilenameBase,
					prefix,
					filenameBit,
					scale,
					shift,
					details,
					prec,
					confirm);
}

void outputSummaryStatMarginals(const subpavings::PiecewiseConstantFunction& pcf,
					int ns,
					int np,
					const std::string& pcfMargFilenameBase,
					const std::string& pcfMargResizedFilenameBase,
					const std::string& prefix,
					const std::string& filenameBit,
					const std::vector < double>& scale,
					const std::vector < double>& shift,
					const std::vector < double>& details,
					int prec,
					bool confirm)
{
	string pcfMargFilenames1D = prefix + filenameBit + "s1D.txt";
	string pcfMargResizedFilenames1D
					= prefix + filenameBit + "Polygons1D.txt";
	bool append1D = false;
	string pcfMargFilenames2D = prefix + filenameBit + "s2D.txt";
	string pcfMargResizedFilenames2D
					= prefix + filenameBit + "Polygons2D.txt";
	bool append2D = false;
				
	for (int s1 = 1; s1 <= ns; ++s1) {
		
		std::vector<int> reqDims1D(1, np+s1);
		std::vector < int > margDims1D(1,s1);
		
		std::vector<double> thisScale1D;
		std::vector<double> thisShift1D;
		
		if (!scale.empty() ) {
			thisScale1D = std::vector<double>(1, scale[np+s1-1]);
			thisShift1D = std::vector<double>(1, shift[np+s1-1]);
		}
	
		PiecewiseConstantFunction pcfMarg1D
								= pcf.makeMarginal(reqDims1D);
								
		outputMarginal(	pcfMarg1D,
						pcfMargFilenameBase,
						pcfMargFilenames1D,
						filenameBit,
						margDims1D,
						details,
						prec,
						confirm,
						append1D);
		
		if (!scale.empty() ) {
			PolygonPiecewiseConstantFunction
				ppcfMarg1D (pcfMarg1D, thisScale1D, thisShift1D);
			
			outputPolygonMarginal(	ppcfMarg1D,
							pcfMargResizedFilenameBase,
							pcfMargResizedFilenames1D,
							filenameBit,
							margDims1D,
							details,
							prec,
							confirm,
							append1D);
		}
		
		if (!append1D) append1D = true;
			
		
		for (int s2 = s1+1; s2 <= ns; ++s2) {
			std::vector<int> reqDims2D = reqDims1D;
			reqDims2D.push_back(np+s2);
			
			//for output filename
			std::vector < int > margDims2D =  margDims1D;
			margDims2D.push_back(s2);
								
			std::vector<double> thisScale2D = thisScale1D;
			std::vector<double> thisShift2D = thisShift1D;
			
			if (!scale.empty() ) {
				thisScale2D.push_back(scale[np+s2-1]);
				thisShift2D.push_back(shift[np+s2-1]);
			}
			
			PiecewiseConstantFunction pcfMarg2D
								= pcf.makeMarginal(reqDims2D);
								
			outputMarginal(	pcfMarg2D,
							pcfMargFilenameBase,
							pcfMargFilenames2D,
							filenameBit,
							margDims2D,
							details,
							prec,
							confirm,
							append2D);
			
			if (!scale.empty() ) {
				PolygonPiecewiseConstantFunction
				ppcfMarg2D (pcfMarg2D, thisScale2D, thisShift2D);
			
				outputPolygonMarginal(	ppcfMarg2D,
								pcfMargResizedFilenameBase,
								pcfMargResizedFilenames2D,
								filenameBit,
								margDims2D,
								details,
								prec,
								confirm,
								append2D);
			}
			
			if (!append2D) append2D = true;
			
		}
		
	}
	
}

void outputSummaryStatMarginals(const subpavings::PiecewiseConstantFunction& pcf,
					int ns,
					int np,
					const std::string& pcfMargFilenameBase,
					const std::string& prefix,
					const std::string& filenameBit,
					int prec,
					bool confirm)
{
	string pcfMargFilenames1D = prefix + filenameBit + "s1D.txt";
	bool append1D = false;
	string pcfMargFilenames2D = prefix + filenameBit + "s2D.txt";
	bool append2D = false;
				
	for (int s1 = 1; s1 <= ns; ++s1) {
		
		std::vector<int> reqDims1D(1, np+s1);
		std::vector < int > margDims1D(1,s1);
	
		PiecewiseConstantFunction pcfMarg1D
								= pcf.makeMarginal(reqDims1D);
								
		outputMarginal(	pcfMarg1D,
						pcfMargFilenameBase,
						pcfMargFilenames1D,
						filenameBit,
						margDims1D,
						prec,
						confirm,
						append1D);
		
		if (!append1D) append1D = true;
		
		for (int s2 = s1+1; s2 <= ns; ++s2) {
			std::vector<int> reqDims2D = reqDims1D;
			reqDims2D.push_back(np+s2);
			
			//for output filename
			std::vector < int > margDims2D =  margDims1D;
			margDims2D.push_back(s2);
								
			PiecewiseConstantFunction pcfMarg2D
								= pcf.makeMarginal(reqDims2D);
								
			outputMarginal(	pcfMarg2D,
							pcfMargFilenameBase,
							pcfMargFilenames2D,
							filenameBit,
							margDims2D,
							prec,
							confirm,
							append2D);
			
			if (!append2D) append2D = true;
			
		}
		
	}
	
}

void outputMarginal(	const subpavings::PiecewiseConstantFunction& pcfMarg,
					const std::string& pcfMargFilenameBase,
					const std::string& pcfMargFilenames,
					const std::string& filenameBit,
					const std::vector < int >& margDims,
					int prec,
					bool confirm,
					bool append)
{
	int setIndex = -1;
	std::vector < double > details;
	outputMarginal(	pcfMarg,
					pcfMargFilenameBase,
					pcfMargFilenames,
					filenameBit,
					margDims,
					details,
					setIndex,
					prec,
					confirm,
					append);
}

void outputMarginal(	const subpavings::PiecewiseConstantFunction& pcfMarg,
					const std::string& pcfMargFilenameBase,
					const std::string& pcfMargFilenames,
					const std::string& filenameBit,
					const std::vector < int >& margDims,
					int setIndex,
					int prec,
					bool confirm,
					bool append)
{
	std::vector < double > details;
	outputMarginal(	pcfMarg,
					pcfMargFilenameBase,
					pcfMargFilenames,
					filenameBit,
					margDims,
					details,
					setIndex,
					prec,
					confirm,
					append);
}

void outputMarginal(	const subpavings::PiecewiseConstantFunction& pcfMarg,
					const std::string& pcfMargFilenameBase,
					const std::string& pcfMargFilenames,
					const std::string& filenameBit,
					const std::vector < int >& margDims,
					const std::vector < double >& details,
					int prec,
					bool confirm,
					bool append)
{
	int setIndex = -1;
	outputMarginal(	pcfMarg,
					pcfMargFilenameBase,
					pcfMargFilenames,
					filenameBit,
					margDims,
					details,
					setIndex,
					prec,
					confirm,
					append);
}

void outputMarginal(	const subpavings::PiecewiseConstantFunction& pcfMarg,
					const std::string& pcfMargFilenameBase,
					const std::string& pcfMargFilenames,
					const std::string& filenameBit,
					const std::vector < int >& margDims,
					const std::vector < double >& details,
					int setIndex,
					int prec,
					bool confirm,
					bool append)
{
	ostringstream oss;
	oss << pcfMargFilenameBase << filenameBit;
	for (size_t i = 0; i < margDims.size(); ++i) {
		oss << "_" << margDims[i];
	} 
	oss << ".txt";
	
	// output marginal
	string pcfMargFilename = oss.str();
	pcfMarg.outputToTxtTabs(pcfMargFilename, prec, confirm);
	
	// output filename to file
	if (details.empty() || (details.size() > 2)) {
		if ((margDims.size() == 2) ||( margDims.size() == 1)) {
			std::vector < int > tmp = margDims;
			if (!(setIndex < 0)) tmp.push_back(setIndex); // arggh 
			outputFileLineAndInts(pcfMargFilenames, stripPath(pcfMargFilename), 
					tmp, append);
		}
		else {
			outputFileLine(pcfMargFilenames, stripPath(pcfMargFilename), append);
		}
	}
	else { // details size 1 or 2
		if ((margDims.size() == 2) || (margDims.size() == 1)) {
			std::vector < int > tmp = margDims;
			if (!(setIndex < 0)) tmp.push_back(setIndex); // arggh 
			outputFileLineAndDblsAndInts(pcfMargFilenames, stripPath(pcfMargFilename), 
					details, tmp, append);
		}
		
		else {
			outputFileLineAndDbls(pcfMargFilenames, 
					stripPath(pcfMargFilename), 
					details, append);
		}
		
	}
	
}

void outputPolygonMarginal(	const subpavings::PolygonPiecewiseConstantFunction& ppcfMarg,
					const std::string& pcfMargResizedFilenameBase,
					const std::string& pcfMargResizedFilenames,
					const std::string& filenameBit,
					const std::vector < int >& margDims,
					int prec,
					bool confirm,
					bool append)
{
	int setIndex = -1;
	std::vector < double > details;
	
	outputPolygonMarginal(ppcfMarg,
					pcfMargResizedFilenameBase,
					pcfMargResizedFilenames,
					filenameBit,
					margDims,
					details,
					setIndex,
					prec,
					confirm,
					append);
}

void outputPolygonMarginal(	const subpavings::PolygonPiecewiseConstantFunction& ppcfMarg,
					const std::string& pcfMargResizedFilenameBase,
					const std::string& pcfMargResizedFilenames,
					const std::string& filenameBit,
					const std::vector < int >& margDims,
					int setIndex,
					int prec,
					bool confirm,
					bool append)
{
	
	std::vector < double > details;
	
	outputPolygonMarginal(ppcfMarg,
					pcfMargResizedFilenameBase,
					pcfMargResizedFilenames,
					filenameBit,
					margDims,
					details,
					setIndex,
					prec,
					confirm,
					append);
}

void outputPolygonMarginal(	const subpavings::PolygonPiecewiseConstantFunction& ppcfMarg,
					const std::string& pcfMargResizedFilenameBase,
					const std::string& pcfMargResizedFilenames,
					const std::string& filenameBit,
					const std::vector < int >& margDims,
					const std::vector < double >& details,
					int prec,
					bool confirm,
					bool append)
{
	int setIndex = -1;
	
	outputPolygonMarginal(ppcfMarg,
					pcfMargResizedFilenameBase,
					pcfMargResizedFilenames,
					filenameBit,
					margDims,
					details,
					setIndex,
					prec,
					confirm,
					append);
}

void outputPolygonMarginal(	const subpavings::PolygonPiecewiseConstantFunction& ppcfMarg,
					const std::string& pcfMargResizedFilenameBase,
					const std::string& pcfMargResizedFilenames,
					const std::string& filenameBit,
					const std::vector < int >& margDims,
					const std::vector < double >& details,
					int setIndex,
					int prec,
					bool confirm,
					bool append)
{
	ostringstream oss;
	oss << pcfMargResizedFilenameBase << filenameBit;
	for (size_t i = 0; i < margDims.size(); ++i) {
		oss << "_" << margDims[i];
	} 
	oss << ".txt";
	
	// output marginal
	string pcfMargResizedFilename = oss.str();
	ppcfMarg.outputToTxtTabs(pcfMargResizedFilename, prec, confirm);
	
	// output filename to file
	if (details.empty() || (details.size() > 2)) {
		if ((margDims.size() == 2) || (margDims.size() == 1)) {
			std::vector < int > tmp = margDims;
			if (!(setIndex < 0)) tmp.push_back(setIndex); // arggh 
			outputFileLineAndInts(pcfMargResizedFilenames, stripPath(pcfMargResizedFilename), 
					tmp, append);
		}
		else {
			outputFileLine(pcfMargResizedFilenames, stripPath(pcfMargResizedFilename), append);
		}
	}
	else { // details size 1 or 2
		if ((margDims.size() == 2) || (margDims.size() == 1)) {
			std::vector < int > tmp = margDims;
			if (!(setIndex < 0)) tmp.push_back(setIndex); // arggh 
			outputFileLineAndDblsAndInts(pcfMargResizedFilenames, stripPath(pcfMargResizedFilename), 
					details, tmp, append);
		}
		else {
			outputFileLineAndDbls(pcfMargResizedFilenames,  
				stripPath(pcfMargResizedFilename), 
				details, append);
		}
		
	}
	
}

std::string setSliceBoxesName(const std::string sliceBoxesFilenameBase,
								const std::string sliceBoxesFilenames,
								const std::vector < double >& slicePts,
								size_t i,
								bool append)
{
	std::ostringstream oss;
		
	oss << sliceBoxesFilenameBase  << (i+1);
	for (size_t j = 0; j < slicePts.size(); ++j) {
			oss << "_" << _double(slicePts[j]);
	}
	oss << ".txt";
	std::string sliceBoxesFilename = oss.str();

	if (slicePts.size() == 2 || slicePts.size() == 1) {
	outputFileLineAndDbls(sliceBoxesFilenames, stripPath(sliceBoxesFilename), 
			slicePts, append);
	}
	else {
		outputFileLine(sliceBoxesFilenames, stripPath(sliceBoxesFilename), append);
	}
	return sliceBoxesFilename; 
}

void outputSlice(	const subpavings::PiecewiseConstantFunction& slice,
					const std::string& pcfSliceFilenameBase,
					const std::string& pcfSliceFilenames,
					size_t i,
					const std::vector < double >& slicePts,
					int prec,
					bool confirm,
					bool append)
{
	ostringstream oss;
	oss << pcfSliceFilenameBase << "_Slice" << (i+1)<< "_";
	ostream_iterator<double> out_it(oss, "_");
	copy ( slicePts.begin(), slicePts.end(), out_it );
	oss << ".txt";
	// output slice
	string pcfSliceFilename = oss.str();
	slice.outputToTxtTabs(pcfSliceFilename, prec, confirm);
	// output filename to file
	if (slicePts.size() == 2 || slicePts.size() == 1) {
		outputFileLineAndDbls(pcfSliceFilenames, stripPath(pcfSliceFilename), 
				slicePts, append);
	}
	else {
		outputFileLine(pcfSliceFilenames, stripPath(pcfSliceFilename), append);
	}
}

void outputPolygonSlice(	
			const subpavings::PolygonPiecewiseConstantFunction& ppcfSlice,
			const std::string& pcfSliceResizedFilenameBase,
			const std::string& pcfSliceResizedFilenames,
			size_t i,
			const std::vector < double >& slicePts,
			int prec,
			bool confirm,
			bool append)
{
	ostringstream oss;
	oss << pcfSliceResizedFilenameBase << "_Slice" << (i+1)<< "_";
	ostream_iterator<double> out_it(oss, "_");
	copy ( slicePts.begin(), slicePts.end(), out_it );
	oss << ".txt";
	// output slice
	string ppcfSliceFilename = oss.str();
	ppcfSlice.outputToTxtTabs(ppcfSliceFilename, prec, confirm);
	// output filename to file
	if (slicePts.size() == 2 || slicePts.size() == 1) {
		outputFileLineAndDbls(pcfSliceResizedFilenames, stripPath(ppcfSliceFilename), 
				slicePts, append);
	}
	else {
		outputFileLine(pcfSliceResizedFilenames, stripPath(ppcfSliceFilename), append);
	}
}


void outputParamSlice(	const subpavings::PiecewiseConstantFunction& slice,
					const std::string& pcfSliceFilenameBase,
					const std::vector < double >& slicePts,
					int prec,
					bool confirm)
{
	ostringstream oss;
	oss << pcfSliceFilenameBase << "_ParamSlice" << "_";
	ostream_iterator<double> out_it(oss, "_");
	copy ( slicePts.begin(), slicePts.end(), out_it );
	oss << ".txt";
	// output slice
	string pcfSliceFilename = oss.str();
	slice.outputToTxtTabs(pcfSliceFilename, prec, confirm);
	
}

void outputPolygonParamSlice(	
			const subpavings::PolygonPiecewiseConstantFunction& ppcfSlice,
			const std::string& pcfSliceResizedFilenameBase,
			const std::vector < double >& slicePts,
			int prec,
			bool confirm)
{
	ostringstream oss;
	oss << pcfSliceResizedFilenameBase << "_ParamSlice" << "_";
	ostream_iterator<double> out_it(oss, "_");
	copy ( slicePts.begin(), slicePts.end(), out_it );
	oss << ".txt";
	// output slice
	string ppcfSliceFilename = oss.str();
	ppcfSlice.outputToTxtTabs(ppcfSliceFilename, prec, confirm);
	
}

void outputPosteriorsAndCoverage(
		const subpavings::PiecewiseConstantFunction& pcfPosteriorNorm,
		const std::string& pcfPosteriorFilenameBase,
		const std::string& pcfPosteriorResizedFilenameBase,
		const std::string& pcfPosteriorFilenames,
		const std::string& pcfPosteriorResizedFilenames,
		const std::string& pcfPosteriorCoverageFilenames,
		const std::string& pcfPosteriorCoverageResizedFilenames,
		const std::vector < double>& scale,
		const std::vector < double>& shift,
		const std::vector < cxsc::real >& covs,
		const std::vector < std::string >& covStrs,
		int i,
		bool appendLik,
		int prec,
		bool confirm)
{
	{
		ostringstream oss;
		oss << pcfPosteriorFilenameBase << "_" << i << ".txt";
		pcfPosteriorNorm.outputToTxtTabs(oss.str(), prec, confirm);
		// output filename to file
		std::vector <int> tmp(1, i);
		outputFileLineAndInts(pcfPosteriorFilenames, 
					stripPath(oss.str()), tmp, appendLik);
	}
	
	if (!scale.empty()) {
		PolygonPiecewiseConstantFunction
			ppcfPosteriorNorm (pcfPosteriorNorm, scale, shift);
		{
			ostringstream oss;
			oss << pcfPosteriorResizedFilenameBase << "_" << i << ".txt";
			ppcfPosteriorNorm.outputToTxtTabs(oss.str(), prec, confirm);
			// output filename to file
			std::vector <int> tmp(1, i);
			outputFileLineAndInts(pcfPosteriorResizedFilenames, 
						stripPath(oss.str()), tmp, appendLik);
		}
	}
	
	string ppcfCoverageFilenameBase("");
	{
		ostringstream oss;
		oss << pcfPosteriorFilenameBase << "_" << i;
		ppcfCoverageFilenameBase = oss.str();
	}
	string ppcfCoverageResizedFilenameBase("");
	{
		ostringstream oss;
		oss << pcfPosteriorResizedFilenameBase << "_" << i;
		ppcfCoverageResizedFilenameBase = oss.str();
	}
	
	std::string filenameBit("");
	
	outputCoverages(pcfPosteriorNorm,
					scale,
					shift,
					covs,
					covStrs,
					ppcfCoverageFilenameBase,
					ppcfCoverageResizedFilenameBase,
					pcfPosteriorCoverageFilenames,
					pcfPosteriorCoverageResizedFilenames,
					i,
					appendLik,
					prec,
					confirm);
}




void processStandardisedTransformedAverage(
			const PiecewiseConstantFunction&  pcf,
			boost::shared_ptr < const mct::AHABCValues  > ahabcValuesPtr,
			const std::string& pcfFilename,
			const std::string& pcfUntransformedFilename,
			const std::string& pcfSliceFilenameBase,
			const std::string& pcfSliceUntransformedFilenameBase,
			const std::vector < int >& sliceDims,
			boost::shared_ptr < mct::SummaryStatisticSet > refset_ptr,
			const std::vector < std::vector < double > >& slicePtsVec,
			int prec,
			bool confirm)
{
	// output pcf
	pcf.outputToTxtTabs(pcfFilename, prec, confirm);
	
	// means and sds - used for reversing transformation
	std::vector < double > scale = ahabcValuesPtr->getSDs();
	std::vector < double > shift = ahabcValuesPtr->getMeans();
	
	// do slices
	std::vector < PiecewiseConstantFunction > slices;
	
	/* For each ref stat, now represented as a slice point
	 * in slicePtsVec, make  slice and output it*/
	size_t nrefstats = slicePtsVec.size();
	
	for (size_t i = 0; i < nrefstats; ++i) {
		
		{
			cout << "This standardised transformed slice point is\t";
			ostream_iterator<double> out_it(std::cout, "\t");
			copy ( slicePtsVec[i].begin(), slicePtsVec[i].end(), out_it );
		}
		std::cout << std::endl;		
		slices.push_back(pcf.makeSlice(sliceDims, slicePtsVec[i]));

		ostringstream oss;
		oss << pcfSliceFilenameBase << "_";
		ostream_iterator<double> out_it(oss, "_");
		copy ( slicePtsVec[i].begin(), slicePtsVec[i].end(), out_it );
		oss << "slice_" << (i+1) << ".txt";
		// output slice
		string pcfSliceFilename = oss.str();
		slices.back().outputToTxtTabs(pcfSliceFilename, prec, confirm);
	}
	// slices holds our slices
	assert(slices.size() == nrefstats);
	
	//untransform the average
	
	std::vector < std::vector < double > > hm
			= ahabcValuesPtr->getHouseholderMatrixStandardised();
	
	PolygonPiecewiseConstantFunction
			ppcfAv (pcf, hm, scale, shift);
	
	ppcfAv.outputToTxtTabs(pcfUntransformedFilename, prec, confirm);

	//untransform the slice
	size_t np = ahabcValuesPtr->getNumberParameters();
	scale.resize(np);
	shift.resize(np);
	hm.resize(np);
	for (std::vector < std::vector < double > >::iterator it = hm.begin();
				it < hm.end(); ++it) {
		it->resize(np);			
					
	}
	
	for (size_t i = 0; i < nrefstats; ++i) {
		
		PolygonPiecewiseConstantFunction
			ppcfSlice (slices[i], hm, scale, shift);
	
		std::vector < double > slicePts = refset_ptr->at(i)->getValues();
	
		ostringstream oss;
		oss << pcfSliceUntransformedFilenameBase << "_";
		ostream_iterator<double> out_it(oss, "_");
		copy ( slicePts.begin(), slicePts.end(), out_it );
		oss << "slice_" << (i+1) << ".txt";
		// output slice
		string ppcfSliceFilename = oss.str();
		ppcfSlice.outputToTxtTabs(ppcfSliceFilename, prec, confirm);

	}
	
	
}

// Method to add a string to a file
// Output goes to file named according to argument s
// adds newline after line
void outputFileLine(const std::string& s, 
				const std::string& line, bool append)
{
	ofstream os;
	if (append) os.open(s.c_str(), ios::app);         // append
	else os.open(s.c_str()); // don't append
	if (os.is_open()) {
		os << line << endl;
		os.close();
	}
	else {
		std::cerr << "Error: could not open file named "
			<< s << std::endl << std::endl;
	}
}

// Method to add a string and integers to a file
// Output goes to file named according to argument s
// adds newline after line
void outputFileLineAndInts(const std::string& s, 
				const std::string& line, 
				const std::vector < int >& ints, 
				bool append)
{
	ofstream os;
	if (append) os.open(s.c_str(), ios::app);         // append
	else os.open(s.c_str()); // don't append
	if (os.is_open()) {
		os << line;
		for (size_t j = 0; j < ints.size(); ++j)
									os << "\t" << ints[j];
		os << endl;
		os.close();
	}
	else {
		std::cerr << "Error: could not open file named "
			<< s << std::endl << std::endl;
	}
}


// Method to add a string and doubles to a file
// Output goes to file named according to argument s
// adds newline after line
void outputFileLineAndDbls(const std::string& s, 
				const std::string& line, 
				const std::vector < double >& dbls,
				bool append)
{
	ofstream os;
	if (append) os.open(s.c_str(), ios::app);         // append
	else os.open(s.c_str()); // don't append
	if (os.is_open()) {
		os << line;
		for (size_t j = 0; j < dbls.size(); ++j)
									os << "\t" << dbls[j];
		os << endl;
		os.close();
	}
	else {
		std::cerr << "Error: could not open file named "
			<< s << std::endl << std::endl;
	}
}


// Method to add a string and doubles 
// and ints to a file
// Output goes to file named according to argument s
// adds newline after line
void outputFileLineAndDblsAndInts(const std::string& s, 
				const std::string& line, 
				const std::vector < double >& dbls,
				const std::vector < int >& ints,
				bool append)
{
	ofstream os;
	if (append) os.open(s.c_str(), ios::app);         // append
	else os.open(s.c_str()); // don't append
	if (os.is_open()) {
		os << line;
		for (size_t j = 0; j < dbls.size(); ++j)
									os << "\t" << dbls[j];
		for (size_t j = 0; j < ints.size(); ++j)
									os << "\t" << ints[j];
		os << endl;
		os.close();
	}
	else {
		std::cerr << "Error: could not open file named "
			<< s << std::endl << std::endl;
	}
}


//strip path from filename
std::string stripPath(const std::string& filename)
{
	/* Search for the last '/' in the scalarsFilenName, break
	 * it there*/

	std::string file = filename;

	size_t found;
	found=filename.find_last_of("/");
	if (found!=string::npos) {
		file = filename.substr(found+1);
	}
	return file;
}

/* output a vector of vectors of doubles to a file, with each inner vector
* as a column*/
void outputToFileVertical(
					const std::vector < std::vector < double > >& dataVec, 
					const std::string& filename,
					int prec)
{
	
	std::ofstream os;
	
	os.open(filename.c_str());         // overwrite 
	

	if (os.is_open()) {
		
		// use manipulators to set the precision for output
		os << fixed << setprecision(prec);
		
		// data
		/* find the longest inner vector, and 
		 * make a vector of iterators to the beginning of each inner vector */
		size_t longest = 0;
		std::vector < std::vector < double >::const_iterator > dataIterators;
		for (std::vector < std::vector < double > >::const_iterator it = dataVec.begin();
				it < dataVec.end(); ++it) {
			dataIterators.push_back( it->begin() );
			if (it->size() > longest) longest = it->size();
		}
		
		size_t n = dataIterators.size();
		
		// use this to output the data
		for (size_t j = 0; j < longest; ++j) {
			os << j+1;
			for (size_t i = 0; i < n; ++i) {
				if ( dataIterators[i] < dataVec[i].end() ) {
					os << "\t" << (*dataIterators[i]);
					++dataIterators[i];
				}
				else {
					os << "\t"; // not sure about this - may not be a good idea?
				}
				
			}
			// if it is not the last row, newline
			if ( j < longest - 1) os << std::endl;
		}
		
		os.close();
	}
	else {
		std::cerr << "Error: could not open file named "
		<< filename << std::endl << std::endl;
	}
	
}

void outputPointwiseExtensions(const std::string& filename,
			const std::vector < double >& trueParams,
			const std::vector < std::vector < double > >& otherParams,
			const std::vector < std::vector < double > >& values)
{
	std::ofstream os;
	
	os.open(filename.c_str());         // overwrite 
	

	if (os.is_open()) {
		
		os << "True parameters";
		for (size_t i = 0; i < trueParams.size(); ++i) 
								os << "\t" << trueParams[i]; 
		os << "\n" << endl;
		
		for (size_t j = 0; j < otherParams.size(); ++j) {
			os << "Other parameters" << (j+1);
			for (size_t i = 0; i < otherParams[j].size(); ++i) 
									os << "\t" << otherParams[j][i]; 
			os << "\n" << endl;
		}
		for (size_t j = 0; j < otherParams.size(); ++j)
				os << "\t" << (j+1); 
		os << endl;

		os << scientific; // formatting

		for (size_t j = 0; j < values.size(); ++j) {
			for (size_t i = 0; i < values[j].size(); ++i) {
				if (i) os << "\t" << values[j][i]; 
				else os << values[j][i]; 
				
			}
			os << endl;
		}
		
		os.close();
	}
	else {
		std::cerr << "Error: could not open file named "
		<< filename << std::endl << std::endl;
	}
	
}

void outputPointwiseExtensions(const std::string& filename,
			const std::vector < double >& values)
{
	std::ofstream os;
	
	os.open(filename.c_str());         // overwrite 
	

	if (os.is_open()) {
		
		ostream_iterator<double> out_it (os,"\n");
		copy ( values.begin(), values.end(), out_it );
	}
	else {
		std::cerr << "Error: could not open file named "
		<< filename << std::endl << std::endl;
	}
	
}




void buildGridPoint(std::vector < std::vector < double > >& gridPoints,
		const std::vector < double >& gridpoint,
		const std::vector < std::vector < double > >& priorGridSegments,
		size_t param)
{
	std::vector < double > copyGP(gridpoint);
	for (size_t i = 0; i < priorGridSegments[param].size(); ++i) {
		/* copy the current pt and send along */
		std::vector < double > copyGP(gridpoint);
		copyGP.push_back(priorGridSegments[param][i]);
		if (param + 1 < priorGridSegments.size()) {
			//send this same one along
			buildGridPoint(gridPoints, copyGP, priorGridSegments,
			param + 1); //pt
		}
		else {
			assert(copyGP.size() == priorGridSegments.size());
			gridPoints.push_back(copyGP);
		}
	}
	
}

double getDistance(
		const std::vector < double >& refs,
		const boost::shared_ptr < const mct::SummaryStatistic >& ss)
{	
	// note no checks for valid pointers here
	std::vector < double > values = ss->getValues();
	
	//  sum of (values-ref)^2
	
	using namespace boost::lambda;
	
	boost::lambda::placeholder1_type v1;
	boost::lambda::placeholder2_type v2;
	
	//size_t n = (stats.front())->size();	// expected size
	transform ( values.begin(), values.end(), refs.begin(), 
						values.begin(), v1 - v2 );
	transform ( values.begin(), values.end(), values.begin(), (v1*v1) );
	
	return std::sqrt(std::accumulate(values.begin(), values.end(), 0.0));
	
}


double getPercentile(const std::vector < double >& distances,
						double percentile) 
{
	std::vector < double > sortedvec(distances);
	sort(sortedvec.begin(), sortedvec.end());
	return ( gsl_stats_quantile_from_sorted_data(&sortedvec[0],
						1, // stride
						sortedvec.size(),
						percentile) );
}
