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

#ifndef _INC_TESTINGTOOLS_H
#define _INC_TESTINGTOOLS_H


#include "polygon_piecewise_function.hpp"
#include "piecewise_constant_function.hpp"

#include "ahabc_values.hpp"


#include <ivector.hpp>

#include <boost/shared_ptr.hpp>


#include <string>
#include <vector>

	void outputAHABCLogFileStart(const std::string& s);
	
	void outputAHABCLogFileString(const std::string& s, 
								const std::string& line, 
								bool append = true);
								
	void printHouseholderMatrix(
			boost::shared_ptr < const mct::AHABCValues  > ahabcValuesPtr);
	
	void printAndOutputMeansAndSDs(const std::vector < double >& means,
				const std::vector < double >& sds,
				const std::string& filename);
	
	void printMeansAndSDs(
			boost::shared_ptr < const mct::AHABCValues  > ahabcValuesPtr);
	
	void printMeansAndSDsAndHouseholderMatrixStandardised(
		boost::shared_ptr < const  mct::AHABCValues  > ahabcValuesPtr);

	void outputStats(const std::string& s,
						const std::vector < double >& stats);

	void printVarCovar(const subpavings::RealVec& varcovar, int dim);
	
	std::vector < std::vector < double > > getSlicePoints(
		boost::shared_ptr < mct::SummaryStatisticSet > refset_ptr);
	
	std::vector < double > standardise(
							const std::vector < double >& invec,
							const std::vector < double >& means,
							const std::vector < double >& sds);
	
	std::vector < double > getStandardisedParamSlicePoint(
		const std::vector < double >& paramPts,
		size_t np,
		size_t ns,
		const std::vector < double >& means,
		const std::vector < double >& sds);
	
	std::vector < double > getStandardisedSlicePoint(
		boost::shared_ptr < const mct::SummaryStatistic  > refStat,
		size_t np,
		const std::vector < double >& means,
		const std::vector < double >& sds);
		
	std::vector < std::vector < double > > getStandardisedSlicePoints(
		boost::shared_ptr < mct::SummaryStatisticSet > refset_ptr,
		size_t np,
		const std::vector < double >& means,
		const std::vector < double >& sds);
		
	std::vector < double > getStandardisedSlicePoint(
		boost::shared_ptr < const mct::SummaryStatistic  > refStat,
		boost::shared_ptr < const mct::AHABCValues  > ahabcValuesPtr);
		
	std::vector < double > getStandardisedParamSlicePoint(
		const std::vector < double >& paramPts,
		boost::shared_ptr < const mct::AHABCValues  > ahabcValuesPtr);
	
	std::vector < std::vector < double > > getStandardisedSlicePoints(
		boost::shared_ptr < mct::SummaryStatisticSet > refset_ptr,
		boost::shared_ptr < const mct::AHABCValues  > ahabcValuesPtr);
		
	std::vector < double > getStandardisedTransformedSlicePoint(
		boost::shared_ptr < const mct::SummaryStatistic  > refStat,
		boost::shared_ptr < const mct::AHABCValues  > ahabcValuesPtr);
	
	std::vector < std::vector < double > > getStandardisedTransformedSlicePoints(
		boost::shared_ptr < mct::SummaryStatisticSet > refset_ptr,
		boost::shared_ptr < const mct::AHABCValues  > ahabcValuesPtr);
	
	cxsc::ivector getBox(
			const std::vector < std::vector < double > >& container,
			const std::vector < std::vector < double > >& slicePtsVec,
			const std::vector < int >& sliceDims,
			cxsc::real padding);
	
			
	void processStandardisedTransformedAverage(
			const subpavings::PiecewiseConstantFunction&  pcf,
			boost::shared_ptr < const mct::AHABCValues  > ahabcValuesPtr,
			const std::string& pcfFilename,
			const std::string& pcfSliceFilename,
			const std::string& pcfUntransformedFilename,
			const std::string& pcfSliceUntransformedFilename,
			const std::vector < int >& sliceDims,
			const std::vector < double >& slicePts,
			int prec,
			bool confirm);
			
	void processAverage(
			const subpavings::PiecewiseConstantFunction&  pcf,
			boost::shared_ptr < const mct::AHABCValues  > ahabcValuesPtr,
			const std::string& prefix,
			const std::string& spec,
			const std::vector < int >& sliceDims,
			const std::vector < std::vector < double > >& slicePtsVec,
			const std::vector < double >& trueParams,
			int prec,
			bool confirm,
			size_t posteriorCaptureInterval = 0);

	std::vector < boost::shared_ptr < subpavings::PiecewiseConstantFunction > > 
		processAverage(
			const subpavings::PiecewiseConstantFunction&  pcf,
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
			size_t posteriorCaptureInterval);

	std::vector < boost::shared_ptr < 
		subpavings::PiecewiseConstantFunction > > 
	processStandardisedAverage(
			const subpavings::PiecewiseConstantFunction&  pcf,
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
			size_t posteriorCaptureInterval = 0);
			
	std::vector < boost::shared_ptr < 
		subpavings::PiecewiseConstantFunction > >
	processStandardisedAverage(
			const subpavings::PiecewiseConstantFunction&  pcf,
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
			const std::vector < std::vector < double > >& otherParams,
			int prec,
			bool confirm,
			size_t posteriorCaptureInterval = 0);
	
	//aabc
	void processStandardisedAverage(
			const subpavings::PiecewiseConstantFunction&  pcf,
			bool setAppend,
			int setIndex,
			const std::vector < double >& sds,
			const std::vector < double >& means,
			const std::string& prefix,
			const std::string& spec,
			int prec,
			bool confirm);

	//aabc
	void processAverage(
			const subpavings::PiecewiseConstantFunction&  pcf,
			bool setAppend,
			int setIndex,
			const std::string& prefix,
			const std::string& spec,
			int prec,
			bool confirm);
			
	void processStandardisedTransformedAverage(
			const subpavings::PiecewiseConstantFunction&  pcf,
			boost::shared_ptr < const mct::AHABCValues  > ahabcValuesPtr,
			const std::string& pcfFilename,
			const std::string& pcfUntransformedFilename,
			const std::string& pcfSliceFilenameBase,
			const std::string& pcfSliceUntransformedFilenameBase,
			const std::vector < int >& sliceDims,
			boost::shared_ptr < mct::SummaryStatisticSet > refset_ptr,
			const std::vector < std::vector < double > >& slicePtsVec,
			int prec,
			bool confirm);

	void processStandardisedPredictivePosteriorAverage(
			const subpavings::PiecewiseConstantFunction&  pcf,
			const std::vector < double >& scale,
			const std::vector < double >& shift,
			int np,
			int ns,
			const std::string& prefix,
			const std::string& spec,
			boost::shared_ptr < mct::SummaryStatisticSet > refset_ptr,
			const std::vector < std::vector < double > >& slicePtsVec,
			int prec,
			bool confirm);
	
	void processPredictivePosteriorAverage(
			const subpavings::PiecewiseConstantFunction&  pcf,
			int np,
			int ns,
			const std::string& prefix,
			const std::string& spec,
			boost::shared_ptr < mct::SummaryStatisticSet > refset_ptr,
			const std::vector < std::vector < double > >& slicePtsVec,
			int prec,
			bool confirm);
	
	void processStandardisedSummaryStatsAverage(
			const subpavings::PiecewiseConstantFunction&  pcf,
			size_t n_reps,
			const std::vector < double >& scale,
			const std::vector < double >& shift,
			int ns,
			const std::string& prefix,
			const std::string& spec,
			const std::vector < std::vector < double > >& slicePtsVec,
			int prec,
			bool confirm);
	
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
					bool confirm);
									
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
					bool confirm);
					
	void outputSummaryStatMarginals(
					const subpavings::PiecewiseConstantFunction& pcf,
					int ns,
					int np,
					const std::string& pcfMargFilenameBase,
					const std::string& pcfMargResizedFilenameBase,
					const std::string& prefix,
					const std::string& filenameBit,
					const std::vector < double>& scale,
					const std::vector < double>& shift,
					int prec,
					bool confirm);
	
	void outputSummaryStatMarginals(
					const subpavings::PiecewiseConstantFunction& pcf,
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
					bool confirm);
	
	void outputSummaryStatMarginals(
					const subpavings::PiecewiseConstantFunction& pcf,
					int ns,
					const std::string& pcfMargFilenameBase,
					const std::string& pcfMargResizedFilenameBase,
					const std::string& prefix,
					const std::string& filenameBit,
					const std::vector < double>& scale,
					const std::vector < double>& shift,
					int prec,
					bool confirm);
	
	void outputSummaryStatMarginals(
					const subpavings::PiecewiseConstantFunction& pcf,
					int ns,
					int np,
					const std::string& pcfMargFilenameBase,
					const std::string& prefix,
					const std::string& filenameBit,
					int prec,
					bool confirm);
	
	std::string setSliceBoxesName(const std::string sliceBoxesFilenameBase,
								const std::string sliceBoxesFilenames,
								const std::vector < double >& slicePts,
								size_t i,
								bool append);
												
	void outputSlice(const subpavings::PiecewiseConstantFunction& slice,
			const std::string& pcfSliceFilenameBase,
			const std::string& pcfSliceFilenames,
			size_t i,
			const std::vector < double >& slicePts,
			int prec,
			bool confirm,
			bool append);
	
	void outputPolygonSlice(	
			const subpavings::PolygonPiecewiseConstantFunction& ppcfSlice,
			const std::string& pcfSliceResizedFilenameBase,
			const std::string& pcfSliceResizedFilenames,
			size_t i,
			const std::vector < double >& slicePts,
			int prec,
			bool confirm,
			bool append);
	
	void outputMarginal(	const subpavings::PiecewiseConstantFunction& pcfMarg,
					const std::string& pcfMargFilenameBase,
					const std::string& pcfMargFilenames,
					const std::string& filenameBit,
					const std::vector < int >& margDims,
					int prec,
					bool confirm,
					bool append);
	
	void outputMarginal(	const subpavings::PiecewiseConstantFunction& pcfMarg,
					const std::string& pcfMargFilenameBase,
					const std::string& pcfMargFilenames,
					const std::string& filenameBit,
					const std::vector < int >& margDims,
					int setIndex,
					int prec,
					bool confirm,
					bool append);
					
	void outputMarginal(	const subpavings::PiecewiseConstantFunction& pcfMarg,
					const std::string& pcfMargFilenameBase,
					const std::string& pcfMargFilenames,
					const std::string& filenameBit,
					const std::vector < int >& margDims,
					const std::vector < double >& details,
					int prec,
					bool confirm,
					bool append);
					
	void outputMarginal(	const subpavings::PiecewiseConstantFunction& pcfMarg,
					const std::string& pcfMargFilenameBase,
					const std::string& pcfMargFilenames,
					const std::string& filenameBit,
					const std::vector < int >& margDims,
					const std::vector < double >& details,
					int setIndex,
					int prec,
					bool confirm,
					bool append);
	
	void outputPolygonMarginal(	const subpavings::PolygonPiecewiseConstantFunction& ppcfMarg,
					const std::string& pcfMargResizedFilenameBase,
					const std::string& pcfMargResizedFilenames,
					const std::string& filenameBit,
					const std::vector < int >& margDims,
					int prec,
					bool confirm,
					bool append);
					
	void outputPolygonMarginal(	const subpavings::PolygonPiecewiseConstantFunction& ppcfMarg,
					const std::string& pcfMargResizedFilenameBase,
					const std::string& pcfMargResizedFilenames,
					const std::string& filenameBit,
					const std::vector < int >& margDims,
					int setIndex,
					int prec,
					bool confirm,
					bool append);
	
	void outputPolygonMarginal(	const subpavings::PolygonPiecewiseConstantFunction& ppcfMarg,
					const std::string& pcfMargResizedFilenameBase,
					const std::string& pcfMargResizedFilenames,
					const std::string& filenameBit,
					const std::vector < int >& margDims,
					const std::vector < double >& details,
					int prec,
					bool confirm,
					bool append);
					
	void outputPolygonMarginal(	const subpavings::PolygonPiecewiseConstantFunction& ppcfMarg,
					const std::string& pcfMargResizedFilenameBase,
					const std::string& pcfMargResizedFilenames,
					const std::string& filenameBit,
					const std::vector < int >& margDims,
					const std::vector < double >& details,
					int setIndex,
					int prec,
					bool confirm,
					bool append);
	
	void outputParamSlice(	const subpavings::PiecewiseConstantFunction& slice,
					const std::string& pcfSliceFilenameBase,
					const std::vector < double >& slicePts,
					int prec,
					bool confirm);
	
	void outputPolygonParamSlice(	
			const subpavings::PolygonPiecewiseConstantFunction& ppcfSlice,
			const std::string& pcfSliceResizedFilenameBase,
			const std::vector < double >& slicePts,
			int prec,
			bool confirm);
			
		
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
		bool confirm);
	
	void outputCoverages(
			const subpavings::PiecewiseConstantFunction&  pcf,
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
			bool confirm);
				
	void outputCoverages(
			const subpavings::PiecewiseConstantFunction&  pcf,
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
			bool confirm);
	
	void outputCoverages(
			const subpavings::PiecewiseConstantFunction&  pcf,
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
			bool confirm);
	
	void outputCoverages(
			const subpavings::PiecewiseConstantFunction&  pcf,
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
			bool confirm);
							
	void outputFileLine(const std::string& s, 
							const std::string& line, bool append);
	
	void outputFileLineAndInts(const std::string& s, 
				const std::string& line, 
				const std::vector < int >& ints, 
				bool append);

	void outputFileLineAndDbls(const std::string& s, 
				const std::string& line, 
				const std::vector < double >& dbls,
				bool append);
	
	void outputFileLineAndDblsAndInts(const std::string& s, 
				const std::string& line, 
				const std::vector < double >& dbls,
				const std::vector < int >& ints,
				bool append);
	
	std::string stripPath(const std::string& filename);
	
	void outputToFileVertical(
					const std::vector < std::vector < double > >& dataVec, 
					const std::string& filename,
					int prec = 5);
	
	void outputPointwiseExtensions(const std::string& filename,
			const std::vector < double >& trueParams,
			const std::vector < std::vector < double > >& otherParams,
			const std::vector < std::vector < double > >& values);
	
	void outputPointwiseExtensions(const std::string& filename,
			const std::vector < double >& values);

	void buildGridPoint(std::vector < std::vector < double > >& gridpoints,
		const std::vector < double >& gridpoint,
		const std::vector < std::vector < double > >& priorGridSegments,
		size_t param);

	double getDistance(
			const std::vector < double >& refs,
			const boost::shared_ptr < const mct::SummaryStatistic >& ss);

	double getPercentile(const std::vector < double >& distances,
						double percentile);


	
#endif
