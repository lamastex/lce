/* MCT - Markov Chains on Trees.

   Copyright (C) 2012 Jenny Harlow

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
\brief AHABCValues definitions.
*/

#include "ahabc_values.hpp"

//#include <gsl/gsl_statistics_double.h>

//#include <cmath>
//#include <limits> 
#include <algorithm>
#include <functional>
#include <numeric>
#include <stdexcept>
//#include <iomanip>
//#include <sstream>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>



//#define MYDEBUG
#define MYDEBUG_EIGEN
#define MYDEBUG_SCATTER
#if defined (MYDEBUG) || defined (MYDEBUG_EIGEN) || defined (MYDEBUG_SCATTER)
	#include <iostream>
	#include <iterator>
#endif

using namespace mct;


AHABCValues::AHABCValues(
			const boost::shared_ptr< const aabc::ParameterParticleSet > p_ptr,
			const boost::shared_ptr < const SummaryStatisticSet > s_ptr)
			: ss_ptr(s_ptr)
{
	try {
		
		if (!s_ptr || s_ptr->empty()) {
			throw std::invalid_argument("s_ptr null pointer or empty set");
		}
		if (!p_ptr || p_ptr->empty()) {
			throw std::invalid_argument("p_ptr null pointer or empty set");
		}
		if (p_ptr->size() != s_ptr->size())
			throw std::invalid_argument("p_ptr and s_ptr sets of different sizes set");
		
		if ( !p_ptr->sizeConsistencyCheck() ) {
				throw std::invalid_argument("inconsistent particle sizes");
		}
		if ( !s_ptr->sizeConsistencyCheck() ) {
				throw std::invalid_argument("inconsistent summary statistic sizes");
		}
		ps_ptr = p_ptr->makeSummaryStatisticSet();
	}
	catch (std::invalid_argument& ia) {
	
		throw std::invalid_argument(
			std::string("AHABCValues::AHABCValues(...) :\n") + ia.what());
	}
}

AHABCValues::AHABCValues(
			const boost::shared_ptr< const aabc::ParameterParticleSet > p_ptr,
			const std::vector <size_t>& reqParams,
			const boost::shared_ptr < const SummaryStatisticSet > s_ptr,
			const std::vector <size_t>& reqStats)
{
	try {
		
		if (!s_ptr || s_ptr->empty()) {
			throw std::invalid_argument("s_ptr null pointer or empty set");
		}
		if (!p_ptr || p_ptr->empty()) {
			throw std::invalid_argument("p_ptr null pointer or empty set");
		}
		if (p_ptr->size() != s_ptr->size())
			throw std::invalid_argument("p_ptr and s_ptr sets of different sizes set");
		
		if ( !p_ptr->sizeConsistencyCheck() ) {
				throw std::invalid_argument("inconsistent particle sizes");
		}
		if ( !s_ptr->sizeConsistencyCheck() ) {
				throw std::invalid_argument("inconsistent summary statistic sizes");
		}
		if ( reqParams.empty() ) {
				throw std::invalid_argument("reqParams.empty()");
		}
		if ( reqStats.empty() ) {
				throw std::invalid_argument("reqStats.empty()");
		}
		
		boost::shared_ptr < const SummaryStatisticSet > 
				ps_ptr_tmp = p_ptr->makeSummaryStatisticSet();
		ps_ptr = ps_ptr_tmp->makeSelectedSummaryStatisticSet(reqParams);
		ss_ptr = s_ptr->makeSelectedSummaryStatisticSet(reqStats);
	}
	catch (std::invalid_argument& ia) {
	
		throw std::invalid_argument(
			std::string("AHABCValues::AHABCValues(...) :\n") + ia.what());
	}
}

AHABCValues::AHABCValues(
			const boost::shared_ptr< const aabc::ParameterParticleSet > p_ptr,
			const boost::shared_ptr < const SummaryStatisticSet > s_ptr,
			const std::vector <size_t>& reqStats)
{
	try {
		
		if (!s_ptr || s_ptr->empty()) {
			throw std::invalid_argument("s_ptr null pointer or empty set");
		}
		if (!p_ptr || p_ptr->empty()) {
			throw std::invalid_argument("p_ptr null pointer or empty set");
		}
		if (p_ptr->size() != s_ptr->size())
			throw std::invalid_argument("p_ptr and s_ptr sets of different sizes set");
		
		if ( !p_ptr->sizeConsistencyCheck() ) {
				throw std::invalid_argument("inconsistent particle sizes");
		}
		if ( !s_ptr->sizeConsistencyCheck() ) {
				throw std::invalid_argument("inconsistent summary statistic sizes");
		}
		if ( reqStats.empty() ) {
				throw std::invalid_argument("reqStats.empty()");
		}
		
		ps_ptr = p_ptr->makeSummaryStatisticSet();
		ss_ptr = s_ptr->makeSelectedSummaryStatisticSet(reqStats);
	}
	catch (std::invalid_argument& ia) {
	
		throw std::invalid_argument(
			std::string("AHABCValues::AHABCValues(...) :\n") + ia.what());
	}
}

AHABCValues::AHABCValues(
			const boost::shared_ptr< const aabc::ParameterParticleSet > p_ptr,
			const std::vector <size_t>& reqParams,
			const boost::shared_ptr < const SummaryStatisticSet > s_ptr)
			: ss_ptr(s_ptr)
{
	try {
		
		if (!s_ptr || s_ptr->empty()) {
			throw std::invalid_argument("s_ptr null pointer or empty set");
		}
		if (!p_ptr || p_ptr->empty()) {
			throw std::invalid_argument("p_ptr null pointer or empty set");
		}
		if (p_ptr->size() != s_ptr->size())
			throw std::invalid_argument("p_ptr and s_ptr sets of different sizes set");
		
		if ( !p_ptr->sizeConsistencyCheck() ) {
				throw std::invalid_argument("inconsistent particle sizes");
		}
		if ( !s_ptr->sizeConsistencyCheck() ) {
				throw std::invalid_argument("inconsistent summary statistic sizes");
		}
		if ( reqParams.empty() ) {
				throw std::invalid_argument("reqParams.empty()");
		}
		
		boost::shared_ptr < const SummaryStatisticSet > 
				ps_ptr_tmp = p_ptr->makeSummaryStatisticSet();
		ps_ptr = ps_ptr_tmp->makeSelectedSummaryStatisticSet(reqParams);
	}
	catch (std::invalid_argument& ia) {
	
		throw std::invalid_argument(
			std::string("AHABCValues::AHABCValues(...) :\n") + ia.what());
	}
}
size_t AHABCValues::size() const
{
	return ps_ptr->size();
}

size_t AHABCValues::getNumberParameters() const
{
	return ps_ptr->at(0)->size();
}

size_t AHABCValues::getNumberSummaryStatistics() const
{
	return ss_ptr->at(0)->size();
}

std::vector < int >  AHABCValues::getParametersDims() const
{
	std::vector < int > pos;
	//positions go in from 1 upwards
	int pars = static_cast<int>( getNumberParameters() );
	for (int i = 1; i <= pars; ++i) pos.push_back(i);
	
	return pos;
}

std::vector < int >  AHABCValues::getSummaryStatisticsDims() const
{
	std::vector < int > pos;
	int pars = static_cast<int>( getNumberParameters() );
	int stats = static_cast<int>( getNumberSummaryStatistics() );
	//positions go in from number of parameters + 1 upwards
	for (int i = 1; i <= stats; ++i) 
			pos.push_back(i+pars);
	
	return pos;
}

void AHABCValues::fillDataContainer(
		std::vector < std::vector < double > >& container) const
{
	boost::shared_ptr < std::vector < std::vector < double > > > dataPtr
	= extractData();
	
	container = *dataPtr;
	
} 

void AHABCValues::fillStandardisedDataContainer(
		std::vector < std::vector < double > >& container) const
{
	boost::shared_ptr < std::vector < std::vector < double > > > dataPtr
	= extractStandardisedData();
	
	container = *dataPtr;
	
} 

void AHABCValues::fillTransformedDataContainer(
		std::vector < std::vector < double > >& container) const
{
	boost::shared_ptr < std::vector < std::vector < double > > > 
				transformedDataPtr = extractTransformedData();
	container = *transformedDataPtr;
	
}

void AHABCValues::fillTransformedStandardisedDataContainer(
		std::vector < std::vector < double > >& container) const
{
	boost::shared_ptr < std::vector < std::vector < double > > > 
				transformedDataPtr = extractTransformedStandarisedData();
	container = *transformedDataPtr;
	
}

std::vector < std::vector < double > > 
				AHABCValues::getHouseholderMatrix() const
{
	createHouseholderMatrix();
	return householderMatrix;
} 

std::vector < std::vector < double > > 
				AHABCValues::getHouseholderMatrixStandardised() const
{
	createHouseholderMatrixStandardised();
	return householderMatrixStandardised;
} 

/* standardise */
std::vector < double > AHABCValues::standardise(
							const std::vector < double >& invec) const
{
	getMeans();
	getSDs();
	
	size_t d = means.size();
	
	if (invec.size() != d) throw std::invalid_argument(
		"AHABCValues::standardised(const std::vector < double >&)");
	
	
	std::vector < double > standardised = invec;
	for (size_t i = 0; i < d; ++i) {
		standardised[i] = (invec[i] - means[i]) / (sds[i]);
	}
	
	
	return standardised;
}


std::vector < double > AHABCValues::householderTransform(
							const std::vector < double >& invec) const
{
	createHouseholderMatrix();
	
	size_t d = householderMatrix.size();
	
	if (invec.size() != d) throw std::invalid_argument(
		"AHABCValues::householderTransform(const std::vector < double >&)");
	
	std::vector < double > transformed(d);
	
	for (size_t j = 0; j < d; ++j) { // row index of new data elements and H
			
		/* inner dot product of this data and jth row of H
		 * matrix is the value at jth element of new data */
		transformed[j] = 
			std::inner_product(invec.begin(), invec.end(),
							householderMatrix[j].begin(), 
							0.0);
			
			
	}
	return transformed;
}

/* standardise and then apply householder transformation for 
standardised data */
std::vector < double > AHABCValues::standardiseAndHouseholderTransform(
							const std::vector < double >& invec) const
{
	createHouseholderMatrixStandardised();
	
	size_t d = householderMatrixStandardised.size();
	
	if (invec.size() != d) throw std::invalid_argument(
		"AHABCValues::standardisedHouseholderTransform(const std::vector < double >&)");
	
	getMeans();
	getSDs();
	
	std::vector < double > st_invec = invec;
	for (size_t i = 0; i < d; ++i) {
		st_invec[i] = (invec[i] - means[i]) / (sds[i]);
	}
	
	
	std::vector < double > transformed(d);
	
	for (size_t j = 0; j < d; ++j) { // row index of new data elements and H
			
		/* inner dot product of this data and jth row of H
		 * matrix is the value at jth element of new data */
		transformed[j] = 
			std::inner_product(st_invec.begin(), st_invec.end(),
							householderMatrixStandardised[j].begin(), 
							0.0);
			
			
	}
	return transformed;
}


boost::shared_ptr < std::vector < std::vector < double > > > 
							AHABCValues::extractData() const
{
	boost::shared_ptr < std::vector < std::vector < double > > >
		dataPtr(new std::vector < std::vector < double > >());
		
	size_t n = size();
	
	dataPtr->reserve(n);
		
	for (size_t i = 0; i < n; ++i) {
		
		dataPtr->push_back( ps_ptr->at(i)->getValues() );
		
		std::vector < double > tmpvec = ss_ptr->at(i)->getValues();
		
		// put all the summary stats in at the end
		dataPtr->at(i).insert(dataPtr->at(i).end(),	tmpvec.begin(), tmpvec.end());
		
	}
	
	return dataPtr;
	
} 

boost::shared_ptr < std::vector < std::vector < double > > > 
				AHABCValues::extractStandardisedData() const
{
	
	std::vector < double > ppMeans;
	std::vector < double > ppSds;
	
	std::vector < double > ssMeans;
	std::vector < double > ssSds;
	
	boost::shared_ptr < const SummaryStatisticSet > ps_std_ptr
		= ps_ptr->makeSummaryStatisticSetStandardised(ppMeans, ppSds);
	
	boost::shared_ptr < const SummaryStatisticSet > ss_std_ptr
		= ss_ptr->makeSummaryStatisticSetStandardised(ssMeans, ssSds);
	
	size_t n = size();
	size_t d = getNumberParameters() + getNumberSummaryStatistics();
	
	means = ppMeans;
		
	means.insert( means.end(), ssMeans.begin(), ssMeans.end() );
	assert(means.size() == d);
	
	sds = ppSds;
		
	sds.insert( sds.end(), ssSds.begin(), ssSds.end() );
	assert(sds.size() == d);
	
	boost::shared_ptr < std::vector < std::vector < double > > >
		dataPtr(new std::vector < std::vector < double > >());
		
	dataPtr->reserve(n);
		
	for (size_t i = 0; i < n; ++i) {
		
		dataPtr->push_back( ps_std_ptr->at(i)->getValues() );
		
		std::vector < double > tmpvec = ss_std_ptr->at(i)->getValues();
		
		// put all the standardised summary stats in at the end
		dataPtr->at(i).insert(dataPtr->at(i).end(),	tmpvec.begin(), tmpvec.end());
		
	}
	
	return dataPtr;
} 

boost::shared_ptr < std::vector < std::vector < double > > > 
				AHABCValues::extractTransformedData() const
{
	boost::shared_ptr < std::vector < std::vector < double > > >
		dataPtr = extractData();
	
	boost::shared_ptr < std::vector < std::vector < double > > >
		transformedDataPtr = transformData(dataPtr);
		
	return transformedDataPtr;
	
} 

boost::shared_ptr < std::vector < std::vector < double > > > 
				AHABCValues::extractTransformedStandarisedData() const
{
	boost::shared_ptr < std::vector < std::vector < double > > >
		dataPtr = extractStandardisedData();
	
	boost::shared_ptr < std::vector < std::vector < double > > >
		transformedDataPtr = transformStandardisedData(dataPtr);
		
	return transformedDataPtr;
	
} 

std::vector < double > AHABCValues::getMeans() const
{

	size_t d = getNumberParameters() + getNumberSummaryStatistics();
	if (means.empty()) {
				
		boost::shared_ptr < const SummaryStatistic > paramMeansSSPtr = 
								ps_ptr->getSummaryStatisticMeans();
		means = paramMeansSSPtr->getValues();
		
		boost::shared_ptr < const SummaryStatistic > ssMeansSSPtr = 
							ss_ptr->getSummaryStatisticMeans();
		std::vector<double> ssMeans = ssMeansSSPtr->getValues();
		
		means.insert( means.end(), ssMeans.begin(), ssMeans.end() );
		
		
	}
	assert(means.size() == d);
	
	#ifdef MYDEBUG
		std::cout << "\nmeans are = " << std::endl;
		{
			std::ostream_iterator<double> out_it (std::cout,"\t");
			copy ( means.begin(), means.end(), out_it );
			std::cout << std::endl;
		}
		std::cout << std::endl;
	#endif
	
	return means;
}

std::vector < double > AHABCValues::getSDs() const
{

	size_t d = getNumberParameters() + getNumberSummaryStatistics();
	if (sds.empty()) {
				
		boost::shared_ptr < const SummaryStatistic > paramSDsSSPtr = 
								ps_ptr->getSummaryStatisticSDs();
		sds = paramSDsSSPtr->getValues();
		
		boost::shared_ptr < const SummaryStatistic > ssSDsSSPtr = 
							ss_ptr->getSummaryStatisticSDs();
		std::vector<double> ssSDs = ssSDsSSPtr->getValues();
		
		sds.insert( sds.end(), ssSDs.begin(), ssSDs.end() );
		
		
	}
	assert(sds.size() == d);
	
	#ifdef MYDEBUG
		std::cout << "\nsds are = " << std::endl;
		{
			std::ostream_iterator<double> out_it (std::cout,"\t");
			copy ( sds.begin(), sds.end(), out_it );
			std::cout << std::endl;
		}
		std::cout << std::endl;
	#endif
	
	return sds;
}

void AHABCValues::createHouseholderMatrix() const
{
	if (householderMatrix.empty()) {
		boost::shared_ptr < std::vector < std::vector < double > > > dataPtr
		= extractData();
		
		createHouseholderMatrix(dataPtr);
	}
}

void AHABCValues::createHouseholderMatrix(
		boost::shared_ptr < std::vector < std::vector < double > > > dataPtr)
														const
{
	if (householderMatrix.empty()) {
		
		std::vector < double > meansForDataPtr = getMeans();
		
		householderMatrix = createHouseholderMatrix(dataPtr, meansForDataPtr);
	}
}


void AHABCValues::createHouseholderMatrixStandardised() const
{
	if (householderMatrixStandardised.empty()) {
		boost::shared_ptr < std::vector < std::vector < double > > > 
			dataPtr = extractStandardisedData();
		
		createHouseholderMatrixStandardised(dataPtr);
	}
}

void AHABCValues::createHouseholderMatrixStandardised(
	boost::shared_ptr < std::vector < std::vector < double > > > dataPtr)
																	const
{
	if (householderMatrixStandardised.empty()) {
		
		size_t d = getNumberParameters() + getNumberSummaryStatistics();
		std::vector < double > meansForDataPtr (d, 0.0);
		
		householderMatrixStandardised 
				= createHouseholderMatrix(dataPtr, meansForDataPtr);
	}
}

std::vector < std::vector < double > > 
		AHABCValues::createHouseholderMatrix(boost::shared_ptr < 
				std::vector < std::vector < double > > > dataPtr,
				const std::vector < double >& meansForDataPtr) const
{
	
	size_t n = size();
	size_t d = getNumberParameters() + getNumberSummaryStatistics();
	
	// vector version
	std::vector<double> scatter; // will be size 1/2 d(d+1)
	for (size_t j = 0; j < d; ++j) {
		for (size_t k = j; k < d; ++k) scatter.push_back(0.0);
	}
	
	for (size_t i = 0; i < n; ++i) {
		std::vector<double> thisData = dataPtr->at(i);
		#ifdef MYDEBUG
			std::cout << "\nData is = " << std::endl;
			{
				std::ostream_iterator<double> out_it (std::cout,"\t");
				copy ( thisData.begin(), thisData.end(), out_it );
				std::cout << std::endl;
			}
			std::cout << std::endl;
		#endif
		
		// subtract means
		std::transform(thisData.begin(), thisData.end(),
						meansForDataPtr.begin(), thisData.begin(), 
						std::minus<double>());
		
		#ifdef MYDEBUG
			std::cout << "\nData centred = " << std::endl;
			{
				std::ostream_iterator<double> out_it (std::cout,"\t");
				copy ( thisData.begin(), thisData.end(), out_it );
				std::cout << std::endl;
			}
			std::cout << std::endl;
		#endif
		size_t index = 0;
		for (size_t j = 0; j < d; ++j) {
			for (size_t k = j; k < d; ++k) {
				// accumulate values
				scatter[index]+=(thisData[j]*thisData[k]);
				++index;
			}
		}
		
	}
	#ifdef MYDEBUG_SCATTER
		std::cout << "\nScatter matrix is = " << std::endl;
		{
			std::ostream_iterator<double> out_it (std::cout,"\t");
			copy ( scatter.begin(), scatter.end(), out_it );
			std::cout << std::endl;
		}
		std::cout << std::endl;
	#endif

	gsl_matrix * m_scatter = gsl_matrix_alloc (d, d);
	gsl_matrix_set_zero(m_scatter);

	{
		size_t index = 0;
		for (size_t i = 0; i < d; ++i) {
			for (size_t j = i; j < d; ++j) {
				gsl_matrix_set (m_scatter, i, j, scatter[index]);
				
				if (j > i) gsl_matrix_set (m_scatter, j, i, scatter[index]);
				
				++index;

		   }
	   }
	}
	
	#ifdef MYDEBUG
	{
		std::cout << "\nm_scatter is = " << std::endl;
		for (size_t i = 0; i < d; ++i) {
			for (size_t j = 0; j < d; ++j) 
				std::cout << gsl_matrix_get (m_scatter, i, j) << "\t";
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
	#endif
	
	
	// Create the vectors to store the eigenvalues and eigenvectors
	gsl_vector* evals = gsl_vector_alloc(d);
	gsl_matrix* evecs = gsl_matrix_alloc(d, d);

	// Allocate an eigenvector / eigenvalue workspace
	gsl_eigen_symmv_workspace* workspace = gsl_eigen_symmv_alloc(d);

	// Compute the eigenvalues and eigenvectors
	gsl_eigen_symmv(m_scatter, evals, evecs, workspace);

	//sort
	gsl_eigen_gensymmv_sort (evals, evecs, GSL_EIGEN_SORT_ABS_DESC);
	
	
	#ifdef MYDEBUG
	
	// Print out the eigenvalues
	std::cout << "\nSorted eigenvalues are = " << std::endl;
	for (size_t i = 0; i < d; i++) {
		std::cout << gsl_vector_get(evals, i) << "\t";
	}
	std::cout << std::endl;
	
	std::cout << "\nSorted eigenvectors are = " << std::endl;
	for (size_t i = 0; i < d; ++i) {
		for (size_t j = 0; j < d; ++j) {
		
			std::cout << gsl_matrix_get(evecs, i, j) << "\t";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
	#endif
	#ifdef MYDEBUG_EIGEN
	// Print out the leading eigenvector
	std::cout << "\nLeading eigenvector (as a row vector) is = " << std::endl;
	for (size_t i = 0; i < d; ++i)
	{
		std::cout << gsl_matrix_get(evecs, i, 0) << "\t";
	}
	std::cout << std::endl;
	#endif

	gsl_vector * u = gsl_vector_alloc(d);
	
	gsl_matrix_get_col(u, evecs, 0);
	
	
	// compute u = e1 - v/||e1 - v||
	gsl_vector_scale(u, -1.0); // -v
	gsl_vector_set(u, 0, gsl_vector_get(u, 0) + 1.0); // -v[0] + 1
	double length = gsl_blas_dnrm2(u);
	#ifdef MYDEBUG
	
		// Print out the unnormalised u
		std::cout << "\nunnormalised u (as a row vector) is = " << std::endl;
		for (size_t i = 0; i < d; ++i)
		{
			std::cout << gsl_vector_get(u, i) << "\t";
		}
		std::cout << std::endl;
		std::cout << "\nlength of u is = " << length << std::endl;
	#endif
	
	gsl_vector_scale(u, (1/length)); // normalise
	
	#ifdef MYDEBUG
	
		// Print out the unnormalised u
		std::cout << "\nnormalised u (as a row vector)is = " << std::endl;
		for (size_t i = 0; i < d; ++i)
		{
			std::cout << gsl_vector_get(u, i) << "\t";
		}
		std::cout << std::endl;
	#endif
	
	//identity matrix
	gsl_matrix * dI = gsl_matrix_alloc (d, d);
	gsl_matrix_set_identity (dI);

	// compute H
	gsl_matrix * H = gsl_matrix_alloc(d, d);
	
	for (size_t i = 0; i < d; ++i) {
		double u_el_scaled = -2.0 * gsl_vector_get(u,i);
		for (size_t j = 0; j < d; ++j) {
	
			gsl_matrix_set(H, i, j, u_el_scaled * gsl_vector_get(u,j));
		}
	}
	

	#ifdef MYDEBUG
	std::cout << "\nH before adding I is = " << std::endl;
	for (size_t i = 0; i < d; ++i) {
		for (size_t j = 0; j < d; ++j) {
		
			std::cout << gsl_matrix_get(H, i, j) << "\t";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
	#endif

	gsl_matrix_add(H, dI);
	
	#ifdef MYDEBUG
	std::cout << "\nH is = " << std::endl;
	for (size_t i = 0; i < d; ++i) {
		for (size_t j = 0; j < d; ++j) {
		
			std::cout << gsl_matrix_get(H, i, j) << "\t";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
	#endif
	
	// use these values to make our householder matrix
	
	std::vector < std::vector < double > > 
				hm(d, std::vector < double> (d, 0.0));
	
	for (size_t i = 0; i < d; ++i) {
		for (size_t j = 0; j < d; ++j) {
		
			hm[i][j] = gsl_matrix_get(H, i, j);
		}
	}
	#ifdef MYDEBUG
		std::cout << "\nfinal Householder matrix now holds " << std::endl;
		for (std::vector < std::vector < double > >::iterator it 
							= householderMatrix.begin();
							it < householderMatrix.end(); ++it) {
			std::ostream_iterator<double> out_it (std::cout,"\t");
			copy ( it->begin(), it->end(), out_it );
			std::cout << std::endl;
		}
		std::cout << std::endl;
	#endif
	
	// Free the gsl memory
	try {
		gsl_eigen_symmv_free(workspace);
		gsl_vector_free(evals);
		gsl_matrix_free(evecs);
		gsl_vector_free(u);
		
		gsl_matrix_free ( m_scatter);
		gsl_matrix_free ( dI);
		gsl_matrix_free ( H);
	}
	catch (...) {} // catch and swallow
	
	return hm;

}

boost::shared_ptr < std::vector < std::vector < double > > >
AHABCValues::transformData(
	boost::shared_ptr < std::vector < std::vector < double > > > dataPtr) const
{

	createHouseholderMatrix(dataPtr);
		
	size_t d = householderMatrix.size();
	
	assert(getNumberParameters() + getNumberSummaryStatistics() == d);
	
	/* premultiply each row of the data (using it as 
	 * a column vector in the calculation) by the Householder matrix */
	
	size_t n = size();
	
	std::vector < double > tmp(d, 0.0);
	boost::shared_ptr < std::vector < std::vector < double > > > 
		transformedDataPtr(new std::vector < std::vector < double > >(n, tmp));
	
	for (size_t i = 0; i < n; ++i ) { // row index for data
	
		std::vector < double > thisData = dataPtr->at(i);
		
		for (size_t j = 0; j < d; ++j) { // row index of new data elements and H
		
			/* add inner dot product of this data and jth row of H
			 * matrix to the 0.0 value at jth element of new data row */
			transformedDataPtr->at(i)[j] = 
				std::inner_product(thisData.begin(), thisData.end(),
									householderMatrix[j].begin(), 
									0.0);
		}
		
	} 

	#ifdef MYDEBUG
		std::cout << "\ntransformed data is " << std::endl;
		for (std::vector < std::vector < double > >::iterator it 
							= transformedDataPtr->begin();
							it < transformedDataPtr->end(); ++it) {
			std::ostream_iterator<double> out_it (std::cout,"\t");
			copy ( it->begin(), it->end(), out_it );
			std::cout << std::endl;
		}
		std::cout << std::endl;
	#endif
	
	return transformedDataPtr;
}


boost::shared_ptr < std::vector < std::vector < double > > >
AHABCValues::transformStandardisedData(
	boost::shared_ptr < std::vector < std::vector < double > > > dataPtr) const
{

	size_t d = getNumberParameters() + getNumberSummaryStatistics();
	
	createHouseholderMatrixStandardised(dataPtr);
		
	assert(householderMatrixStandardised.size() == d);
	
	/* premultiply each row of the data (using it as 
	 * a column vector in the calculation) by the Householder matrix */
	
	size_t n = size();
	
	std::vector < double > tmp(d, 0.0);
	boost::shared_ptr < std::vector < std::vector < double > > > 
		transformedDataPtr(new std::vector < std::vector < double > >(n, tmp));
	
	for (size_t i = 0; i < n; ++i ) { // row index for data
	
		std::vector < double > thisData = dataPtr->at(i);
		
		for (size_t j = 0; j < d; ++j) { // row index of new data elements and H
		
			/* add inner dot product of this data and jth row of H
			 * matrix to the 0.0 value at jth element of new data row */
			transformedDataPtr->at(i)[j] = 
				std::inner_product(thisData.begin(), thisData.end(),
									householderMatrixStandardised[j].begin(), 
									0.0);
		}
		
	} 

	#ifdef MYDEBUG
		std::cout << "\ntransformed standardised data is " << std::endl;
		for (std::vector < std::vector < double > >::iterator it 
							= transformedDataPtr->begin();
							it < transformedDataPtr->end(); ++it) {
			std::ostream_iterator<double> out_it (std::cout,"\t");
			copy ( it->begin(), it->end(), out_it );
			std::cout << std::endl;
		}
		std::cout << std::endl;
	#endif
	
	return transformedDataPtr;
}


void AHABCValues::dataSVD() const
{
	
	boost::shared_ptr < std::vector < std::vector < double > > > dataPtr
		= extractData();

	size_t n = size();
	size_t d = getNumberParameters() + getNumberSummaryStatistics();
	
	
	
	boost::shared_ptr < const SummaryStatistic > paramMeansSSPtr = 
						ps_ptr->getSummaryStatisticMeans();
	std::vector<double> means = paramMeansSSPtr->getValues();
	
	boost::shared_ptr < const SummaryStatistic > ssMeansSSPtr = 
						ss_ptr->getSummaryStatisticMeans();
	std::vector<double> ssMeans = ssMeansSSPtr->getValues();
	
	means.insert( means.end(), ssMeans.begin(), ssMeans.end() );
	assert(means.size() == d);
	
	#ifdef MYDEBUG
		std::cout << "\nmeans are = " << std::endl;
		{
			std::ostream_iterator<double> out_it (std::cout,"\t");
			copy ( means.begin(), means.end(), out_it );
			std::cout << std::endl;
		}
		std::cout << std::endl;
	#endif

	// data matrix has n rows and d columns, ie each x as a row vec
	gsl_matrix * data_matrix = gsl_matrix_alloc(n, d);
	
	/* could maybe do this more quickly
	 * - dump data to an array [] and use gsl_matrix_view? but... */
	for (size_t i = 0; i < n; ++i) {
		
		for (size_t j = 0; j < d; ++j) {
			gsl_matrix_set(data_matrix, i, j, 
				dataPtr->at(i)[j] - means[j]);
		}
	}
	
	gsl_matrix * v = gsl_matrix_alloc(d, d);
	gsl_matrix * workspace_x = gsl_matrix_alloc(d, d);
	gsl_vector * work = gsl_vector_alloc(d);
	gsl_vector * sigma = gsl_vector_alloc(d);
	
	gsl_linalg_SV_decomp_mod (data_matrix, workspace_x, v, sigma, work);

	//This function computes the SVD using the modified Golub-Reinsch algorithm, which is faster for M>>N. It requires the vector work of length N and the N-by-N matrix X as additional working space.  

	std::cout << "sigma is : " << std::endl; 
	{
		for (size_t i = 0; i < d; ++i) {
			std::cout << gsl_vector_get(sigma, i) << "\t";
	   }
	   std::cout << std::endl;
	}
	
	
	{
		std::cout << "\nv is = " << std::endl;
		for (size_t i = 0; i < d; ++i) {
			for (size_t j = 0; j < d; ++j) 
				std::cout << gsl_matrix_get (v, i, j) << "\t";
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}

}
