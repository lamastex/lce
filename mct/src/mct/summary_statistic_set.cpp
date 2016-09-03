/* MCT - Markov Chains on Trees.

   Copyright (C) 2011 Jenny Harlow

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
\brief SummaryStatisticSet definitions.
*/

#include "summary_statistic_set.hpp"

#include <gsl/gsl_math.h> //gsl_isnan()

#include <boost/lambda/bind.hpp>
#include <boost/lambda/lambda.hpp>

#include <cmath>
#include <limits>
#include <algorithm> // transform, swap
#include <iterator> // distance
#include <sstream>  // to be able to manipulate strings as streams
#include <iostream>
#include <fstream>
#include <stdexcept>


using namespace mct;


SummaryStatisticSet::SummaryStatisticSet() {}

SummaryStatisticSet::SummaryStatisticSet(size_t reservedSpace)
{
		container.reserve(reservedSpace); 
}

SummaryStatisticSet::SummaryStatisticSet(
				const SummaryStatisticSet& other)
		: container(other.container) {}

SummaryStatisticSet& SummaryStatisticSet::operator=(SummaryStatisticSet tmp)
{
	std::swap( container, tmp.container );

	return *this;
}

SummaryStatisticSet::~SummaryStatisticSet()
{}
			
const boost::shared_ptr < const SummaryStatistic >& SummaryStatisticSet::at(
										const size_t index) const
{
	return container.at(index);
}

void SummaryStatisticSet::add(
			const boost::shared_ptr < const SummaryStatistic >& sptr)
{
	container.push_back(sptr);
}

void SummaryStatisticSet::add(
			const boost::shared_ptr < SummaryStatisticSet >& sptr)
{
	std::vector< boost::shared_ptr < const SummaryStatistic > > otherContainer
		= sptr->_getContainer();
	container.insert( container.end(),
					otherContainer.begin(), otherContainer.end() );
	
}

size_t SummaryStatisticSet::size() const
{
	return container.size();
}

bool SummaryStatisticSet::empty() const
{
	return container.empty();
}

bool SummaryStatisticSet::sizeConsistencyCheck() const 
{
	bool result = true;
	if (!empty()) {
		
		SummaryStatisticSet::const_iterator it = begin();
		size_t n = (*it)->size();
		++it;				
		while ( result && it < end()) {
			
			if ((*it)->size() != n) {
				result = false;
			}
			++it;
		}
	}
	return result;
}

boost::shared_ptr < SummaryStatisticSet >
		SummaryStatisticSet::makeSummaryStatisticSetStandardised() const
{
	std::vector < double > means;
	std::vector < double > sds;
	
	return makeSummaryStatisticSetStandardised(means, sds);

}


boost::shared_ptr < SummaryStatisticSet >
		SummaryStatisticSet::makeSummaryStatisticSetStandardised(
				std::vector < double >& means,
				std::vector < double >& sds) const
{
	//get all the standardised values as vec of vecs
	boost::shared_ptr < std::vector < std::vector < double > > > 
				stanvals = getAllValuesStandardised(means, sds);
	
	size_t nss = size();
		
	// empty set with enough reserved space
	boost::shared_ptr < SummaryStatisticSet > ss_set( new SummaryStatisticSet(nss) );
	
	// construct and add summary statistics to the set
	if (!empty()) {
		
		size_t nstats = stanvals->size();
		
		for (size_t i = 0; i < nss; ++i) {
			
			std::vector < double > tmp;
			tmp.reserve(nstats);
			
			for (size_t j = 0; j < nstats; ++j) {
				tmp.push_back( (stanvals->at(j)).at(i) );
			}
			
			boost::shared_ptr < SummaryStatistic > sptr( new SummaryStatistic(tmp) );
			ss_set->add(sptr);
		}
	}
	
	return ss_set;
}


boost::shared_ptr < SummaryStatisticSet >
		SummaryStatisticSet::makeSelectedSummaryStatisticSet(
				std::vector < size_t > reqStats) const
{
	
	if (reqStats.empty()) 
		throw std::invalid_argument(
		"SummaryStatisticSet::makeSelectedSummaryStatisticSet(const std::vector < size_t >&)");

	if (*std::min_element(reqStats.begin(), reqStats.end()) < 1)
		throw std::invalid_argument(
		"SummaryStatisticSet::makeSelectedSummaryStatisticSet(const std::vector < size_t >&)");

	if (!empty() && ( this->at(0)->size() < *std::max_element(reqStats.begin(), reqStats.end())) ) 
		throw std::invalid_argument(
		"SummaryStatisticSet::makeSelectedSummaryStatisticSet(const std::vector < size_t >&)");

	

	//get all the values as vec of vecs
	boost::shared_ptr < std::vector < std::vector < double > > > 
				vals = getAllValues();
	
	size_t nss = size();
	
	// empty set with enough reserved space
	boost::shared_ptr < SummaryStatisticSet > ss_set( new SummaryStatisticSet(nss) );
	
	// construct and add summary statistics to the set
	if (!empty()) {
		
		std::sort(reqStats.begin(), reqStats.end());
		
		size_t nstats = reqStats.size();
		
		for (size_t i = 0; i < nss; ++i) {
			
			std::vector < double > tmp;
			tmp.reserve(nstats);
			
			for (size_t k = 0; k < nstats; ++k) {
				
				//reqStats indexing is 1, 2, ...
				size_t j = reqStats[k] - 1;
				
				tmp.push_back( (vals->at(j)).at(i) );
			}
			
			boost::shared_ptr < SummaryStatistic > sptr( new SummaryStatistic(tmp) );
			ss_set->add(sptr);
		}
	}
	
	return ss_set;
}

boost::shared_ptr < const SummaryStatistic > 
					SummaryStatisticSet::getSummaryStatistic() const
{
	
	return getSummaryStatisticMeans();
}
			

boost::shared_ptr < std::vector < std::vector < double > > >
				SummaryStatisticSet::getAllValues() const
{
	if (empty() ) {
		return boost::shared_ptr < std::vector < std::vector < double > > >
		( new std::vector < std::vector <double> >() ); 
	}
	else {
		size_t n = (*begin())->size();
		
		std::vector < double > tmp;
		boost::shared_ptr < std::vector < std::vector < double > > >
		results ( new std::vector < std::vector < double > > (n, tmp));
		size_t set_size = size();
		
		using namespace boost::lambda;
		
		for_each(results->begin(), results->end(), 
			bind(&std::vector<double>::reserve, _1, set_size));
		
		for (SummaryStatisticSet::const_iterator it = begin();
				it< end();
				++it) {
				
				std::vector < double > v = (*it)->getValues();	
				size_t this_n = v.size();
			
				size_t extra = this_n > n ? this_n - n : 0;
				if (extra) n = this_n;
			
				while (extra > 0) {
					results->push_back( std::vector < double > 
							( distance( begin(), it ),
							std::numeric_limits<double>::quiet_NaN() ) );
					(results->back()).reserve( set_size );	
					
					extra--;
				}
					
			for (size_t i = 0; i < this_n; ++i) {
				(results->at(i)).push_back( (*it)->at(i) );
			}
			//stuff the rest with nans
			for (size_t i = this_n; i < n; ++i) {
				(results->at(i)).push_back(std::numeric_limits<double>::quiet_NaN());
			}
					
		}
		
		return results;
	}
}

boost::shared_ptr < std::vector < std::vector < double > > >
				SummaryStatisticSet::getAllValuesSumStatLayout() const
{
	if (empty() ) {
		return boost::shared_ptr < std::vector < std::vector < double > > >
		( new std::vector < std::vector <double> >() ); 
	}
	else {
		size_t set_size = size();
		
		std::vector < double > tmp;
		boost::shared_ptr < std::vector < std::vector < double > > >
		results ( new std::vector < std::vector < double > > ());
		
		results->reserve(set_size);
		
		for (SummaryStatisticSet::const_iterator it = begin();
				it< end();
				++it) {
				
				results->push_back( (*it)->getValues() );	
					
		}
		
		return results;
	}
}

boost::shared_ptr < std::vector < std::vector < double > > >
				SummaryStatisticSet::getAllValuesStandardised() const
{
	std::vector < double > means;
	std::vector < double > sds;
	
	return getAllValuesStandardised(means, sds);
}

boost::shared_ptr < std::vector < std::vector < double > > >
				SummaryStatisticSet::getAllValuesStandardised(
				std::vector < double >& means,
				std::vector < double >& sds) const
{
	try {
		if (empty() ) {
			return boost::shared_ptr < std::vector < std::vector < double > > >
			( new std::vector < std::vector <double> >() ); 
		}
		
		if ( !sizeConsistencyCheck() ) {
			throw std::runtime_error(":\n Inconsistent SummaryStatistic sizes in set");
		}
		
		/* we know this is not empty, so now should
		 * only have a problem if size is 1 in which case we'll try to
		 * divide by the sample sd 0, which will be a bad plan */
		if (size() == 1 ) { 
			throw std::runtime_error(": size() == 1");
		}
			
		boost::shared_ptr < std::vector < std::vector < double > > >
				allValuesPtr = getAllValues();
		
		/*nstat inner containers */
		
		means = getSummaryStatisticsMeanVec(means);
		
		sds = getSummaryStatisticsVarVec(sds);
		using namespace boost::lambda;
		transform(sds.begin(), sds.end(), sds.begin(),
		bind<double>(sqrt, _1 ) );
	
		
		for (size_t i = 0; i < allValuesPtr->size(); ++i) {
			
			double mean = means.at(i);
			double sd = sds.at(i);
			
			//using namespace boost::lambda;
			transform(allValuesPtr->at(i).begin(),allValuesPtr->at(i).end(), 
				allValuesPtr->at(i).begin(),
				bind(&SummaryStatisticSet::standardise, _1 , mean, sd) );
			
		}

		return allValuesPtr;	
	}
	catch (std::runtime_error& re) {
		throw std::runtime_error(
				std::string("SummaryStatisticSet::getAllValuesStandardised()")
				+ re.what());
	}
}

boost::shared_ptr < std::vector < std::vector < double > > >
	SummaryStatisticSet::getAllValuesStandardisedSumStatLayout() const
{
	std::vector < double > means;
	std::vector < double > sds;
	
	return getAllValuesStandardisedSumStatLayout(means, sds);
}

boost::shared_ptr < std::vector < std::vector < double > > >
	SummaryStatisticSet::getAllValuesStandardisedSumStatLayout(
				std::vector < double >& means,
				std::vector < double >& sds) const
{
	try {
		if (empty() ) {
			return boost::shared_ptr < std::vector < std::vector < double > > >
			( new std::vector < std::vector <double> >() ); 
		}
		
		if ( !sizeConsistencyCheck() ) {
			throw std::runtime_error(":\n Inconsistent SummaryStatistic sizes in set");
		}
		
		/* we know this is not empty, so now should
		 * only have a problem if size is 1 in which case we'll try to
		 * divide by the sample sd 0, which will be a bad plan */
		if (size() == 1 ) { 
			throw std::runtime_error(": size() == 1");
		}
			
		boost::shared_ptr < std::vector < std::vector < double > > >
				allValuesPtr = getAllValuesSumStatLayout();
		
		size_t set_size = allValuesPtr->size();
		
		means = getSummaryStatisticsMeanVec(means);
		
		sds = getSummaryStatisticsVarVec(sds);
		using namespace boost::lambda;
		transform(sds.begin(), sds.end(), sds.begin(),
		bind<double>(sqrt, _1 ) );
	
		
		for (size_t i = 0; i < set_size; ++i) {
			
			//using namespace boost::lambda;
			std::transform(allValuesPtr->at(i).begin(),allValuesPtr->at(i).end(), 
				means.begin(), allValuesPtr->at(i).begin(),
				_1 - _2 );
			
			//using namespace boost::lambda;
			std::transform(allValuesPtr->at(i).begin(),allValuesPtr->at(i).end(), 
				sds.begin(), allValuesPtr->at(i).begin(),
				_1/_2 );
			
		}

		return allValuesPtr;	
	}
	catch (std::runtime_error& re) {
		throw std::runtime_error(
				std::string(
				"SummaryStatisticSet::getAllValuesStandardisedSumStatLayout()")
				+ re.what());
	}
}

			
boost::shared_ptr < const SummaryStatistic > 
				SummaryStatisticSet::getSummaryStatisticMeans() const
{
	std::vector < double > result;
	
	result = getSummaryStatisticsMeanVec(result);
	
	return boost::shared_ptr < const SummaryStatistic > ( 
						new SummaryStatistic(result) );
}

						
boost::shared_ptr < const SummaryStatistic > 
				SummaryStatisticSet::getSummaryStatisticVars() const
{
	std::vector < double > result;
	
	result = getSummaryStatisticsVarVec(result);
	
	return boost::shared_ptr < const SummaryStatistic > ( 
						new SummaryStatistic(result) );
	
}


boost::shared_ptr < const SummaryStatistic > 
				SummaryStatisticSet::getSummaryStatisticSDs() const
{
	std::vector < double > result;
	
	result = getSummaryStatisticsVarVec(result);
	
	using namespace boost::lambda;
	transform(result.begin(), result.end(), result.begin(),
	bind<double>(sqrt, _1 ) );
	
	return boost::shared_ptr < const SummaryStatistic > ( 
						new SummaryStatistic(result) );
}	



void  SummaryStatisticSet::outputToFile(	const std::string& s,
											bool append) const

{
	 std::ofstream os;
        
	if (append) os.open(s.c_str(), std::ios::app);         // append
	else os.open(s.c_str()); // overwrite
        
    if (os.is_open()) {

        os << plainString();
        os.close();
    }
    else {
        std::cerr << "Error: could not open file named "
            << s << std::endl << std::endl;
    }
}


SummaryStatisticSet::const_iterator SummaryStatisticSet::begin() const
{
	return container.begin();
}

SummaryStatisticSet::const_iterator SummaryStatisticSet::end() const
{
	return container.end();
}
	/*		
double SummaryStatisticSet::getDistance(
					const boost::shared_ptr < SummaryStatistic >& ref)
					const
{	
	
	try {
		//std::vector < double > Ds = getTajimasDMeans();
		std::vector < double > Hs = getHeterozygosityMeans();
		std::vector < double > Ss = getSegSitesMeans();
		
		//double refD = ref->getMeanTajimasD();
		double refH = ref->getMeanHeterozygosity();
		double refS = ref->getMeanSegSites();
		
		//double distanceD = normalisedDistance(refD, Ds);
		double distanceH = normalisedDistance(refH, Hs);
		double distanceS = normalisedDistance(refS, Ss);
		double distance = sqrt(+ distanceH*distanceH
								+ distanceS*distanceS);
		
		#ifdef DEBUG
			
			
			
		#endif
		
						
		return distance;
	}
	catch (std::exception& e) {
		throw AABCexception(
				"Error in SummaryStatisticSet::getDistance :\n"
						+ std::string(e.what()) );
		
	}	
	
	
	
}
* */

std::string SummaryStatisticSet::toString() const
{
	int prec = 3;
	return toString(prec);

}

std::string SummaryStatisticSet::toString(int prec) const
{
	std::ostringstream stm;
	
	for (std::vector < 
				boost::shared_ptr < const SummaryStatistic > 
				>::const_iterator it = container.begin();
		it < container.end();
		++it) {
		stm << (*it)->toString(prec) << "\n";
	}
	
	return stm.str();

}


std::string SummaryStatisticSet::plainString(int prec) const
{
	std::ostringstream stm;
	
	
	for (std::vector < 
				boost::shared_ptr < const SummaryStatistic > 
				>::const_iterator it = container.begin();
		it < container.end();
		++it) {
		stm << (*it)->plainString(prec) << "\n";
	}
	
	return stm.str();

}


// protected

std::vector<double>& SummaryStatisticSet::getSummaryStatisticsMeanVec(
								std::vector < double >& result) const
{
	if (empty() ) {
		throw std::runtime_error(
		std::string("SummaryStatisticSet::getSummaryStatisticsMeanVec(std::vector < double >&)")
		+ std::string(":\n Set is empty") );
	}
	if ( !sizeConsistencyCheck() ) {
		throw std::runtime_error(
		std::string("SummaryStatisticSet::getSummaryStatisticsMeanVec(std::vector < double >&)")
		+ std::string(":\n Inconsistent SummaryStatistic sizes in set") );
	}
	
	//size_t n = (stats.front())->size();	// expected size
	result = (container.front())->getValues();
			
	size_t n = size();
	
	using namespace boost::lambda;
	
	for (size_t i = 1; i < n; ++i) {
		std::vector < double > tmp = (container.at(i))->getValues();
		
		transform ( result.begin(), result.end(), tmp.begin(), 
						result.begin(), _1 + _2 );
	}
	
	transform ( result.begin(), result.end(), 
					result.begin(), _1 / n );
	
	
	return result; // return by ref
}

std::vector<double>& SummaryStatisticSet::getSummaryStatisticsVarVec(
									std::vector<double>& result) const
{
	if (empty() ) {
		throw std::runtime_error(
		std::string("SummaryStatisticSet::getSummaryStatisticsVarVec(std::vector < double >&)")
		+ std::string(":\n Set is empty") );
	}
	
		
	if (size() == 1 ) { // return a container of 0.0s
		result = std::vector < double> 
					((container.front())->size(), 0.0);
		
	}
			
	else { // we know size() > 1	 
		if ( !sizeConsistencyCheck() ) {
			throw std::runtime_error(
			std::string("SummaryStatisticSet::getSummaryStatisticsVarVec(std::vector < double >&)")
			+ std::string(":\n Inconsistent SummaryStatistic sizes in set") );
		}
		
		using namespace boost::lambda;
		
		//size_t n = (stats.front())->size();	// expected size
		result = (container.front())->getValues();
		transform ( result.begin(), result.end(),  
							result.begin(), _1 * _1 );
							
		size_t n = size();
		for (size_t i = 1; i < n; ++i) {
			std::vector < double > tmp = (container.at(i))->getValues();
			transform ( result.begin(), result.end(), tmp.begin(), 
							result.begin(), _1 + (_2 * _2) );
		}
		
		std::vector < double > means;
	
		means = getSummaryStatisticsMeanVec(means);
		
		transform ( result.begin(), result.end(), means.begin(), 
							result.begin(), (_1 - (n * (_2 * _2)))/(n-1) );
	}			
	return result;
}

double SummaryStatisticSet::standardise(double value, double mean, double sd)
{
	if (!gsl_isnan(value)) {
		value = (value - mean)/sd;
	}
	
	return value;
	 
}


const std::vector< boost::shared_ptr < const SummaryStatistic > >&
				SummaryStatisticSet::_getContainer() const
{
	return container;
}

// non member functions
std::ostream& mct::operator<< (std::ostream &out, const SummaryStatisticSet &s)
{
	out << s.toString();
	return out;
}


