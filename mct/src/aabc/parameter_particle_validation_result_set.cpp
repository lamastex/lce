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
\brief ParameterParticleValidationResultSet definitions.
*/

#include "parameter_particle_validation_result_set.hpp"

#include <boost/function.hpp>
#include <boost/bind.hpp>

#include <gsl/gsl_statistics_double.h>

#include <cmath> 
#include <limits> 
#include <algorithm> // for partition. partial_sort and for_each
#include <sstream>  // to be able to manipulate strings as streams
#include <iterator> // for ostream_iterator
#include <stdexcept>



//#define MYDEBUG
#ifdef MYDEBUG
	#include <iostream>
#endif

using namespace aabc;

ParameterParticleValidationResultSet::ParameterParticleValidationResultSet(){}

ParameterParticleValidationResultSet::ParameterParticleValidationResultSet(
			size_t reservedSpace)
{
	validation_results.reserve(reservedSpace);
	
}

ParameterParticleValidationResultSet::~ParameterParticleValidationResultSet(){}

void ParameterParticleValidationResultSet::add(
		const boost::shared_ptr< ParameterParticleValidationResult >& vr)
{
	validation_results.push_back(vr);
}

ParameterParticleValidationResultSet::iterator ParameterParticleValidationResultSet::begin()
{
	return validation_results.begin();
}
			
ParameterParticleValidationResultSet::iterator ParameterParticleValidationResultSet::end()
{
	return validation_results.end();
}

size_t ParameterParticleValidationResultSet::size() const
{
	return validation_results.size();
}
			
bool ParameterParticleValidationResultSet::empty() const
{
	return validation_results.empty();
}

void ParameterParticleValidationResultSet::reserve(
									size_t reservedSpace)
{
	validation_results.reserve(reservedSpace);
}

boost::shared_ptr< ParameterParticleSet > 
		ParameterParticleValidationResultSet::getParameterParticleSet()
{
	boost::shared_ptr< ParameterParticleSet >  
			pset(new ParameterParticleSet( size() ));
		
	boost::function< void (boost::shared_ptr < const ParameterParticle >) >
				// need cast to specify exactly which add to use - args are not automatically used to resolve function address
				f1 = boost::bind(static_cast<void (ParameterParticleSet::*)( const boost::shared_ptr < const ParameterParticle >& )>(&ParameterParticleSet::add),
				pset, _1);
		
	boost::function<boost::shared_ptr < const ParameterParticle > 
				(boost::shared_ptr < ParameterParticleValidationResult >)>
					f2 = boost::bind(&ParameterParticleValidationResult::particle, _1);
		
		
		std::for_each( begin(), end(),
			boost::bind(f1, boost::bind(f2, _1) ) );
			 
	 
	return pset;
}


boost::shared_ptr< mct::SummaryStatisticSet >  
		ParameterParticleValidationResultSet::getSummaryStatisticSet()
{
	boost::shared_ptr< mct::SummaryStatisticSet >  
		ss_set(new mct::SummaryStatisticSet( size() ));
	
	boost::function< void (boost::shared_ptr < const mct::SummaryStatistic >) >
				// need cast to specify exactly which add to use - args are not automatically used to resolve function address
				f1 = boost::bind(static_cast<void (mct::SummaryStatisticSet::*)( const boost::shared_ptr < const mct::SummaryStatistic >& )>(&mct::SummaryStatisticSet::add),
							ss_set, _1);
		
	boost::function<boost::shared_ptr < const mct::SummaryStatistic > 
				(boost::shared_ptr < ParameterParticleValidationResult >)>
					f2 = boost::bind(&ParameterParticleValidationResult::summary_statistic, _1);
		
		
	std::for_each( begin(), end(),
		boost::bind(f1, boost::bind(f2, _1) ) );
		 
	 
	return ss_set;
}

/* make sure we discard  only > e, not >= e
* otherwise we will discard too many if e falls in the middle of a
* a stretch of consecutive distances all == e */
void ParameterParticleValidationResultSet::discardByDistance(double e)
{
	ParameterParticleValidationResultSet::iterator bound
				= std::partition( 	begin(), 
									end(),
		(boost::bind(&ParameterParticleValidationResult::distance, _1) <= e) );
		/* and take out all the results that fail, ie where distance > e
		* we don't use copy and swap because we actually want to retain
		* all capacity that is currently there.*/
		validation_results.erase(bound, end());
		//std::vector< boost::shared_ptr < ParameterParticleValidationResult > >
			//( validation_results.begin(), bound ).swap(
				//						validation_results);
}


void ParameterParticleValidationResultSet::discardN(size_t n)
{
	if (n >= size()) {
		throw std::invalid_argument(
			"ParameterParticleValidationResultSet::discardN(size_t)");
	}
	
	/* do an n-th element sort, supplying a custom sorter that puts things in
	 * order largest first.  n-th element works on positions 0, 1, 2 etc
	 * so to discard 1, say, we'd use n-th element 1 and keep from
	 * that element onwards inclusive*/
	ParameterParticleValidationResultSet::iterator bound = begin()+n;
	
	nth_element (begin(), bound, end(), 
	((boost::bind(&ParameterParticleValidationResult::distance, _1) > 
	(boost::bind(&ParameterParticleValidationResult::distance, _2))))
	);
	
	// use copy and swap to shrink-to-fit the memory allocated
	std::vector < boost::shared_ptr< ParameterParticleValidationResult > >
		( bound, end() ).swap(validation_results);
	
	#if(0)
		// reverse around the unwanted stuff
		ParameterParticleValidationResultSet::iterator first = begin();
		ParameterParticleValidationResultSet::iterator last = end();
		ParameterParticleValidationResultSet::iterator new_bound = end()-n;
		while (( first != bound )&&( first != new_bound )&&( first != --last )) {
			#ifdef MYDEBUG
			{
				std::cout << "Swapping " << **first << " and " << **last << std::endl;
			}
			#endif
			swap (*first++,*last);
		}
		

		validation_results.erase(new_bound, end());
	#endif
}




double ParameterParticleValidationResultSet::getDistancePercentile(
												double percentile)
{
	if (empty()) {
		throw std::runtime_error(
			std::string("ParameterParticleValidationResultSet::getDistancePercentile(double)")
				+ std::string(":\nSet is empty)") );
	}
	if (percentile <= 0 || percentile > 1) {
		throw std::invalid_argument(
			"ParameterParticleValidationResultSet::getDistancePercentile(double)");
	}
	
	double result = std::numeric_limits<double>::quiet_NaN();
	size_t n = size();
	double completepos = (n-1) * percentile; 
	double integralpos = std::floor(completepos); 
	double delta = completepos - integralpos;
	size_t pos = static_cast <size_t> (integralpos);
	
	#ifdef MYDEBUG
		std::cout << "In getDistancePercentile, percentile = " << percentile << ", size = " << size() << std::endl;
		std::cout << "completepos = " << completepos << "integralpos = " << integralpos<<  ", pos = " << pos 
					<< " delta = " << delta << std::endl;
		
	#endif
	
	if (pos >= (n-1)) {
		result = ( * max_element( begin(), end(), 
		((boost::bind(&ParameterParticleValidationResult::distance, _1) < 
		(boost::bind(&ParameterParticleValidationResult::distance, _2))))
		) )->distance * ( 1.0 + std::numeric_limits<double>::epsilon());
		/*
		#ifdef MYDEBUG
		{
			std::cout << std::scientific;
			std::cout.precision(17);
			std::cout << "max element is " << (* max_element( begin(), end(), 
				((boost::bind(&ParameterParticleValidationResult::distance, _1) < 
				(boost::bind(&ParameterParticleValidationResult::distance, _2))))
				) )->distance << " and epsilon is " << result << std::endl; 
		}
		#endif
		*/
	}
	
	else if ( fabs(delta - 0.0) <= std::numeric_limits<double>::epsilon() ) {
		ParameterParticleValidationResultSet::iterator bound = begin()+pos;
		nth_element (begin(), bound, end(),
		((boost::bind(&ParameterParticleValidationResult::distance, _1) < 
		(boost::bind(&ParameterParticleValidationResult::distance, _2))))
		);
		
		#ifdef MYDEBUG
			std::cout << "[" << pos <<"] = " << (validation_results.at(pos))->distance << std::endl;
		
		#endif
		result = (*bound)->distance;
	}
	else {
		
		/* use nth element to get elements in positions [pos] and 
		 [pos+1] correct, using sorting by distance. 
		 the only way that pos can be n-1 is if percentile = 1, and then
		 delta will be 0*/
		
		ParameterParticleValidationResultSet::iterator bound1 = begin()+pos+1;
		nth_element (begin(), bound1, end(),
		((boost::bind(&ParameterParticleValidationResult::distance, _1) < 
		(boost::bind(&ParameterParticleValidationResult::distance, _2))))
		);
		ParameterParticleValidationResultSet::iterator bound2 = begin()+pos;
					nth_element (begin(), bound2, bound1,
		((boost::bind(&ParameterParticleValidationResult::distance, _1) < 
		(boost::bind(&ParameterParticleValidationResult::distance, _2))))
		);
		/* 
		partial_sort( begin(), begin() + pos + 1, end(),
		((boost::bind(&ParameterParticleValidationResult::distance, _1) < 
		(boost::bind(&ParameterParticleValidationResult::distance, _2))))
		);
		*/
		#ifdef MYDEBUG
			std::cout << "[" << pos <<"] = " << (validation_results.at(pos))->distance << std::endl;
			std::cout << "[" << pos+1 <<"] = " << (validation_results.at(pos+1))->distance << std::endl;
		
		#endif
		
		
		result = ( (1.0 - delta)*(validation_results.at(pos))->distance 
				+ delta*(validation_results.at(pos+1))->distance );
	}
	
	
	#ifdef MYDEBUG
		std::cout << "result = " << result << std::endl;
		
	#endif
	return result;
}

