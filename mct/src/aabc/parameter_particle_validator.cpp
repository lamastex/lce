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
\brief ParameterParticleValidator definitions.
*/
#include "parameter_particle_validator.hpp"

#include <boost/lambda/bind.hpp>
#include <boost/lambda/lambda.hpp>

#include <boost/function.hpp>
#include <boost/bind.hpp>

#include <cmath> //sqrt

#include <algorithm> // transform, for_each
#include <numeric> // accumulate
#include <stdexcept>

//#define MYDEBUG
#ifdef MYDEBUG
	#include <iostream>
#endif

using namespace aabc;

ParameterParticleValidator::ParameterParticleValidator(
			const boost::shared_ptr < const mct::SummaryStatistic >& _r)
				: 	referenceStats(_r) 
{
	normaliser = boost::shared_ptr < mct::SummaryStatistic >
		( new mct::SummaryStatistic (
	std::vector < double > (referenceStats->size(), 1.0) ) );
		
}

		
ParameterParticleValidator::~ParameterParticleValidator() 
{
		
		
}


boost::shared_ptr< ParameterParticleValidationResult >
		ParameterParticleValidator::getValidationResult(
		const boost::shared_ptr< const ParameterParticle>& p) const

{
	boost::shared_ptr< ParameterParticleValidationResult  > 
					pr (new ParameterParticleValidationResult() );
	
	pr->particle = p;
	pr->summary_statistic = getSummaryStatistic(p);
	pr->distance = getDistance(pr->summary_statistic);
	
	return pr;
	
}


boost::shared_ptr< ParameterParticleValidationResultSet >
		ParameterParticleValidator::getValidationResultSet(
		const boost::shared_ptr< const ParameterParticleSet>& pset) const

{
	boost::shared_ptr< ParameterParticleValidationResultSet > 
			results( new ParameterParticleValidationResultSet(
			pset->size()) );
	for (size_t i = 0; i < pset->size(); ++i) {
		results->add(getValidationResult(pset->at(i)));
	}
	return results;
}


boost::shared_ptr< ParameterParticleValidationResultSet >
	ParameterParticleValidator::setNormaliserAndGetValidationResultSet(
	const boost::shared_ptr< const ParameterParticleSet>& pset)
													const
{
	boost::shared_ptr< ParameterParticleValidationResultSet > 
			results( new ParameterParticleValidationResultSet() );
	
	boost::shared_ptr< mct::SummaryStatisticSet > 
		ss_set = getSummaryStatisticSet(pset);
	
	_setNormaliser(ss_set);
	
	results = _updateValidationResultSet(results, pset, ss_set);
	
	return results;
}


boost::shared_ptr< ParameterParticleValidationResultSet >&
	ParameterParticleValidator::updateValidationResultSet(
	boost::shared_ptr< ParameterParticleValidationResultSet >& results,
	const boost::shared_ptr< const ParameterParticleSet>& pset)
													const
{
	boost::shared_ptr< mct::SummaryStatisticSet > 
		ss_set = getSummaryStatisticSet(pset);
	
	results = _updateValidationResultSet(results, pset, ss_set);
	return results;
}


boost::shared_ptr < mct::SummaryStatistic >  
			ParameterParticleValidator::setNormaliser(
			const boost::shared_ptr< const ParameterParticleSet>& pset)
																const
{
	boost::shared_ptr < mct::SummaryStatisticSet > 
		ss_set = getSummaryStatisticSet(pset);
	return ParameterParticleValidator::_setNormaliser(ss_set);
	
}

void ParameterParticleValidator::setNormaliser(
				const boost::shared_ptr < const mct::SummaryStatistic >& s)
																const
{
	if (s->size() != referenceStats->size() ) {
		throw std::invalid_argument(
		"ParameterParticleValidator::setNormaliser(const boost::shared_ptr < const mct::SummaryStatistic >&");
	}
		
	// set normaliser to a copy of the given SummaryStatistic
	normaliser = boost::shared_ptr < mct::SummaryStatistic >
							( new mct::SummaryStatistic(*s) );
}


boost::shared_ptr < mct::SummaryStatistic >  
					ParameterParticleValidator::getNormaliser() const
{
	// return a copy of the normaliser
	return boost::shared_ptr < mct::SummaryStatistic >
						( new mct::SummaryStatistic(*normaliser) );
	
}


boost::shared_ptr< ParameterParticleValidationResultSet >&
	ParameterParticleValidator::_updateValidationResultSet(
	boost::shared_ptr< ParameterParticleValidationResultSet >& results,
	const boost::shared_ptr< const ParameterParticleSet>& pset,
	const boost::shared_ptr< mct::SummaryStatisticSet >& ss_set)
													const
{
	//make sure we have enough space for additions 
	results->reserve(results->size() + pset->size());
	
	for (size_t i = 0; i < pset->size(); ++i) {
		boost::shared_ptr< const mct::SummaryStatistic > ss = ss_set->at(i);
		
		#ifdef MYDEBUG
			std::cout << "i = " << i << std::endl;
			std::cout << "particle = " << (pset->at(i))->toString() << std::endl;
			std::cout << "summary stats = " << ss->toString() << std::endl;
		#endif
		
		double distance = getDistance(ss);
		
		#ifdef MYDEBUG
			std::cout << "distance = " << distance << std::endl;
		#endif
	
		boost::shared_ptr< ParameterParticleValidationResult > new_r
				(new ParameterParticleValidationResult(
							pset->at(i), ss, distance) );
		results->add(new_r);
	}
	return results;
}


boost::shared_ptr < mct::SummaryStatistic > 
			ParameterParticleValidator::_setNormaliser(
			const boost::shared_ptr< const mct::SummaryStatisticSet>& ss_set)
																const
{
	boost::shared_ptr < const mct::SummaryStatistic > 
				n = ss_set->getSummaryStatisticSDs();
	normaliser = n;
	#ifdef MYDEBUG
		std::cout << "\nnormaliser now is" << std::endl;
		std::cout << normaliser->toString() << std::endl;
	#endif
	return getNormaliser();
}


		
double ParameterParticleValidator::getDistance(
		const boost::shared_ptr < const mct::SummaryStatistic >& ss)
		const
{	
	// note no checks for valid pointers here
	std::vector < double > values = ss->getValues();
	std::vector < double > refs = referenceStats->getValues();
	std::vector < double > normals = normaliser->getValues();
	
	//  sum of ((values-ref)/normals)^2
	
	using namespace boost::lambda;
	
	boost::lambda::placeholder1_type v1;
	boost::lambda::placeholder2_type v2;
	
	//size_t n = (stats.front())->size();	// expected size
	transform ( values.begin(), values.end(), refs.begin(), 
						values.begin(), v1 - v2 );
	#ifdef MYDEBUG
		std::cout << "\t refs =  " << referenceStats->toString() << std::endl;
		std::cout << "\tunnormalisedDistances = ";
		for (std::vector<double>::iterator it = values.begin(); it < values.end(); ++it) {
			std::cout << "\t"<< (*it);
		}
		std::cout << std::endl;
	#endif
	transform ( values.begin(), values.end(), normals.begin(), 
						values.begin(), (v1*v1)/(v2*v2) );
	#ifdef MYDEBUG
		std::cout << "\tnormalisedSqDistances = ";
		for (std::vector<double>::iterator it = values.begin(); it < values.end(); ++it) {
			std::cout << "\t"<< (*it);
		}
		std::cout << std::endl;
	#endif
	return std::sqrt(std::accumulate(values.begin(), values.end(), 0.0));
	
}
	

