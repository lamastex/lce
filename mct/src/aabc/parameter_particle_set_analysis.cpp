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
\brief ParameterParticleSetAnalysis definitions.
*/

#include "parameter_particle_set_analysis.hpp"

#include <boost/lambda/bind.hpp>
#include <boost/lambda/lambda.hpp>

#include <gsl/gsl_math.h> //gsl_isnan()

#include <cmath>

#include <algorithm>
#include <stdexcept>
#include <numeric>
#include <stdexcept>

//#define MYDEBUG
#ifdef MYDEBUG
	#include <iostream>
#endif

using namespace aabc;

ParameterParticleSetAnalysis::ParameterParticleSetAnalysis(
		const boost::shared_ptr< const ParameterParticleSet >& pset)
		: n(0)
{
	if (!pset || pset->empty()) {
		throw std::invalid_argument(
			std::string("ParameterParticleSetAnalysis::ParameterParticleSetAnalysis(")
			+ std::string("const boost::shared_ptr< const ParameterParticleSet >&)") );
	}
	if ( !pset->sizeConsistencyCheck() ) {
			throw std::invalid_argument(
			std::string("ParameterParticleSetAnalysis::ParameterParticleSetAnalysis(")
			+ std::string("const boost::shared_ptr< const ParameterParticleSet >&)")
			+ std::string(" :\ninconsistent particle sizes") );
	}
	
	parameters =  pset->getAllValues();
	
	size_t nParams = parameters->size();
	
	#ifdef MYDEBUG
			std::cout << "nParams = " << nParams << std::endl;
		#endif
		
	if (nParams > 0) {
		n = (parameters->begin())->size();
		
		#ifdef MYDEBUG
			std::cout << "n = " << n << std::endl;
		#endif
		
	
		//boost::shared_ptr < std::vector< std::vector< double > > > 
		//( new std::vector < std::vector < double> >(nParams, tmp) );
		sums.reserve(nParams);
		sumsqs.reserve(nParams);
	}
	
	for (std::vector < std::vector < double > >::iterator it = parameters->begin();
			it < parameters->end();
			++it)
	{
		sums.push_back( accumulate( it->begin(), it->end(), 0.0 ) );
		using namespace boost::lambda;
		sumsqs.push_back( 
			accumulate( it->begin(), it->end(), 0.0, _1 + (_2 * _2) ) );
	}
}

			
ParameterParticleSetAnalysis::~ParameterParticleSetAnalysis() 
{
	#if(0)
		#ifdef DEBUG
			std::cout << "ParameterParticleSetAnalysis destructor called" << std::endl;
		#endif
	#endif	
}

size_t ParameterParticleSetAnalysis::size() const
{
	return parameters->size();
}

bool ParameterParticleSetAnalysis::empty() const
{
	return parameters->empty();
}

size_t ParameterParticleSetAnalysis::nParticles() const
{
	return n;
}

boost::shared_ptr < mct::SummaryStatisticSet >
		ParameterParticleSetAnalysis::makeSummaryStatisticSetStandardisedParameterValues() const
{
	//get all the standardised values as vec of vecs
	boost::shared_ptr < std::vector < std::vector < double > > > 
				stanvals = getAllParametersStandardised();
	
	// empty set with enough reserved space
	boost::shared_ptr < mct::SummaryStatisticSet > 
						ss_set( new mct::SummaryStatisticSet(n) );
	
	// construct and add summary statistics to the set
	if (!empty()) {
		
		size_t nstats = stanvals->size();
		
		for (size_t i = 0; i < n; ++i) {
			
			std::vector < double > tmp;
			tmp.reserve(nstats);
			
			for (size_t j = 0; j < nstats; ++j) {
				tmp.push_back( (stanvals->at(j)).at(i) );
			}
			
			boost::shared_ptr < mct::SummaryStatistic > sptr( new mct::SummaryStatistic(tmp) );
			ss_set->add(sptr);
		}
	}
	
	return ss_set;
}


boost::shared_ptr < std::vector< std::vector< double > > >
			ParameterParticleSetAnalysis::getAllParameters() const
{
	// return a copy of the contents of parameters
	
	boost::shared_ptr < std::vector< std::vector< double > > > result
		(new std::vector< std::vector< double > > (*parameters));
	return result;
}
						

boost::shared_ptr < std::vector < std::vector < double > > >
				ParameterParticleSetAnalysis::getAllParametersStandardised() const
{
	try {
		/* we know this is not empty, so now should
		 * only have a problem if size is 1 in which case we'll try to
		 * divide by the sample sd 0, which will be a bad plan */
		if (n == 1 ) { 
			throw std::runtime_error(": n == 1");
		}
			
		/*nstat inner containers */
		
		std::vector < double > means;
		means = vec_sampleMeans(means);
		
		std::vector < double > sds;
		sds = vec_sampleVars(sds);

		boost::lambda::placeholder1_type _s;
		
		//sds as sqrt vars
		transform(sds.begin(), sds.end(), sds.begin(),
		boost::lambda::bind<double>(sqrt, _s ) );
	
		boost::shared_ptr < std::vector< std::vector< double > > > results
			( new std::vector< std::vector < double > >(*parameters) );
		
		for (size_t i = 0; i < results->size(); ++i) {
			
			double mean = means.at(i);
			double sd = sds.at(i);
			
			boost::lambda::placeholder1_type _v;
			
			transform(results->at(i).begin(),results->at(i).end(), 
				results->at(i).begin(),
				boost::lambda::bind(&ParameterParticleSetAnalysis::standardise, _v , mean, sd) );
			
		}

		return results;	
	}
	catch (std::runtime_error& re) {
		throw std::runtime_error(
				std::string(
				"ParameterParticleSetAnalysis:::getAllParametersStandardised()")
				+ re.what());
	}
}

boost::shared_ptr < const mct::SummaryStatistic > 
				ParameterParticleSetAnalysis::getMeans() const
{
	std::vector < double > means;
	means = vec_sampleMeans(means);
	
	return boost::shared_ptr < const mct::SummaryStatistic >
		( new mct::SummaryStatistic(means) );
	
}

		
		
boost::shared_ptr < const mct::SummaryStatistic > 
				ParameterParticleSetAnalysis::getSampleVars() const
{
	std::vector < double > vars;
	vars = vec_sampleVars(vars);
	
	return boost::shared_ptr < const mct::SummaryStatistic >
		( new mct::SummaryStatistic(vars) );
	
}

			
boost::shared_ptr < const mct::SummaryStatistic > 
				ParameterParticleSetAnalysis::getSampleSDs() const
{
	std::vector < double > sds;
	sds = vec_sampleVars(sds);
	
	transform(sds.begin(), sds.end(), sds.begin(),
								static_cast < double (*)(double)> (&std::sqrt));
	
	return boost::shared_ptr < const mct::SummaryStatistic >
		( new mct::SummaryStatistic(sds) );
	
}

			
boost::shared_ptr < const mct::SummaryStatistic > 
				ParameterParticleSetAnalysis::getRMSEs(
				const boost::shared_ptr < const ParameterParticle >& p) const
{
	std::vector < double > result = sumsqs;
	std::vector < double > tmp = sums;
	std::vector < double > tru = p->get_vector();
	
	if (tru.size() != sums.size()) {
		throw std::invalid_argument(
		"ParameterParticleSetAnalysis::getRMSEs(const boost::shared_ptr < ParameterParticle >&");
	}
	
	using namespace boost::lambda;
	// -2*sum(p_i)/n * true_pi
	transform(tmp.begin(), tmp.end(), tru.begin(), tmp.begin(),
		(-2 * _1*_2/n));
	// sum(theta_i ^2)/n + true_pi^2
	transform(result.begin(), result.end(), tru.begin(), result.begin(),
		(_1/n + (_2 * _2)));
	// sum(theta_i ^2)/n + true_pi^2 - 2*sum(p_i)/n * true_pi
	transform(result.begin(), result.end(), tmp.begin(), result.begin(),
		(_1 + _2));

	
	return boost::shared_ptr < const mct::SummaryStatistic >
		( new mct::SummaryStatistic(result) );
	
}

std::vector < double >& 
				ParameterParticleSetAnalysis::vec_sampleMeans(
						std::vector < double >& vec) const
{
	
	vec = sums;
	using namespace boost::lambda;
	
	transform(vec.begin(), vec.end(), vec.begin(), _1/n);
	
	return vec;
	
}
			
std::vector < double >& 
				ParameterParticleSetAnalysis::vec_sampleVars(
						std::vector < double >& vec) const
{
	vec = sumsqs;
		
	using namespace boost::lambda;
	
	if (n == 1) {
		transform(vec.begin(), vec.end(), vec.begin(),
								(_1 = 0.0) );
	}
	else { 
		
			transform(vec.begin(), vec.end(), sums.begin(), vec.begin(),
									(_1 - (_2 * _2)/n)/(n-1) );
	}
	return vec;
			
}

double ParameterParticleSetAnalysis::standardise(
							double value, double mean, double sd)
{
	if (!gsl_isnan(value)) {
		value = (value - mean)/sd;
	}
	
	return value;
	 
}

