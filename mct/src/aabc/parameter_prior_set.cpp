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
\brief ParameterPriorSet definitions.
*/

#include "parameter_prior_set.hpp"

#include <stdexcept>
#include <sstream>  // to be able to manipulate strings as streams

using namespace aabc;

ParameterPriorSet::ParameterPriorSet(size_t _n)
	: n_params(_n), 
	priors(std::vector< boost::shared_ptr< ParameterPrior > >(_n))
	{}

ParameterPriorSet::ParameterPriorSet(const ParameterPriorSet& other)
	: n_params(other.n_params), 
	priors(other.priors)
	{}

ParameterPriorSet::~ParameterPriorSet() 
{
	//std::cout << "ParameterPriorSet destructor called" << std::endl;
}
			
const boost::shared_ptr< ParameterPrior >& ParameterPriorSet::at(
											const size_t index) const
{
	return priors.at(index);
}

	
boost::shared_ptr< ParameterPrior >& ParameterPriorSet::at(
											const size_t index)
{
	return priors.at(index);
}

// will throw an exception if particle is shorter than prior
double ParameterPriorSet::pdf(
				const boost::shared_ptr < const ParameterParticle >& p) const
{
	double result = 1.0;
	size_t i = 0;
	while (i < n_params && result > 0.0) {
		result *= pdfFrom(i, p->at(i));
		++i;
	}
	return result;
}


double ParameterPriorSet::cdfFrom(const size_t index, 
									const double x) const
{
	return (priors.at(index))->cdf(x);
}


double ParameterPriorSet::pdfFrom(const size_t index, 
									const double x) const
{
	return (priors.at(index))->pdf(x);
}
	
boost::shared_ptr < const ParameterParticle >
		ParameterPriorSet::generateParameterParticle() const
					
{
	boost::shared_ptr < ParameterParticle > 
							p ( new ParameterParticle(n_params) );
	for (size_t j = 0; j < n_params; ++j)
	{
		p->at(j) = randFrom(j);
	}
	return p;
	
}

boost::shared_ptr < ParameterParticleSet >
	ParameterPriorSet::generateParameterParticleSet(size_t n_reps) const
					
{
	boost::shared_ptr< ParameterParticleSet > 
						pset( new ParameterParticleSet(n_reps) );
	
	for (size_t i = 0; i < n_reps; ++i) {
		#ifdef DEBUG
			if ( (n_reps <= 10) || ( i % (n_reps/10) == 0) ) {
				std::cout << ".";
				std::cout.flush();
			}
		#endif
		
		pset->add( generateParameterParticle() );
	}
	return pset; 
	
}

double ParameterPriorSet::randFrom(const size_t index) const
{
	return (priors.at(index))->rand();
}


void ParameterPriorSet::testFrom(const size_t index) const
{
	(priors.at(index))->test();
}
			
void ParameterPriorSet::testAll() const
{
	for (size_t i = 0; i < n_params; ++i)
	{
		(priors.at(i))->test();
	}
}


size_t ParameterPriorSet::size() const
{
	return n_params;
}

std::string ParameterPriorSet::toString() const
{
	std::ostringstream stm;
	
	for (std::vector < boost::shared_ptr< ParameterPrior > >::const_iterator 
			it = priors.begin();
			it < priors.end();
			++it) {
		stm << (*it)->toString() << "  ";		
				
	}
		
	return stm.str();
}
