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
\brief ParameterParticleValidationResult definitions.
*/

#include "parameter_particle_validation_result.hpp"

//#include <boost/function.hpp>
//#include <boost/bind.hpp>


#include <sstream>  // to be able to manipulate strings as streams
//#include <iterator> // for ostream_iterator
#include <iostream>

using namespace aabc;

ParameterParticleValidationResult::ParameterParticleValidationResult()
				: distance(0.0) {}
		
ParameterParticleValidationResult::ParameterParticleValidationResult(
				const boost::shared_ptr< const ParameterParticle >& _p,
				const boost::shared_ptr< const mct::SummaryStatistic >& _s,
				double _d)
				: particle(_p), summary_statistic(_s), distance(_d) {}

ParameterParticleValidationResult::~ParameterParticleValidationResult()
{
	//std::cout << "ParameterParticleValidationResult destructor called, distance = " << distance << std::endl;
}


void ParameterParticleValidationResult::setDistance(double d)
{
	distance = d;
}
		
std::string ParameterParticleValidationResult::toString() const
{
	std::ostringstream stm;
	stm << "Particle: " << particle->toString() << "\n" 
		<< "SummaryStatistic: " << summary_statistic->toString() << "\n" 
		<< "Distance: " << distance << std::endl;
	return stm.str();
}


// non member functions
std::ostream& operator<< (std::ostream &out, 
						const ParameterParticleValidationResult &pr)
{
	out << pr.toString();
	return out;
}

