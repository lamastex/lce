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
\brief ParameterParticleValidationResult declarations.
*/

#ifndef _INC_PARAMETERPARTICLE_VALIDATION_RESULT_H
#define _INC_PARAMETERPARTICLE_VALIDATION_RESULT_H

#include "parameter_particle_set.hpp"
#include "parameter_particle.hpp"
#include "summary_statistic_set.hpp"
#include "summary_statistic.hpp"
#

#include <boost/smart_ptr.hpp>

#include <vector>
#include <string>


namespace aabc {

	/** @brief A struct for the results of a 
	 * validating a ParameterParticle.
	 * 
	 * The struct contains:
	 * <ul>
	 * <li>The original ParameterParticle</li>
	 * <li>The mct::SummmaryStatistic generated for that ParameterParticle</li>
	 * <li>The validation distance for that SummaryStatistic</li>
	 * </ul>
	 * 
	 * A ParameterParticleValidationResult is contingent on a particular
	 * validation method: this determines how the summary statistic
	 * is generated and how the validation distance is calculated.
	 * 
	 * ParameterParticleValidationResult overloads operator< so that 
	 * results can be sorted and ordered.  
	  */
	
	
	struct ParameterParticleValidationResult {
		
		/*! \brief The Parameter Particle. */
		boost::shared_ptr< const ParameterParticle > particle;
		
		/*! \brief The Summary Statistic associated with the 
		 * particle by whatever validation method was used. */
		boost::shared_ptr< const mct::SummaryStatistic > summary_statistic;
		
		/*! \brief The validation distance.  This will depend
		 on the validation method and validation reference object. */
		double distance;
		
		/*! \brief No-args constructor. */
		ParameterParticleValidationResult();
		
		/*! \brief Constructor with all args. */
		ParameterParticleValidationResult(const boost::shared_ptr< const ParameterParticle >& _p,
										const boost::shared_ptr< const mct::SummaryStatistic >& _s,
										double _d);
		
		~ParameterParticleValidationResult();
		
		
		/*! \brief Less than operator. 
		 * 
		 * a < b iff a.distance < b.distance.*/
		bool operator <(const ParameterParticleValidationResult& rhs) const;
		
		void setDistance(double d);
		
		/*! \brief Return a string representation of the result. */
		std::string toString() const;
	};
	
	std::ostream& operator<< (std::ostream &out, 
						const ParameterParticleValidationResult &pr);

} // end namespace aabc

#endif 

