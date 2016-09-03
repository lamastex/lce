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
\brief GenericDataObjectSampler declarations.
*/

#ifndef _GENERIC_DATA_OBJECT_SAMPLER_H
#define _GENERIC_DATA_OBJECT_SAMPLER_H


#include "parameter_particle_set.hpp"
#include "parameter_particle.hpp"
#include "generic_data_object_set.hpp"

#include <boost/smart_ptr.hpp>


namespace aabc {
	
	
	/** @brief An abstract type capable of 
	 generating a GenericDataObjectSet tests.
	 
	 */

	class GenericDataObjectSampler {
		public:
		
			
			virtual ~GenericDataObjectSampler();
			
			/*! \brief  Generate a collection of data using 
			 a given ParameterParticle.
			  
			\param n_reps The number of replications to generate.
			\param p The ParameterParticle to use to generate the
			simulations.
			\return A pointer to the GenericDataObjectSet 
			holding the simulations.			*/
			virtual boost::shared_ptr < GenericDataObjectSet >
				sample(const size_t n_reps,
				const boost::shared_ptr< const ParameterParticle> p) const;
		
			/*! \brief  Generate a collection of data using 
			 a given ParameterParticleSet.
			  
			\param n_reps_per_particle The number of replications to 
			generate for each particle in the ParameterParticleSet.
			\param pset The ParameterParticleSet to use to generate the
			simulations.
			\return A pointer to the GenericDataObjectSet 
			holding the simulations.			*/
			virtual boost::shared_ptr < GenericDataObjectSet > 
				sample(const size_t n_reps_per_particle, 
				const boost::shared_ptr< ParameterParticleSet>& pset) 
															const;
			
		
		protected:
		
			virtual boost::shared_ptr < GenericDataObjectSet >& 
				_sample(
					boost::shared_ptr < GenericDataObjectSet >& 
							gsetptr,
					const size_t n_reps, 
					const boost::shared_ptr< const ParameterParticle>& p) const = 0;

		
			
			
		
	};

} // end namespace aabc

#endif /* _SAMPLER_TEST_H */

