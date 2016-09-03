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
\brief ParameterParticleValidatorGenericDataObject declarations.
*/

#ifndef _PARAMETER_PARTICLE_VALIDATOR_GENERIC_DATA_OBJECT_H
#define _PARAMETER_PARTICLE_VALIDATOR_GENERIC_DATA_OBJECT_H


#include "parameter_particle_validator.hpp"



namespace aabc {
	
	//forward declarations
	class GenericDataObjectSampler;	
	class DataObject;	
		
	/** @brief A class capable of 
	 validating a ParameterParticle for a DataObject tests.
	 
	 See the interface ParameterParticleValidator. 
	 	 
	 A ParameterParticleValidatorGenericDataObject
	 will have access to a reference MixtureTest and
	 will be able to take a ParameterParticle and use it to calculate a
	 one-number summary of the distance between sample(s) 
	 of MixtureTests generated
	 by the given ParameterParticle and the reference %DataObject.*/

	class ParameterParticleValidatorGenericDataObject 
								: public ParameterParticleValidator {
		public:
		
			/*! \brief Constructor. 
			 
			 \param _rs The reference object to use to check  
			 samples simulated.
			 \param _s A sampler.
			 \param _n_validation_reps The number of replications to 
			 produce using the sampler for validation against the
			 reference object.
			 */
			ParameterParticleValidatorGenericDataObject(
						const boost::shared_ptr < const DataObject >& _rs,
						const boost::shared_ptr < GenericDataObjectSampler >& _s,
						const size_t _n_validation_reps);
		
		
			virtual ~ParameterParticleValidatorGenericDataObject();
			
			/*! Get a mct::SummaryStatistic for a given ParameterParticle.
			 
			 This method must be implemented by 
			 concrete ParameterParticleValidators.  
			 
			 \param p The ParameterParticle to get the SummaryStatistic for.
			 \return A SummaryStatistic for the
			 ParameterParticle \a p.  */
			virtual boost::shared_ptr< const mct::SummaryStatistic >
						getSummaryStatistic(
				const boost::shared_ptr< const ParameterParticle>& p) const;
			
			/*! Get a mct::SummaryStatisticSet for a given 
			 ParameterParticleSet.
			 
			 \param pset The ParameterParticleSet to get
			  the SummaryStatisticSet for.
			 \return A SummaryStatisticSet for the
			 ParameterParticleSet \a pset.  
			 \post the SummaryStatisticSet returned contains
			 one SummaryStatistic every every ParameterParticle
			 in \a pset, in the same order as those ParameterParticles 
			 are in \a pset, ie the SummaryStatistic associated with
			 the i<sup>th</sup> %ParameterParticle is at the 
			 i<sup>th</sup> position in the SummaryStatisticSet.*/
			virtual boost::shared_ptr< mct::SummaryStatisticSet >
					getSummaryStatisticSet(
					const boost::shared_ptr< const ParameterParticleSet>& pset)
															const;
			
		
		protected:
		
			/*! \brief The reference object. */
			const boost::shared_ptr < const DataObject > referenceObj;
		
			/*! \brief The sampler. */
			const boost::shared_ptr < GenericDataObjectSampler > sampler;
		
		private:
		
			const size_t n_validation_reps;
		
	};

} // end namespace aabc

#endif /* _PARAMETER_PARTICLE_VALIDATOR_TEST_H */

