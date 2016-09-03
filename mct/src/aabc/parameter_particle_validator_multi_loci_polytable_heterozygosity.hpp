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
\brief ParameterParticleValidatorMultiLociPolyTableHeterozygosity declarations.
*/

#ifndef _PARAMETER_PARTICLE_VALIDATOR_MULTI_LOCI_POLYTABLE__HETEROZYGOSITY_H
#define _PARAMETER_PARTICLE_VALIDATOR_MULTI_LOCI_POLYTABLE__HETEROZYGOSITY_H


#include "parameter_particle_validator.hpp"
#include "multi_loci_polytable_sampler.hpp"
#include "multi_loci_polytable.hpp"


namespace aabc {
	
	/** @brief A class capable of 
	 validating a ParameterParticle for mct::MultiLociPolyTable.
	 
	 See the interface ParameterParticleValidator. 
	 	 
	 A ParameterParticleValidatorMultiLociPolyTableHeterozygosity
	 will have access to a reference mct::MultiLociPolyTable and
	 will be able to take a ParameterParticle and use it to a
	 calculate a mct::SummaryStatistic based on Heterozygosity
	 for that particle.
	 
	 A ParameterParticleValidatorMultiLociPolyTableHeterozygosity is 
	 configured with a 
	 mct::MultiLociPolyTable reference object and a
	 mct::MultiLociPolyTableSampler which it uses with the
	 given ParameterParticle to generate the
	 mct::MultiLociPolyTableSet used to calculate
	 a mct::SummaryStatistic using heterozygosity for the particle.
	 
	 The number of repetitions used to calculate 
	 the SummaryStatistic to be associated with %ParameterParticle 
	 is configured when the validator is constructed.  Usually this is
	 just one, but the summary statistic could also 
	 be calculated as the average over summary statistics from
	 a number of replications from the sampler.
	 
	 */

	class ParameterParticleValidatorMultiLociPolyTableHeterozygosity 
								: public ParameterParticleValidator {
		public:
		
			/*! \brief Constructor. 
			 
			 \param _rs A pointer to the reference object to use to check  
			 samples simulated.
			 \param _s A pointer to a mct::MultiLociPolyTableSampler.
			 \param _n_validation_reps The number of repetitions
			 to be used for calculating each SummaryStatistic.			 */
			ParameterParticleValidatorMultiLociPolyTableHeterozygosity(
				const boost::shared_ptr < const mct::MultiLociPolyTable >& _rs,
				const boost::shared_ptr < mct::MultiLociPolyTableSampler >& _s,
				const size_t _n_validation_reps);
		
		
			virtual ~ParameterParticleValidatorMultiLociPolyTableHeterozygosity();
			
			/*! Get a mct::SummaryStatistic for a given ParameterParticle.
			 
			 This method must be implemented by 
			 concrete ParameterParticleValidators.  
			 
			 This ParameterParticleValidator calculates the summary 
			 statistic using just heterozygosity.
			 
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
			 \return A SummaryStatisticSet using 
			 heterozygosity for the
			 ParameterParticleSet \a pset.  
			 \post the SummaryStatisticSet returned contains
			 one SummaryStatistic (containing only heterozygosity)
			 for every every ParameterParticle
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
			const boost::shared_ptr < const mct::MultiLociPolyTable > referenceObj;
			
			/*! \brief The sampler. */
			const boost::shared_ptr < mct::MultiLociPolyTableSampler > sampler;
		
		
		private:
		
			/*! \brief The number of repetitions to use to calculate
			 a SummaryStatistic for a ParameterParticle. 
			 
			 Configured when the validator is constructed and 
			 cannot be changed.*/
			const size_t n_validation_reps;
		
	};

} // end namespace aabc

#endif 

