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
\brief ParameterParticleValidator declarations.
*/

#ifndef _PARAMETER_PARTICLE_VALIDATOR_H
#define _PARAMETER_PARTICLE_VALIDATOR_H

#include "parameter_particle_validation_result_set.hpp"
#include "parameter_particle_validation_result.hpp"
#include "parameter_particle_set.hpp"
#include "parameter_particle.hpp"
#include "summary_statistic_set.hpp"
#include "summary_statistic.hpp"

#include <boost/smart_ptr.hpp>



namespace aabc {
	
		
	/** @brief An abstract type capable of 
	 validating a ParameterParticle.
	 
	 The type is configured with a reference mct::SummaryStatistic
	 and is capable of calculating the distance between
	 that reference and a SummaryStatistic associated with 
	 a ParameterParticle to be validated.
	 
	 The distance calculation can include a normalisation factor.  The
	 normalisation used is contained in the 
	 \link ParameterParticleValidator::normaliser normaliser\endlink.  
	 
	 Validation results are summarised in a 
	 ParameterParticleValidationResult structure.
	 
	 A concrete %ParameterParticleValidator 
	 will be able to take a ParameterParticle and use it to 
	 generate a mct::SummaryStatistic for that %ParameterParticle.
	 
	 \note A concrete %ParameterParticleValidator may use some
	 form of pseudo-random process to generate a SummaryStatistic 
	 for a given %ParameterParticle: there is no guarantee that
	 the same SummaryStatistic will generated on
	 repeated validations of the same or identical ParameterParticles.
	 */

	class ParameterParticleValidator {
		public:
		
			/*! Constructor.
			
			 \note The constructor should set the 
			 \link aabc::ParameterParticleValidator::normaliser normaliser\endlink
			 to some default value (e.g., all ones).
			 
			 \param _r A pointer to a mct::SummaryStatistic reference to calculate 
			 distances against.
			 \pre _r must be a non-null pointer.*/
			ParameterParticleValidator(
				const boost::shared_ptr < const mct::SummaryStatistic >& _r);
		
			virtual ~ParameterParticleValidator();
			
			/*! Get a validation result for a ParameterParticle.
			 
			 \param p The ParameterParticle to validate.
			 \return A ParameterParticleValidationResult for the
			 ParameterParticle \a p.	 */
			virtual boost::shared_ptr< ParameterParticleValidationResult >
				getValidationResult(
				const boost::shared_ptr< const ParameterParticle>& p)
																const;
			
			
			/*! Get a set of validation results for a 
			 ParameterParticleSet.
			 
			 \param pset The ParameterParticleSet to validate.
			 \return A ParameterParticleValidationResultSet for the
			 ParameterParticleSet \a pset.  */
			virtual boost::shared_ptr< ParameterParticleValidationResultSet >
				getValidationResultSet(
				const boost::shared_ptr< const ParameterParticleSet>& pset)
																const;

			/*! Use a ParameterParticleSet to set the 
			 \link ParameterParticleValidator::normaliser normaliser\endlink
			 and also then return the validation results set, using 
			 that normaliser, for the ParameterParticleSet.
			 
			 \param pset The ParameterParticleSet to use to set the
			 normaliser and then to validate.
			 \return A ParameterParticleValidationResultSet for the
			 ParameterParticleSet \a pset.  */
			virtual boost::shared_ptr< ParameterParticleValidationResultSet >
				setNormaliserAndGetValidationResultSet(
				const boost::shared_ptr< const ParameterParticleSet>& pset)
																	const;
			
			/*! Add the results from a ParameterParticleSet to the given
			  ParameterParticleValidationResultSet.
			 
			 This method does not update the 
			 \link ParameterParticleValidator::normaliser normaliser\endlink,
			 i.e. distances are normalised using the existing normaliser
			 or default normaliser.
			 
			 \param results A reference to a 
			 ParameterParticleValidationResultSet 
			 to add the new results to.
			\param pset The ParameterParticleSet to use generate
			new results to be added to \a results.
			 \return A reference to \a results.  */
			virtual boost::shared_ptr< ParameterParticleValidationResultSet >&
				updateValidationResultSet(
				boost::shared_ptr< ParameterParticleValidationResultSet >& results,
				const boost::shared_ptr< const ParameterParticleSet>& pset)
													const;

			/*! Get a mct::SummaryStatistic for a given ParameterParticle.
			 
			 This method must be implemented by 
			 concrete ParameterParticleValidators.  
			 
			 \param p The %ParameterParticle to get the SummaryStatistic for.
			 \return A SummaryStatistic for the
			 %ParameterParticle \a p.  */
			virtual boost::shared_ptr< const mct::SummaryStatistic >
					getSummaryStatistic(
					const boost::shared_ptr< const ParameterParticle>& p)
															const = 0;

			/*! Get a mct::SummaryStatisticSet for a given 
			 ParameterParticleSet.
			 
			 This method must be implemented by 
			 concrete ParameterParticleValidators.  
			 
			 \param pset The %ParameterParticleSet to get
			  the SummaryStatisticSet for.
			 \return A SummaryStatisticSet for the
			 %ParameterParticleSet \a pset.  
			 \post the SummaryStatisticSet returned contains
			 one SummaryStatistic every every ParameterParticle
			 in \a pset, in the same order as those ParameterParticles 
			 are in \a pset, ie the SummaryStatistic associated with
			 the i<sup>th</sup> %ParameterParticle is at the 
			 i<sup>th</sup> position in the SummaryStatisticSet.*/
			virtual boost::shared_ptr< mct::SummaryStatisticSet >
					getSummaryStatisticSet(
					const boost::shared_ptr< const ParameterParticleSet>& pset)
															const = 0;
			
			/*! Use a ParameterParticleSet to set the 
			 \link ParameterParticleValidator::normaliser normaliser\endlink, and 
			 return a copy of that new normaliser to the caller.
			 
			 This method sets the normalising values to be the sample
			 standard deviations of the 
			 \link mct::SummaryStatistic mct::SummaryStatistics\endlink
			 generated for the %ParameterParticles in the 
			 %ParameterParticleSet. 
			 			 
			 \param pset The ParameterParticleSet to use to set the
			 normaliser.
			 \return A pointer to a copy of the normaliser.  */
			virtual boost::shared_ptr < mct::SummaryStatistic > 
				setNormaliser(
				const boost::shared_ptr< const ParameterParticleSet>& pset)
																const;

			/*! Set the 
			 \link ParameterParticleValidator::normaliser normaliser\endlink.
			 
			 \param s The mct::SummaryStatistic to use for the normaliser.*/
			virtual void 
				setNormaliser(
				const boost::shared_ptr < const mct::SummaryStatistic >& s)
																const;


			/*! Get a pointer to a copy of the 
			 \link ParameterParticleValidator::normaliser normaliser\endlink.
			 
			 \return A pointer to a copy of the normaliser.  */
			virtual boost::shared_ptr < mct::SummaryStatistic > 
							getNormaliser() const;
			
		
				
				
		protected:
			
			/*! \brief The reference statistics. */
			const boost::shared_ptr < const mct::SummaryStatistic > 
														referenceStats;
														
			/*! \brief The normalising statistics. 
			 
			 The normaliser is used in
			 calculation of validation distances.
			 
			 \note The normaliser should be initialised to some 
			 default in the constructor.*/
			mutable boost::shared_ptr < const mct::SummaryStatistic > 
														normaliser;
		
			virtual boost::shared_ptr< ParameterParticleValidationResultSet >&
				_updateValidationResultSet(
				boost::shared_ptr< ParameterParticleValidationResultSet >& results,
				const boost::shared_ptr< const ParameterParticleSet>& pset,
				const boost::shared_ptr< mct::SummaryStatisticSet >& ss_set)
															const;

			/*! \brief Get distance between the reference 
			 mct::SummaryStatistic 
			 \link ParameterParticleValidator::referenceStats referenceStats\endlink
			 and a given SummaryStatistic.
			 
			 If each SummaryStatistic has \f$ nstat \f$ dimensions, and 
			 \f$ s^{ref}_j \f$ is the j<sup>th</sup> value in the
			 reference SummaryStatistic, and \f$ s^{\prime}_j \f$ 
			 is the j<sup>th</sup> value in the
			 supplied SummaryStatistic \a ss, then the distance 
			 is calculated as:
			  
			 \f[ distance = \sqrt{\sum^{nstat}_{j=1} \left( \frac{  s^{ref}_j - s^{\prime}_j }{n_j} \right)^2} \f]
			 
			 where \f$ n_j \f$ is the normalisation value for the j<sup>th</sup> statistic.
			  
			 \pre \link ParameterParticleValidator::normaliser normaliser\endlink
			 should be a non-null pointer.
			 \param ss The SummaryStatistic to calculate the distance for.
			 \return The distance between \a ss and the 
			 \link ParameterParticleValidator::referenceStats referenceStats\endlink.*/
			double getDistance(
				const boost::shared_ptr < const mct::SummaryStatistic >& ss)
				const;
		
														
			virtual boost::shared_ptr < mct::SummaryStatistic > 
				_setNormaliser(
				const boost::shared_ptr< const mct::SummaryStatisticSet>& ss_set)
																const;

			
	};

} // end namespace aabc

#endif 

