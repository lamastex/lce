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
\brief ParameterParticleValidationResultSet declarations.
*/

#ifndef _INC_PARAMETERPARTICLE_VALIDATION_RESULT_SET_H
#define _INC_PARAMETERPARTICLE_VALIDATION_RESULT_SET_H

#include "parameter_particle_validation_result.hpp"
#include "parameter_particle_set.hpp"
#include "parameter_particle.hpp"
#include "summary_statistic_set.hpp"
#include "summary_statistic.hpp"

#include <boost/smart_ptr.hpp>

#include <vector>
#include <string>


namespace aabc {

	
	/** @brief A class for an ordered collection of 
	 * \link aabc::ParameterParticleValidationResult ParameterParticleValidationResults\endlink.
	 */
	 class ParameterParticleValidationResultSet {
		
		public :
		
			typedef std::vector < boost::shared_ptr< ParameterParticleValidationResult > >::iterator 
					iterator;
		
			/*! \brief No-args constructor.*/
			ParameterParticleValidationResultSet();
			
			/*! \brief Constructor.
			
			Constructs an empty structure in which to store 
			\link aabc::ParameterParticleValidationResult ParameterParticleValidationResults\endlink,
			but reserves enough memory space for
			\a reservedSpace elements.  This is useful if the 
			set is to be used to store large numbers of
			objects because it can avoid or reduce memory
			reallocations when add is called.
			
			\param reservedSpace The number of elements to reserve
			space for in memory.
			\note The reserved space is not a limit on the capacity
			of the set nor does it stuff the container with empty
			elements, it just ensures that the set's memory 
			allocation is sufficient for that number of elements
			to be added without reallocation.  If more elements are 
			then added, memory will be reallocated.
			*/
			explicit ParameterParticleValidationResultSet(size_t reservedSpace);
			
			
			virtual ~ParameterParticleValidationResultSet();
			
			
			/*! \brief Add a ParameterParticleValidationResult to the collection.
			
			\param vr The  ParameterParticleValidationResult to add.*/
			void add(const boost::shared_ptr< ParameterParticleValidationResult >& vr);
			
			/*! \brief Get iterator to the element at the beginning of 
			 the collection.*/
			ParameterParticleValidationResultSet::iterator begin();
			
			/*! \brief Get iterator to the element past the end of 
			 * the collection.*/
			ParameterParticleValidationResultSet::iterator end();
			
			boost::shared_ptr< ParameterParticleSet >  getParameterParticleSet();
			
			boost::shared_ptr< mct::SummaryStatisticSet >  getSummaryStatisticSet();
			
			/*! \brief Get the number of elements in the collection.
			 
			 \return The number of elements in the collection.*/
			size_t size() const;
			
			/*! \brief Get whether the collection is empty.
			  
			 \return True if the collection is empty, false otherwise.*/
			bool empty() const;

			/*! \brief Ensure at least \a reservedSpace capacity 
			 in the collection.
			 
			 This is not a limit on the collection size but the number
			 of elements the collection can hold without reallocating
			 memory
			 
			 \param reservedSpace The number of elements the collection
			 can hold without memory reallocation being necessary.*/
			void reserve(size_t reservedSpace);
			
			/*! \brief Discard all elements of the collection with 
			 \link aabc::ParameterParticleValidationResult::distance distance\endlink
			 >= \a e.
			  
			 \param e The criteria used for discards:  all elements with
			 \link aabc::ParameterParticleValidationResult::distance distance\endlink >= \a e will be dropped from the collection.
			 \post The collection will only contain elements with 
			 \link aabc::ParameterParticleValidationResult::distance distance\endlink
			  > e.  The original ordering of the elements may have been altered.  
			  The capacity of the collection will not have been changed.*/
			virtual void discardByDistance(double e);
			
			/*! \brief Discard the largest \n elements of the collection by 
			 \link aabc::ParameterParticleValidationResult::distance distance\endlink
			 and shrink the capacity to fit the new number of elements
			  
			 \param n The number of elements to discard.
			 \pre n < size()
			 \post The n elements with the largest
			 \link aabc::ParameterParticleValidationResult::distance distance\endlink
			 will have been discarded.  The original ordering of the elements
			  may have been altered.  
			  The capacity of the collection will have been shrunk 
			  to fit the new number of elements.*/
			virtual void discardN(size_t n);
		
			/*! \brief Get the distance corresponding to the 
			 \a percentile of the distances of all the elements
			 in the collection.
			 
			 \note The method used to calculate the percentile 
			 uses interpolation and is supposed to be equivalent
			 to the method used for gsl_stats_quantile_from_sorted_data
			 (see the GNU Scientific Library).
			 
			 if \a percentile = \f$ p \f$, \f$ 0 \leq p \leq 1 \f$,
			 
			 then
			 
			 result = \f$ (1-\delta) x_{[i]} + \delta x_{[i+1]} \f$
			 
			 where 
			 <ul>
			 <li>\f$ x_{[i]} \f$ is the value in position \f$ i \f$ of 
			 the ordered distances (\f$ i = 0, 1, \ldots, n-1 \f$, 
			 \f$ n = \f$ size of data set),</li>
			 <li>\f$ i = \lfloor{ (n-1)p }\rfloor \f$</li>
			 <li>\f$ \delta = (n-1)p - i \f$</li>
			 </ul>
			 
			 \param percentile The percentile distance to find.
			 \return the value of the \a percentile of the distances
			 of the elements in the collection.
			 \post The original ordering of the elements
			  may have been altered.  The capacity of the collection
			  will be unchanged.*/
			
			virtual double getDistancePercentile(double percentile);
		
		
		protected:
		
			explicit ParameterParticleValidationResultSet(const ParameterParticleValidationResultSet& other);
			
			ParameterParticleValidationResultSet& operator=(ParameterParticleValidationResultSet tmp);
			
			/*! \brief The ParameterParticleValidationResults. */
			std::vector < boost::shared_ptr< ParameterParticleValidationResult > > validation_results;
		

		private:
			
			
	};
	
	//std::ostream& operator<< (std::ostream &out, 
						//const ParameterParticleValidationResultSet &vrs);

} // end namespace aabc

#endif 

