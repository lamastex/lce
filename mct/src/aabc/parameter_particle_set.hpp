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
\brief ParameterParticleSet declarations.
*/

#ifndef _INC_PARAMETER_PARTICLE_SET_H
#define _INC_PARAMETER_PARTICLE_SET_H

#include "parameter_particle.hpp"
#include "summary_statistic_set.hpp"


#include <boost/function.hpp>
#include <boost/smart_ptr.hpp>

#include <gsl/gsl_matrix.h>

#include <vector>
#include <string>


namespace aabc {
	
	/** @brief A class for a set of
	\link aabc::ParameterParticle ParameterParticles\endlink for 
	the population Monte Carlo approximate Bayesian computation routine.
	
	A %ParameterParticleSet manages an ordered collection 
	of pointers to 
	\link aabc::ParameterParticle ParameterParticles\endlink.
	
	*/
	 

	class ParameterParticleSet {
		public:
		
			/*! Iterator into the ordered collection of pointers to 
			\link aabc::ParameterParticle ParameterParticles\endlink.*/
			typedef 
			std::vector< 
					boost::shared_ptr< const ParameterParticle > >
					::const_iterator
					const_iterator;
		
			/*! \brief No-argument constructor. */
			ParameterParticleSet();
			
			/*! \brief Constructor.
			
			Constructs an empty structure in which to store 
			\link aabc::ParameterParticle ParameterParticles\endlink,
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
			explicit ParameterParticleSet(size_t reservedSpace);
			
			/*! \brief Constructor with collection of ParameterParticle pointers. */
			explicit ParameterParticleSet(
				const std::vector< boost::shared_ptr< const ParameterParticle > >& vec);
			
			/*! \brief Copy constructor. */
			ParameterParticleSet(const ParameterParticleSet& other);
			
			/*! \brief Equals operator. */
			ParameterParticleSet& operator=(ParameterParticleSet tmp);
						
			/*! \brief Destructor. */			
			virtual ~ParameterParticleSet();
									
			/*! \brief Get const reference to the index<sup>th</sup>
			ParameterParticle pointer in the collection.*/
			virtual const boost::shared_ptr< const ParameterParticle >& at(const size_t index) const;
			
			/*! \brief Get reference to the index<sup>th</sup>
			ParameterParticle pointer in the collection.*/
			virtual boost::shared_ptr< const ParameterParticle >& at(const size_t index);
			
			/*! \brief Add a ParameterParticle pointer to the collection.
			
			\param toAdd Pointer to the ParameterParticle to add.*/
			virtual void add(const boost::shared_ptr< const ParameterParticle >& toAdd);
			
			/*! \brief Add all the  ParameterParticle pointers
			from another ParameterParticleSet to the collection.
			
			\param toAdd Pointer to the ParameterParticleSet 
			whose contents are to be added.*/
			virtual void add(const boost::shared_ptr < ParameterParticleSet >& toAdd);
			
			/*! \brief Get the number of 
			\link aabc::ParameterParticle ParameterParticles\endlink in 
			the collection.*/
			size_t size() const;
			
			/*! \brief Get whether the container is empty. */
			bool empty() const;
			
			
			/*! \brief Return whether the sizes of the 
			 \link aabc::ParameterParticle ParameterParticles\endlink
			 in the collection are consistent.
			 
			 \return true if the collection is empty or if the sizes of
			 of the contained 
			 ParameterParticles
			 are all the same, false otherwise. */
			bool sizeConsistencyCheck() const;
			
			/*! \brief Make and return a SummaryStatisticSet
			holding the values from the particles in this set.
			* 
			The i<sup>th</sup> SummaryStatistic in the set returned 
			holds the values
			from the i<sup>th</sup> ParameterParticle in this. 
			
			\return a pointer to a %SummaryStatisticSet holding	
			values from the ParameterParticles in this.*/				
			virtual boost::shared_ptr < mct::SummaryStatisticSet >
				makeSummaryStatisticSet() const;
			
			/*! \brief Get all the values in the 
			\link ParameterParticle ParameterParticles\endlink
			held.
			
			 If this contains n
			 \link aabc::ParameterParticle ParameterParticles\endlink
			 with at most nparam elements in each, then the outer container in
			 the returned structure contains nparam inner containers,
			 and each inner container contains n values, and the
			 inner containers are in the same order as the 
			 elements are in the ParameterParticles.  ie, the 
			 i<sup>th</sup> value in the j<sup>th</sup> 
			 %ParameterParticle in the set will be the j<sup>th</sup>
			 element of the i<sup>th</sup> inner container. If any
			 %ParameterParticle in the set contains less than
			 nparam elements then the first nparam inner containers
			 will contain the values from that %ParameterParticle and
			 the remaining inner containers will have a NaN at
			 that position.
			\return a pointer to a container of containers of	
			values in the ParameterParticles in this.*/		
			boost::shared_ptr < std::vector < std::vector < double > > >
				getAllValues() const;
			
			
			
			/*! \brief Get the variance-covariance matrix for the N
			\link aabc::ParameterParticle ParameterParticles\endlink 
			in the collection.
			
			Each ParameterParticle is an observation; 
			the variance-covariance matrix is n x n, where n is the number
			of parameter elements in each ParameterParticle.
			 
			The element in row \f$ i \f$ and column \f$ j \f$ of the 
			variance-covariance
			matrix is \f$ q_{ij} \f$ where 
			
			\f[
			q_{ij} = \frac{1}{N}\sum_{k=1}^N(x_{ik}-\bar{x}_i)(x_{jk}-\bar{x}_j)
			\f]
			 
			\note The variance-covariance matrix is calculated using the MLE 
			for the variance (covariance), i.e. division by N not N-1
			where N is the number of ParameterParticles in the container. 
			*/
			virtual gsl_matrix* covariance() const;	
			
			/*! \brief Get the weighted variance-covariance matrix for 
			the 
			\link aabc::ParameterParticle ParameterParticles\endlink 
			in the collection.
			
			Each ParameterParticle is an observation; 
			the weighted variance-covariance matrix is n x n, 
			where n is the number
			of elements in each ParameterParticle.
			 
			The element in row \f$ i \f$ and column \f$ j \f$ of the 
			weighted variance-covariance
			matrix is \f$ q_{ij} \f$ where 
			
			\f[
			q_{ij} = \sum_{k=1}^Nw_k(x_{ki}-\bar{x}_i)(x_{kj}-\bar{x}_j)
			\f]
			
			and \f$ w_k \f$ is the weight of the k<sup>th</sup>
			ParameterParticle in the collection.
			 
			\note The weighted variance-covariance matrix is calculated 
			using the MLE 
			for the variance (covariance), i.e. division by N not N-1
			where N is the number of ParameterParticles in the container. 
			 
			\note Because of the multiple stages of the calculations 
			using doubles for the weighted variance-covariance matrix,
			accuracy in the lower decimal places will be lost. In
			particular, elements which should be 0.0 may end up 
			as very very small non-zero values.  
			*/
			virtual gsl_matrix* weightedCovariance() const;	
			
			/*! \brief Get the mean for  
			the 
			\link aabc::ParameterParticle ParameterParticles\endlink 
			in the collection. */
			virtual gsl_vector* mean() const;
			
			/*! \brief Get the weighted mean for  
			the 
			\link aabc::ParameterParticle ParameterParticles\endlink 
			in the collection. */
			virtual gsl_vector* weightedMean() const;
			
			
			/*! \brief Return a const iterator to the first 
			element in the set.*/
			const_iterator begin() const;
			
			/*! \brief Return a const iterator to the
			'just past the end' element in the set.
			
			\note The idiom for iterating through a 
			%ParameterParticleSet pset is 
			for (ParameterParticleSet::const_iterator it = pet.begin();
						it < pset.end();
						++it) {...}.
			*/
			const_iterator end() const;
			
			/*! \brief Randomly pick a ParameterParticle index from the 
			container with probability proportional relative to
			ParameterParticle weights.
			
			\param rand_num_func Function which generates a uniform 
			distribution in the range (0, 1].
			\return An index identifying one of the container of
			\link aabc::ParameterParticle ParameterParticles\endlink.
			 */
			virtual size_t pickOne(
					boost::function<double ()> rand_num_func) const;
			
			/*! \brief Get a fancy string representation of this.
			
			One line per ParameterParticle.
			Uses ParameterParticle::toString() to format each line.
			*/
			virtual std::string toString() const;
			
			/*! \brief Get a plain string representation of this.
			
			One line per ParameterParticle.
			\see ParameterParticle::plainString().
			*/
			virtual std::string plainString() const;
			
			/*! \brief Get a plain string representation of this, rescaled.
			
			For each %ParameterParticle in the set, the string gives
			the value of the \f$ i^{th} \f$ element
			in the particle is multiplied by the \f$ i^{th} \f$ element
			in \a mults.  ie, \a mults must contain as many elements as
			there are elements in the particles in the set.
			
			One line per ParameterParticle.
			
			\param mults the collection of rescaling multipliers to 
			be applied to the values of the elements in the particles in this
			to give the string result. 
			\see ParameterParticle::plainString(const std::vector < double >& mults).
			\pre mults.size() is the same as the size of each %ParameterParticle
			in this.
			*/
			virtual std::string plainString(const std::vector < double >& mults) const;
			
			/*! \brief Output this to a file.
			
			One line in the output per ParameterParticle.
			\see ParameterParticle::plainString() for format.
			\param s The filename for the output (which can include a path).
			\param append An indicator for whether output is appended to
			the file (if it exists already) or overwrites existing file content.
			\post A file called \a s exists containing output from this.*/
			virtual void outputToFile(const std::string& s,
								bool append = true) const;

			/*! \brief Output this to a file, rescaled.
			
			For each %ParameterParticle in the set, the output gives the
			value of the \f$ i^{th} \f$ element
			in the particle is multiplied by the \f$ i^{th} \f$ element
			in \a mults.  ie, \a mults must contain as many elements as
			there are elements in the particles in the set.
			
			One line in the output per ParameterParticle.
			
			\see ParameterParticle::plainString() for format.
			\param mults the collection of rescaling multipliers to 
			be applied to the values of the elements in the particles in this
			to give the output. 
			\param s The filename for the output (which can include a path).
			\param append An indicator for whether output is appended to
			the file (if it exists already) or overwrites existing file content.
			\post A file called \a s exists containing output from this.
			\pre mults.size() is the same as the size of each %ParameterParticle
			in this.			*/
			virtual void outputToFile(const std::vector < double >& mults,
								const std::string& s,
								bool append) const;

					
		protected:
			
			
			/*accumulate all the weights for all the 
			ParameterParticles in the container. */
			virtual double accumulateWeights() const;
			
			virtual const std::vector< boost::shared_ptr < const ParameterParticle > >&
				_getContainer() const;
				
					
		private:
		
			/* \brief A container of pointers to 
			 \link aabc::ParameterParticle ParameterParticles\endlink. */
			//std::vector< ParameterParticle > container;
			std::vector< boost::shared_ptr< const ParameterParticle > > 
												container;
			
			
			
	};
	
	std::ostream& operator<< (std::ostream& out, 
							const ParameterParticleSet& p);

} // end namespace aabc

#endif /* _INC_PARAMETERPARTICLEGENERATION_H */

