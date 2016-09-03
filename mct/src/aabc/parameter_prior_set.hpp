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
\brief ParameterPriorSet declarations.
*/

#ifndef _INC_PARAMETERPRIORSET_H
#define _INC_PARAMETERPRIORSET_H

#include "parameter_prior.hpp"
#include "parameter_particle_set.hpp"
#include "parameter_particle.hpp"

#include <boost/smart_ptr.hpp>

#include <vector>
#include <string>


namespace aabc {
	
	/** @brief A class for an ordered collection of 
	 \link aabc::ParameterPrior ParameterPriors\endlink for
	 a parameter vector for a model.
	 */
	 

	class ParameterPriorSet {
		public:
		
			/*! \brief Constructor. 
			
			\param _n The number of priors to be held.*/
			explicit ParameterPriorSet(size_t _n);
			
			/*! \brief Copy constructor. 
			
			\param other The ParameterPriorSet to copy.
			*/
			ParameterPriorSet(const ParameterPriorSet& other);
						
			virtual ~ParameterPriorSet();
			
			/*! \brief Get const reference to the index<sup>th</sup>
			ParameterPrior pointer in the collection.*/
			const boost::shared_ptr< ParameterPrior >& at(const size_t index) const;
			
			/*! \brief Get reference to the index<sup>th</sup>
			ParamaterPrior pointer in the collection.*/
			boost::shared_ptr< ParameterPrior >& at(const size_t index);
			
			/*! \brief Get the density of a ParameterParticle \a p.*/
			virtual double pdf(
				const boost::shared_ptr < const ParameterParticle >& p) const;
			
			/*! \brief Get the cdf of x using the 
			index<sup>th</sup> prior.*/
			virtual double cdfFrom(const size_t index,
									const double x) const;
			
			/*! \brief Get the density of x using the 
			index<sup>th</sup> prior.*/
			virtual double pdfFrom(const size_t index,
									const double x) const;
			
			/*! \brief Generate a single parameter particle.*/
			virtual boost::shared_ptr < const ParameterParticle >
					generateParameterParticle() const;
			
			/*! \brief Generate a set of parameter particles.*/
			virtual boost::shared_ptr < ParameterParticleSet >
					generateParameterParticleSet(size_t n_reps) const;
			
			/*! \brief Generate a random value using the 
			index<sup>th</sup> prior.*/
			virtual double randFrom(const size_t index) const;
			
			/*! \brief Test the index<sup>th</sup> prior.*/
			virtual void testFrom(const size_t index) const;
			
			/*! \brief Test all the priors.*/
			virtual void testAll() const;
			
			/*! \brief Get the number of prior elements.*/
			virtual size_t size() const;	
			
			/*! \brief Get a string representation of this prior set.*/
			virtual std::string toString() const;	
						
						
		protected:
			ParameterPriorSet() : n_params(0) {}
			
			ParameterPriorSet& operator=(ParameterPriorSet tmp);
			
			
		private:
		
			/*! \brief Number of priors. */
			const size_t n_params;
			
			/*! \brief A container of pointers to ParameterPriors. */
			std::vector< boost::shared_ptr< ParameterPrior > > priors;
			
			
			
	};

} // end namespace aabc

#endif /* _INC_PARAMETERPRIORSET_H */

