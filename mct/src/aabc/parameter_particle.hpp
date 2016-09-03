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
\brief ParameterParticle declarations.
*/

#ifndef _INC_PARAMETER_PARTICLE_H
#define _INC_PARAMETER_PARTICLE_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include <boost/smart_ptr.hpp>
#include <boost/function.hpp>

#include <vector>
#include <string>


namespace aabc {
	
	

	/** @brief A class for a set of parameters for the model
	 * wrapped together with the weight given to that set of parameters
	 * in an importance sampling context.
	 * 
	 * ParameterParticle at the moment is basically a data object, 
	 * holding the set of parameter values and the weight with very 
	 * little behaviour of its own.
	 * 
	 * For example, in a population genetics contact, a 
	 * ParameterParticle mightt includes elements for migration rates,
	 * recombination rates, etc.
	  */
	 

	class ParameterParticle {
		public:
		
			/*! Iterator into the ordered collection of parameter values.*/
			typedef std::vector< double >::const_iterator
											const_iterator;
		
			/*! \brief Constructor.
			
			\param values The values for the parameters.*/
			ParameterParticle(const std::vector < double >& values);
			
			/*! \brief Constructor.
			
			\param a[] The start of an array of doubles.
			\param _n The number of elements \a a[].
			*/
			ParameterParticle(const double a[],
										const size_t _n);
			
			/*! \brief Constructor.
			
			\param _n The number of parameter elements.
			*/
			explicit ParameterParticle(const size_t _n);
			
			
			/*! \brief Constructor.
			
			\param v A gsl::vector of values to use for the particle.
			*/
			explicit ParameterParticle(const gsl_vector * const v);



			/*! \brief Copy constructor.
			
			\param other The ParameterParticle to copy.
			*/
			ParameterParticle(const ParameterParticle& other);
			
			/*! \brief Equals operator.
			
			\param tmp The ParameterParticle to set this equal to.
			*/
			ParameterParticle& operator=(ParameterParticle tmp);
			
			virtual ~ParameterParticle();
			
			/*! \brief Get a gsl_vector representation of the particle.*/
			gsl_vector * get_gsl_vector() const;
			
			/*! \brief Get a std::vector representation of the particle.*/
			std::vector < double > get_vector() const;
			
			/*! \brief Set the element at index.*/
			double& at(const size_t index);
			
			/*! \brief Get the element at index.*/
			const double& at(const size_t index) const;
			
			
			
			/*! \brief Get the weight.*/
			virtual double getWeight() const;
			
			/*! \brief Set the weight.*/
			virtual void setWeight(const double w) const;
			
			/*! \brief Get the number of parameter elements.*/
			size_t size() const;	
			
			/*! \brief Return a const iterator to the first 
			element in the particle.*/
			const_iterator begin() const;
			
			/*! \brief Return a const iterator to the
			'just past the end' element in the particle.
			
			\note The idiom for iterating through a 
			particle p is 
			for (ParameterParticle::const_iterator it = p.begin();
						it < p.end();
						++it) {...}.
			*/
			const_iterator end() const;
			
			/*! \brief Subtraction operator. */
			ParameterParticle operator-(
					const ParameterParticle& b) const;
			
			/*! \brief Get a fancy string representation of this.
			
			String is formatted as:
			([param1]\\t[param2]\\t...\\t[paramn]) weight: [weight]
			*/
			std::string toString() const;
			
			/*! \brief Get a plain string representation of this
			 String is formatted as:
			 [param1]\\t[param2]\\t...\\t[paramn] 
			 This format is suitable as a line for a file
			*/
			std::string plainString() const;
			
			/*! \brief Get a plain string representation of this, rescaled.
			
			Form as for plainString(), but the \f$ i^{th} \f$ element
			in thsi is multiplied by the \f$ i^{th} \f$ element
			in \a mults.  
			
			One line per ParameterParticle.
			
			\param mults the collection of rescaling multipliers to 
			be applied to the values of the elements in this
			to give the string result. 
			
			\pre mults.size() is the same as the size of this.
			*/
			std::string plainString(const std::vector < double >& mults) const;
						
		protected:
			ParameterParticle();

				
		private:
		
			/*! \brief Number of parameters. */
			size_t n; 
			
			/*! \brief A container of doubles. */
			std::vector< double > parameters;
			
			/*! \brief Weight for this set of parameters. */
			mutable double weight; 
			
			
	};
	
	std::ostream& operator<< (std::ostream &out, const ParameterParticle &p);

} // end namespace aabc

#endif /* _INC_PARAMETERPARTICLE_H */

