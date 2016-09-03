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
\brief ParameterPriorNormal declarations.
*/

#ifndef _INC_PARAMETERPRIORNORMAL_H
#define _INC_PARAMETERPRIORNORMAL_H

#include "parameter_prior.hpp"
#include "prng.hpp"

#include <boost/function.hpp>

#include <string>



namespace aabc {
	
	

	/** @brief A class for normal (gaussian) priors for
	 a single (one-dimensional) parameter for a model.
	 
	 See ParameterPrior interface.
	 */
	 

	class ParameterPriorNormal : public ParameterPrior{
		public:
		
			
			/*! \brief Constructor. 
			
			\param _m The mean \f$ \mu \f$.
			\param _s The standard deviation \f$\sigma \f$.
			\param _r_ptr A pointer to a random number generator.
			*/
			ParameterPriorNormal(
					double _m, double _s,
					const boost::shared_ptr < mct::PRNGen >& _r_ptr);
			
			virtual ~ParameterPriorNormal();
			
			/*! \brief Get the cdf at \a x.*/
			virtual double cdf(const double x) const;
			
			/*! \brief Get the density at \a x.*/
			virtual double pdf(const double x) const;
			
			/*! \brief Generate a random value.*/
			virtual double rand() const;
			
			/*! \brief String representation of this.*/
			virtual std::string toString() const;
						
						
		protected:
		
			ParameterPriorNormal();
			
			ParameterPriorNormal(const ParameterPriorNormal& other);
			
			ParameterPriorNormal& operator=(ParameterPriorNormal tmp);
			
			
		private:
		
			/*! \brief The mean mu. */
			const double m;
			
			/*! \brief The standard deviation sigma. */
			const double s;
		
		
			/*! \brief The cdf function. */
			const boost::function<double (double x)> cdfFunc;
				
			/*! \brief The pdf function. */
			const boost::function<double (double x)> pdfFunc;
			
			/*! \brief The prng. */
			const boost::shared_ptr < mct::PRNGen > r_rand_ptr;
	};

} // end namespace aabc

#endif /* _INC_PARAMETERPRIORNORMAL_H */

