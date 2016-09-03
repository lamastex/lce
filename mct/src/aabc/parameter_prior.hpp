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
\brief ParameterPrior declarations.
*/

#ifndef _INC_PARAMETERPRIOR_H
#define _INC_PARAMETERPRIOR_H

#include <string>

namespace aabc {
	
	

	/** @brief An abstract class for a priors for
	 a single (one-dimensional) parameter for a model.
	 
	 A ParameterPrior can:
	 <ul>
	 <li>generate a random value according to the density,</li>
	 <li>and calculate the density at a given value,</li>
	 <li>and calculate the cdf at a given value.</li>
	 </ul>
	 A ParameterPrior can also send test output to standard output.
	 
	 */
	 

	class ParameterPrior {
		public:
			
			virtual ~ParameterPrior() {};
			
			/*! \brief Get the cdf at \a x.*/
			virtual double cdf(const double x) const = 0;
			
			/*! \brief Get the density at \a x.*/
			virtual double pdf(const double x) const = 0;
			
			/*! \brief Generate a random value.*/
			virtual double rand() const = 0;
			
			/*! \brief Output some test results to standard output.*/
			virtual void test() const;
			
			/*! \brief String representation of this.*/
			virtual std::string toString() const = 0;
			
		
	};

} // end namespace aabc

#endif /* _INC_PARAMETERPRIOR_H */

