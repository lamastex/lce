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
\brief PRNGenGSL declarations.
*/

#ifndef _INC_PRNG_GSL_H
#define _INC_PRNG_GSL H

#include "prng.hpp"

#include <gsl/gsl_rng.h>

namespace mct {

	/** @brief A class for a pseudo-random
	 * number generator using the gsl prng functions.
	 * 
	  */
	 

	class PRNGenGSL : public PRNGen{
		public:
		
			PRNGenGSL(int seed);
		
			virtual ~PRNGenGSL();
			
			/*! \brief Generate a pseudo-random integer in [0,max).*/
			int uniform_int(unsigned int max) const;
			
			/*! \brief Generate a pseudo-random uniform double in [0,1).*/
			double uniform01_double() const;
			
			/*! \brief Generate a pseudo-random uniform double in [a, b].*/
			double uniform_double(double a, double b) const;
			
			/*! \brief Generate a pseudo-random exponential with rate a.*/
			double exponential_double(double a) const;
			
			/*! \brief Generate a pseudo-random poisson with rate a.*/
			virtual int poisson_int(double a) const;
						
			/*! \brief Generate a pseudo-random standard gaussian.*/
			double ugaussian_double() const;
			
			/*! \brief Generate a pseudo-random gaussian with mean m and sigma s.*/
			double gaussian_double(double m, double s) const;
			
			/*! \brief Get the maximum number this can generate.*/
			int getMax() const;
			
			/*! \brief Spawn another one like this, wrapped in a shared pointer.
			 
			  Uses this to generate the seed for the new one.*/
			virtual boost::shared_ptr < PRNGen > spawnAnother() const;
			
			/*! \brief Spawn another one like this, wrapped in a shared pointer.
			 
			 \param seed The seed for the new prng.*/
			virtual boost::shared_ptr < PRNGen > spawnAnother(int seed) const;
			
			private:
			
			gsl_rng* rr;
			
			PRNGenGSL();
			
			explicit PRNGenGSL(const PRNGenGSL& other);
			
			PRNGenGSL& operator=(PRNGenGSL tmp);
			
	};
	

} 

#endif 

