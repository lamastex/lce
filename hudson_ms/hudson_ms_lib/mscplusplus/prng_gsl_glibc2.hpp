/* mscplusplus - ms for c++.

   Copyright (C) 2011, 2012 Jenny Harlow
   
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
\brief HudsonMSPRNGenGLS_glibc2 declarations.
*/

#ifndef _INC_HUDSON_MS_PRNG_RAND_H
#define _INC_HUDSON_MS_PRNG_RAND H

#include <gsl/gsl_rng.h>

#include "prng.hpp"

namespace hudson_ms {

	/** @brief A class for a pseudo-random
	 * number generator for hudson_ms using the rand() function.
	 * 
	 * \note This class can set the seed for rand() but cannot 
	 * guarantee that rand() is not used between uses of this
	 * object, or that the seed has not been reset.
	 * 
	  */
	 

	class HudsonMSPRNGenGLS_glibc2 : public HudsonMSPRNGen {
		public:
		
			HudsonMSPRNGenGLS_glibc2(int seed);
		
			virtual ~HudsonMSPRNGenGLS_glibc2();
			
			/*! \brief Generate a pseudo-random double in [0,1).*/
			virtual double genRand() const;
			
			/*! \brief Generate a pseudo-random double in (0,1).*/
			virtual double genRandGtZero() const;
			
			
			private:
			
			gsl_rng* rr;
			
			explicit HudsonMSPRNGenGLS_glibc2();
			
			explicit HudsonMSPRNGenGLS_glibc2(const HudsonMSPRNGenGLS_glibc2& other);
			
			HudsonMSPRNGenGLS_glibc2& operator=(HudsonMSPRNGenGLS_glibc2 tmp);
			
	};
	

} 

#endif 

