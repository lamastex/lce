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
\brief HudsonMSPRNGenGSL_mt declarations.
*/

#ifndef _INC_HUDSON_MS_PRNG_GSL_H
#define _INC_HUDSON_MS_PRNG_GSL H

#include "prng.hpp"

#include <gsl/gsl_rng.h>

namespace hudson_ms {

	/** @brief A class for a pseudo-random
	 * number generator for hudson_ms using the gsl prng functions.
	 * 
	  */
	 

	class HudsonMSPRNGenGSL_mt : public HudsonMSPRNGen {
		public:
		
			HudsonMSPRNGenGSL_mt(int seed);
		
			virtual ~HudsonMSPRNGenGSL_mt();
			
			/*! \brief Generate a pseudo-random double in [0,1).*/
			virtual double genRand() const;
			
			/*! \brief Generate a pseudo-random double in (0,1).*/
			virtual double genRandGtZero() const;
			
			
			private:
			
			gsl_rng* rr;
			
			explicit HudsonMSPRNGenGSL_mt();
			
			explicit HudsonMSPRNGenGSL_mt(const HudsonMSPRNGenGSL_mt& other);
			
			HudsonMSPRNGenGSL_mt& operator=(HudsonMSPRNGenGSL_mt tmp);
			
	};
	

} 

#endif 

