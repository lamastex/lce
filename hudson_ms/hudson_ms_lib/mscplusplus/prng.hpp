/* mscplusplus - ms for c++.

   Copyright (C) 2011, 2012 Jenny Harlow
   
   his program is free software; you can redistribute it and/or modify
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
\brief HudsonMSPRNGen declarations.
*/

#ifndef _INC_HUDSON_MS_PRNG_H
#define _INC_HUDSON_MS_PRNG_H


namespace hudson_ms {

	/** @brief An abstract class for a pseudo-random
	 * number generator for hudson_ms work.
	 * 
	 * Base type for concrete sub-types which act as wrappers
	 * for particular prngs.
	  */
	 

	class HudsonMSPRNGen {
		public:
		
			virtual ~HudsonMSPRNGen() {};
			
			/*! \brief Generate a pseudo-random double in [0,1).*/
			virtual double genRand() const = 0;
			
			/*! \brief Generate a pseudo-random double in (0,1).*/
			virtual double genRandGtZero() const = 0;
	};
	

} 

#endif 

