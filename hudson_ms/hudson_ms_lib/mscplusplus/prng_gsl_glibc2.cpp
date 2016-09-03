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
\brief HudsonMSPRNGenGLS_glibc2 definitions.
*/
#include "prng_gsl_glibc2.hpp"

#include <cstdlib>


using namespace hudson_ms;

HudsonMSPRNGenGLS_glibc2::HudsonMSPRNGenGLS_glibc2(int seed)
{
	try {
		rr = NULL;
		rr = gsl_rng_alloc(gsl_rng_random_glibc2);
		
		gsl_rng_set(rr, seed);
	}
	catch (...) {
		try {
			if (rr != NULL) gsl_rng_free(rr);
		}
		catch (...) {}
		throw;
	}
	
}
		
HudsonMSPRNGenGLS_glibc2::~HudsonMSPRNGenGLS_glibc2()
{
	try {
		gsl_rng_free(rr);
	}
	catch (...) {} // catch and swallow
}

// Generate a pseudo-random integer in [0,1).
double HudsonMSPRNGenGLS_glibc2::genRand() const
{
	return gsl_rng_uniform(rr);
}

// Generate a pseudo-random uniform double in (0,1).
double HudsonMSPRNGenGLS_glibc2::genRandGtZero() const
{
	return gsl_rng_uniform_pos(rr);
}

