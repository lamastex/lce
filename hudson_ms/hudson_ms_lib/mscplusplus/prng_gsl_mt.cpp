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
\brief HudsonMSPRNGenGSL_mt definitions.
*/
#include "prng_gsl_mt.hpp"

#include <gsl/gsl_randist.h>

using namespace hudson_ms;

HudsonMSPRNGenGSL_mt::HudsonMSPRNGenGSL_mt(int seed)
{

	try {
		rr = NULL;
		rr = gsl_rng_alloc(gsl_rng_mt19937);
		
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
		
HudsonMSPRNGenGSL_mt::~HudsonMSPRNGenGSL_mt()
{
	try {
		gsl_rng_free(rr);
	}
	catch (...) {} // catch and swallow
}

// Generate a pseudo-random integer in [0,1).
double HudsonMSPRNGenGSL_mt::genRand() const
{
	return gsl_rng_uniform(rr);
}

// Generate a pseudo-random uniform double in (0,1).
double HudsonMSPRNGenGSL_mt::genRandGtZero() const
{
	return gsl_rng_uniform_pos(rr);
}

