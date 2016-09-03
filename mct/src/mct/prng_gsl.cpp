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
\brief PRNGenGSL definitions.
*/
#include "prng_gsl.hpp"

#include <gsl/gsl_randist.h>

using namespace mct;

PRNGenGSL::PRNGenGSL(int seed)
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
		
PRNGenGSL::~PRNGenGSL()
{
	try {
		gsl_rng_free(rr);
	}
	catch (...) {} // catch and swallow
}

// Generate a pseudo-random integer in [0,max).
int PRNGenGSL::uniform_int(unsigned int max) const
{
	return gsl_rng_uniform_int(rr, max);
}

// Generate a pseudo-random uniform double in [0,1).
double PRNGenGSL::uniform01_double() const
{
	return gsl_rng_uniform(rr);
}

// Generate a pseudo-random uniform double in [a, b].
double PRNGenGSL::uniform_double(double a, double b) const
{
	return gsl_ran_flat(rr, a, b);
}

// Generate a pseudo-random exponential with rate a.
double PRNGenGSL::exponential_double(double a) const
{
	return gsl_ran_exponential(rr, a);
}

// Generate a pseudo-random poisson with rate a.
int PRNGenGSL::poisson_int(double a) const
{
	return gsl_ran_poisson(rr, a);
}

// Generate a pseudo-random standard gaussian.
double PRNGenGSL::ugaussian_double() const
{
	return gsl_ran_ugaussian(rr);
}

// Generate a pseudo-random gaussian with mean m and sigma s.
double PRNGenGSL::gaussian_double(double m, double s) const
{
	return m + gsl_ran_gaussian(rr, s);
}

// what is the max
int PRNGenGSL::getMax() const
{
	return gsl_rng_max(rr);
}

// Spawn another one like this, wrapped in a shared pointer.
boost::shared_ptr < PRNGen > PRNGenGSL::spawnAnother() const
{
	int seed = 1 + uniform_int( getMax() );	
			
	boost::shared_ptr < PRNGen > ret ( new PRNGenGSL(seed) );
	
	return ret;
}

// Spawn another one like this, wrapped in a shared pointer.
boost::shared_ptr < PRNGen > PRNGenGSL::spawnAnother(int seed) const
{
	boost::shared_ptr < PRNGen > ret ( new PRNGenGSL(seed) );
	
	return ret;
}
