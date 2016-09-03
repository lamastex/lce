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
\brief MixtureTest2DSampler definitions.
*/


#include "mixture_test_2d_sampler.hpp"
#include "parameter_particle.hpp"
#include "data_object.hpp"
#include "generic_data_object_set.hpp"

#include "rmv.h" // multivariate normal density and random values

#include <boost/bind.hpp>
#include <boost/function.hpp>

using namespace aabc;

MixtureTest2DSampler::MixtureTest2DSampler(const boost::shared_ptr < mct::PRNGen >& _r_ptr)
{
	initialisePRNGs(_r_ptr);
}



		
MixtureTest2DSampler::~MixtureTest2DSampler() 
{}

boost::shared_ptr < GenericDataObjectSet >& 
		MixtureTest2DSampler::_sample(
			boost::shared_ptr < GenericDataObjectSet >& gsetptr,
			const size_t n_reps, 
			const boost::shared_ptr< const ParameterParticle>& p) const
{
	for (size_t i = 0; i < n_reps; ++i) {
		gsetptr->add( genNew(p) );
		
	}
	//return by ref
	return gsetptr;
	
}

// this is just a dummy at present
boost::shared_ptr < MixtureTest2D > MixtureTest2DSampler::genNew(
					const boost::shared_ptr < const ParameterParticle >& p)
					const
{
	gsl_vector* result = NULL;
	gsl_vector* mean = NULL;
	gsl_matrix* cov = NULL;
	
	try {
		boost::function<double ()> ugaussian_func =
		boost::bind<double>(&mct::PRNGen::ugaussian_double, rand);
	    //boost::bind<double>(gsl_ran_ugaussian, rand);
			
		result = gsl_vector_calloc(2);
		
		mean = p->get_gsl_vector();
		//mean is (p->at(0), p->at(1))
		
		if ( rand->uniform01_double() < 0.5 ) {
			
			cov = gsl_matrix_calloc(2,2);
			gsl_matrix_set(cov, 0, 0, 1.0);
			gsl_matrix_set(cov, 0, 1, 0.0);
			gsl_matrix_set(cov, 1, 0, 0.0);
			gsl_matrix_set(cov, 1, 1, 1.0);
			
			rmvnorm(ugaussian_func, mean, cov, result); // fill in result
			gsl_matrix_free(cov);
		}
		
		else {
			gsl_matrix* cov = gsl_matrix_calloc(2,2);
			gsl_matrix_set(cov, 0, 0, 0.01);
			gsl_matrix_set(cov, 0, 1, 0.0);
			gsl_matrix_set(cov, 1, 0, 0.0);
			gsl_matrix_set(cov, 1, 1, 0.01);
			
			rmvnorm(ugaussian_func, mean, cov, result); // fill in result
			gsl_matrix_free(cov);
			cov = NULL;
		}
		
		double v0 = gsl_vector_get(result, 0);
		double v1 = gsl_vector_get(result, 1);
		
		gsl_vector_free(mean);
		mean = NULL;
		gsl_vector_free(result);
		result = NULL;
		
		return boost::shared_ptr < MixtureTest2D > ( new MixtureTest2D(v0, v1) );
	}
	catch (...) {
		try {
			if (mean != NULL) gsl_vector_free(mean);
		}
		catch (...) {}// catch and swallow
		try {
			if (result != NULL) gsl_vector_free(result);
		}
		catch (...) {}// catch and swallow
		try {
			if (cov != NULL) gsl_matrix_free(cov);
		}
		catch (...) {}// catch and swallow
		
		throw;
	}
		
}

void MixtureTest2DSampler::initialisePRNGs(
				const boost::shared_ptr < mct::PRNGen >& _r_ptr)
{
	rand = _r_ptr->spawnAnother();
	/*
	rand = gsl_rng_alloc(gsl_rng_mt19937);

	//use the current random number generator to get a seed
	// seed will be between 1 and gsl_rng_max(r) inclusive
	size_t seed = 1 + 
				gsl_rng_uniform_int(r, gsl_rng_max(rand));		
	gsl_rng_set(rand, seed);
	*/
}
