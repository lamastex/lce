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
\brief NormalSDTestSampler definitions.
*/

#include "normalsd_test_sampler.hpp"
#include "parameter_particle.hpp"
#include "data_object.hpp"
#include "generic_data_object_set.hpp"

#include <boost/bind.hpp>
#include <boost/function.hpp>

#include <gsl/gsl_randist.h>


using namespace aabc;

NormalSDTestSampler::NormalSDTestSampler(
			const boost::shared_ptr < mct::PRNGen >& _r_ptr)
{
	initialisePRNGs(_r_ptr);
}

		
NormalSDTestSampler::~NormalSDTestSampler() 
{}

boost::shared_ptr < GenericDataObjectSet >& 
		NormalSDTestSampler::_sample(
			boost::shared_ptr < GenericDataObjectSet >& gsetptr,
			const size_t n_reps, 
			const boost::shared_ptr< const ParameterParticle>& p) const
{
	for (size_t i = 0; i < n_reps; ++i) {
		gsetptr->add( genNew(p) );
	
	}
	return gsetptr;
}	



boost::shared_ptr < NormalSDTest > NormalSDTestSampler::genNew(
					const boost::shared_ptr < const ParameterParticle >& p)
					const
{
	
	double value = rand->gaussian_double(0.0, p->at(0));
	
	return boost::shared_ptr < NormalSDTest > ( new NormalSDTest(value) );
	
}

void NormalSDTestSampler::initialisePRNGs(
				const boost::shared_ptr < mct::PRNGen >& _r_ptr)
{
	rand = _r_ptr->spawnAnother();
}
