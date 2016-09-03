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
\brief GenericDataObjectSampler definitions.
*/

#include "generic_data_object_sampler.hpp"
#include "data_object.hpp"

using namespace aabc;

	
GenericDataObjectSampler::~GenericDataObjectSampler() {}

boost::shared_ptr < GenericDataObjectSet >
		GenericDataObjectSampler::sample(
		const size_t n_reps,
		const boost::shared_ptr< const ParameterParticle> p) const
{
	boost::shared_ptr < GenericDataObjectSet > gsetptr 
			(new GenericDataObjectSet ());
	
	return _sample(
				gsetptr,
				n_reps,
				p);
}		

boost::shared_ptr < GenericDataObjectSet > 
		GenericDataObjectSampler::sample(
		const size_t n_reps_per_particle, 
		const boost::shared_ptr< ParameterParticleSet>& pset) const
{
	boost::shared_ptr <GenericDataObjectSet > gsetptr 
			(new GenericDataObjectSet ());

	//for each particle in the set, generate its data set and add to the dataset
	// could do this more efficiently with iterators on the set types
	
	for (size_t i = 0; i < pset->size(); ++i) {
		_sample( gsetptr,
				n_reps_per_particle,
				pset->at(i) );
	}
	
	return gsetptr;	
}
			
