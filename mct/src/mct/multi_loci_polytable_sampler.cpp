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
\brief MultiLociPolyTableSampler definitions.
*/

#include "multi_loci_polytable_sampler.hpp"
#include "parameter_particle.hpp"
#include "parameter_particle_set.hpp"
#include "multi_loci_polytable_set.hpp"
#include "multi_loci_polytable.hpp"

#include <stdexcept>

//#define MYDEBUG_COUNT_REPS
#if defined (MYDEBUG_COUNT_REPS)
	#include <iostream>
	
#endif

using namespace mct;

boost::shared_ptr < MultiLociPolyTableSet > 
		MultiLociPolyTableSampler::sample(const size_t n_reps, 
				const MultiLociPolyTable& example,
				const boost::shared_ptr< const aabc::ParameterParticle>& p) const
{
	return sample(n_reps,
					example.populationStructure(),
					example.nsites(), 
					p);
}

boost::shared_ptr < MultiLociPolyTableSet > 
		MultiLociPolyTableSampler::sample(const size_t n_reps, 
				const boost::shared_ptr< const hudson_ms::PopulationStructure >& pop,
				const std::vector < size_t >& nsites,
				const boost::shared_ptr< const aabc::ParameterParticle>& p) const
{
	try {	
		boost::shared_ptr < MultiLociPolyTableSet > mlsetptr 
					(new MultiLociPolyTableSet ());
		
		
		
		return sample(
					mlsetptr,
					n_reps,
					pop,
					nsites, 
					p);
		
	}
	catch (std::invalid_argument& ia) {
		
		std::string msg = std::string("MultiLociPolyTableSampler::sample(")
		+ std::string("const size_t, ")
		+ std::string("const boost::shared_ptr< hudson_ms::PopulationStructure >&, ")
		+ std::string("const std::vector < size_t >&, ")
		+ std::string("const boost::shared_ptr< aabc::ParameterParticle>&) :\n");
	
		throw (std::invalid_argument(msg + ia.what()));
	}	
}


boost::shared_ptr < MultiLociPolyTableSet > 
	MultiLociPolyTableSampler::sample(
				const size_t n_reps_per_particle, 
				const MultiLociPolyTable& example,
				const boost::shared_ptr< aabc::ParameterParticleSet>& pset) 
															const
{
	return sample(n_reps_per_particle,
					example.populationStructure(),
					example.nsites(), 
					pset);
}



boost::shared_ptr < MultiLociPolyTableSet > 
	MultiLociPolyTableSampler::sample(
				const size_t n_reps_per_particle,
				const boost::shared_ptr< const hudson_ms::PopulationStructure >& pop,
				const std::vector < size_t >& nsites,
				const boost::shared_ptr< aabc::ParameterParticleSet>& pset) const
{
	try {	
		
		if (!(pset)) { // p ptr is empty
			throw std::invalid_argument(": pset");
		}
		
		boost::shared_ptr < MultiLociPolyTableSet > mlsetptr 
				(new MultiLociPolyTableSet ());

		//for each particle in the set, generate its data set and add to the dataset
		// could do this more efficiently with iterators on the set types
		
		for (size_t i = 0; i < pset->size(); ++i) {
			
			#ifdef MYDEBUG_COUNT_REPS
				if ((i+1)%10000 == 0) 
					std::cout << "\nRep number " << (i+1) << std::endl;
			#endif
	
			
			sample( mlsetptr,
					n_reps_per_particle,
					pop,
					nsites,
					pset->at(i) );
		}
		
		return mlsetptr;	
	}
	catch (std::invalid_argument& ia) {
		
		std::string msg = std::string("MultiLociPolyTableSampler::sample(")
		+ std::string("const size_t, ")
		+ std::string("const boost::shared_ptr< hudson_ms::PopulationStructure >&, ")
		+ std::string("const std::vector < size_t >&, ")
		+ std::string("const boost::shared_ptr< aabc::ParameterParticleSet>&) :\n");
	
		throw (std::invalid_argument(msg + ia.what()));
	}

	
}
