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
\brief MultiLociMicrosatSampler definitions.
*/

#include "multi_loci_microsat_sampler.hpp"
#include "mct_arg_factory.hpp"
#include "partitioned_arg_sequence_generator.hpp"
#include "libseq.hpp"
#include "population_labeled_arg.hpp"
#include "mct_sample_sequence.hpp"

#include <Sequence/Coalescent/SimTypes.hpp>
#include <Sequence/PolySites.hpp>

#include <stdexcept>

//#define MYDEBUG
//#define MYDEBUGTMP
//#define XTRAMYDEBUG

#if defined (MYDEBUG) || defined (XTRAMYDEBUG) || defined (MYDEBUGTMP)
	#include <iostream>
	
#endif

using namespace mct;

MultiLociMicrosatSampler::MultiLociMicrosatSampler(
				const boost::shared_ptr < MultiLociModelSetSet >& _m_ptr,
				const boost::shared_ptr < PRNGen >& _r_ptr)
	: models_ptr(_m_ptr)
{
	initialisePRNGs(_r_ptr);
}

		
MultiLociMicrosatSampler::~MultiLociMicrosatSampler() 
{
		
}


boost::shared_ptr < MultiLociMicrosatSet > 
		MultiLociMicrosatSampler::sample(const size_t n_reps, 
				const MultiLociMicrosat& example,
				const boost::shared_ptr< const aabc::ParameterParticle>& p) const
{
	return sample(n_reps,
					example.populationStructure(),
					example.nsites(), 
					p);
}


boost::shared_ptr < MultiLociMicrosatSet > 
		MultiLociMicrosatSampler::sample(const size_t n_reps, 
				const boost::shared_ptr< const hudson_ms::PopulationStructure >& pop,
				const std::vector < size_t >& nsites, 
				const boost::shared_ptr< const aabc::ParameterParticle>& p) const
{
	try {
		boost::shared_ptr < MultiLociMicrosatSet > mlsetptr 
				(new MultiLociMicrosatSet ());
		
		return sample(
					mlsetptr,
					n_reps,
					pop,
					nsites, 
					p);
	}
	catch (std::invalid_argument& ia) {
		
		std::string msg = std::string("MultiLociMicrosatSampler::sample(")
		+ std::string("const size_t, ")
		+ std::string("const boost::shared_ptr< hudson_ms::PopulationStructure >&, ")
		+ std::string("const std::vector < size_t >&, ")
		+ std::string("const boost::shared_ptr< aabc::ParameterParticle>&)");
	
		throw (std::invalid_argument(msg + ia.what()));
	}	
}


boost::shared_ptr < MultiLociMicrosatSet > 
	MultiLociMicrosatSampler::sample(
				const size_t n_reps_per_particle, 
				const MultiLociMicrosat& example,
				const boost::shared_ptr< aabc::ParameterParticleSet>& pset) 
															const
{
	return sample(n_reps_per_particle,
					example.populationStructure(),
					example.nsites(), 
					pset);
}



boost::shared_ptr < MultiLociMicrosatSet > 
	MultiLociMicrosatSampler::sample(
				const size_t n_reps_per_particle,
				const boost::shared_ptr< const hudson_ms::PopulationStructure >& pop,
				const std::vector < size_t >& nsites, 
				const boost::shared_ptr< aabc::ParameterParticleSet>& pset) const
{
	try {
		
		if (!(pset)) { // p ptr is empty
			throw std::invalid_argument(": pset");
		}	
		
		boost::shared_ptr < MultiLociMicrosatSet > mlsetptr 
				(new MultiLociMicrosatSet ());

		//for each particle in the set, generate its data set and add to the dataset
		// could do this more efficiently with iterators on the set types
		
		for (size_t i = 0; i < pset->size(); ++i) {
			sample( mlsetptr,
					n_reps_per_particle,
					pop,
					nsites,
					pset->at(i) );
		}
	
		return mlsetptr;
	}
	catch (std::invalid_argument& ia) {
		
		std::string msg = std::string("MultiLociMicrosatSampler::sample(")
		+ std::string("const size_t, ")
		+ std::string("const boost::shared_ptr< hudson_ms::PopulationStructure >&, ")
		+ std::string("const std::vector < size_t >&, ")
		+ std::string("const boost::shared_ptr< aabc::ParameterParticleSet>&)");
	
		throw (std::invalid_argument(msg + ia.what()));
	}
	
}



boost::shared_ptr < MultiLociMicrosatSet >& 
		MultiLociMicrosatSampler::sample(
			boost::shared_ptr < MultiLociMicrosatSet >& mlsetptr,
			const size_t n_reps, 
			const boost::shared_ptr< const hudson_ms::PopulationStructure >& pop,
			const std::vector < size_t >& nsites, 
			const boost::shared_ptr< const aabc::ParameterParticle>& p) const
	
{
	if (!(pop) || pop->empty()) { // pop ptr is empty or has no subpops
		throw std::invalid_argument(": pop");
	}
	
	if (!(p)) { // p ptr is empty
		throw std::invalid_argument(": p");
	}
	
	if (nsites.size() != models_ptr->size() ) {
		throw std::invalid_argument(
			": nsites.size() != loci in models");
	}
	
	// copy the population structure
	boost::shared_ptr< hudson_ms::PopulationStructure > _pop(
									new hudson_ms::PopulationStructure( *pop) );
	
	for (size_t i = 0; i < n_reps; ++i) {
		
		#ifdef MYDEBUG
			std::cout << "\nRep number " << i << std::endl;
		#endif
		
		mlsetptr->add( genNew(_pop, nsites, p) );
		
	}
	return mlsetptr;
	
	#ifdef MYDEBUG
		std::cout << std::endl;
	#endif
		
}



/* return pointer to new MultiLociMicrosat
 * If this method throws an exception 
 * the original set will not have changed.
 */
boost::shared_ptr < MultiLociMicrosat > 
		MultiLociMicrosatSampler::genNew(
					const boost::shared_ptr< hudson_ms::PopulationStructure >& pop,
					const std::vector < size_t >& nsites, 
					const boost::shared_ptr < const aabc::ParameterParticle >& p)
					const
{
	// get parameters from the particle
	double theta = p->at(0);
	if (p->size() > 1) {	
		double growth = p->at(1);
	
		//growth is fixed for all sub-pops here:  we may want to be more flexible than this
		//eg use the pop's addToEventList method
		pop->setAllGrowthRates(growth);
	}
	
	boost::shared_ptr < MultiLociMicrosat > mlptr(
						new MultiLociMicrosat(pop));
	
	//reset the thetas in the models
	// theta is in terms of theta per site already
	resetModelThetas(theta);
		
	/*
	 * use this prng to make another.
	 * This ensures that if we increase the number of loci from n to 
	 * n+1, the first n loci will be the same
	 * */
	
	boost::shared_ptr < const PRNGen > r_loci_args 
						= r_args_ptr->spawnAnother();
	
	size_t nloci = nsites.size();

	for (size_t j = 0; j < nloci; ++j) {
		
		// Use a different ARG for each locus
		// ie assume infinite recombination between loci
		double rho = 0.0;
		
		// nsites for this locus taken from the vector
		size_t this_nsites = nsites.at(j);
		
		boost::shared_ptr < PopulationLabeledARG > labeled_hist = 
			ARGFactory::makeHistFromMS(pop, this_nsites, rho, 
						r_loci_args);
				
		PartitionedARGSequenceGenerator partitioned_arg(this_nsites, labeled_hist);
		
		#if  defined (XTRAMYDEBUG)  || defined (MYDEBUG)
			std::cout << "j = " << j << " and arg is: " << std::endl;
			std::cout << partitioned_arg.newick_string() << std::endl;
		#endif
		
		/* extract sequences from the sample population
		 * using the ModelSet at j in the MultiLociModelSetSet 
		 * pointed to by models_ptr
		 * the total sites covered by the ModelSet (*(models_ptr->at(j)) 
		 * must be equal to this_nsites or we get an exception. */
		boost::shared_ptr <MCTSampleSequence > sample_seqs_ptr
					= partitioned_arg.createSampleSeqs(models_ptr->at(j));

		#ifdef MYDEBUG
			/* print the sequences of the sample in full form */
			std::cout << "sequence in full format is" << std::endl;
			std::cout << *sample_seqs_ptr << std::endl;
			std::cout << "nsites = " << this_nsites << std::endl;
		#endif
		
		mlptr->add( sample_seqs_ptr->microsatRepeatsFormat() ) ;
		
	}
	
	return mlptr;	
}



void MultiLociMicrosatSampler::resetModelThetas(double new_theta) const
{
	models_ptr->resetThetas(new_theta);
}

void MultiLociMicrosatSampler::initialisePRNGs(
					const boost::shared_ptr < PRNGen > _r_ptr)
{
	r_args_ptr = _r_ptr->spawnAnother();
}





