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
\brief MultiLociPolySitesSampler definitions.
*/

#include "multi_loci_polysites_sampler.hpp"
#include "mct_arg_factory.hpp"
#include "partitioned_arg_sequence_generator.hpp"
#include "libseq.hpp"
#include "population_labeled_arg.hpp"
#include "mct_sample_sequence.hpp"

#include <Sequence/Coalescent/SimTypes.hpp>
#include <Sequence/PolySites.hpp>


//#define CONVERT // no don't use this - just for experimentation
//#define INF_SITES_MODEL // no don't use this - just for experimentation

//#define MYDEBUG
//#define MYDEBUG_SEQ_PRINT
//#define MYDEBUG_INF // no don't use this - just for experimentation
#if defined (MYDEBUG) || defined (MYDEBUG_INF) || defined (MYDEBUG_SEQ_PRINT)
	#include <iostream>
	
#endif

using namespace mct;

MultiLociPolySitesSampler::MultiLociPolySitesSampler(
				const boost::shared_ptr < MultiLociModelSetSet >& _m_ptr,
				const boost::shared_ptr < PRNGen >& _r_ptr)
	: models_ptr(_m_ptr)
{
	initialisePRNGs(_r_ptr);
}

		
MultiLociPolySitesSampler::~MultiLociPolySitesSampler() 
{
		
}



boost::shared_ptr < MultiLociPolyTableSet >& 
		MultiLociPolySitesSampler::sample(
			boost::shared_ptr < MultiLociPolyTableSet >& mlsetptr,
			const size_t n_reps, 
			const boost::shared_ptr< const hudson_ms::PopulationStructure >& pop,
			const std::vector < size_t >& nsites,
			const boost::shared_ptr< const aabc::ParameterParticle>& p) const
	
{
	if (!(pop) || pop->empty()) { // pop ptr is empty or has no subpops
		throw std::invalid_argument("pop");
	}
	
	if (!pop->totalNsam()) { // no samples
		throw std::invalid_argument("pop contains no samples");
	}
	
	if (!(p)) { // p ptr is empty
		throw std::invalid_argument("p");
	}
	
	if (nsites.size() != models_ptr->size() ) {
		throw std::invalid_argument(
			"nsites.size() != loci in models");
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
	
	#ifdef MYDEBUG
		std::cout << std::endl;
	#endif
	
	return mlsetptr;
	
		
}


/* return pointer to new MultiLociPolyTable
 * If this method throws an exception 
 * the original set will not have changed.
 */
boost::shared_ptr < MultiLociPolyTable > 
		MultiLociPolySitesSampler::genNew(
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

	boost::shared_ptr < MultiLociPolyTable > mlptr(new MultiLociPolyTable(pop));

	//reset the thetas in the models
	// theta is in terms of theta per site already
	resetModelThetas(theta);
		
	/*
	 * use this prng to spawn another.
	 * This ensures that if we increase the number of loci from n to 
	 * n+1, the first n loci will be the same
	 * */
	
	boost::shared_ptr < PRNGen > r_loci_args 
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
				
		// now do the partitioned_arg
		PartitionedARGSequenceGenerator partitioned_arg(
											this_nsites, labeled_hist);
		
		#ifdef MYDEBUG
			std::cout << "j = " << j << " and arg is: " << std::endl;
			std::cout << partitioned_arg.newick_string() << std::endl;
		#endif
		
		/* extract sequences from the sample population
		 * using the ModelSet at j in the MultiLociModelSetSet 
		 * pointed to by models_ptr
		 * the total sites covered by the ModelSet (*(models_ptr->at(j)) 
		 * must be equal to this_nsites or we get an exception. */
		boost::shared_ptr <MCTSampleSequence > sample_seqs_ptr
					= partitioned_arg.createSampleSeqs(
					models_ptr->at(j));

		#if defined (MYDEBUG_INF) || defined (MYDEBUG) || defined (MYDEBUG_SEQ_PRINT)
			/* print the sequences of the sample in full form */
			std::cout << "sequence in full format is" << std::endl;
			std::cout << *sample_seqs_ptr << std::endl;
		#endif
		
		#if defined (MYDEBUG) || defined (MYDEBUG_SEQ_PRINT)
			/* print the sequences of the sample population in the phylip format */
			std::cout << "sequence in stringPhylipFormat() is" << std::endl;
			std::cout << sample_seqs_ptr->stringPhylipFormat();
		#endif
		
		
		boost::shared_ptr < Sequence::PolySites > new_ps
				= sample_seqs_ptr->polyNucTableFormatPtr();
		
		#if defined (MYDEBUG) || defined (MYDEBUG_SEQ_PRINT)
			std::cout << "polytable format is: " << std::endl;
			std::cout << LibSeq::polytableToString( *new_ps )<< std::endl;
		#endif	
		
		#ifdef CONVERT
			// recode the site data as binary 0,1
			
			#ifndef INF_SITES_MODEL
				LibSeq::startConversion (*new_ps,
					false, 0, false, // do not enforce strict inf sites model
					this_nsites);
			#endif
			#ifdef INF_SITES_MODEL
				LibSeq::startConversion (*new_ps,
					false, 0, true, // enforce strict inf sites model
					this_nsites);
			#endif
			
			
			#ifdef MYDEBUG_INF
				std::cout << "after start of conversion format is: " << std::endl;
				std::cout << LibSeq::polytableToString( *new_ps )<< std::endl;
			#endif
			
			boost::shared_ptr < Sequence::SimData > new_sd(new Sequence::SimData(
						new_ps->GetPositions(), new_ps->GetData() ));
			
			mlptr->add(new_sd, this_nsites);
			
		#endif
		
		#ifndef CONVERT				
			mlptr->add(new_ps, this_nsites);
		#endif
		
		#if defined (MYDEBUG) || defined (MYDEBUG_SEQ_PRINT)
			std::cout << "and seq is: " << std::endl;
			std::cout << LibSeq::polytableToString( *(mlptr->at(j)) )<< std::endl;
			std::cout << "nsites = " << this_nsites << std::endl;
		#endif		
	}
	
	return mlptr;	
}



void MultiLociPolySitesSampler::resetModelThetas(double new_theta) const
{
	models_ptr->resetThetas(new_theta);
	
}

void MultiLociPolySitesSampler::initialisePRNGs(
						const boost::shared_ptr < PRNGen >& _r_ptr)
{
	r_args_ptr = _r_ptr->spawnAnother();
	
}
	




