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
\brief MultiLociSimDataSampler definitions.
*/

#include "multi_loci_simdata_sampler.hpp"
#include "parameter_particle.hpp"
#include "multi_loci_polytable_set.hpp"
#include "multi_loci_polytable.hpp"
#include "population_labeled_arg.hpp"
#include "libseq.hpp"
#include "mct_arg_factory.hpp"

#include <Sequence/Coalescent/SimTypes.hpp>

//#define MYDEBUG
#if defined (MYDEBUG)
	#include <iostream>
	
#endif

using namespace mct;

MultiLociSimDataSampler::MultiLociSimDataSampler(
					const boost::shared_ptr < PRNGen >& _r_ptr)
{
	initialisePRNGs(_r_ptr);
}


		
MultiLociSimDataSampler::~MultiLociSimDataSampler() 
{
		
}

//protected 

boost::shared_ptr < MultiLociPolyTableSet >& 
			MultiLociSimDataSampler::sample(
				boost::shared_ptr < MultiLociPolyTableSet >& 
						mlsetptr,
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
	
	// copy the population structure
	boost::shared_ptr< hudson_ms::PopulationStructure > _pop(
									new hudson_ms::PopulationStructure( *pop) );
	
	for (size_t i = 0; i < n_reps; ++i) {
		
		#ifdef MYDEBUG
			std::cout << "\n*********Rep number********** " << i << std::endl;
		#endif
		
		mlsetptr->add( genNew(_pop, nsites, p) );
		#ifdef MYDEBUG
			std::cout << "rep " << i << " added" << std::endl;
		#endif
		
	}
	
	#ifdef MYDEBUG
		std::cout << std::endl;
	#endif
	
	return mlsetptr;
	
}

/* return pointer to new MultiLociPolyTable
 * If this method throws an exception 
 * the original set will not have changed.
 * 
 * I tried replacing this with a method that took the set by reference
 * and added directly to the set, to save passing pointers around and
 * it made no difference to the timings.
 */
boost::shared_ptr < MultiLociPolyTable > 
		MultiLociSimDataSampler::genNew(
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
	
	/*
	 * use our rand to seed 2 more random number generators.
	 * This ensures that if we increase the number of loci from n to 
	 * n+1, the first n loci will be the same, and if we have
	 * have different methods of creating args, the args will be 
	 * the same.  Do the args one first so that it makes no
	 * difference if we don't want one for sims at all.
	 * */
	boost::shared_ptr < PRNGen > r_loci_args 
				= r_args_ptr->spawnAnother();
				
	boost::shared_ptr < PRNGen > r_loci_sims 
				= r_sims_ptr->spawnAnother();
	
	size_t nloci = nsites.size();
	
	for (size_t j = 0; j < nloci; ++j) {
		
		// Use a different ARG for each locus
		// ie assume infinite recombination between loci
		double rho = 0.0;
		
		// nsites for this locus taken from the vector
		size_t this_nsites = nsites.at(j);
		// theta is converted to the per locus rate libseq wants
		double this_theta = theta * this_nsites;  
		
		// ARG generation depending on the population structure
		boost::shared_ptr < PopulationLabeledARG > labeled_hist;
		
							
		labeled_hist = 
			ARGFactory::makeHistFromMS(pop, this_nsites, rho, r_loci_args);
			
							
		#ifdef MYDEBUG
			std::cout << "j = " << j << " and arg is: " << std::endl;
			std::cout << LibSeq::argNewickString(*(labeled_hist->getARG())) << std::endl;
		#endif
		
		Sequence::SimData tmp = 
						LibSeq::infinite_sites_sim_data(
								this_nsites,
								r_loci_sims,
								*( labeled_hist->getARG() ),
								this_theta);
		
		#ifdef MYDEBUG
			std::cout << " sizeof (tmp) = " << sizeof(tmp) << std::endl;
		#endif
		
		boost::shared_ptr < Sequence::SimData > 
				new_sd(new Sequence::SimData( tmp ) );
		
		mlptr->add(new_sd, this_nsites) ;
		
		#ifdef MYDEBUG								
			std::cout << "and seq is: " << std::endl;
			std::cout << LibSeq::polytableToString( *(mlptr->at(j)) )<< std::endl;
			std::cout << "nsites = " << this_nsites << std::endl;
		#endif
	}
				

	return mlptr;	

}

void MultiLociSimDataSampler::initialisePRNGs(
						const boost::shared_ptr < PRNGen >& _r_ptr)
{
	r_args_ptr = _r_ptr->spawnAnother();
	r_sims_ptr = _r_ptr->spawnAnother();
	
}
		
