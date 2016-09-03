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
\brief MultiLociPolyTableSampler declarations.
*/

#ifndef _MULTI_LOCI_POLYTABLE_SAMPLER_H
#define _MULTI_LOCI_POLYTABLE_SAMPLER_H

#include "multi_loci_polytable_set.hpp"
#include "multi_loci_polytable.hpp"
#include "parameter_particle_set.hpp"
#include "parameter_particle.hpp"
#include <mscplusplus/population_structure.hpp>

#include <boost/smart_ptr.hpp>


namespace mct {
	
	/** @brief An abstract type capable of 
	 generating a collection of MultiLociPolyTable given a 
	 aabc::ParameterParticle or a aabc::ParameterParticleSet.
	 
	 A MultiLociPolyTableSampler should be able to produce a 
	 MultiLociPolyTableSet given the number of replications 
	 required and the required configuration of each.  
	  
	 @todo the samplers should be configured with a map telling them where
	 to look for each parameter they need in the parameter particles.
	 */

	class MultiLociPolyTableSampler {
		public:
		
			virtual ~MultiLociPolyTableSampler() {};
			
			/*! \brief  Generate a collection of the data based on an
			 example.
			  
			 The example is used to get the configuration needed 
			 (e.g. number of samples, loci, sites, etc).
			
			\param n_reps The number of replications to generate.
			\param example the MultiLociPolyTable to get configurations
			from.
			\param p The aabc::ParameterParticle to use to generate the
			simulations.
			\return A pointer to the MultiLociPolyTableSet 
			holding the simulations.			*/
			virtual boost::shared_ptr < MultiLociPolyTableSet > 
				sample(const size_t n_reps, 
				const MultiLociPolyTable& example,
				const boost::shared_ptr< const aabc::ParameterParticle>& p) const;
			
			/*! \brief  Generate a collection of the data.
			 
			 The configuration of the required replications is given 
			 with individual parameter arguments.
			
			\param n_reps The number of replications to generate.
			\param pop The 
			\link hudson_ms::PopulationStructure PopulationStructure\endlink
			for each replicant.
			\param nsites An ordered collection giving the number
			of sites in each locus (the size of \a nsites is the number of loci).
			\param p The aabc::ParameterParticle to use to generate the
			simulations.
			\return A pointer to the MultiLociPolyTableSet 
			holding the simulations.
			\pre \a pop is not an empty pointer.
			\pre \a pop contains at least one subpopulation (one
			subpopulation is taken as the whole population).*/
			virtual boost::shared_ptr < MultiLociPolyTableSet > 
				sample(const size_t n_reps, 
				const boost::shared_ptr< const hudson_ms::PopulationStructure >& pop,
				const std::vector < size_t >& nsites,
				const boost::shared_ptr< const aabc::ParameterParticle>& p) const;
		
			/*! \brief  Generate a collection of the data based on an
			 example.
			  
			 The example is used to get the configuration needed 
			 (e.g. number of samples, loci, sites, etc).
			
			\param n_reps_per_particle The number of replications to 
			generate for each particle in the aabc::ParameterParticleSet.
			\param example the MultiLociPolyTable to get configurations
			from.
			\param pset The aabc::ParameterParticleSet to use to generate the
			simulations.
			\return A pointer to the MultiLociPolyTableSet 
			holding the simulations.			*/
			virtual boost::shared_ptr < MultiLociPolyTableSet > 
				sample(const size_t n_reps_per_particle, 
				const MultiLociPolyTable& example,
				const boost::shared_ptr< aabc::ParameterParticleSet>& pset) 
															const;
			
			/*! \brief  Generate a collection of the data.
			 
			 The configuration of the required replications is given 
			 with individual parameter arguments.
			
			\param n_reps_per_particle The number of replications to 
			generate for each particle in the aabc::ParameterParticleSet.
			\param pop The 
			\link hudson_ms::PopulationStructure PopulationStructure\endlink
			for each replicant.
			\param nsites An ordered collection giving the number
			of sites in each locus (the size of \a nsites is the number of loci).
			\param pset The aabc::ParameterParticleSet to use to generate the
			simulations.
			\return A pointer to the MultiLociPolyTableSet 
			holding the simulations.
			\pre \a pop is not an empty pointer.
			\pre \a pop contains at least one subpopulation (one
			subpopulation is taken as the whole population).
			\pre The size of \a nsites is equal to \a nloci.*/
			virtual boost::shared_ptr < MultiLociPolyTableSet > 
				sample(const size_t n_reps_per_particle, 
				const boost::shared_ptr< const hudson_ms::PopulationStructure >& pop,
				const std::vector < size_t >& nsites,
				const boost::shared_ptr< aabc::ParameterParticleSet>& pset) 
															const;
		
		protected :
			virtual boost::shared_ptr < MultiLociPolyTableSet >& 
				sample(
					boost::shared_ptr < MultiLociPolyTableSet >& 
							mlsetptr,
					const size_t n_reps, 
					const boost::shared_ptr< const hudson_ms::PopulationStructure >& pop,
					const std::vector < size_t >& nsites, 
					const boost::shared_ptr< const aabc::ParameterParticle>& p) 
															const = 0;
	};

} 

#endif 

