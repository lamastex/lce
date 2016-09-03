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
\brief MultiLociPolySitesSampler declarations.
*/

#ifndef _MULTI_LOCI_POLYSITES_SAMPLER_H
#define _MULTI_LOCI_POLYSITES_SAMPLER_H

#include "multi_loci_polytable_sampler.hpp"
#include "multi_loci_modelset_set.hpp"

#include "prng.hpp"


namespace mct {
	
		
	/** @brief A class capable of 
	 generating a collection of MultiLociPolyTable given a 
	 aabc::ParameterParticle.
	 
	 See the MultiLociPolyTableSampler interface.
	 	 
	 This class is used to generate 
	 \link MultiLociPolyTable MultiLociPolyTables\endlink
	 using simulations of polymorphic site data under 
	<b>finite sites models</b>.  
	 
	 
	 
	 \internal
	 At present the implementation uses mct models
	 to simulate of DNA sequence data.  The DNA sequence
	 data is then analysed into polymorphic site data and turned
	 into a PolyTable. The mct models are <b>finite sites models</b>. 
	 */

	
	class MultiLociPolySitesSampler : public MultiLociPolyTableSampler {
		public:
		
			/*! \brief Constructor. 
			
			\param  Pointer to a MultiLociModelSetSet, used to simulate
			the loci data. 
			\param _r_ptr A pointer to a random number generator.
			 */
			MultiLociPolySitesSampler(
				const boost::shared_ptr < MultiLociModelSetSet >& _m_ptr,
				const boost::shared_ptr < PRNGen >& _r_ptr);
		
		
			virtual ~MultiLociPolySitesSampler();
			
			
		
				
		protected:
		
			virtual boost::shared_ptr < MultiLociPolyTableSet >& 
				sample(
					boost::shared_ptr < MultiLociPolyTableSet >& 
							mlsetptr,
					const size_t n_reps, 
					const boost::shared_ptr< const hudson_ms::PopulationStructure >& pop,
					const std::vector < size_t >& nsites, 
					const boost::shared_ptr< const aabc::ParameterParticle>& p) const;

			virtual void initialisePRNGs(const boost::shared_ptr < PRNGen >& _r_ptr);
		
			virtual boost::shared_ptr < MultiLociPolyTable >
					genNew(
					const boost::shared_ptr< hudson_ms::PopulationStructure >& pop,
					const std::vector < size_t >& nsites,
					const boost::shared_ptr < const aabc::ParameterParticle >& p)
					const;
			
			virtual void resetModelThetas(double new_theta) const;
			
		
		private:
			
			boost::shared_ptr < MultiLociModelSetSet > models_ptr;
			
			boost::shared_ptr < PRNGen > r_args_ptr;
			
			
		
	};

} 

#endif 

