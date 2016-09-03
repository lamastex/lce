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
\brief MultiLociSimDataSampler declarations.
*/

#ifndef _MULTI_LOCI_SIMDATA_SAMPLER_H
#define _MULTI_LOCI_SIMDATA_SAMPLER_H

#include "multi_loci_polytable_sampler.hpp"

#include "prng.hpp"



namespace mct {
	
	/** @brief A class capable of 
	 generating a collection of MultiLociPolyTable given a 
	 aabc::ParameterParticle.
	  
	 See the MultiLociPolyTableSampler interface.
	 
	 This class is used to generate 
	 \link MultiLociPolyTable MultiLociPolyTables\endlink
	 using simulations of polymorphic site data under an
	<b>infinite sites model</b>.  
	 * 
	
	\internal
	This implementation uses Libsequence to simulate snp data (simdata objects)
	but my mscplusplus replication of ms to make the args.  The main
	advantage of this is that ms can deal with
	populations of subpopulations. In the one-population case, 
	using the ms route means
	that args will be the same as if using the polysites sampler
	(if same arg seed is used) and the only difference in the output
	is the differences in the ways the polymorphic site data data is simulated. 

	
	 \internal
	 At present the implementation uses Libsequence to directly
	 simulate the snp data given the arg.  Libsequence uses 
	 an <b>infinite sites model</b>. Now that mscplusplus
	 can also simulate segregated sites, we could use this 
	 instead of Libsequence (and then use whatever method
	 Libsequence offers to turn mscplusplus results into
	 a Sequence::SimData object to use here.  I am not sure 
	 if there would be any advantage in the mscplusplus route
	 but it is an option for later.	 */

	
	class MultiLociSimDataSampler  : public MultiLociPolyTableSampler {
		public:
		
			/*! \brief Constructor. 
			 
			\param _r_ptr A pointer to a random number generator.
			 */
			explicit MultiLociSimDataSampler(
						const boost::shared_ptr < PRNGen >& _r_ptr);
		
		
			virtual ~MultiLociSimDataSampler();
			
				
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
		
		private:
		
			boost::shared_ptr < PRNGen > r_args_ptr;
			boost::shared_ptr < PRNGen > r_sims_ptr;
			
	};

} 

#endif 

