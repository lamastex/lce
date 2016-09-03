/* MCT - Markov Chains on Trees.


   Copyright (C) 2011, 2012 Jennifer Harlow

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
\brief Tools for doing various simdata routines

*/

#include "testing_toolsSimdata.hpp"


#include "config.h" // mct

#include "multi_loci_simdata_sampler.hpp"
#include "parameter_particle.hpp"
#include "prng.hpp"
#include "prng_gsl.hpp"

#include <iostream>
#include <cassert>


using namespace std;
using namespace mct;

// make some reference statistics
boost::shared_ptr < mct::MultiLociPolyTableSet > makeSimDataRefSet(
		size_t n_obs,
		size_t nloci,
		size_t nsites,
		boost::shared_ptr< hudson_ms::PopulationStructure > pop,
		const std::vector < double >& trueParams,
		int trueseed)
{
	
	//vector of nloci copies of nsites
	std::vector < size_t > nloci_nsites(nloci, nsites);
		
	boost::shared_ptr < mct::PRNGen > rr_sampler ( new mct::PRNGenGSL(trueseed) );
	
	// make a simdata (infinite sites model) sampler with its prng
	boost::shared_ptr < MultiLociPolyTableSampler > 
					sampler_ptr ( new MultiLociSimDataSampler(rr_sampler) );
	
	
	boost::shared_ptr < aabc::ParameterParticle > 
			p( new aabc::ParameterParticle(trueParams) );

	cout << "sampling 'true' values " << ( p->toString() )<< endl;
	

	// get a set of values from the sampler, using the particles
	boost::shared_ptr < MultiLociPolyTableSet > refValues_ptr = 
			sampler_ptr->sample(n_obs, 
				pop,
				nloci_nsites,
				p);
	
	return refValues_ptr;
	
}


// make some reference statistics, specifying required stats
boost::shared_ptr < mct::SummaryStatisticSet > makeAHABCRefSet(
		size_t n_obs,
		size_t nloci,
		size_t nsites,
		boost::shared_ptr< hudson_ms::PopulationStructure > pop,
		const std::vector < double >& trueParams,
		const std::vector < size_t >& reqStats,
		int trueseed)
{
	
	boost::shared_ptr < MultiLociPolyTableSet > refValues_ptr
							= makeSimDataRefSet(
										n_obs,
										nloci,
										nsites,
										pop,
										trueParams,
										trueseed);

	
	// and get a set of reference summary statistics from the values
	boost::shared_ptr < mct::SummaryStatisticSet > refset_ptr_tmp =
					refValues_ptr->getSummaryStatisticSet();
	boost::shared_ptr < mct::SummaryStatisticSet > refset_ptr;
	if(!reqStats.empty()) {

		refset_ptr =
				refset_ptr_tmp->makeSelectedSummaryStatisticSet(reqStats);
	}
	else {

		refset_ptr = refset_ptr_tmp;
	}
	
	return refset_ptr;
}
