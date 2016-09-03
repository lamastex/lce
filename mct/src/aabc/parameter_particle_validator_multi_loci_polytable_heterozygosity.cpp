/* MCT - Markov Chains on Trees.

   Copyright (C) 2012 Jenny Harlow

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
\brief ParameterParticleValidatorMultiLociPolyTableHeterozygosity definitions.
*/

#include "parameter_particle_validator_multi_loci_polytable_heterozygosity.hpp"
#include "multi_loci_polytable_set.hpp"

#include <stdexcept>

using namespace aabc;

ParameterParticleValidatorMultiLociPolyTableHeterozygosity::ParameterParticleValidatorMultiLociPolyTableHeterozygosity(
			const boost::shared_ptr < const mct::MultiLociPolyTable >& _rs,
			const boost::shared_ptr < mct::MultiLociPolyTableSampler >& _s,
			const size_t _n_validation_reps)
				:	ParameterParticleValidator( _rs->getMeanHeterozygosity() ),	
					referenceObj(_rs),
					sampler(_s),
					n_validation_reps(_n_validation_reps) {}

		
ParameterParticleValidatorMultiLociPolyTableHeterozygosity::~ParameterParticleValidatorMultiLociPolyTableHeterozygosity() 
{}
	
	

boost::shared_ptr< const mct::SummaryStatistic >
		ParameterParticleValidatorMultiLociPolyTableHeterozygosity::getSummaryStatistic(
				const boost::shared_ptr< const ParameterParticle>& p) const

{
	return (sampler->sample(n_validation_reps,
				*referenceObj, p))->getHeterozygositySummaryStatistic();
	/*
	return (sampler->sample(n_validation_reps,
								referenceObj->populationStructure(),
								referenceObj->nloci(),
								referenceObj->nsites(),
								p))->getSummaryStatistic();
	*/
}


boost::shared_ptr< mct::SummaryStatisticSet >
	ParameterParticleValidatorMultiLociPolyTableHeterozygosity::getSummaryStatisticSet(
			const boost::shared_ptr< const ParameterParticleSet>& pset)
															const
{
	boost::shared_ptr< mct::SummaryStatisticSet > 
			ss_set( new mct::SummaryStatisticSet( pset->size() ) );
	
	for (std::size_t i = 0; i < pset->size(); ++i) {
		
		( ss_set->add( getSummaryStatistic( pset->at(i) ) ) );
	}
	return ss_set;
	
}


