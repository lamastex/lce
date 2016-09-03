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
\brief ParameterParticleValidatorGenericDataObject definitions.
*/

#include "parameter_particle_validator_generic_data_object.hpp"
#include "generic_data_object_sampler.hpp"
#include "generic_data_object_set.hpp"
#include "data_object.hpp"
#include "parameter_particle_set.hpp"
#include "summary_statistic_set.hpp"

#include <stdexcept>

using namespace aabc;

ParameterParticleValidatorGenericDataObject::ParameterParticleValidatorGenericDataObject(
						const boost::shared_ptr < const DataObject >& _rs,
						const boost::shared_ptr < GenericDataObjectSampler >& _s,
						const size_t _n_validation_reps)
				:	ParameterParticleValidator( _rs->getSummaryStatistic() ),	
					referenceObj(_rs),
					sampler(_s),
					n_validation_reps(_n_validation_reps) {}


		
ParameterParticleValidatorGenericDataObject::~ParameterParticleValidatorGenericDataObject() 
{
		//std::cout << "ParameterParticleValidatorGenericDataObject destructor called" << std::endl;
		
}

boost::shared_ptr< const mct::SummaryStatistic >
	ParameterParticleValidatorGenericDataObject::getSummaryStatistic(
			const boost::shared_ptr< const ParameterParticle>& p) const
{
	return (sampler->sample(n_validation_reps,
								p))->getSummaryStatistic();
	
}


boost::shared_ptr< mct::SummaryStatisticSet >
	ParameterParticleValidatorGenericDataObject::getSummaryStatisticSet(
			const boost::shared_ptr< const ParameterParticleSet>& pset) const
{
	boost::shared_ptr< mct::SummaryStatisticSet > 
			ss_set( new mct::SummaryStatisticSet( pset->size() ) );
	
	for (std::size_t i = 0; i < pset->size(); ++i) {
		( ss_set->add( getSummaryStatistic( pset->at(i) ) ) );
	}
	return ss_set;
}		

