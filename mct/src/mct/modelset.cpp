/* MCT - Markov Chains on Trees.

   Copyright (C) 2009 Brendan Bycroft <brb44@student.canterbury.ac.nz>

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
\brief ModelSet definitions.
*/

#include "modelset.hpp"
#include "modelfactory.hpp"

#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <gsl/gsl_randist.h>

#include <sstream>
#include <stdexcept>
#include <cassert>

//#define MYDEBUG
#if defined (MYDEBUG)
	#include <iostream>
	
#endif

using namespace mct;

ModelSet::ModelSet(const std::vector<ModelData>& _mdata_list,
					const boost::shared_ptr < PRNGen >& _r_ptr)
					: mdata_list(_mdata_list)
{
	locusMap = boost::shared_ptr < MCTLocusMap >( new MCTLocusMap() );

	r_rand_ptr = _r_ptr->spawnAnother();
	
	boost::function<double ()> func =
		boost::bind<double>(&mct::PRNGen::uniform01_double, r_rand_ptr);
		//boost::bind<double>(gsl_ran_flat, r, 0.0, 1.0);
	
	/*  construct the appropriate models using the mdata_list.*/
	size_t from = 0;

	for (unsigned int i=0; i<mdata_list.size(); ++i) {
		ModelData *mdata = &mdata_list[i];
		
		/* create the models */
		try {
			if (mdata->repeats < 1) {
				throw std::invalid_argument(
					std::string("ModelSet::ModelSet(")
					+ std::string("const std::vector<ModelData>&, ")
					+ std::string("const boost::shared_ptr < PRNGen >&)")
					+ std::string(":\nmdata repeats < 1") );
			}
			boost::shared_ptr < Model > m 
					= ModelFactory::create_model(mdata->name, func,
						   mdata->ss, mdata->theta, mdata->params);
			// each model is using the same prng function func               
			m_set.push_back(m);
			for (size_t i=0; i<mdata->repeats; ++i) {
				m_repeated_set.push_back(m);
			}
			size_t to = from + mdata->repeats - 1;
			
			// add to the map
			locusMap->add(from, to, m->modelType() );
			
			from = to + 1;
			
		}
		catch (std::exception &e) {
			std::ostringstream stm;
			stm << "Error creating model " << (i+1) << ":\n" << e.what() << std::endl;
			throw;
		}
		catch (...) {
			std::ostringstream stm;
			stm << "Error creating model " << (i+1) << ":\n" << std::endl;
			throw;
		}
	}
			
	assert(from == m_repeated_set.size());
}
	
	#ifdef MYDEBUG
		size_t all_models_i = 0;
		for (size_t model_i = 0; model_i < m_set.size(); model_i++) {
				std::cout << "this model_i = " << model_i << " and address = " << &(*(m_set[model_i])) << std::endl;
				ModelData *mdata = &mdata_list[model_i];
				for (size_t ami = 0 ; ami < mdata->repeats; ami++, all_models_i++) {
					std::cout << "this all_model_i = " << all_models_i << " and address = " << &(*(m_repeated_set[all_models_i])) << std::endl;
				}
		}
	#endif

ModelSet::~ModelSet() {}

std::string ModelSet::toString() const
{
	std::ostringstream stm;
	for (unsigned int i=0; i<mdata_list.size(); ++i) {
        //const ModelData *mdata = &mdata_list[i];
	    
		//stm << "Model " << mdata->name << ": ss = " << mdata->ss;
		stm << (m_set[i])->toString();
		stm << "\t" << mdata_list[i].repeats << " repeats\n\n";
		/*
		std::map<std::string, double>::const_iterator it;
		for (it=mdata->params.begin(); it!=mdata->params.end(); ++it) {
			stm << ", " << it->first << " = " << it->second;
		}
		*/
		
	}

	return stm.str();
}

size_t ModelSet::getNsites() const
{
	return m_repeated_set.size();
}

size_t ModelSet::size() const
{
	return m_set.size();
}

int ModelSet::step_at(const size_t index,
						const int current_state,
						const double len) const
{
	return (m_repeated_set.at(index))->step(current_state, len);
}

int ModelSet::ancestor_at(const size_t index) const
{
	return (m_repeated_set.at(index))->ancestor();
}

void ModelSet::resetThetas(const std::vector < double >& new_thetas)
{
	if( new_thetas.size()!= size() ) {
		throw std::invalid_argument(
		"ModelSet::resetThetas(const std::vector < double >&)");
	}
	std::vector< double >::const_iterator cit = new_thetas.begin();
	for (model_set::iterator it = m_set.begin();
			it < m_set.end();
			++it, ++cit) {
				(*it)->resetTheta(*cit);
	}
}

void ModelSet::resetThetas(double new_theta)
{
	for (model_set::iterator it = m_set.begin();
			it < m_set.end();
			++it) {
				(*it)->resetTheta(new_theta);
	}
}

const boost::shared_ptr <MCTLocusMap> ModelSet::getLocusMap() const
{
	return locusMap;
}




