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
\brief ModelSetConfigBuilder definitions.
*/

#include "model_set_config_builder.hpp"
#include "model_options_parser.hpp" 


using namespace mct;

ModelSetConfigBuilder::ModelSetConfigBuilder() {}

ModelSetConfigBuilder::~ModelSetConfigBuilder() {}

void ModelSetConfigBuilder::parseModelConfigsFromString(
									const std::string& stringModels)
{
	mdata_list = ModelOptionsParser::parse(mdata_list, stringModels);
}

										
void ModelSetConfigBuilder::addModelConfig(const std::string& _name,
										size_t _repeats,
										double _theta,
										ModelData::ParamD _params,
										int _ss) 
{
	std::vector < ModelData>::iterator it;
	size_t oldsize = mdata_list.size();
	
	try {
		it = mdata_list.insert( mdata_list.end(), ModelData() );
		
		mdata_list.back().name = _name;
		mdata_list.back().ss = _ss;
		mdata_list.back().theta = _theta;
		mdata_list.back().repeats = _repeats;
		mdata_list.back().params = _params;
			
	}
	
	catch (...) {
		
		//take out the model data that failed
		try {
			if ( mdata_list.size() > oldsize ) mdata_list.erase(it);
		}
		catch (...) {} // catch and swallow
		
		throw;
	}
}

void ModelSetConfigBuilder::addModelConfig(const std::string& name,
										size_t repeats,
										double th,
										int ss) 
{
	ModelData::ParamD params;
	addModelConfig(name, repeats, th, params, ss);
}


void ModelSetConfigBuilder::addModelConfig(const std::string& name,
										size_t repeats,
										int ss) 
{
	double theta = 0.0;
	ModelData::ParamD params;
	addModelConfig(name, repeats, theta, params, ss);
}
										
void ModelSetConfigBuilder::addModelConfig(const std::string& name,
										size_t repeats,
										ModelData::ParamD params,										
										int ss) 
{
	double theta = 0.0;
	addModelConfig(name, repeats, theta, params, ss);
}
	
const std::vector< ModelData >& ModelSetConfigBuilder::getModelSetConfigs()
																const
{
	return mdata_list;
}
