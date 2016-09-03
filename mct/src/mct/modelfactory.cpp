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
\brief ModelFactory definitions.
*/

#include "modelfactory.hpp"

#include <boost/lambda/lambda.hpp>
#include <boost/lambda/construct.hpp>
#include <boost/lambda/bind.hpp>

#include <sstream>
//#include <exception>

using namespace mct;


// static members
ModelFactory::ModelInfoMap ModelFactory::model_info;

bool ModelFactory::_dummy = load_model_info(model_info);

/* Model information.

   createFunc type parameter (eg ModelJC69) is a constructor matching
   specification in the ModelInfo structure; boost::lambda::new_ptr uses
   this constructor to create a pointer to a newed object of that type.
*/
    
bool ModelFactory::load_model_info(ModelInfoMap& m_info) {
    m_info.clear();

    ModelInfo m;

    m.ss = 4;
	m.model_name = "Jukes and Cantor, 1969";
    m.DNA_model = DNA;
	m.lower_bound = 0;
	m.theta = 0.0;
	m.createFunc = boost::lambda::new_ptr<ModelJC69>();
    m.opts.clear();
    m_info["JC69"] = m;
    
    m.ss = 40;
	m.model_name = "Stepwise mutation model";
    m.DNA_model = MICROSAT;
	m.lower_bound = 10;
	m.theta = 0.0;
	m.createFunc = boost::lambda::new_ptr<ModelSMM>();
    m.opts.clear();
    m_info["SMM"] = m;

    m.ss = 40;
	m.model_name = "REV";
    m.DNA_model = MICROSAT;
	m.lower_bound = 10;
	m.theta = 0.0;
	m.createFunc = boost::lambda::new_ptr<ModelREV>();
    m.opts.clear();
    m_info["REV"] = m;

    m.ss = 40;
	m.model_name = "AQB";
    m.DNA_model = MICROSAT;
	m.lower_bound = 10;
	m.theta = 0.0;
	m.createFunc = boost::lambda::new_ptr<ModelAQB>();
    m.opts.clear();
    m.opts["threskappa"] = 0.0;
    m.opts["u0"] = 0.5;
    m.opts["d0"] = 0.5;
    m.opts["u1"] = 0.005;
    m.opts["d1"] = 0.005;
    m.opts["u2"] = 0.03;
    m.opts["d2"] = 0.03;
    m.opts["avgmr"] = 1.0;
    m_info["AQB"] = m;

    m.ss = 40;
	m.model_name = "EXP";
    m.DNA_model = MICROSAT;
	m.lower_bound = 10;
	m.theta = 0.0;
	m.createFunc = boost::lambda::new_ptr<ModelEXP>();
    m.opts.clear();
    m.opts["gammaU"] = 3.1e-6;
    m.opts["gammaD"] = 4e-7;
    m.opts["alphaU"] = 0.2;
    m.opts["alphaD"] = 0.302;
    m.opts["lambdaU"] = 1.06;
    m.opts["lambdaD"] = 1.06;
    m.opts["avgmr"] = 1.0;
    m_info["EXP"] = m;

    m.ss = 40;
	m.model_name = "PL";
    m.DNA_model = MICROSAT;
	m.lower_bound = 10;
	m.theta = 0.0;
	m.createFunc = boost::lambda::new_ptr<ModelPL>();
    m.opts.clear();
    m.opts["u"] = 0.6246;
    m.opts["v"] = 0.01542;
    m.opts["m"] = 1.0;
    m.opts["s"] = 0.8752;
    m.opts["p"] = 0.0;
    m.opts["bdryPLo"] = 0.5;
    m.opts["bdryPHi"] = 0.5;
    m.opts["avgmr"] = 1.0;
    m_info["PL1"] = m;

    m.ss = 40;
	m.model_name = "PL";
    m.DNA_model = MICROSAT;
	m.lower_bound = 10;
	m.theta = 0.0;
	m.createFunc = boost::lambda::new_ptr<ModelPL>();
    m.opts.clear();
    m.opts["u"] = 0.8158;
    m.opts["v"] = 0.03947;
    m.opts["m"] = 0.5475;
    m.opts["s"] = 0.7638;
    m.opts["p"] = 0.0;
    m.opts["bdryPLo"] = 0.5;
    m.opts["bdryPHi"] = 0.5;
    m.opts["avgmr"] = 0.761;
    m_info["PL2"] = m;
	
	return true;
}

/*
boost::shared_ptr < Model >
ModelFactory::create_model(std::string name, boost::function<double ()> rand) {
    return ModelFactory::create_model(name, rand, 0, ModelData::ParamD());
}
*/

boost::shared_ptr < Model >
ModelFactory::create_model(
		std::string name, boost::function<double ()> rand,
		int ss,
		double th,
		ModelData::ParamD p) {
	
	/* load_model_info should have been called. */
	if (model_info.empty())
		load_model_info(model_info);

	ModelData::ParamD::iterator it_a, it_b;

	if (model_info.find(name) == model_info.end()) {
		std::ostringstream stm;
			stm << ":\nUnknown model: " << name << std::endl;
			throw std::invalid_argument(
			std::string("ModelFactory::create_model(std::string, boost::function<double ()>, ")
			+ std::string("int, double, ModelData::ParamD)") 
			+ stm.str() );
	}
		
	/* pad params with default values */
	if (ss <= 0) {
		ss = model_info[name].ss;
	}
	if (th <= 0.0) {
		th = model_info[name].theta;
	}

	ModelData::ParamD p_new = model_info[name].opts;

	for (it_a = p.begin(); it_a != p.end(); ++it_a) {
		it_b = p_new.find(it_a->first);

		if (it_b == p_new.end()) {
			std::ostringstream stm;
			stm << ":\nParameter " << it_a->first << 
			" does not exist for the " << name << " model" << std::endl;
			throw std::invalid_argument(
			std::string("ModelFactory::create_model(std::string, boost::function<double ()>, ")
			+ std::string("int, double, ModelData::ParamD)") 
			+ stm.str() );
		}
		it_b->second = it_a->second;
	}

	return  boost::shared_ptr <Model> ((Model *) (model_info[name].createFunc) (rand,
									ss, 
									model_info[name].model_name,
									model_info[name].DNA_model,
									model_info[name].lower_bound,
									th,
									p_new));
	
}



/** Print all model parameters to cout. */
void ModelFactory::printAvailableModels() {
    /* load_model_info should have been called. */
	if (model_info.empty())
		load_model_info(model_info);

	std::cout << "Models and default parameters:\n\n";
	
	ModelInfoMap::iterator it;
	for (it=model_info.begin(); it!=model_info.end(); ++it)
	{
		ModelInfo* info = &it->second;
		std::cout << "Model " << it->first << ": State Space (ss) = " << info->ss << "\n";
		std::cout << "Theta (t) = " << info->theta << "\n";
		
		ModelData::ParamD::iterator jt;
		for (jt=info->opts.begin(); jt!=info->opts.end(); ++jt) {
			std::cout << "\t" << jt->first << " = " << jt->second << "\n";
		}
		std::cout << std::endl;
	}
}


