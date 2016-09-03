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
\brief ModelSetConfigBuilder declarations.
*/

#ifndef _INC_MODEL_SET_CONFIG_BUILDER_H
#define _INC_MODEL_SET_CONFIG_BUILDER_H

#include "model_data.hpp"

#include <string>
#include <vector>


namespace mct {
	
	/** @brief A class for building the  
	ModelData configuration options 
	for a ModelSet.
	
	A ModelSetConfigBuilder makes it easy to build up the configuration
	or specification (see ModelData) for a ModelSet  within 
	a program.
	
	The ModelSetConfigBuilder can either parse a model sites (model
	configuration) file (see parseModelConfigsFromString()) or
	the configuration can be built up by adding the name of a Model
	and the number of sites it covers (see addModelConfig()).  
	
	The configuration can be retrieved (for example to be passed 
	to the ModelSet constructor) using getModelSetConfigs().
	  */
	 

	class ModelSetConfigBuilder {
		public:
		
			/*! \brief Constructor.*/
			ModelSetConfigBuilder();
			
			virtual ~ModelSetConfigBuilder();
			
			/*! \brief Build the configuration options by parsing a string.
			 
			\param stringModels The string of model configurations to 
			parse. */
			virtual void parseModelConfigsFromString(
										const std::string& stringModels);
			
			/*! \brief Add a model configuration using given parameters.
			 
			\param _name The name of the model to add.  Must match an
			available model name.
			\param _repeats The number of consecutive sites covered
			by this model.
			\param _theta The per site coalescent-scaled mutation rate
			for the model, ie \f$ 4 N_0 \mu_{site} \f$ where
			\f$ \mu_{site} \f$ is the neutral per site migration rate.
			\param _params Extra parameters in the form of a 
			ModelData::ParamD object.
			\param _ss The state space, defaulting to -1.  */
			virtual void addModelConfig(const std::string& _name,
										size_t _repeats,
										double _theta,
										ModelData::ParamD _params,
										int _ss = -1);
										
			virtual void addModelConfig(const std::string& _name,
										size_t _repeats,
										ModelData::ParamD _params,
										int _ss = -1);
										
			virtual void addModelConfig(const std::string& _name,
										size_t _repeats,
										int _ss = -1);
										
			virtual void addModelConfig(const std::string& _name,
										size_t _repeats,
										double _theta,
										int _ss = -1);
			
			/*! \brief Return a collection of Model configurations.
			 
			\return A copy of the currrent collection of 
			configurations. */			
			const std::vector< ModelData >& getModelSetConfigs() const; 
			
			

		protected:

		
		private:
			
			/*! \brief Model configuration options for one or
			more \link mct::Model Models\endlink.*/
			std::vector< ModelData > mdata_list; 

			
	};

} // end namespace mct

#endif 

