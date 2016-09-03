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
\brief ModelFactory declarations.
*/

#ifndef _INC_MODELFACTORY_H
#define _INC_MODELFACTORY_H

#include "model.hpp"
#include "mct_model_data_types.hpp"

#include <string>
#include <map>
//#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/smart_ptr.hpp>

namespace mct {

	/*! \brief Class to create Markov chain \link mct::Model Models\endlink.

	\note ModelFactory is a static class, i.e. a ModelFactory object 
	cannot be instantiated and the class only contains static members.
	\note The available models, their model-specific parameter names, 
	and their default parameter values, are specified
	in the protected and private members of this class.   
	 */
	class ModelFactory {
		
		
		public:
			
			/*! \brief Create a model using given state space and parameters.
			
			\param name The name of the Model to create.
			\param rand A boost::function returning a value from a
			uniform(0,1) distribution.
			\param ss The state space for the Model.
			\param th The per-site mutation rate \f$ \theta \f$
			 for the Model.  \f$ \theta = 4N_0 \mu_{site} \f$
			\param p Parameters to override the defaults for the Model.
			\exception MCTexception error creating Model.
			\return A pointer to an instance of the specified Model with
			given state space and parameters.
			*/
			static boost::shared_ptr < Model > create_model(
					std::string name, boost::function<double ()> rand,
					int ss,
					double th,
					ModelData::ParamD p);
										
			/*! \brief Print available models and default parameters
			to standard output.
			*/
			static void printAvailableModels();

		protected:

			/* A structure to hold model information.

			   Model constructors are specified in the ModelInfo as the createFunc.
			*/
			struct ModelInfo {
				int ss;
				double theta;
				std::string model_name;
				DATATYPE DNA_model;
				int lower_bound;
				boost::function<Model* (boost::function<double ()> rand,
										int ss,
										std::string& _name,
										DATATYPE _dna,
										int _lb,
										double th, 
										ModelData::ParamD& p)> createFunc;
				ModelData::ParamD opts;
				
			};

			typedef std::map<std::string, ModelInfo> ModelInfoMap;

			static ModelInfoMap model_info;
			
			static bool _dummy;
					
			static bool load_model_info(ModelInfoMap& m_info);
		
		private:
		
		ModelFactory(); // can't instantiate class

	};

} // end namespace mct

#endif /* _INC_MODELFACTORY_H */

