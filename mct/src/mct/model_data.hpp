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

/*! \file     model_data.hpp
\brief ModelData struct to hold data about models.
*/

#ifndef _INC_MODEL_DATA_H
#define _INC_MODEL_DATA_H

#include <string>
#include <map>

namespace mct {

	 
	 
	
	/*! \brief Struct to hold data about a Model. */
	struct ModelData {
		
		/*! \brief A mapping from string identifiers to double values.
		
		Used for model-specific parameter names and values.*/ 
		typedef std::map<std::string, double> ParamD;

	
		/*! \brief The name of the model (eg "JC69"). */
		std::string name; 
		
		/*! \brief The state space of the model.  
		
		See \ref models_statespace "model state space parameter"
		and \ref models_sites "model sites". */
		int ss; 
		
		/*! \brief The per-site mutation rate \f$ \theta \f$
		for the model.*/
		double theta; 
		
		/*! \brief The number of contiguous sites to which the model applies. */
		size_t repeats; 
		
		/*! \brief \link mct::Model Model-specific\endlink 
							parameter names and values. */
		ParamD params; 
	};

} // end namespace mct


#endif
