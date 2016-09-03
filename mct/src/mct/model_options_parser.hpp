/* MCT - Markov Chains on Trees.

   Copyright (C) 2010 Brendan Bycroft <brb44@student.canterbury.ac.nz>

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

/*! \file     model_options_parser.hpp
\brief Declarations for a class to parse model options from a string.
*/


#ifndef _INC_MODEL_OPTIONS_PARSER_H
#define _INC_MODEL_OPTIONS_PARSER_H

#include "model_data.hpp"

#include <regex.h>

#include <vector>
#include <string>

namespace mct {

	/*! \brief Class to parse Model options from a string.
	 \note ModelOptionsParser is a static class, 
	i.e. a ModelOptionsParser object cannot
	be instantiated and the class only contains static members.
	 */
	class ModelOptionsParser
	{
		public:
		
			/*! \brief Parse Model options from a string.
			\param mdata_list A reference to the ModelData to be filled in
			and returned by the parsing operation.
			\param input The string to be parsed.
			\return A reference to the ModelData \a mdata_list supplied to be filled in
			and returned by the parsing operation.
			\post \a mdata_list contains the ModelData parsed from the 
			input string.
			*/
			static std::vector<ModelData>& //return by reference
					parse(std::vector<ModelData>& mdata_list,
									const std::string& input);
		
		protected:
		
			static std::vector<std::string> split_at(const std::string& s,
													char c);

			static std::string match_str(const std::string &s, 
									regmatch_t *m);
		
		private:
			ModelOptionsParser(); // class cannot be instantiated
	};


} // end namespace mct

#endif
