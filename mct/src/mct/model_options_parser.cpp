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

/*! \file     model_options_parser.cpp
\brief Definitions for a class to parse model options from a string.
*/

#include "model_options_parser.hpp"

#include <map>
#include <iostream> 
#include <cstdlib> // atoi
#include <cstdio> // sscanf etc

#define N_MATCH 4

using namespace mct;

std::vector<ModelData>& 
	ModelOptionsParser::parse(std::vector<ModelData>& mdata_list,
									const std::string& input)
{
    regex_t reg_model, reg_ws, reg_param;
	
	//take a backup of the given mdata_list in case of problems
	std::vector<ModelData>& mdata_list_backup = mdata_list;
	
	try {
		regmatch_t reg_match[N_MATCH];
		
		/* matches " AcD4 [123] blah" */
		const char *re1 = "([a-zA-Z0-9]+)[ \t\v\f]*"
						  "\\[([0-9]+)\\]"
						  "(.*)";
		
		/* matches "ABc3_ = 234.34 " */
		const char *re2 = "([a-zA-Z0-9_]+)[ \t]*=[ \t]*([^ \t]+)[ \t]*";
		
		regcomp(&reg_ws, "^[ \t\v\f]*$", REG_EXTENDED);
		regcomp(&reg_model, re1, REG_EXTENDED);
		regcomp(&reg_param, re2, REG_EXTENDED);
		
		std::vector<std::string> lines = split_at(input, '\n');
		
		for (size_t i=0; i<lines.size(); ++i) {
			
			
			// remove comments
			size_t comment = lines[i].find_first_of('#');
			if (comment != std::string::npos) lines[i].erase(comment);
			
			// skip empty lines
			if (!regexec(&reg_ws, lines[i].c_str(), 0, NULL, 0))
				continue;
			
			ModelData model;
					
			if (!regexec(&reg_model, lines[i].c_str(), N_MATCH, reg_match, 0)) {
				std::vector<std::string> matches;
				for (int j=0; j<N_MATCH; ++j) {
					matches.push_back(match_str(lines[i], &reg_match[j]));
				}
				
				// save data
				model.name = matches[1];
				model.ss = -1;
				model.theta = 0.0;
				model.repeats = atoi(matches[2].c_str());
				
				std::string param_str = matches[3];
				
				std::map<std::string, std::string> tmp_params;
				
				// extract parameter info
				while(!regexec(&reg_param, param_str.c_str(), 3, reg_match, 0)) {
					std::string fst, snd;
					fst = match_str(param_str, &reg_match[1]);
					snd = match_str(param_str, &reg_match[2]);
					tmp_params[fst] = snd;
					param_str = param_str.substr(reg_match[0].rm_eo);
				}
				
				if (param_str.size() > 0) {
					std::cerr << "error on line " << i 
							<< ". Ignoring." 
											<< std::endl;
					continue;
				}
				
				// copying params to structure, and reading ints/floats
				std::map<std::string, std::string>::iterator iter;
				for (iter = tmp_params.begin(); iter != tmp_params.end(); ++iter) {
					
					// special case of state space, which is an integer
					if (iter->first == "ss") {
						int val2;
						if (sscanf(iter->second.c_str(), "%d", &val2) != 1) {
							std::cerr << "error: not an integer" << std::endl;
							continue;
						}
						model.ss = val2;
						continue;  // jump to next one in for loop
					}
					// special case of theta
					if ((iter->first == "t") || (iter->first == "theta")) {
						double val1;
						if (sscanf(iter->second.c_str(), "%lf", &val1) != 1) {
							std::cerr << "error: not a floating point number" << std::endl;
							continue;
						}
						model.theta = val1;
						continue;  // jump to next one in for loop
					}
					
					double val;
					if (sscanf(iter->second.c_str(), "%lf", &val) != 1) {
						std::cerr << "error: not a floating point number" << std::endl;
						continue;
					}
					model.params[iter->first] = val;
				}
			}
			mdata_list.push_back(model);
			
		}
		
		regfree(&reg_ws);
		regfree(&reg_model);
		regfree(&reg_param);
		
		return mdata_list; // return by reference
	}
	catch (...) {
		
		//restore the old state of mdata_list
		mdata_list = mdata_list_backup;
		
		try {
			regfree(&reg_ws);
		}
		catch (...) {}
		try {
			regfree(&reg_model);
		}
		catch (...) {}
		try {
			regfree(&reg_param);
		}
		catch (...) {}
		
		throw;
	}
}

/* splits a std::string into a std::vector of std::strings, delimited by c */
std::vector<std::string> 
			ModelOptionsParser::split_at(const std::string& s, char c)
{
    std::vector<std::string> result;
    size_t i, prev = 0;
    for (i=0; i<s.size(); ++i) {
        if (s[i] == c) {
            result.push_back(s.substr(prev, i - prev));
            prev = i+1;
        }
    }
    result.push_back(s.substr(prev, i - prev));
    
    return result; 
}

std::string 
	ModelOptionsParser::match_str(const std::string &s, regmatch_t *m)
{
    if (m->rm_so == -1 || m->rm_eo - m->rm_so == 0) return "";
    return s.substr(m->rm_so, m->rm_eo - m->rm_so);
}



