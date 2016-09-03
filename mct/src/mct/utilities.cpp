/* MCT - Markov Chains on Trees.


   Copyright (C) 2012 Jennifer Harlow

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
\brief Utilities for testing and using mct

*  
*/

#include "config.h"

#include "utilities.hpp"


#include <fstream>
#include <iostream>
#include <iterator>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <cstdlib>


namespace mct_utilities {


	void outputToFile(const std::string& line, 
				const std::string& s,
				bool append)
	{
		std::ofstream os;
			
		if (append) os.open(s.c_str(), std::ios::app);         // append
		else os.open(s.c_str()); // overwrite
		
		if (os.is_open()) {
			os << line << "\n";
			os.close();
		}
		
		else {
			std::cerr << "Error: could not open file named "
				<< s << std::endl << std::endl;
		}
	}

	void outputVectorVectorToFile(const std::vector < std::vector < double > >& vecvec,
						const std::string& s,
						bool append)
	{
		std::ofstream os;
			
		if (append) os.open(s.c_str(), std::ios::app);         // append
		else os.open(s.c_str()); // overwrite
		
		if (os.is_open()) {
			
			os.precision(5);
			os <<   std::fixed;
			
			if (!vecvec.empty()) {		
				for (std::vector< std::vector< double  > >::const_iterator it = vecvec.begin();
						it < vecvec.end(); ++it) {
					
					if (!(it->empty())) {
						std::vector< double >:: const_iterator e_it = it->end();
						std::advance(e_it,-1); // move it back one	
						copy(it->begin(), e_it, std::ostream_iterator<double>(os, "\t"));
						os << it->back(); // last one
						os << "\n";			// new line
					}
				}
			}
			os.close();
			
		}
		else {
			std::cerr << "Error: could not open file named "
				<< s << std::endl << std::endl;
		}
	}


	std::vector< std::vector< double  > >& rotate(std::vector< std::vector< double  > >& vec)
	{
		std::vector< std::vector< double  > > tmp;
		tmp.swap(vec);
		
		if ( !tmp.empty() && !tmp[0].empty() ) {
			for (size_t i = 0; i < tmp[0].size(); ++i) {
		
					
				vec.push_back( std::vector< double  >() );
			
				for (size_t j = 0; j < tmp.size(); ++j) {
					
					vec.back().push_back( tmp[j][i] );
				}
			}
		}
		
		
		return vec;
	}

	int fitPrecision(double data)
	{
		// fit the precision to the magnitudes of the data
		int result = 1;
		while(data != 0.0 && data < 0.1) {
			data*= 10;
			result++;
		}
		return result;
	}
	
	/* Open a file and return the data as a string */
	std::string string_from_file(const std::string name)
	{
		std::ifstream is(name.c_str());
		if (!is.is_open())
			throw std::runtime_error("unable to open file '" + name + "'");
		
		return std::string((std::istreambuf_iterator<char>(is)),
					 std::istreambuf_iterator<char>());
	}
	
	std::string makeDir(const std::string& baseOutputDir, 
				const std::string& thisDir)
	{
		std::string dStr; // string identifying this directory
		{	
			std::ostringstream oss;
			oss << baseOutputDir << thisDir;
			dStr = oss.str(); 
		}
		
		/*command make the output directory (-p switch makes intermediates
		 *  as well, okay if already exists */
		std::string  mkdirCmd;
		{
			std::ostringstream oss;
			oss << "mkdir -p " 
				<< dStr;
			 mkdirCmd = oss.str(); 
		}	
		// and system call to make it
		system( mkdirCmd.c_str() );
		
		return dStr;
	}



}	// end namespace
