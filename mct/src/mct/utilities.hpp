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

#ifndef __UTILITIES_HPP__
#define __UTILITIES_HPP__

#include <vector>
#include <string>

namespace mct_utilities {

	void outputToFile(const std::string& line, 
				const std::string& s,
				bool append = true);
				
	void outputVectorsToFile(const std::vector < double >& first,
						const std::vector < double >& second,
						const std::string& s,
						bool append = true);
						
	void outputVectorToFile(const std::vector < double >& first,
						const std::string& s,
						bool append = true);

	void outputVectorVectorToFile(const std::vector < std::vector < double > >& vecvec,
						const std::string& s,
						bool append = true);

	std::vector< std::vector< double  > >& rotate(std::vector< std::vector< double  > >& vec);

	int fitPrecision(double data);
	
	/*! \brief Open a file and return the data as a string.
	 
	 \param name The name of the file to open.
	 \return a string version of the file.
	  */
	std::string string_from_file(const std::string name);
	
	std::string makeDir(const std::string& baseOutputDir, 
				const std::string& thisDir);

}
#endif

