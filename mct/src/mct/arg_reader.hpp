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
\brief Declarations for ARGReader class to read ARGs from a file.
*/

#ifndef __ARG_READER_HPP__
#define __ARG_READER_HPP__

#include "arg_parser.hpp"
#include "nodes_and_trees.hpp"

#include <Sequence/Coalescent/SimTypes.hpp>

#include <boost/smart_ptr.hpp>

#include <string>
#include <vector>


namespace mct {
	
	/*! \brief Class to read an ARG from a file.
	 
	 Allows potential for different readers to be configured with
	 different parsing strategies, but for the moment just have a
	 single constructor which sets up a default parser.*/
	class ARGReader {
		
		public:
		
			/*! \brief No-argument constructor.
			
			Will set up a default parser.*/
			ARGReader();
			
			/*! \brief Destructor. */
			virtual ~ARGReader();

			/*! \brief Reads a string into a 
			 * \ref libsequence_arg "Sequence::arg"
			 format and returns a pointer to the new arg.
			 
			 The ARGReader uses its ARGParser parserPtr to read a string
			 into a \ref libsequence_arg "Sequence::arg".
			 
			 \param stringARG A string containing an ARG.
			 \post hist contains the ARG read in from the specified file.
			 \param nsites The number of sites this arg represents.
			 \return A pointer to the new Sequence::arg.*/
			boost::shared_ptr < Sequence::arg > 
								makeTreefileIntoArg(
								const std::string& stringARG,
								const size_t nsites);
		
		protected:
		
		/*! \brief A pointer to an ARGParser parsing strategy. */
		boost::shared_ptr < ARGParser > parserPtr; 
		
		/*! \brief A collection of \link mct::Tree Trees\endlink. */
		std::vector<Tree> tree_list; 
		
		private:
		
	};

} //end namespace mct

#endif
