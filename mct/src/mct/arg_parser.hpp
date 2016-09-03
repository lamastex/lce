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

/*! \file     arg_parser.hpp
\brief Declarations for abstract ARGParser class and subclasses.
* 
* Parser classes parse an ARG given in a string.
*/


#ifndef INC_ARGPARSER_H
#define INC_ARGPARSER_H

#include "nodes_and_trees.hpp"

#include <string>
#include <vector>

namespace mct {

	/*! @brief Interface for a parser which can parse an Ancestral
	Recombination Graph (ARG) */
	class ARGParser {
		
		public:
			/*! \brief Parses a string into a collection of 
			 \link mct::Tree Trees\endlink. */
			virtual std::vector<Tree> parse(const std::string& str) = 0;
			
			/*! \brief Destructor. */
			virtual ~ARGParser() {}

	};


	/*! @brief A recursive-descent parser which can parse an Ancestral
	Recombination Graph (ARG) given in Newick format.
	* 
	
	By default each line of the file is assumed to describe the marginal
	tree apply to one locus, unless the number of loci to which the 
	marginal tree applies is specified with an integer in square 
	brackets preceeding the tree description.  For example, 
	\code
	(2:1.059,((1:0.042,4:0.042):0.313,(3:0.058,5:0.058):0.297):0.704);
	\endcode
	is a marginal tree applying to one locus, as is
	\code
	[1](2:1.059,((1:0.042,4:0.042):0.311,(3:0.058,5:0.058):0.296):0.706);
	\endcode
	whereas 
	\code
	[5]((1:0.042,4:0.042):1.017,(2:0.258,(3:0.058,5:0.058):0.200):0.802);
	\endcode
	specifies a marginal tree for 5 contiguous loci.
	
	*/
	class NewickParser : public ARGParser {
		
		public:
			/*! \brief Parses a string employing Newick format into a
			 collection of \link mct::Tree Trees\endlink. 
			 
			 \param str The string to be parsed.
			 \return A collection of \link Tree Trees\endlink. */
			virtual std::vector<Tree> parse(const std::string& str);
		
		private:
		
			enum types {
				STRING, 
				SYMBOL 
			};
			
			struct Token {
				std::string val;
				types type;
				int line;
				int col;
				Token(std::string _val, types _type, int _line, int _col)
					: val(_val), type(_type), line(_line), col(_col) {};
			};
			
			std::vector<Token> tokens;
			
			size_t pos;
			
			bool accept(std::string c, bool required=false);
			
			bool accept_name(std::string &d, bool required=false);
			
			bool accept_double(double &d, bool required=false);
			
			bool accept_integer(int &d, bool required=false);
			
			bool read_node(Node &node);
			
			void set_node_ids(Tree &tree);
			
			bool read_tree(Tree &tree);
			
			bool is_symbol(char c);
			
			/** converts a string into a vector of tokens */
			void get_tokens(const std::string& str);
			
		
	};

} // end namespace mct

#endif /* INC_ARGPARSER_H */



