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

/*! \file     arg_parser.cpp
\brief Definitions for parser subclasses.
*/

#include "arg_parser.hpp"

#include <sstream>
#include <algorithm>
#include <stdexcept>
#include <cassert>

//#define MYDEBUG
#if defined (MYDEBUG)
	#include <iostream>
	
#endif

using namespace mct;

//public:
std::vector<Tree> NewickParser::parse(const std::string& str) {
	
	std::vector<Tree> arg;
	get_tokens(str); // convert the std::string to tokens
	pos = 0;
	while (pos < tokens.size()) {
		Tree tree;
		read_tree(tree); // read the tokens into the tree.
		arg.push_back(tree);
	}
	return arg;
}

    
bool NewickParser::accept(std::string c, bool required) {
	if (tokens[pos].val != c) {
		if (required) {
			std::ostringstream stm;
			stm << "NewickParser::accept(std::string, bool) :\nexpected '" << c << "'";
			throw std::logic_error(stm.str());
		}
		else
			return false;
	}
	pos += 1;
	return true;
}

bool NewickParser::accept_name(std::string &d, bool required) {
	if (tokens[pos].type != STRING) {
		if (required) {
			throw std::logic_error(
			"NewickParser::accept_name(std::string, bool) :\nexpected a name");
		}
		else
			return false;
	}
	
	d = tokens[pos].val;
	pos += 1;
	return true;
}

bool NewickParser::accept_double(double &d, bool required) {
	std::istringstream ss(tokens[pos].val);
	ss >> d;
	
	#ifdef MYDEBUG
		std::cout << "accept_double: tokens[pos].val = " << tokens[pos].val << std::endl;
		std::cout << "d = " << d << std::endl;
	#endif
	
	if (!ss || !ss.eof()) {
		if (required)
			throw std::logic_error(
			"NewickParser::accept_double(double, bool) :\nexpected a floating point number");
		else
			return false;
	}
	pos += 1;
	return true;
}

bool NewickParser::accept_integer(int &d, bool required) {
	std::istringstream ss(tokens[pos].val);
	
	ss >> d;
	if (!ss || !ss.eof()) {
		if (required)
			throw std::logic_error(
			"NewickParser::accept_integer(int, bool) :\nexpected an integer\n");
		else
			return false;
	}
	
	
	
	pos += 1;
	return true;
}

bool NewickParser::read_node(Node &node) {
	std::string name;
		
	if (accept("(")) {
		while (true) {
			Node child;
			read_node(child);
			node.children.push_back(child);
			if (!accept(",")) break;
		}
		accept(")", true);
	}
	
	node.name = "";
	accept_name(node.name);
	
	if (accept(":"))
		accept_double(node.length, true);
	else
		node.length = 0.0;
	
	if (node.name == "" && node.children.size() == 0) {
		throw std::logic_error("NewickParser::read_node(Node &) :\nexpected a node");
	}
	
	return true;
}

bool NewickParser::read_tree(Tree &tree) {
	bool haveRepeats = accept("[", false);
	if (haveRepeats) {
		int temp = 0;
		accept_integer(temp, true);
		if (temp < 0 ) {
			throw std::invalid_argument(
			"NewickParser::read_tree(Tree &) :\nrepeats < 0");
		}
		tree.repeats = static_cast<size_t> (temp);
		accept("]", true);
	}
	else {
		tree.repeats = 1;
	}
	
	read_node(tree.root);
	accept(";", true);
	
	tree.set_heights_and_ids_and_parents();
					
	return true;
}
/*
bool NewickParser::read_tree(Tree &tree) {
	accept("[", true);
	accept_integer(tree.repeats, true);
	accept("]", true);
	read_node(tree.root);
	accept(";", true);
	
	tree.set_heights_and_ids_and_parents();
					
	return true;
}
*/
bool NewickParser::is_symbol(char c) {
	return c == '(' || c == ')' || c == '[' || c == ']' || c == ';' ||
		   c == ',' || c == ':';
}

void NewickParser::get_tokens(const std::string& str) {
	tokens.clear();
	char c;
	size_t cur_col = 0, cur_line = 0, i=0, i_prev=0;
	while (i < str.size()) {
		
		assert(i >= i_prev);
		
		cur_col += i - i_prev;
		i_prev = i;
		c = str[i++];
			
		if (is_symbol(c)) {
			Token tok(std::string(&c, 1), SYMBOL, cur_line, cur_col);
			tokens.push_back(tok);
		} else if (isspace(c)) {
			if (c == '\n' || c == '\r') {
				cur_col = 0;
				cur_line += 1;
			}
		} else {
			Token tok(std::string(&c, 1), STRING, cur_line, cur_col);
			while (i < str.size()) {
				c = str[i++];
				if (isspace(c) || is_symbol(c)) {
					i -= 1;
					break;
				}
				tok.val += c;
			}
			tokens.push_back(tok);
		}
	}
}




