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
\brief Definitions for a class to read ARGs from a file.
*/

#include "arg_reader.hpp"
#include "nodes_and_trees.hpp"

#include "libseq.hpp" // only need for debugging output of marginal tree

#include <sstream>
#include <algorithm>
#include <stdexcept>

//#define MYDEBUG
#if defined (MYDEBUG)
	#include <iostream>
	
#endif

using namespace mct;

ARGReader::ARGReader() 
{
	boost::shared_ptr < ARGParser > p(new NewickParser);
	parserPtr = p;
}

ARGReader::~ARGReader()
{}  // parserPtr will take care of itself


boost::shared_ptr < Sequence::arg > ARGReader::makeTreefileIntoArg(
							const std::string& stringARG,
							const size_t nsites)
{
	boost::shared_ptr < Sequence::arg > hist(new Sequence::arg );
	
	tree_list.clear();
	size_t part_start = 0;
	
	// parse the string into a vector of Trees.
	tree_list = parserPtr->parse(stringARG);
	
	/*if there is only one tree in the tree list it may not have
	 * the correct number of repeats because eg MS does not
	 * put repeats (number of sites covered) unless there
	 * is recombination 
	*/
	
	if (tree_list.size() == 1) {
		tree_list.begin()->repeats = nsites;
	}
	
	for (size_t i=0; i<tree_list.size(); ++i) {
		Tree *tree = &tree_list[i];
		std::vector<Node*> node_list;
		std::vector<Sequence::node> seq_node_list;
		
		node_list = tree->root.get_node_list();
		
		//nodes must go in in height order
		std::sort(node_list.begin(), node_list.end(), node_comp_height);
		
		// converting from a Node list to a Sequence::node list
		for (size_t j=0; j<node_list.size(); ++j) {
			
			#ifdef MYDEBUG
				std::cout << "Node " << j << " with id " << node_list[j]->id << " and node_list[j].length is " << node_list[j]->length 
					<< " and height is " << node_list[j]->height<< std::endl;
			#endif
			
			Node *n = node_list[j];
			//int parent_id = n->parent ? n->parent->id : -1;
			seq_node_list.push_back(Sequence::node(n->height, n->parent_id));
		}
	 
		
		int n_nodes = node_list.size();
		int n_sample = (n_nodes + 1)/2;
		
		//debug
		//std::cout << "i = " << i << " and tree is " << std::endl;
		//std::cout << tree->toString();
		//std::cout << "n_nodes = " << n_nodes << " and n_sample = " << n_sample << std::endl;
		//std::cout << std::endl;
		
		// creating the marginal tree and adding it to the ARG
		Sequence::marginal marg(part_start, n_sample, n_nodes, seq_node_list);
		
		#ifdef MYDEBUG
			std::cout << "Marginal is\n" <<LibSeq::margToString(marg) << std::endl;
		#endif
		
		
		hist->push_back(marg);
		
		part_start += tree->repeats; // final part_start is number of sites in total
		
	}
	if (nsites != part_start) {
		
			#ifdef MYDEBUG
				std::cout << "nsites is " << nsites <<  std::endl;
				std::cout << "part_start is "  << static_cast <size_t> (part_start) << std::endl;
			#endif
			throw std::invalid_argument(
			std::string("ARGReader::makeTreefileIntoArg(const std::string&, const size_t)")
			+ std::string(":\nnsites != nsites in tree from file") );
	}
	
	return hist;
}




