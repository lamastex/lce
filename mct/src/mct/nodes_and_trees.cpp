/* MCT - Markhov Chains on Trees.

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

/*! \file     nodes_and_trees.cpp
\brief Definitions for structures to help parsing files into ARGs.
*/

#include "nodes_and_trees.hpp"

#include <iostream>
#include <sstream>
#include <algorithm> // for sort

//Node

using namespace mct;    
    
void Node::print(int depth) {
	for (int i=0; i<depth; ++i)
		std::cout << "   ";
	std::cout << "`-- " << length << " " << name << "\n";
	for (size_t i=0; i<children.size(); ++i)
		children[i].print(depth + 1);
}

/*
 * Note that using this method, the height of a parent is effectively calculated
 * from the length and height of the last of his children, in order, to be processed
 * This means that rounding in the lengths given in the file being processed
 * will affect heights, especially since the last child processed
 * may have the most layers below him and hence most potential 
 * compounding of roundings.  
 * 
 * eg, looking at the results of the process that this call is all part of,
 * (ARGFactory: after ms call - the tree we got back is
	((1:0.034,(2:0.014,5:0.014):0.021):0.577,(7:0.313,(3:0.273,(8:0.087,(4:0.021,6:0.021):0.065):0.186):0.040):0.298);

	The tree is
	((1:0.035,(2:0.014,5:0.014):0.021):0.577,(7:0.313,(3:0.273,(8:0.087,(4:0.021,6:0.021):0.066):0.186):0.04):0.299);
 
 In this case the bottom branch with length 0.034 did not actually determine
 * the height of his parent:  the sum of 0.021 and 0.014 did and hence it was 0.035
 * and this in turn became his length when that was recalculated from his height (0) 
 * and his parent's height (0.035)
 */ 

double Node::set_height(void) {
	
	//std::cout << "In set height, length is " << this->length << std::endl;
	
	this->height = 0.0;
	for (size_t i=0; i<children.size(); ++i) {
		//std::cout << "doing child i = " << i << std::endl;
		double child_height = children[i].set_height();
		//std::cout << "got back = " << child_height << std::endl;
		if (child_height > this->height)
			this->height = child_height;
		//std::cout << "set height = " << this->height << std::endl;
	}
	//std::cout << "returning = " << this->height + this->length << std::endl;
	return this->height + this->length;
}

void Node::set_parent(void) {
	this->parent_id = -1;
	for (size_t i=0; i<children.size(); ++i) {
		children[i].set_parent();
		children[i].parent_id = this->id;
	}
}

std::vector<Node*> Node::get_node_list(void) {
	std::vector<Node*> nodes(1, this);
	for (size_t i=0; i<children.size(); ++i) {
		std::vector<Node*> child_nodes = children[i].get_node_list();
		nodes.insert(nodes.end(), child_nodes.begin(), child_nodes.end());
	}
	
	return nodes;
}

	
std::string Node::toString(size_t level) const
{
	std::ostringstream stm;
	for (size_t i = 0; i < level; ++i) {
		stm << '\t';
	}
	stm << "Name: " << name;
	stm << ", Id: " << id;
	stm << ", Length: " << length;
	stm << ", Height: " << height;
	stm << ", Parent_id: " << parent_id;
	stm << '\n';
	std::string retStr = stm.str();
	for (std::vector<Node>::const_iterator cit = children.begin();
			cit < children.end();
			++cit) {
		retStr += cit->toString(level+1);
	}
	return retStr;
}
// end of Node definitions



bool mct::node_comp_height(Node *n1, Node *n2) {
    if (n1->height == n2->height) return n1->name < n2->name;
    return n1->height < n2->height;
}


// Tree
    
std::string Tree::toString() const
{
	size_t starting_level = 0;
	return root.toString(starting_level) + "\n";
}

void Tree::set_heights_and_ids_and_parents()
{
	root.set_height();
	std::vector<Node*> node_list = root.get_node_list();
	std::sort(node_list.begin(), node_list.end(), node_comp_height);
	// setting id's
	for (size_t j=0; j<node_list.size(); ++j)
		node_list[j]->id = j;
	root.set_parent();	

}





