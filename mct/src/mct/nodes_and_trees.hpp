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

/*! \file     nodes_and_trees.hpp
\brief Declarations for structures to help parsing files into ARGs.
*/

#ifndef INC_NODES_AND_TREES_H
#define INC_NODES_AND_TREES_H

#include <string>
#include <vector>

namespace mct {

	/*! @brief A structure for a Node in a Tree
	 * 
	 * Nodes are used as part of a local representation of a marginal Tree. 
	*/
	struct Node {
		
		/*! \brief Name of this Node (leaf nodes only, identifies individual in sample). */
		std::string name; 
		
		 /*! \brief ID number for this Node (ultimate root is -1). */
		int id;
		
		/*! \brief Time from parent Node to this Node. */
		double length; 
		
		/*! \brief Time from time 0 to this Node. */
		double height; 
		
		/*! \brief id of the parent of this Node. */
		int parent_id; 
		
		/*! \brief Collection of the children of this Node. */
		std::vector<Node> children; 
		
		/*! \brief Print to standard output a representation of the Nodes
		 in the Tree rooted at this Node.
		 
		 \param depth The depth of the Node in the Tree (root depth = 0).*/   
		void print(int depth=0);
		
		/*! \brief Set the height of Nodes in the Tree rooted at this Node.
		 
		 Works recursively, setting heights of parents from heights 
		 already worked out for children.
		 
		 \return The height of the parent of this Node. */   
		double set_height(void);
		
		/*! @brief Set the parents for the Nodes descended from this Node.*/   
		void set_parent(void);
		
		/*! @brief Returns a collection of the Nodes in
					 the tree rooted at this node.*/   
		std::vector<Node*> get_node_list(void);
		
		/*! @brief Returns a string representation of the Node.*/   
		std::string toString(size_t level) const;
	};


	/*! @brief Comparison function for \link mct::Node Nodes\endlink,
		to order Nodes in a list.
	 
	 Gives the result of n1 < n2 by ordering Nodes by height.
	 
	 \param n1 The left hand side operand Node.
	 \param n2 The right hand side operand Node.
	 \return n1.height < n2.height (strict weak ordering). */   
	bool node_comp_height(Node *n1, Node *n2);

	/*! @brief A structure for a marginal Tree.
	 * 
	 * Trees are used as a local representation of a marginal tree.
	 * 
	 * This local representation is used as an intermediary between a 
	 * tree read in from a file and the creation of a sequence::arg.
	*/
	struct Tree {
		
		/*! \brief The root Node of the tree. */
		Node root; 
		
		/*! \brief The number of repeats of this tree in the ARG.
		 * 
		 * The number of consecutive sites to which the Tree applies. */ 
		size_t repeats; 
		
		/*! \brief Returns a string representation of the Tree.*/   
		std::string toString() const;
		
		/*! \brief Set heights, ids and parents for the 
		 \link mct::Node Nodes\endlink in the Tree.*/   
		void set_heights_and_ids_and_parents();
	};


} // end namespace mct

#endif /* INC_NODES_AND_TREES_H */



