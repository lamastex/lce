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

/*! \file     partition.cpp
\brief Definitions for a class holding data for one marginal tree
* and additional site information.
*/

#include "partition.hpp"

#include <algorithm>
#include <sstream>

//#define MYDEBUG
#if defined (MYDEBUG)
	#include <iostream>
	
#endif

using namespace mct;

/* Get outstream representation of a Partition. */
std::ostream& mct::operator<<(std::ostream & s, const Partition &p)
{
    s << "Partition(" << p.site_start << ":" << p.site_end << " "
      << p.num_sites << ')';
    return s;
}


/* Get outstream representation of a list of site_seq. */
std::ostream & mct::operator<<(std::ostream & s, 
						const std::vector<mct::site_seq> &p)
{
    for (unsigned i=0; i<p.size(); ++i) {
        mct::site_seq p2 = p[i];
        s << "(";
        for (unsigned j=0; j<p2.size(); ++j) {
	    if (j != 0) s << ".";
	    s.fill('0');
	    s.width(2);
            s << p2[j];
	}
        s << ")\n";
    }
    return s;
}



std::string Partition::toString() const
{
	std::ostringstream stm;
	stm << (*this) << ", seqs are\n" << this->seqs;
	return stm.str();
}

std::string Partition::newick_string() const
{
	std::ostringstream stm;
	stm << "[" << num_sites << "] "
		 << Sequence::newick_stream_marginal_tree(*tree);
		stm << std::endl;

	return stm.str();
}

/*
 * Evolve a tree, and write sequence data to the partition.
 */
void Partition::evolve_tree(
			const boost::shared_ptr < ModelSet >& models)
{
	
	initialise_sequences();
	
	//Sequence::marginal *margtree = this->tree;

    /* keeps track of the site wrt the whole DNA
     * goes through each site in the partition (tree) */
    int site_id = site_start;
    /* run through each site covered by this partition (tree) */
    for (int i=0; i<num_sites; ++i, ++site_id)
    {
		
	    /* run through each node in the tree (top to bottom). */
        for (int i_node = tree->tree.size() - 1; i_node >= 0; i_node--)
        {
			#ifdef MYDEBUG
				std::cout << "i_node = " << i_node << std::endl;
			#endif
			
			Sequence::node *node = &tree->tree[i_node];

	        if (node->abv == -1)
            {
	            /* assign the ancestral node */
				(seqs[i_node])[i] = models->ancestor_at(site_id);
				
				#ifdef MYDEBUG
					std::cout << "Ancestor assignment is " << ((seqs[i_node])[i]) << std::endl;
				#endif
            
            }
            else
            {
                /* the inner loop of the program.
                 * apply markov models to individual sites.
                 */
                Sequence::node *parent_node = &tree->tree[node->abv];

                double length = parent_node->time - node->time;
                
                #ifdef MYDEBUG
					std::cout << "My above is  " << node->abv << " and length is " << length << std::endl;
				#endif
		
                (seqs[i_node])[i] = models->step_at(site_id, 
													(seqs[node->abv])[i],
													 length);
            
            }
        }
    }
}

void Partition::initialise_sequences()
{
	//swap temporary and seqs around - efficient way to clear and resize
	(std::vector< mct::site_seq >( tree->tree.size(), mct::site_seq(num_sites) )).swap(seqs);
	
	/* note that tree->nnodes is NOT a reliable source for nodes in the marg.tree
	 * see libsequence documentation for Sequence::marginal::nnodes */
	
	
}

