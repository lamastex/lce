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
\brief Definitions for ARGmsConverter class to convert ms args to ours.
*/

#include "arg_ms_converter.hpp"

#include "libseq.hpp" // for node height comparison


#include <mscplusplus/arg_factory.hpp>
#include <mscplusplus/recomb_params.hpp>
#include <mscplusplus/prng_gsl_glibc2.hpp>


#include <numeric> // for accumulate
#include <algorithm> // for stable sort
#include <vector>
#include <cassert>


//#define MYDEBUG_OUT
//#define MYDEBUG

#if defined (MYDEBUG) || defined (MYDEBUG_OUT)
	#include <iostream>
	
#endif

using namespace mct;

boost::shared_ptr < Sequence::arg > ARGmsConverter::makeArg(
							size_t nsites,
							double rho,
							const boost::shared_ptr< const hudson_ms::PopulationStructure > pop,
							int seed)
{
	boost::shared_ptr < Sequence::arg > hist(new Sequence::arg );
	
	//size_t nsam = accumulate(popnsams.begin(), popnsams.end(), 0);
	
	boost::shared_ptr < hudson_ms::HudsonMSPRNGen > 
				rr ( new hudson_ms::HudsonMSPRNGenGLS_glibc2(seed) );

	hudson_ms::HudsonMSRecombParams recombParams(rho, nsites);
	
	hudson_ms::HudsonMSArgFactory factory(recombParams,
											pop,
											rr);
	
	boost::shared_ptr < hudson_ms::HudsonMSArg > arg_ptr = factory.makeArg();

	
	#ifdef MYDEBUG
		std::cout << "nsegs = " << arg_ptr->nsegs() << std::endl;
		
	#endif
	
	/*
	if (tree_list.size() == 1) {
		tree_list.begin()->repeats = nsites;
	}
	*/
	
	int segmnt=0; 
	size_t totlen = 0;
	int ns = arg_ptr->nsegs();
	for(int k=0; k < ns; ++k) { 
		
		#ifdef MYDEBUG
			std::cout << "k = " << k << " segmnt = " << segmnt << std::endl;
			
		#endif
		
		//int end = ( k<nsegs-1 ? seglist[seglist[segmnt].next].beg -1 : nsites-1 );
		//int start = seglist[segmnt].beg ;
		
		int len = arg_ptr->at(segmnt).len();
		totlen += len;
		
		//hudson_ms Tree *tree = &tree_list[i];
		
		std::vector<Sequence::node> seq_node_list;
	
		//(seglist[segmnt].ptree + 2*nsam-2) -> abv = -1; // root

		for (hudson_ms::HudsonMSArg::const_tree_iterator it = arg_ptr->at(segmnt).tree_begin();
			it < arg_ptr->at(segmnt).tree_end();
			++it) {
			
			#ifdef MYDEBUG
				std::cout << "Transferring node with parent" << it->abv 
					<< " and height is " << it->time << std::endl;
			#endif
			
			//int parent_id = n->parent ? n->parent->id : -1;
			seq_node_list.push_back(Sequence::node(it->time, it->abv));
		}
	 
		
		//nodes must go in in height order
		std::stable_sort(seq_node_list.begin(), 
						seq_node_list.end(), 
						seq_node_comp_height);
		
		
		int n_nodes = seq_node_list.size();
		int n_sample = (n_nodes + 1)/2;
		
		// creating the marginal tree and adding it to the arg_ptr
		Sequence::marginal marg(arg_ptr->at(segmnt).beg, n_sample, n_nodes, seq_node_list);
		
		#ifdef MYDEBUG
			std::cout << "Marginal is\n" <<LibSeq::margToString(marg) << std::endl;
		#endif
		
		
		hist->push_back(marg);
		
		segmnt = arg_ptr->at(segmnt).next; 
		
	}
	
	#ifdef MYDEBUG
		std::cout << "nsites is " << nsites << std::endl;
		std::cout << "totlen is " << totlen << std::endl;
	#endif
	
	if (nsites > 1) assert(nsites == totlen);
	else assert(totlen == 2);

	#ifdef MYDEBUG_OUT
		std::cout << "arg_ptr is " << std::endl;
		std::cout << LibSeq::argNewickString(*hist) << std::endl;
	#endif
	
	return hist;
}
	

	

