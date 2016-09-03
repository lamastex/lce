/* MCT - Markov Chains on Trees.

   Copyright (C) 2011, 2012 Jenny Harlow

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
\brief PolyDataSequenceMaker definitions.
*/


#include "poly_data_sequence_maker.hpp"

//#define MYDEBUG
//#define MYDEBUGDETAIL

#if defined (MYDEBUG) || defined (MYDEBUGDETAIL)
	#include <iostream>
	#include <iterator>
	
#endif

using namespace mct;


std::vector< site_seq >& PolyDataSequenceMaker::reduceToPolyData(
						std::vector< site_seq >& data)
{
	#ifdef MYDEBUG
	{
		std::cout << "In PolyDataSequenceMaker::reduceToPolyData, data = " << std::endl;
		for (std::vector < site_seq >::const_iterator it = data.begin();
				it < data.end();
				++it) {
			std::ostream_iterator<int> out_it (std::cout," ");
			copy ( it->begin(), it->end(), out_it );
			std::cout << std::endl;
		}
	}
	#endif
	if (!data.empty()) {
		size_t n = data.size();
		
		if (n > 1) {
			
			size_t nsites = data.front().size();
			
			if (nsites > 0) {
				
				std::vector< site_seq > tmp(n, site_seq());
		
				/* get the 'polymorphic' alleles, 
				 * ie where there is a difference
				 * size of difference does not matter */
				for (size_t s = 0; s < nsites; ++s) {
					bool poly = false;
					
					#ifdef MYDEBUGDETAIL
						std::cout << "s = " << s << std::endl;
					#endif
					
					/* stop the loops when we have found one position not the same */
					for (size_t i = 0; (i < n - 1) && !poly; ++i) {
						#ifdef MYDEBUGDETAIL
							std::cout << "i = " << i << std::endl;
						#endif
					
						for (size_t j = i + 1; (j < n) && !poly; ++j) {
							#ifdef MYDEBUGDETAIL
								std::cout << "j = " << j << std::endl;
							#endif
							
							if ( data[i][s] != data[j][s] ) {
								poly = true;
								
								#ifdef MYDEBUGDETAIL
									std::cout << "data[i][s] " << data[i][s] << "!=" 
											<< " data[j][s] = " << data[j][s] << std::endl;
								#endif
								
								// and fill in the relevant position into tmp
								for (size_t k = 0; k < n; ++k) {
									tmp[k].push_back(data[k][s]);
								}
								
							}
							
						}
					}
				}
				tmp.swap(data);
			}
		}
		// vector of empty vectors if there are no sites or if nsam = 1
	}
	// empty vector if no samples at all
	#ifdef MYDEBUG
	{
		std::cout << "end of PolyDataSequenceMaker::reduceToPolyData, data now = " << std::endl;
		for (std::vector < site_seq >::const_iterator it = data.begin();
				it < data.end();
				++it) {
			std::ostream_iterator<int> out_it (std::cout," ");
			copy ( it->begin(), it->end(), out_it );
			std::cout << std::endl;
		}
	}
	#endif
	return data; // return by reference
}

