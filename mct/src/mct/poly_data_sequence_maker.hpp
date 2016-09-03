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
\brief PolyDataSequenceMaker declarations.
*/

#ifndef _INC_POLY_SITE_SEQUENCE_MAKER_H
#define _INC_POLY_SITE_SEQUENCE_MAKER_H


#include "site_seq.hpp"

#include <vector>

namespace mct {
	
	
	/** @brief A static class for 
	 reducing a sample of site_seq data into a sample containing just
	 the polymorphic sites in the sequences .
	
	*/
	 
	

	class PolyDataSequenceMaker {
		public:
			
			/*! \brief Turn a sample of site_seq data into a sample
			containing just the polymorphic sites. 
			
			A 'polymorphic site' means any site position where there
			is a difference between two or more members of the sample
			in the element at that position. 
			
			If \a data contains no site_seqs, the returned sample will
			be empty. Otherwise the returned sample will contain
			as many site_seqs as there are site_seqs in the \a data
			passed in.  If each site_seq in the sample \a data has
			length 0, then \a data will be unchanged after the operation.
			If each site_seq in the sample \a data is identical to 
			every other site_seq in the sample \a data then there will
			be no polymorphic sites and after the operation \a data 
			will contain empty site_seqs (as many as there were site_seq
			in \a data passed in).  
			* 
			For example, if the sample passed in \a data contains
			three site_seqs and can be described as follows:
			< < 1 2 2 > < 1 1 2 > < 1 2 1> >
			* 
			then the second and third sites are polymorphic 
			(but every element at the first position in the site_seqs
			in the sample is the same so the first position is not
			polymorphic) and after the operation \a data will be  
			< < 2 2 > < 1 2 > < 2 1> >
			 
			\param data A reference to the sample to be reduced
			to its polymorphic sites.
			\post The sample \a data is reduced to 
			its polymorphic site constituents.
			\pre each site_seq in the sample \a data passed in is the 
			same length (which can be 0).  The results are not defined 
			if any site_seq in the sample \a data has a different length
			to any other site_seq in \a data. 		   		 */
			static std::vector< site_seq >& reduceToPolyData(
						std::vector< site_seq >& data);
			
					
		private:
		
					
			PolyDataSequenceMaker();
			
			explicit PolyDataSequenceMaker(
			 						const PolyDataSequenceMaker& other);
			
			PolyDataSequenceMaker& operator=(PolyDataSequenceMaker tmp);
			
	};
	
	
} // end namespace mct

#endif 

