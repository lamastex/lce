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
\brief Definition for a sequence of genetic data over a number of sites.
*/

#ifndef _INC_SITE_SEQ_H
#define _INC_SITE_SEQ_H


#include <vector>

namespace mct {

	/*! \brief A simple structure holding a coded sequence of sites on
	a chromosone for one individual. The state at each site is given by an integer.
	The integer could represent a nucleotide state or a 
	microsatellite state:  the %site_seq itself does not 
	know what the codes represent.
	* 
	\internal
	The interpretation of the 
	integer depends on the Model applying to that site, 
	which the %site_sequence is unaware of (see 
	IndividualSiteSequence and MCTSampleSequence
	to see how %site_seq build into labelled sample
	sequences in nucleotide or microsatellite format. */
	typedef std::vector< int > site_seq;

}

#endif /* _INC_SITE_SEQ_H */
