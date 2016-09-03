/* MCT - Markov Chains on Trees.

   Copyright (C) 2011 Jenny Harlow

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
\brief IndividualSiteSequence declarations.
*/

#ifndef _INC_INDIVIDUAL_SITE_SEQUENCE_H
#define _INC_INDIVIDUAL_SITE_SEQUENCE_H

#include "site_seq.hpp"

namespace mct {

	/** @brief A struct for an individual site sequence.
	 * 
	 * The struct contains:
	 * <ul>
	 * <li>The index of the population the individual belongs to</li>
	 * <li>The site_seq for that individual</li>
	 * </ul>
	 * 
	 */
	
	
	struct IndividualSiteSequence {
		
		/*! \brief No-args constructor.
		 
		 \a pop will default to 0.
		 */
		IndividualSiteSequence();
		
		/*! \brief Constructor.
		 
		 \param s The site sequence.
		 \param pl The population label.
		 */
		IndividualSiteSequence(size_t pl, const site_seq& s);
		
		~IndividualSiteSequence();
		
		/*! \brief Get the population index. */
		size_t getPop() const;
		
		/*! \brief Get the site_sequence. */
		site_seq getSeq() const;
		
		/*! \brief The population index. */
		size_t pop;
		
		/*! \brief The site_sequence. */
		site_seq seq;
		
		
		
	};
	
	
} // end namespace mct

#endif 

