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
\brief IndividualSiteSequence definitions.
*/

#include "ind_site_sequence.hpp"

using namespace mct;

IndividualSiteSequence::IndividualSiteSequence()
				: pop(0) {}

IndividualSiteSequence::IndividualSiteSequence(size_t pl,
												const site_seq& s)
				: pop(pl), seq(s) {}
		
IndividualSiteSequence::~IndividualSiteSequence()
{
	
}

size_t IndividualSiteSequence::getPop() const
{
	return pop;
}

		
site_seq IndividualSiteSequence::getSeq() const
{
	return seq;
}
