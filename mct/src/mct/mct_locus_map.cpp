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
\brief MCTLocusMap definitions.
*/

#include "mct_locus_map.hpp"

#include <string>
#include <stdexcept>

using namespace mct;

MCTLocusMap::MCTLocusMap() {}

MCTLocusMap::~MCTLocusMap() 
{
	//std::cout << "locus map destructor called " << std::endl;
}

MCTLocusMap::Block::Block(size_t _from, 
								size_t _to, 
								DATATYPE _type)
		: from(_from), to(_to), type(_type) {}

void MCTLocusMap::add(size_t _from, size_t _to, DATATYPE _type)
{
	if (!modelmap.empty() && (_from != (modelmap.back())->to + 1)) {
		throw std::invalid_argument(
			std::string("MCTLocusMap::add(size_t, size_t, DATATYPE)")
			+ std::string(":\nnew block inconsistent with existing blocks") );
	}
			
	if ( (modelmap.empty()) || 
		(!modelmap.empty() && _type != (modelmap.back())->type) ){
		// need to add new block
		modelmap.push_back(
			boost::shared_ptr < MCTLocusMap::Block >(
				new MCTLocusMap::Block(_from, _to, _type) ) );
	}
	
	else  { 
		// same type as last block, just extend last block
		(modelmap.back())->to = _to;
	}
		
}
			
const boost::shared_ptr < MCTLocusMap::Block >& MCTLocusMap::at(
						size_t index) const
{
	return modelmap.at(index);
}

			
size_t MCTLocusMap::size() const
{
	return modelmap.size();
}
			
bool MCTLocusMap::empty() const
{
	return modelmap.empty();
}			


