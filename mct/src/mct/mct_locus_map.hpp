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
\brief MCTLocusMap declarations.
*/

#ifndef _INC_MCTLOCUSMAP_H
#define _INC_MCTLOCUSMAP_H

#include "mct_model_data_types.hpp"

#include <boost/smart_ptr.hpp>

#include <vector>



namespace mct {
	
	/** @brief A class for holding a map of the type of data 
	 in locus represented in a MCTSequence.
	 
	 Data can be DNA data or microsatellite data.  The basic sequence
	 is just integers and needs to know how to code itself for
	 analysis.
	  
	 */
	 

	class MCTLocusMap {
		
		public:
			
			virtual ~MCTLocusMap();
			
		protected:
		
			friend class MCTSampleSequence;
			friend class ModelSet;
		
			/*! \brief a description of a block of sites in the sequence.
			
			A block is a contiguous section of sites within a 
			the sequence
			with the same datatype,
			described by the index of the site where
			the block starts from, the index where the block goes to,
			and the type of data in the block.
			
			Site indexing runs from 0 to nsites-1 where nsites 
			is the number of sites in the sequence.
			*/
			struct Block {
				
				Block(size_t _from, size_t _to, DATATYPE _type);
				
				size_t from;
				size_t to;
				DATATYPE type;
			};
		
			/*! \brief Constructor.*/
			MCTLocusMap();
			
			
			/*! \brief Add new information to the map.
			 * 
			\pre The information to be added must be consistent with
			the existing map.  Specifically, either the map must be 
			empty and \a _from == 0, or _from == 1 + \c last_to where 
			\c last_to is the the site index of the end of the last
			block already in the map.
			
			If the new information describes a section of the sequence
			which has the same datatype as the last existing block
			in the map, then that block is extended to include
			the new information. If the new information describes 
			a section of the sequence
			which has a different datatype to the last existing block
			in the map, then a new block is added to the map.
			
			\param _from The index of the site within the sequence
			where the new information starts from.
			\param _to The index of the site within the sequence
			where the new information goes to.
			\param _type The datatype of the sites 
			\a _from-_to inclusive.
			
			\exception MCTException "new information inconsistent with
			existing map."*/
			void add(size_t _from, size_t _to, DATATYPE _type);
			
			/*! \brief Get the index<sup>th</sup> block.
			\param index The index of the Block to get. */
			const boost::shared_ptr < Block >& at(size_t index)
																const;
			
			/*! \brief Get the number of blocks of different kinds of data in 
			the map.*/
			size_t size() const;
			
			/*! \brief Get if the map is empty.*/
			bool empty() const;
			
			explicit MCTLocusMap(const MCTLocusMap& other);
			
			MCTLocusMap& operator=(MCTLocusMap tmp);
						
			
			
		private:
		
			std::vector < boost::shared_ptr <Block> > modelmap;
			
			
	};
	
	
} // end namespace mct

#endif 

