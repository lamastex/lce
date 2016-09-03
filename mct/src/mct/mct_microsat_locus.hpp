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
\brief MCTMicrosatLocus declarations.
*/

#ifndef _INC_MCT_MICROSAT_LOCUS_H
#define _INC_MCT_MICROSAT_LOCUS_H

#include "site_seq.hpp"

#include <vector>
#include <string>


namespace mct {

	
	/** @brief A class for an ordered collection of 
	 microsatellite site sequence data for one or more sites 
	 in one locus for a sample.
	 
	 The locus can cover one or more sites, i.e. have data for more than
	 one microsat repeats count, but would usually only have one.
	  
	 */
	 

	class MCTMicrosatLocus { 
				
		public:
		
			/*! Iterator into the ordered collection of site_seq.*/
			typedef 
			std::vector< site_seq >::const_iterator
					const_iterator;
				
			/*! \brief Constructor.
			
			Constructs an empty structure.
			*/
			MCTMicrosatLocus();
			
			/*! \brief Constructor.
			
			Constructs a MCTMicrosatLocus from a collection of
			\link mct::site_seq site_seqs\endlink.
			
			\param seqs the collection to copy into this.
			*/
			explicit MCTMicrosatLocus( 
						const std::vector < site_seq >& seqs);
			
			virtual ~MCTMicrosatLocus();
			
			/*! \brief Get a copy of an indexed site_seq.
			
			\param index The index of the %site_seq to get.
			\return A const_reference to the site_seq at \a index. */
			const site_seq& at(const size_t index) const;
			
			/*! \brief Add a site_seq to the collection.
			
			\param seq The %site_seq to add.*/
			void add(const site_seq& seq);
			
			/*! \brief Get the number of 
			 \link mct::site_seq site_seqs\endlink
			 in the collection.	*/
			size_t size() const;
			
			/*! \brief Return whether this contains no data.	*/
			bool empty() const;
			
			/*! \brief Get the number of sites covered by each
			 \link mct::site_seq site_seqs\endlink
			 in the collection.	*/
			size_t nsites() const;
			
			/*! \brief Get a fancy string representation of this.*/
			virtual std::string toString() const;
			
			/*! \brief Get a plain string representation of this.*/
			virtual std::string plainString() const;
			
			/*! \brief Return a const iterator to the first 
			element in the collection.*/
			const_iterator begin() const;
			
			/*! \brief Return a const iterator to the
			'just past the end' element in the set.
			
			\note The idiom for iterating through a 
			%MCTMicrosatLocus mtab is 
			for (MCTMicrosatLocus::const_iterator it = mtab.begin();
						it < mtab.end();
						++it) {...}.
			*/
			const_iterator end() const;
									
		protected:
		
			explicit MCTMicrosatLocus(
								const MCTMicrosatLocus& other);
			
			MCTMicrosatLocus& operator=(MCTMicrosatLocus tmp);
			
			
				
		private:
		
			/*! \brief The collection of microsat repeats data
			 held.*/
			std::vector< site_seq > container;
			
			size_t _nsites;
			
	};
	
	// does not need to be friend
	std::ostream& operator<< (std::ostream &out, 
								const MCTMicrosatLocus &m);

} 

#endif 

