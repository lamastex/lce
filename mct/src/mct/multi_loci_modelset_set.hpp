/* MCT - Markov Chains on Trees.

   Copyright (C) 2012 Jenny Harlow

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
\brief MultiLociModelSetSet declarations.
*/

#ifndef _INC_MULTI_LOCI_MODELSET_SET_H
#define _INC_MULTI_LOCI_MODELSET_SET_H

#include "modelset.hpp"

#include <boost/smart_ptr.hpp>

#include <vector>


namespace mct {
	
	/** @brief A class for an ordered collection of 
	 * \link ModelSet ModelSets\endlink.
	  
	There are no checks on what is added to the collection.
	 
	 */
	 

	class MultiLociModelSetSet {
				
		public:
		
			/*! \brief Constructor.
			
			Constructs an empty structure in which to store 
			\link ModelSet ModelSets\endlink.
			*/
			MultiLociModelSetSet();
			
			/*! \brief Constructor.
			
			Constructs a structure storing \a nloci copies of \a ms.
			
			\param ms A pointer to the ModelSet to use for all \a nloci loci.
			\param nloci The number of loci for which holds represents models.
			*/
			MultiLociModelSetSet(
					const boost::shared_ptr < ModelSet >& ms,
					size_t nloci);
			
			virtual ~MultiLociModelSetSet();
			
			/*! \brief Get a pointer to an indexed 
			 ModelSet in the set.
			
			\param index The index of the ModelSet to get.
			\return A pointer to the ModelSet at \a index. */
			const boost::shared_ptr < ModelSet >&
					at(const size_t index) const;
			
			/*! \brief Add a ModelSet to the 
			end of the collection.
			
			\param mlptr The ModelSet to add.*/
			void add(
				const boost::shared_ptr < ModelSet >& ms);
			
			
			/*! \brief Get the number of 
			 \link ModelSet ModelSets\endlink.
			 in the collection.	*/
			size_t size() const;
			
			/*! \brief Return whether this contains no data.	*/
			bool empty() const;
			
							
			/*! \brief Get a string representation of this.*/
			virtual std::string toString() const;
			
			/*! \brief Set the thetas for each  
			ModelSet in the collection.
			
			\param new_thetas The new thetas to set the models to. 
			\pre thetas.size() == this.size().*/
			virtual void resetThetas(
				const std::vector < std::vector < double > >& new_thetas);

			/*! \brief Set the thetas for every
			ModelSet in the collection to the same value.
			
			\param new_theta The new theta to set all the models to. */
			virtual void resetThetas(double new_theta);

			
						
		protected:
		
			explicit MultiLociModelSetSet(const MultiLociModelSetSet& other);
			
			MultiLociModelSetSet& operator=(MultiLociModelSetSet tmp);
					
						
					
		private:
		
			/*! \brief The collection of pointers to 
			(non const) ModelSets held.*/
			std::vector< boost::shared_ptr < ModelSet > > container;
			
			
	};
	
	// does not need to be friend
	std::ostream& operator<< (std::ostream &out, 
								const MultiLociModelSetSet &m);

} // end namespace mct

#endif 

