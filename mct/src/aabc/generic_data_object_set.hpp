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
\brief GenericDataObjectSet declarations.
*/

#ifndef _INC_GENERIC_DATA_OBJECT_SET_H
#define _INC_GENERIC_DATA_OBJECT_SET_H

#include "data_object_set.hpp"
#include "data_object.hpp"


#include <boost/smart_ptr.hpp>


#include <vector>
//#include <string>


namespace aabc {
	
	/** @brief A class for an ordered collection of 
	 DataObject types.
	 * 
	 This class acts as a generic container for DataObject types. 
	 It implements the methods of the interface 
	 DataObjectSet but provides no other functionality.	 */
	 

	class GenericDataObjectSet : public aabc::DataObjectSet {
				
		public:
		
			/*! \brief Constructor.
			
			Constructs an empty structure in which to store 
			DataObject.
			*/
			GenericDataObjectSet();
			
			virtual ~GenericDataObjectSet();
			
			/*! \brief Get a pointer to an indexed 
			 DataObject in the set.
			
			\param index The index of the DataObject to get.
			\return A pointer to the DataObject at \a index. */
			virtual const boost::shared_ptr < const DataObject >&
					at(const size_t index) const;
			
			/*! \brief Add a DataObject to the 
			end of the collection.
			
			\param doptr The DataObject to add.*/
			virtual void add(
				const boost::shared_ptr < const DataObject >& doptr);
			
			/*! \brief Add the contents of another 
			GenericDataObjectSet to the end of collection.
			
			\param sptr The GenericDataObjectSet to add.*/
			virtual void add(
				const boost::shared_ptr < const GenericDataObjectSet >& sptr);
			
			/*! \brief Get the number of DataObject
			 in the collection.	*/
			virtual size_t size() const;
			
			/*! \brief Return whether this contains no data.	*/
			virtual bool empty() const;
			
			/*! \brief Get the mct::SummaryStatisticSet for the set.*/
			virtual boost::shared_ptr < mct::SummaryStatisticSet > 
						getSummaryStatisticSet() const;
			
			/*! \brief Get the mct::SummaryStatistic for the set.*/
			virtual boost::shared_ptr < const mct::SummaryStatistic >
											getSummaryStatistic() const;
			
			/*! \brief Get a string representation of this.*/
			virtual std::string toString() const;
			
						
		protected:
		
			explicit GenericDataObjectSet(const GenericDataObjectSet& other);
			
			GenericDataObjectSet& operator=(GenericDataObjectSet tmp);
					
			virtual const std::vector< boost::shared_ptr < const DataObject > >&
				_getContainer() const;
			
					
		private:
		
			/*! \brief The collection of pointers to 
			DataObject held.*/
			std::vector< boost::shared_ptr < const DataObject > > container;
			
			
	};
	
	// does not need to be friend
	std::ostream& operator<< (std::ostream &out, 
								const GenericDataObjectSet &m);

} // end namespace aabc

#endif 

