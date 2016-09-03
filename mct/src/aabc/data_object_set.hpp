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
\brief DataObjectSet declarations.
*/

#ifndef _INC_DATA_OBJECT_SET_H
#define _INC_DATA_OBJECT_SET_H

#include "data_object.hpp"
#include "summary_statistic_set.hpp"
#include "summary_statistic.hpp"

#include <boost/smart_ptr.hpp>

//#include <vector>
#include <string>

namespace aabc {
	
	/** @brief A interface for a type for an ordered collection of 
	 DataObject types.
	  
	 There are no checks on what is added to the collection.
	 
	 A DataObjectSet is a \link aabc::DataObject DataObject\endlink
	type (i.e uses the composite design pattern, allowing 
	clients to treat the composition or set of DataObjects 
	as a DataObject itself). */
	 
	class DataObjectSet : public DataObject {
				
		public:
		
			~DataObjectSet() {}
		
			/*! \brief Get the number of 
			 \link DataObject DataObjects\endlink
			 in the collection.	*/
			virtual size_t size() const = 0;
			
			/*! \brief Return whether this contains no data.	*/
			virtual bool empty() const = 0;
			
			/*! \brief Get the mct::SummaryStatistic for the set.*/
			virtual boost::shared_ptr < const mct::SummaryStatistic >
						getSummaryStatistic() const = 0;
			
			/*! \brief Get the mct::SummaryStatisticSet for the set.*/
			virtual boost::shared_ptr < mct::SummaryStatisticSet > 
						getSummaryStatisticSet() const = 0;
			
			/*! \brief Get a mct::SummaryStatistic containing 
			the means of the 
			\link SummaryStatistic SummaryStatistics\endlink
			from the collection of
			\link DataObject DataObjects\endlink.*/				
			virtual boost::shared_ptr < const mct::SummaryStatistic >
						getSummaryStatisticMeans() const;
						
			/*! \brief Get a mct::SummaryStatistic containing 
			the sample variances of the 
			\link SummaryStatistic SummaryStatistics\endlink
			from the collection of
			\link DataObject DataObjects\endlink.*/					
			virtual boost::shared_ptr < const mct::SummaryStatistic  > 
						getSummaryStatisticVars() const;
			
			
									
		protected:
		
					
		private:
		
			
			
	};
	
	// does not need to be friend
	std::ostream& operator<< (std::ostream &out, 
								const DataObjectSet &d);

} // end namespace aabc

#endif 

