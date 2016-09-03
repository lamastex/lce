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
\brief DataObject declarations.
*/

#ifndef _INC_DATA_OBJECT_H
#define _INC_DATA_OBJECT_H


#include "summary_statistic.hpp"

#include <boost/smart_ptr.hpp>

namespace aabc {
	
	/** @brief An interface for types that can behave as data objects
	 for aabc.
	 
	 A DataObject must be able to return a mct::SummaryStatistic and provide
	 a string representation of itself. 
	 */
	 

	class DataObject {
		public:
		
			virtual ~DataObject() {};
			
			/*! \brief Get the mct::SummaryStatistic.*/
			virtual boost::shared_ptr < const mct::SummaryStatistic  > 
										getSummaryStatistic() const = 0;
			
			/*! \brief Get a string representation of this.*/
			virtual std::string toString() const = 0;
			
			
	};
	
} // end namespace aabc

#endif 

