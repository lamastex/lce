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
\brief MixtureTest2D declarations.
*/

#ifndef _INC_MIXTURETEST_2D_H
#define _INC_MIXTURETEST_2D_H

#include "data_object.hpp"
#include "summary_statistic.hpp"

#include <gsl/gsl_vector.h>

#include <boost/smart_ptr.hpp>

//#include <vector>
#include <string>


namespace aabc {
	
	/** @brief A class for testing a pop gen model with 2-d data.
	 */
	 

	class MixtureTest2D : public DataObject {
		public:
		
			/*! \brief Constructor.*/
			explicit MixtureTest2D(double v0, double v1);
			
			virtual ~MixtureTest2D();
			
			/*! \brief Get the SummaryStatistic.*/
			virtual boost::shared_ptr < const mct::SummaryStatistic  > 
										getSummaryStatistic() const;
			
			/*! \brief Get a string representation of this.*/
			virtual std::string toString() const;
			
			// /*! \brief Get the value.*/
			//virtual const gsl_vector* getValue() const;
			
						
		protected:
		
			MixtureTest2D();
			
			explicit MixtureTest2D(const MixtureTest2D& other);
			
			MixtureTest2D& operator=(MixtureTest2D tmp);
						
			
		private:
		
			gsl_vector* v;
			
			
	};

} // end namespace aabc

#endif 

