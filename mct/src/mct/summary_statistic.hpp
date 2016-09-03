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
\brief SummaryStatistic declarations.
*/

#ifndef _INC_SUMMARYSTATISTIC_H
#define _INC_SUMMARYSTATISTIC_H



#include <vector>
#include <string>


namespace mct {
	
	

	/** @brief A class for an ordered collection of summary statistics
	 * from a data object.
	 * 
	 * SummaryStatistic at the moment is basically a data object, wrapping
	 * the collection of summary stats for one data object together.
	 * 
	 * The SummaryStatistic has no knowledge of what the values it holds
	 * represent or where they came from
	 * (means, variances, totals, medians, Tajimas's D's, 
	 * segregating sites, elephants ... it is all the same provided
	 * that each value can be represented as a double.
     */
	 

	class SummaryStatistic {
		public:
		
			/*! \brief No-args constructor.*/
			SummaryStatistic();
		
			/*! \brief Constructor.
			
			\param v A collection of values of summary statistics.
			*/
			explicit SummaryStatistic(const std::vector<double>& v);

			/*! \brief Constructor.
			
			\param a An array of values of summary statistics.
			\param n Number of elements in \a a.
			*/
			SummaryStatistic(const double a[],
									const size_t n);
									
			/*! \brief Copy constructor.
			
			\param other The SummaryStatistic to copy.
			*/
			SummaryStatistic(const SummaryStatistic& other);
			
			/*! \brief Equals operator.
			
			\param tmp The SummaryStatistic to set this equal to.
			*/
			SummaryStatistic& operator=(SummaryStatistic tmp);
			
			virtual ~SummaryStatistic();
			
			/*! \brief Get a vector representation of the values.*/
			const std::vector<double>& getValues() const;
			
			/*! \brief Get the element at index.*/
			const double& at(const size_t index) const;
			
			/*! \brief Add an element to the back of the collection.*/
			virtual void pushBack(double d);
			
			/*! \brief Get whether the collection is empty.*/
			virtual bool empty() const;	
			
			/*! \brief Get the number of values in the collection.*/
			virtual size_t size() const;	
			
			/*! \brief Get a fancy string representation of this.
			
			String is formatted as:
			([value1]\\t[value2]\\t...\\t[valuen])
			\param prec The precision to which to show the output. 
			\return A formatted string representation of this.*/
			std::string toString(int prec = 3) const;
			
			/*! \brief Get a plain string representation of this.
			 
			 Numbers are printed using fixed formatting.
			\todo Make this more flexible 
			
			 String is formatted as:
			 [value1]\\t[value2]\\t...\\t[valuen] 
			 This format is suitable as a line for a file
			\param prec The precision to which to show the output.
			\return A plain string representation of this.*/
			std::string plainString(int prec = 3) const;
			
			/*! \brief Get the Euclidean distance between this
			 and another %SummaryStatistic.
			 
			 If the elements of this are \f$ s_1, s_2, \ldots s_n \f$
			 and the elements of \a other are 
			 \f$ s^{\prime}_1, s^{\prime}_2, \ldots, s^{\prime}_n \f$
			 then the Euclidean distance between this and \a other is:
			 \f[ \sqrt{\sum_{i=1}^n \left(s_i - s^{\prime}_i\right)^2} \f]
			 
			 \param other The %SummmaryStatistic to calculate the distance
			 against.
			 \return The Euclidean distance between this and \a other.
			 \pre This and \a other should be the same size.  The behaviour
			 of this method is undefined if this condition is not met.
			*/
			virtual double getEuclideanDistance(
				const SummaryStatistic& other)
				const;
		
					
		protected:
			
				
		private:
		
			/*! \brief A container of values of statistics. */
			std::vector< double > values;
			
			
			
	};
	
	std::ostream& operator<< (std::ostream &out, 
								const SummaryStatistic &ss);
	
	
} // end namespace mct

#endif 

