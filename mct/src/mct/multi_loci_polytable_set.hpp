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
\brief MultiLociPolyTableSet declarations.
*/

#ifndef _INC_MULTI_LOCI_POLYTABLE_SET_H
#define _INC_MULTI_LOCI_POLYTABLE_SET_H

#include "data_object_set.hpp"
#include "multi_loci_polytable.hpp"

#include <boost/smart_ptr.hpp>

#include <vector>


namespace mct {
	
	/** @brief A class for an ordered collection of 
	 * \link MultiLociPolyTable MultiLociPolyTables\endlink.
	  
	 The MultiLociPolyTable stored can have any number of loci or 
	 samples; there are no checks on what is added to the collection.
	 
	 A MultiLociPolyTableSet is a 
	\link aabc::DataObjectSet DataObjectSet\endlink
	type.
	 
	 */
	 

	class MultiLociPolyTableSet : public aabc::DataObjectSet {
				
		public:
		
			/*! \brief Constructor.
			
			Constructs an empty structure in which to store 
			\link MultiLociPolyTable MultiLociPolyTables\endlink.
			*/
			MultiLociPolyTableSet();
			
			virtual ~MultiLociPolyTableSet();
			
			/*! \brief Get a pointer to an indexed 
			 MultiLociPolyTable in the set.
			
			\param index The index of the MultiLociPolyTable to get.
			\return A pointer to the MultiLociPolyTable at \a index. */
			const boost::shared_ptr < const MultiLociPolyTable >&
					at(const size_t index) const;
			
			/*! \brief Add a MultiLociPolyTable to the 
			end of the collection.
			
			\param mlptr The MultiLociPolyTable to add.*/
			void add(
				const boost::shared_ptr < const MultiLociPolyTable >& mlptr);
			
			/*! \brief Add the contents of another 
			MultiLociPolyTableSet to the end of collection.
			
			\param sptr The MultiLociPolyTableSet to add.*/
			void add(
				const boost::shared_ptr < const MultiLociPolyTableSet >& sptr);
			
			/*! \brief Get the number of 
			 \link MultiLociPolyTable MultiLociPolyTables\endlink.
			 in the collection.	*/
			size_t size() const;
			
			/*! \brief Return whether this contains no data.	*/
			bool empty() const;
			
			/*! \brief Get the SummaryStatisticSet for the set.*/
			virtual boost::shared_ptr < SummaryStatisticSet > 
						getSummaryStatisticSet() const;
			
			/*! \brief Get the set of all available statistics 
			for the set.*/
			virtual boost::shared_ptr < SummaryStatisticSet > 
						getAllAvailableStatisticsSet() const;
						
			/*! \brief Get the SummaryStatistic for the set.
			
			The %SummaryStatistic for the set is the mean
			%SummaryStatistic for the MultiLociPolyTables in this 
			set.  See MultiLociPolyTable::getSummaryStatistic().*/
			virtual boost::shared_ptr < const SummaryStatistic >
											getSummaryStatistic() const;
			
			/*! \brief Get a Tajima's D %SummaryStatistic 
			for the set.
			
			The %SummaryStatistic for the set is the mean
			Tajima's D for the MultiLociPolyTables in this 
			set.*/
			virtual boost::shared_ptr < const SummaryStatistic >
							getTajimasDSummaryStatistic() const;
			
			/*! \brief Get a heterozygosity %SummaryStatistic 
			for the set.
			
			The %SummaryStatistic for the set is the mean
			heterozygosity for the MultiLociPolyTables in this 
			set.*/
			virtual boost::shared_ptr < const SummaryStatistic >
							getHeterozygositySummaryStatistic() const;
								
			/*! \brief Get a segregating sites %SummaryStatistic 
			for the set.
			
			The %SummaryStatistic for the set is the mean
			segregating sites for the MultiLociPolyTables in this 
			set.*/
			virtual boost::shared_ptr < const SummaryStatistic >
							getSegSitesSummaryStatistic() const;
			
			/*! \brief Get the means of the Tajimas D's 
			for the collection of MultiLociPolyTables
			as a SummaryStatisticSet
			
			Each SummaryStatistic in the set contains the value
			of the mean of Tajimas D's for the corresponding 
			MultiLociPolyTable in this.	*/
			virtual boost::shared_ptr < SummaryStatisticSet > 
										getTajimasDMeans() const;
			
			/*! \brief Get the means of the heterozygosities 
			for the collection of MultiLociPolyTables
			as a SummaryStatisticSet
			
			Each SummaryStatistic in the set contains the value
			of the mean of the heterozygosity statistic for the corresponding 
			MultiLociPolyTable in this.	*/
			virtual boost::shared_ptr < SummaryStatisticSet > 
										getHeterozygosityMeans() const;
			
			/*! \brief Get the means of the seg sites 
			for the collection of MultiLociPolyTables
			as a SummaryStatisticSet
			
			Each SummaryStatistic in the set contains the value
			of the mean of the number of segregating sites for the corresponding 
			MultiLociPolyTable in this.	*/
			virtual boost::shared_ptr < SummaryStatisticSet > 
										getSegSitesMeans() const;
		
					
			/*! \brief Get a string representation of this.*/
			virtual std::string toString() const;
			
						
		protected:
		
			explicit MultiLociPolyTableSet(const MultiLociPolyTableSet& other);
			
			MultiLociPolyTableSet& operator=(MultiLociPolyTableSet tmp);
					
			/* \brief Get the Tajimas D's for the collection.	
			 
			 The inner vectors are the statistics for each 
			 MultiLociPolyTable in the collection. The outer vector 
			 packages these together.  To get the statistic for the 
			 j<sup>th</sup> element within the i<sup>th</sup> 
			 MultiLociPolyTable, one would index
			 return value by [i][j].  
			 */
			virtual std::vector < std::vector < double > > 
						getTajimasD() const;
			
			/* \brief Get the heterozygosity for the collection 
			of Sequence::SimData.	
			The inner vectors are the statistics for each 
			 MultiLociPolyTable in the collection. The outer vector 
			 packages these together.  To get the statistic for the 
			 j<sup>th</sup> element within the i<sup>th</sup> 
			 MultiLociPolyTable, one would index
			 return value by [i][j].*/
			virtual std::vector < std::vector < double > > 
						getHeterozygosity() const;
			
			/* \brief Get the number of segregating sites 
			for the collection. 
			The inner vectors are the statistics for each 
			 MultiLociPolyTable in the collection. The outer vector 
			 packages these together.  To get the statistic for the 
			 j<sup>th</sup> element within the i<sup>th</sup> 
			 MultiLociPolyTable, one would index
			 return value by [i][j].*/
			virtual std::vector < std::vector < size_t > > 
						getSegSites() const;
			
					
			
			const std::vector< boost::shared_ptr < const MultiLociPolyTable > >&
				_getContainer() const;
			
				
			
					
		private:
		
			/*! \brief The collection of pointers to 
			MultiLociPolyTable held.*/
			std::vector< boost::shared_ptr < const MultiLociPolyTable > > container;
			
			
	};
	
	// does not need to be friend
	std::ostream& operator<< (std::ostream &out, 
								const MultiLociPolyTableSet &m);

} // end namespace mct

#endif 

