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
\brief SummaryStatisticSet declarations.
*/

#ifndef _INC_SUMMARY_STATISTIC_SET_H
#define _INC_SUMMARY_STATISTIC_SET_H

#include "data_object.hpp"
#include "summary_statistic.hpp"

#include <boost/smart_ptr.hpp>

#include <vector>
#include <string>


namespace mct{

	//forward declaration
	//class SummaryStatistic;
	
	/** @brief A class for an ordered collection of 
	 \link SummaryStatistic SummaryStatistics\endlink.
	  
	 There are no checks on what is added to the collection.
	 
	 */
	 

	class SummaryStatisticSet : public aabc::DataObject { 
				
		public:
		
			/*! Iterator into the ordered collection of pointers to 
			\link mct::SummaryStatistic SummaryStatistics\endlink.*/
			typedef 
			std::vector< 
					boost::shared_ptr< const SummaryStatistic > >
					::const_iterator
					const_iterator;
				
			/*! \brief Constructor.
			
			Constructs an empty structure in which to store 
			SummaryStatistic.
			*/
			SummaryStatisticSet();
			
			/*! \brief Constructor.
			
			Constructs an empty structure in which to store 
			\link mct::SummaryStatistic SummaryStatistics\endlink,
			but reserves enough memory space for
			\a reservedSpace elements.  This is useful if the 
			set is to be used to store large numbers of
			objects because it can avoid or reduce memory
			reallocations when add is called.
			
			\param reservedSpace The number of elements to reserve
			space for in memory.
			\note The reserved space is not a limit on the capacity
			of the set nor does it stuff the container with empty
			elements, it just ensures that the set's memory 
			allocation is sufficient for that number of elements
			to be added without reallocation.  If more elements are 
			then added, memory will be reallocated.
			*/
			explicit SummaryStatisticSet(std::size_t reservedSpace);
			
			explicit SummaryStatisticSet(
								const SummaryStatisticSet& other);
			
			
			virtual ~SummaryStatisticSet();
			
			SummaryStatisticSet& operator=(SummaryStatisticSet tmp);
			
			/*! \brief Get a pointer to an indexed 
			 SummaryStatistic in the set.
			
			\param index The index of the %SummaryStatistic to get.
			\return A pointer to the SummaryStatistic at \a index. */
			virtual const boost::shared_ptr < const SummaryStatistic >&
					at(const size_t index) const;
			
			/*! \brief Add a SummaryStatistic to the collection.
			
			\param sptr The %SummaryStatistic to add.*/
			virtual void add(
				const boost::shared_ptr < const SummaryStatistic >& sptr);
			
			/*! \brief Add a SummaryStatisticSet to the collection.
			\param setptr The %SummaryStatisticSet to add.*/
			virtual void add(
				const boost::shared_ptr < SummaryStatisticSet >& setptr);
			
			/*! \brief Get the number of 
			 \link mct::SummaryStatistic SummaryStatistics\endlink
			 in the collection.	*/
			size_t size() const;
			
			/*! \brief Return whether this contains no data.	*/
			bool empty() const;
			
			/*! \brief Get the SummaryStatistic for this.
			 
			 \return A pointer to a %SummaryStatistic containing
			 the mean values of the SummaryStatistics in this.*/
			virtual boost::shared_ptr < const SummaryStatistic  > 
										getSummaryStatistic() const;
			
			
			/*! \brief Make and return a %SummaryStatisticSet
			holding the values in this, standardised.
			 
			Standardising means centering each summary statistic in this on its
			mean (ie subtracting the summary statistic mean from each
			summary statistic held) and scaling each summary centred
			summary statistic by its sample standard deviation (ie for each
			centred summary statistic, dividing
			the each element by the
			corresponding element of the sample standard deviations of
			the summary statistics.
			
			\return a pointer to a %SummaryStatisticSet holding	
			standardised values in the SummaryStatistics in this.
			\pre If this is not empty, there should be at least two
			elements in the collection
			and the sizes of all SummaryStatistics in the collection 
			should be consistent.*/		
			virtual 
				boost::shared_ptr < SummaryStatisticSet >
					makeSummaryStatisticSetStandardised() const;
			
			
			/*! \brief Make and return a %SummaryStatisticSet
			holding the values in this, standardised, and also
			fill in values for means and sds in supplied containers.
			 
			Standardising means centering each summary statistic in this on its
			mean (ie subtracting the summary statistic mean from each
			summary statistic held) and scaling each summary centred
			summary statistic by its sample standard deviation (ie for each
			centred summary statistic, dividing
			the each element by the
			corresponding element of the sample standard deviations of
			the summary statistics.
			
			\param means A container to be filled with the means of 
			the SummaryStatistics in this.
			\param sds A container to be filled with the 
			sample standard deviations of 
			the SummaryStatistics in this.
			\return a pointer to a %SummaryStatisticSet holding	
			standardised values in the SummaryStatistics in this.
			\pre If this is not empty, there should be at least two
			elements in the collection
			and the sizes of all SummaryStatistics in the collection 
			should be consistent.
			\post \a means contains the means of the SummaryStatistics
			in this, \a sds contains the sample standard deviations
			of the SummaryStatistics in this*/				
			virtual 
				boost::shared_ptr < SummaryStatisticSet >
					makeSummaryStatisticSetStandardised(
						std::vector < double >& means,
						std::vector < double >& sds) const;
						
			/*! \brief Make and return a %SummaryStatisticSet
			holding a selection of the values in this.
			 
			\param reqStats specifies the indices of the required 
			SummaryStatistics for the new set.  The indices 
			can take values from 1 upwards.  The maximum index
			is the number of statistics in each %SummaryStatistic
			in this.
			\return a pointer to a %SummaryStatisticSet holding	
			the required values from this.
			\pre \a reqStats is not empty, and the minimum value in
			\a reqStats is 1, and the maximum is the number of 
			statistics in each SummaryStatistic in this.
			\pre The sizes of all SummaryStatistics in the collection 
			should be consistent.*/		
			virtual boost::shared_ptr < SummaryStatisticSet >
				makeSelectedSummaryStatisticSet(
				std::vector < size_t > reqStats) const;
			
			/*! \brief Get all the values in the 
			\link SummaryStatistic SummaryStatistics\endlink
			held.
			
			 If this contains n
			 \link mct::SummaryStatistic SummaryStatistics\endlink
			 with nstats elements in each, then the outer container in
			 the returned structure contains nstats inner containers,
			 and each inner container contains n values, and the
			 inner containers are in the same order as the 
			 elements are in the SummaryStatistics.  ie, the 
			 i<sup>th</sup> value in the j<sup>th</sup> 
			 %SummaryStatistic in the set will be the j<sup>th</sup>
			 element of the i<sup>th</sup> inner container.
			 If any
			 %SummaryStatistic in the set contains less than
			 nstats elements then the first nstats inner containers
			 will contain the values from that %SummaryStatistic and
			 the remaining inner containers will have a NaN at
			 that position.
			\return a pointer to a container of containers of	
			values in the SummaryStatistics in this.*/		
			virtual 
				boost::shared_ptr < std::vector < std::vector < double > > >
					getAllValues() const;
					
			/*! \brief Get all the values in the 
			\link SummaryStatistic SummaryStatistics\endlink
			held, laid out like summary statistics.
			
			 If this contains n
			 \link mct::SummaryStatistic SummaryStatistics\endlink,
			 then the outer container in
			 the returned structure contains n inner containers,
			 and the i<sup>th</sup> inner container will contain the 
			 values from the i<sup>th</sup> %SummaryStatistic in this.
			
			\return a pointer to a container of containers of	
			values in the SummaryStatistics in this.*/		
			virtual 
				boost::shared_ptr < std::vector < std::vector < double > > >
					getAllValuesSumStatLayout() const;
			
					
			/*! @name Get all the values in the 
			\link SummaryStatistic SummaryStatistics\endlink
			held, standardised.
			 
			Standardising means centering each summary statistic on its
			mean (ie subtracting the summary statistic mean from each
			summary statistic held) and scaling each centred
			summary statistic by its sample standard deviation (ie for each
			centred summary statistic, dividing
			the each element by the
			corresponding element of the sample standard deviations of
			the summary statistics.
			
			 If this contains n
			 \link mct::SummaryStatistic SummaryStatistics\endlink
			 with nstats elements in each, then the outer container in
			 the returned structure contains nstats inner containers,
			 and each inner container contains n values, and the
			 inner containers are in the same order as the 
			 elements are in the SummaryStatistics.  ie, the 
			 i<sup>th</sup> value in the j<sup>th</sup> 
			 %SummaryStatistic in the set will be the result of
			 standardising the j<sup>th</sup>
			 element of the i<sup>th</sup> inner container
			 using mean and standard deviation of the j<sup>th</sup>
			 elements over all inner containers.
			 \return a pointer to a container of containers of	
			standardised values in the SummaryStatistics in this.
			\pre If this is not empty, there should be at least two
			elements in the collection
			and the sizes of all SummaryStatistics in the collection 
			should be consistent.*/	
		
			//@{
			/*! \brief Version which takes references to two containers
			 *  as parameters; the operation fills these containers with
			 * the means and sample standard deviations of the data.*/ 
			virtual 
				boost::shared_ptr < std::vector < std::vector < double > > >
					getAllValuesStandardised(
						std::vector < double >& means,
						std::vector < double >& sds) const;
			
			virtual 
				boost::shared_ptr < std::vector < std::vector < double > > >
					getAllValuesStandardised() const;
			
			//@}
			
			/*! \brief Get all the standardised values in the 
			\link SummaryStatistic SummaryStatistics\endlink
			held, laid out like summary statistics.
			
			Standardising means centering each summary statistic on its
			mean (ie subtracting the summary statistic mean from each
			summary statistic held) and scaling each centred
			summary statistic by its sample standard deviation (ie for each
			centred summary statistic, dividing
			the each element by the
			corresponding element of the sample standard deviations of
			the summary statistics.
						
			 If this contains n
			 \link mct::SummaryStatistic SummaryStatistics\endlink,
			 then the outer container in
			 the returned structure contains n inner containers,
			 and the i<sup>th</sup> inner container will contain the 
			 standardised values from 
			 the i<sup>th</sup> %SummaryStatistic in this.
			
			\return a pointer to a container of containers of	
			standardised values in the SummaryStatistics in this.*/		
			
			//@{
			/*! \brief Version which takes references to two containers
			 *  as parameters; the operation fills these containers with
			 * the means and sample standard deviations of the data.*/ 
			virtual 
				boost::shared_ptr < std::vector < std::vector < double > > >
					getAllValuesStandardisedSumStatLayout(
						std::vector < double >& means,
						std::vector < double >& sds) const;
			
			virtual 
				boost::shared_ptr < std::vector < std::vector < double > > >
					getAllValuesStandardisedSumStatLayout() const;
			//@}
			/*! \brief Get a SummaryStatistic containing 
			the means of the 
			\link SummaryStatistic SummaryStatistics\endlink
			from the collection of
			\link SummaryStatistic SummaryStatistics\endlink.
			
			\return A pointer to a SummaryStatistic
			whose values are the means of the values in 
			the SummaryStatistics in the set, i.e. the i<sup>th</sup>
			value in the returned %SummaryStatistic is the mean
			of the i<sup>th</sup> values of all the 
			SampleStatistics in the set. 
			\pre There should be at least one element in the collection
			and the sizes of all SummaryStatistics in the collection 
			should be consistent.*/				
			virtual boost::shared_ptr < const SummaryStatistic >
						getSummaryStatisticMeans() const;
						
			/*! \brief Get a SummaryStatistic containing 
			the sample variances of the 
			\link SummaryStatistic SummaryStatistics\endlink
			from the collection of
			\link SummaryStatistic SummaryStatistics\endlink.
			
			\return A pointer to a SummaryStatistic
			whose values are the sample variances of the values in 
			the SummaryStatistics in the set, i.e. the i<sup>th</sup>
			value in the returned %SummaryStatistic is the sample
			variance of the i<sup>th</sup> values of all the 
			SampleStatistics in the set. If there is only one
			%SampleStatistic in the set the sample variances
			in the returned %SummaryStatistic will all be 0.0.					
			\pre There should be at least one element in the collection
			and the sizes of all SummaryStatistics in the collection 
			should be consistent.*/				
			virtual boost::shared_ptr < const SummaryStatistic  > 
						getSummaryStatisticVars() const;
			
			/*! \brief Get a SummaryStatistic containing 
			the sample standard deviations of the 
			\link SummaryStatistic SummaryStatistics\endlink
			from the collection of
			\link SummaryStatistic SummaryStatistics\endlink.					
			
			\return A pointer to a SummaryStatistic
			whose values are the sample standard deviations of the values in 
			the SummaryStatistics in the set, i.e. the i<sup>th</sup>
			value in the returned %SummaryStatistic is the sample
			standard deviation of the i<sup>th</sup> values of all the 
			SampleStatistics in the set. If there is only one
			%SampleStatistic in the set the sample standard deviations
			in the returned %SummaryStatistic will all be 0.0.					
			\pre There should be at least one element in the collection
			and the sizes of all SummaryStatistics in the collection 
			should be consistent.*/
			virtual boost::shared_ptr < const SummaryStatistic  > 
						getSummaryStatisticSDs() const;
			
			/*! \brief Get a fancy string representation of this.
			
			\return A formatted string representation of this.*/
			virtual std::string toString() const;
			
			/*! \brief Get a fancy string representation of this.
			
			\param prec The precision to which to show the output. 
			\return A formatted string representation of this.*/
			virtual std::string toString(int prec) const;
			
			/*! \brief Get a plain string representation of this.
			
			Numbers are printed using fixed formatting.
			\todo Make this more flexible 
			
			\param prec The precision to which to show the output.
			\return A plain string representation of this.*/
			virtual std::string plainString(int prec = 3) const;
			
			/*! \brief Output to a file, same format as for plainString().
			
			\param s The filename for the output (which can include a path).
			\param append An indicator for whether output is appended to
			the file (if it exists already) or overwrites existing file content.
			\post A file called \a s exists containing output from this.*/
			virtual void outputToFile(	const std::string& s,
									bool append = true) const;

			
			/*! \brief Return whether the sizes of the 
			 \link mct::SummaryStatistic SummaryStatistics\endlink
			 in the collection are consistent.
			 
			 \return true if the collection is empty or if the sizes of
			 of the contained 
			 SummaryStatistics
			 are all the same, false otherwise.
			 */
			virtual bool sizeConsistencyCheck() const;
			
			/*! \brief Return a const iterator to the first 
			element in the set.*/
			const_iterator begin() const;
			
			/*! \brief Return a const iterator to the
			'just past the end' element in the set.
			
			\note The idiom for iterating through a 
			%SummaryStatisticSet ss_set is 
			for (SummaryStatisticSet::const_iterator it = ss_set.begin();
						it < ss_set.end();
						++it) {...}.
			*/
			const_iterator end() const;
									
		protected:
		
			std::vector<double>& getSummaryStatisticsMeanVec(
								std::vector < double >& result) const;
			
			std::vector<double>& getSummaryStatisticsVarVec(
									std::vector<double>& result) const;
			
			virtual const std::vector< boost::shared_ptr < const SummaryStatistic > >&
				_getContainer() const;
									
			static double standardise(double value, double mean, double sd);
				
		private:
		
			/*! \brief The collection of pointers to 
			SummaryStatistic held.*/
			std::vector< boost::shared_ptr < const SummaryStatistic > >
														container;
			
	};
	
	// does not need to be friend
	std::ostream& operator<< (std::ostream &out, 
								const SummaryStatisticSet &s);
	
	

} 

#endif 

