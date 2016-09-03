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
\brief DescriptiveStats static class declarations.
*/

#ifndef _DESCRIPTIVE_STATS_H
#define _DESCRIPTIVE_STATS_H

#include "summary_statistic_set.hpp"

#include "boost/tuple/tuple_io.hpp"

#include <vector>
#include <string>

namespace mct {
	
	/*! \brief A structure for holding the descriptive 
	 statistics.
	 
	 The ordering is:
	 <ul>
	 <li>Mean</li>
	 <li>Sample standard deviation</li>
	 <li>Minimum</li>
	 <li>Lower critical value for an inter-percentile 
	 range</li> 
	 <li>Lower quartile</li>
	 <li>Median</li>
	 <li>Upper quartile</li>
	 <li>Upper critical value for an inter-percentile 
	 range</li>
	 <li>Maximum</li>
	 <li>Collection size n</li>
	 </ul>
	 
	 The inter-percentile range for which to calculate the upper and 
	 lower critical values is configured in the constructor and defaults
	 to 0.9 (90\%).
	 For example, if the default value of 0.90 (90\%) is given,
	 the lower critical value will be the 
	 0.05<sup>th</sup> percentile in the collection and 
	 the upper critical value will be the 
	 0.95<sup>th</sup> percentile in the collection.
	 
	 */
	
		
	class DescriptiveStats {
		public:
			/*! \brief A structure for holding the descriptive 
			 statistics.
			 
			 \note The medians, quartiles and percentiles are 
			 are calculated using interpolation equivalent
			 to the method used in gsl_stats_quantile_from_sorted_data.
			 
			 The ordering is:
			 <ul>
			 <li>Mean</li>
			 <li>Sample standard deviation</li>
			 <li>Minimum</li>
			 <li>Lower critical value for an inter-percentile 
			 range</li> 
			 <li>Lower quartile</li>
			 <li>Median</li>
			 <li>Upper quartile</li>
			 <li>Upper critical value for an inter-percentile 
			 range</li>
			 <li>Maximum</li>
			 <li>Collection size n</li>
			 </ul>
			 * 
			 * @todo add ability to take a vector of vectors to this -
			 *  maybe with indicator for which way to collate stats.
			 */
			typedef boost::tuple< double, double, double, double, 
						double, double, double, double, double, 
						size_t > DStats;

			
			/*! \brief Constructor.
			 
			 \param vec The vector of doubles to calculate stats for.
			 \param c The inter-percentile range for which to calculate 
			 the upper and lower critical values for. Defaults to 0.9
			 (90\%).
			 \pre 0.0 < c <= 1.0.
			 */
			DescriptiveStats(const std::vector < double >& vec,
								double c = 0.90);
			
			/*! \brief Constructor.
			 
			 This constructor assumes it will get a vector of inner vectors
			 where each inner vector represents a sample of values
			 for one random variable.  The constructor makes a DescriptiveStats
			 object describing a sample of all the elements in the 
			 all the inner vector with index \a index.
			 
			 \param vecvec The vector of vector of doubles to calculate stats from.
			 \param index The index of the inner vector in \a vecvec for
			 which to calculate the stats.
			 \param c The inter-percentile range for which to calculate 
			 the upper and lower critical values for. Defaults to 0.9
			 (90\%).
			 \pre 0 <= \a index < vecvec.size().
			 \pre 0.0 < c <= 1.0.
			 
			 */
			DescriptiveStats(
				const std::vector < std::vector < double > >& vecvec,
				int index, 
				double c = 0.90);
			
			/*! \brief Constructor.
			 
			 Constructor to calculate descriptive statistics for
			 a sample of all 
			 the values at an indexed position in the SummaryStatistics in a 
			 given SummaryStatisticsSet.
			 
			 \param ss_set The %SummaryStatisticSet to calculate stats from.
			 \param index The index of the elements in the 
			 SummaryStatistics in \a ss_set for
			 which to calculate the stats.
			 \param c The inter-percentile range for which to calculate 
			 the upper and lower critical values for. Defaults to 0.9
			 (90\%).
			 \pre 0 <= \a index < maximum number of values 
			 in any %SummaryStatistic in \a ss_set.
			 \pre 0.0 < c <= 1.0.
			 */
			DescriptiveStats(
				const SummaryStatisticSet& ss_set,
				int index, 
				double c = 0.90);
			
			/*! \brief Constructor.
			 
			 Constructor to calculate descriptive statistics for
			 a sample of all 
			 the values at an indexed position in the SummaryStatistics in a 
			 SummaryStatisticsSet given with a pointer.
			 
			 \param ss_ptr Pointer to the %SummaryStatisticSet to calculate stats from.
			 \param index The index of the elements in the 
			 SummaryStatistics in \a ss_set for
			 which to calculate the stats.
			 \param c The inter-percentile range for which to calculate 
			 the upper and lower critical values for. Defaults to 0.9
			 (90\%).
			 \pre 0 <= \a index < maximum number of values 
			 in any %SummaryStatistic in \a ss_set.
			 \pre 0.0 < c <= 1.0.
			 */
			DescriptiveStats(
				const boost::shared_ptr < const SummaryStatisticSet > ss_ptr,
				int index, 
				double c = 0.90);
				
			virtual ~DescriptiveStats();
			
			/*! \brief Return the mean value \f$ \bar{x} \f$.*/	
			double mean() const;
			
			/*! \brief Return the sample standard deviation value.
			 
			 \note this calculation returns
			 \f$ \frac{\sum (x_i - \bar{x})^2}{(n-1)} \f$.*/
			double sampleSD() const;
			
			
			/*! \brief Return the minimum value.*/	
			double min() const;
			
			/*! \brief Return the lower end of the critical 
			 inter-percentile range.*/	
			double lowerCriticalValue() const;
			
			/*! \brief Return the lower quartile 
			(25th percentile) value.*/	
			double LQ() const;
			
			/*! \brief Return the median value.*/	
			double median() const;
			
			/*! \brief Return the upper quartile 
			(75th percentile) value.*/	
			double UQ() const;
			
			/*! \brief Return the upper end of the critical 
			 inter-percentile range.*/	
			double upperCriticalValue() const;
			
			/*! \brief Return the maximum value.*/	
			double max() const;
			
			/*! \brief Return the collection size n.*/	
			size_t n() const;
			
			/*! \brief Return all the statistics 
			 in a DStats structure.*/	
			DStats all() const;
			
			/*! \brief Return the inter-percentile range being used
			 for the upper and lower critical values.*/	
			double getCriticalPercentileRange() const;
			
			/*! \brief Get a string representation of this.
			
			Numbers are printed using fixed formatting.
			\todo Make this more flexible 
			
			\param prec The precision to which to show the output.
			\return A string representation of this.*/
			std::string toString(int prec = 3) const;
			
			
		protected:
		
			typedef boost::tuple<double, double> expStats;
			
			typedef boost::tuple<double, double,
				double, double, double, double, double, 
				size_t > rngStats;
		
			DescriptiveStats();
			
			virtual expStats calcMeanAndSD(
							const std::vector < double >& vec);
		
			virtual rngStats calcRangeStats(
							const std::vector < double >& vec);
								
			static double getPercentile(
							const std::vector < double >& sortedvec,
							double percentile);
		
			double criticalPercentileRange;
		
			DStats stats;
			
		private:
		
			explicit DescriptiveStats(const DescriptiveStats& other);
			 
			DescriptiveStats& operator=(DescriptiveStats rhs);
			
			void swap (DescriptiveStats& ds);
		
	};

} // end namespace mct

#endif 

