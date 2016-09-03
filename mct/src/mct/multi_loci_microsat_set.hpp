/* MCT - Markov Chains on Trees.

   Copyright (C) 2011, 2012 Jenny Harlow

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
\brief MultiLociMicrosatSet declarations.
*/

#ifndef _INC_MULTI_LOCI_MICROSAT_SET_H
#define _INC_MULTI_LOCI_MICROSAT_SET_H

#include "data_object_set.hpp"
#include "multi_loci_microsat.hpp"

#include <boost/smart_ptr.hpp>


#include <vector>


namespace mct {
	
	/** @brief A class for an ordered collection of 
	 \link MultiLociMicrosat MultiLociMicrosats\endlink.
	  
	 The MultiLociMicrosats stored can have any number of loci or 
	 samples; there are no checks on what is added to the collection.
	 
	 A %MultiLociMicrosatSet is a 
	\link aabc::DataObjectSet DataObjectSet\endlink
	type.
	 
	 */
	 

	class MultiLociMicrosatSet : public aabc::DataObjectSet {
				
		public:
		
			/*! \brief Constructor.
			
			Constructs an empty structure in which to store 
			\link MultiLociMicrosat MultiLociMicrosats\endlink.
			*/
			MultiLociMicrosatSet();
			
			virtual ~MultiLociMicrosatSet();
			
			/*! \brief Get a pointer to an indexed 
			 MultiLociMicrosat in the set.
			
			\param index The index of the MultiLociMicrosat to get.
			\return A pointer to the MultiLociMicrosat at \a index. */
			const boost::shared_ptr < const MultiLociMicrosat >&
					at(const size_t index) const;
			
			/*! \brief Add a MultiLociMicrosat to the 
			end of the collection.
			
			\param mlptr The MultiLociMicrosat to add.*/
			void add(
				const boost::shared_ptr < const MultiLociMicrosat >& mlptr);
			
			/*! \brief Add the contents of another 
			MultiLociMicrosatSet to the end of collection.
			
			\param sptr The MultiLociMicrosatSet to add.*/
			void add(
				const boost::shared_ptr < const MultiLociMicrosatSet >& sptr);
			
			/*! \brief Get the number of 
			\link MultiLociMicrosat MultiLociMicrosats\endlink.
			 in the collection.	*/
			size_t size() const;
			
			/*! \brief Return whether this contains no data.	*/
			bool empty() const;
			
			/*! \brief Get the SummaryStatisticSet for the set.*/
			virtual boost::shared_ptr < SummaryStatisticSet > 
						getSummaryStatisticSet() const;
			
			/*! \brief Get the SummaryStatistic for the set.*/
			virtual boost::shared_ptr < const SummaryStatistic >
											getSummaryStatistic() const;
			
			
			/*! \brief Get a SummaryStatisticSet holding
			the averages of the repeats 
			for the collection of MultiLociMicrosats.	*/
			virtual boost::shared_ptr < SummaryStatisticSet > 
										getAverageRepeats() const;
			
			/*! \brief Get a SummaryStatisticSet holding
			the variance-covariance
			representations 
			for the collection of MultiLociMicrosats.
			
			Each SummaryStatistic in the set holds the top diagonal 
			half of the
			inter-individual variance covariance matrix in row
			major order.
			*/
			virtual boost::shared_ptr < SummaryStatisticSet > 
										getCovariances() const;
			
			/** @name Population comparison statistics.
			
			These functions return values for statistics 
			 used to compare and evaluate 'distances' between the
			 subpopulations with the populations represented
			 by this.  
			  
			Where between-subpopulation comparison statistics are given 
			 and there are \f$ P > 0 \f$ subpopulations,
			  indexed 
			 \f$ j = 0, 1, 2, \ldots, P-1 \f$, and \f$ stat_{(j,k)} \f$
			 , \f$ j \neq k \f$ is the \f$ stat \f$ is the 
			 appropriate statistic calculated between
			  subpopulations indexed \f$ j \f$ and \f$ k \f$, 
			  the the values in each summmary statistic will be ordered
			  \f$ stat_{(0,1)},\ldots, stat_{(0,P-1)},stat_{(1,2)}, \ldots, stat_{(1,P-1)}, \ldots, stat_{(P-2,P-1)} \f$.
			.  If there is just one subpopulation in the population 
			represented by this, any collection of between-population comparison
			statistics will be empty. 
			
			 Where pairwise comparison statistics are are given 
			 and there are \f$ P > 0 \f$ subpopulations,
			  indexed 
			 \f$ j = 0, 1, 2, \ldots, P-1 \f$, and \f$ stat_{(j,k)} \f$
			 , \f$ j \neq k \f$ is the \f$ stat \f$ is the 
			 appropriate statistic calculated between
			  subpopulations indexed \f$ j \f$ and \f$ k \f$, 
			  the the values in each summmary statistic will be ordered
			  \f$ stat_{(0,0)},stat_{(0,1)},\ldots, stat_{(0,P-1)},stat_{(1,1)},stat_{(1,2)}, \ldots, stat_{(1,P-1)}, \ldots, stat_{(P-2,P-2)},stat_{(P-2,P-1)},stat_{(P-1,P-1)} \f$.
			.  If there is just one subpopulation in the population 
			represented by this, any collection of pairwise comparison
			statistics will be contain just one element. 
			
			 */
			//@{
			
			/*! \brief Get a SummaryStatisticSet holding
			the between-population
			Goldstein Delta-Mu Squared distances for the collection
			of MultiLociMicrosats.			*/
			virtual boost::shared_ptr < SummaryStatisticSet > 
								getDistanceDeltaMuSqBetweenPop() const;
			
			/*! \brief Get a SummaryStatisticSet holding
			the pairwise inter-population
			average sum squared difference in repeat sizes
			for the collection of MultiLociMicrosats.			*/
			virtual boost::shared_ptr < SummaryStatisticSet > 
								getPairwiseAverageSumSqDifferences() const;
			
			/*! \brief Get a SummaryStatisticSet holding
			the pairwise inter-population
			Slatkin Rst distances for the collection of MultiLociMicrosats.*/
			virtual boost::shared_ptr < SummaryStatisticSet >
								getDistanceSlatkinRstBetweenPop() const;
			
			/*! \brief Get a SummaryStatisticSet holding
			the between-population \f$ Fst \f$ values for
			for the collection of MultiLociMicrosats.			*/					
			virtual boost::shared_ptr < SummaryStatisticSet > 
								getBetweenPopFst() const;
			
			/*! \brief Get a SummaryStatisticSet holding
			the \f$ Rst \f$ values for
			for the collection of MultiLociMicrosats.			*/
			virtual boost::shared_ptr < SummaryStatisticSet > 
								getRst() const;
			
			/*! \brief Get a SummaryStatisticSet holding
			the \f$ Fst \f$ values for
			for the collection of MultiLociMicrosats.			*/					
			virtual boost::shared_ptr < SummaryStatisticSet > 
								getFst() const;
			
			
			//@}
			
			/*! \brief Get a string representation of this.*/
			virtual std::string toString() const;
			
						
		protected:
		
			explicit MultiLociMicrosatSet(const MultiLociMicrosatSet& other);
			
			MultiLociMicrosatSet& operator=(MultiLociMicrosatSet tmp);
					
			const std::vector< boost::shared_ptr < const MultiLociMicrosat > >&
				_getContainer() const;
			
				
			
					
		private:
		
			/*! \brief The collection of pointers to 
			MultiLociMicrosat held.*/
			std::vector< boost::shared_ptr < const MultiLociMicrosat > > container;
			
			
	};
	
	// does not need to be friend
	std::ostream& operator<< (std::ostream &out, 
								const MultiLociMicrosatSet &m);

} // end namespace mct

#endif 

