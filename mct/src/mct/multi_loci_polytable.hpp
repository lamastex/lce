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
\brief MultiLociPolyTable declarations.
*/

#ifndef _INC_MULTI_LOCI_POLYTABLE_H
#define _INC_MULTI_LOCI_POLYTABLE_H


#include "genetic_pop_structure_analysable.hpp"
#include <mscplusplus/population_structure.hpp>

#include <Sequence/PolyTable.hpp>

#include <boost/smart_ptr.hpp>

#include <vector>
#include <map>


namespace mct {
	
	
	/** @brief A class for an ordered collection of Sequence::PolyTable.
	
	A MultiLociPolyTable is a \link aabc::DataObject DataObject\endlink
	type.
	 
	 Each Sequence::PolyTable represents the sample sequence data 
	 encoded as polymorphic site data at one locus.
	 
	 Loci are independent. 
	 
	 Each locus can cover a different number of sites.   
	 
	 */
	 

	class MultiLociPolyTable : public mct::GeneticPopStructureAnalysable {
		public:
		
			/*! \brief Constructor.
			
			\param _p Pointer to the population structure for this.*/
			explicit MultiLociPolyTable(
				const boost::shared_ptr< const hudson_ms::PopulationStructure > _p);
			
			virtual ~MultiLociPolyTable();
			
			/*! \brief Get the Sequence::PolyTable for an
			indexed locus in the set.
			
			\param index The index of the locus to get the PolyTable for.
			\return A pointer to the PolyTable at \a index. */
			const boost::shared_ptr< Sequence::PolyTable >&
			at(const size_t index) const;
			
			/*! \brief Get the sites represented by an indexed 
			locus.
			
			\param index The index of the locus to get the sites for.
			\return The number of sites represented by the locus at \a index. */
			size_t nsitesAt(const size_t index) const;
			
			/*! \brief Add a Sequence::PolyTable to the set.
			
			\param toAdd The PolyTable to add.
			\param _sites The number of sites represented by \toAdd. */
			void add(
				const boost::shared_ptr <Sequence::PolyTable >& toAdd,
				size_t _sites);
			
			/*! \brief Get the population structure.	*/
			const boost::shared_ptr< const hudson_ms::PopulationStructure > 
								populationStructure() const;
			
			/*! \brief Get the number of samples in the 
			 Sequence::PolyTable in the collection.	*/
			size_t nsam() const;
			
			/*! \brief Get the number of Sequence::PolyTable
			 in the collection.	*/
			size_t nloci() const;
			
			/*! \brief Get the number of sites in each locus .
			
			\return an ordered collection showing the number of
			sites in each locus in this.*/
			std::vector < size_t > nsites() const;
			
			/*! \brief Return whether this contains no data.	*/
			bool empty() const;
			
			/*! \brief Get the 
			 \link SummaryStatistic SummaryStatistic\endlink
			 for this.
			 
			The summary statistic is (meanSegSite, meanHeterozygosity) */
			virtual boost::shared_ptr < const SummaryStatistic  > 
									getSummaryStatistic() const;
			
			/*! \brief Get all the available statistics for this
			 in the same form as a 
			 \link SummaryStatistic SummaryStatistic\endlink.
			  
			The available stats are
			(meanSegSite, meanHeterozygosity, meanTajimasD) */
			virtual boost::shared_ptr < const SummaryStatistic  > 
									getAllAvailableStatistic() const;
			
			/*! \brief Get the mean of the Tajimas Ds over the collection 
			of Sequence::PolyTable in a SummaryStatistic.	*/
			virtual boost::shared_ptr < const SummaryStatistic  >
										getMeanTajimasD() const;
			
			/*! \brief Get the mean heterozygosity over the collection 
			of Sequence::PolyTable  in a SummaryStatistic.	*/
			virtual boost::shared_ptr < const SummaryStatistic  >
										getMeanHeterozygosity() const;
			
			/*! \brief Get the mean number of segregating sites 
			over the collection of Sequence::PolyTable
			in a SummaryStatistic. */
			virtual boost::shared_ptr < const SummaryStatistic  >
										getMeanSegSites() const;
			
			/*! \brief Get the value of the mean of the Tajimas Ds 
			over the collection of Sequence::PolyTable.	*/
			virtual double getMeanTajimasDValue() const;
			
			/*! \brief Get the value of the mean heterozygosity
			over the collection of Sequence::PolyTable.	*/
			virtual double getMeanHeterozygosityValue() const;
			
			/*! \brief Get the value of the mean number of segregating sites 
			over the collection of Sequence::PolyTable. */
			virtual double getMeanSegSitesValue() const;
			
			
			/*! \brief Get a string representation of this.*/
			virtual std::string toString() const;
			
			/*! \brief Get the Tajimas Ds as a collection, one 
			element for every locus in this.	*/
			virtual std::vector < double > getTajimasD() const;
			
			/*! \brief Get the number of segregating sites as a collection, one 
			element for every locus in this.	*/
			virtual std::vector < size_t > getSegSites() const;
			
			/*! \brief Get the heterozygosity as a collection, one 
			element for every locus in this.	*/
			virtual std::vector < double > getHeterozygosity() const;
			
			/*! \brief Print the polysites information to an output
			stream, separating the subpopulation data. 
			* 
			Each line is information over all loci for one individual in
			the sample.	*/
			virtual std::ostream& printSubpopData(std::ostream& out) const;
			
									
		protected:
		
			typedef std::map < size_t, std::vector < std::string > >  popStringsMap;
			typedef popStringsMap::iterator popStringsMapIt;
			typedef popStringsMap::const_iterator popStringsMapConstIt;
		
			MultiLociPolyTable();
			
			explicit MultiLociPolyTable(
			 						const MultiLociPolyTable& other);
			
			MultiLociPolyTable& operator=(MultiLociPolyTable tmp);
			
			static double getD (const boost::shared_ptr < Sequence::PolyTable >& simptr);
			
			static double getThetaPi (const boost::shared_ptr < Sequence::PolyTable >& simptr);
			
			// return true if cached differences are considered to be up to date
			bool checkCachedDifferences() const;
			
			void recalcCachedDifferences() const;
			
			popStringsMap& fillPopStringsMap(popStringsMap& mp) const;
			
			/* I originally thought that differences were site by site,
			 * but now I think that they are calculated looking at
			 * an entire string, based on comparing my results with
			 * Arlequin */
			static size_t countDifferences(const std::string& str, 
					const std::vector < std::string >& strs );
			
			void countAllDiffsBetweenPops(
					const MultiLociPolyTable::popStringsMap& mp) const;
							
			static size_t countDiffsBetweenPops(
					const MultiLociPolyTable::popStringsMap& mp,
					size_t pop1, size_t pop2);
			
			void countAllDiffsWithinPops(
					const MultiLociPolyTable::popStringsMap& mp) const;
							
			static size_t countDiffsWithinPop(
					const MultiLociPolyTable::popStringsMap& mp,
					size_t pop);
			
			// implement from base class
			virtual std::vector < int > getPopnsams() const;
			
			virtual size_t getNpops() const;
			
			virtual std::vector < size_t > getCountAllDiffsBetweenPops() const;
			
			virtual std::vector < size_t > getCountAllDiffsWithinPops() const;
			
						
		private:
		
			/*! \brief A pointer to the population structure.*/
			const boost::shared_ptr< const hudson_ms::PopulationStructure > _pop;
			
			/*! \brief The collection of the number of sites in each locus 
			represented by this.*/
			std::vector < size_t > _nsites;
			
			/*! \brief The collection of pointers to 
			Sequence::PolyTable held.*/
			std::vector< boost::shared_ptr <Sequence::PolyTable > > container;
			
			//cache the number of loci last time we did population comparison calcs
			mutable size_t cachedLociInPopCompCalcs;
			
			//cache counts of differences between populations
			mutable std::vector < size_t > cachedCountAllDiffsBetweenPops;
			
			//cache counts of difference within populations
			mutable std::vector < size_t > cachedCountAllDiffsWithinPops;
			
			
			
	};
	
	// does not need to be friend
	std::ostream& operator<< (std::ostream &out, const MultiLociPolyTable &m);

} // end namespace mct

#endif 

