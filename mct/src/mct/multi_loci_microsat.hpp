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
\brief MultiLociMicrosat declarations.
*/

#ifndef _INC_MULTI_LOCI_MICROSAT_H
#define _INC_MULTI_LOCI_MICROSAT_H


#include "genetic_pop_structure_analysable.hpp"
#include "mct_microsat_locus.hpp"

#include <mscplusplus/population_structure.hpp>

#include <boost/smart_ptr.hpp>

#include <vector>
#include <map>
#include <utility> // for pairs

namespace mct {
	
	
	/** @brief A class for an ordered collection of microsat repeats data.
	
	A MultiLociMicrosat is a \link aabc::DataObject DataObject\endlink
	type.
	 
	 A copy of data at each locus can be accessed as a pointer 
	 to a MCTMicroSatTable representing the sample sequence data at one 
	 locus.
	 
	 Loci are independent.
	 
	 Different loci can have different numbers of sites in them, but
	 we would normally expect one site per locus. 	 */
	 
	

	class MultiLociMicrosat : public mct::GeneticPopStructureAnalysable {
		public:
		
			/*! \brief Constructor.
			
			\param _p Pointer to the population structure for this.
			*/
			explicit MultiLociMicrosat(
				const boost::shared_ptr< const hudson_ms::PopulationStructure > _p);
			
			
			virtual ~MultiLociMicrosat();
			
			/*! \brief Get a pointer to an MCTMicrosatLocus representing
			data for a particular locus in the set.
			
			\param index The index of the locus to get 
			MCTMicrosatLocus data for.
			\return A pointer to the MCTMicrosatLocus 
			representing data at \a index. */
			boost::shared_ptr< MCTMicrosatLocus >
			at(const size_t index) const;
			
			/*! \brief Add a MCTMicrosatLocus pointer to the set.
			
			\param toAdd The MCTMicrosatLocus pointer to add. */
			void add(
				const boost::shared_ptr <MCTMicrosatLocus >& toAdd);
			
			/*! \brief Add a collection of site_seq to the set.
			
			\param toAdd The collection of site_seq to add.
			\pre The site_seq contains finite integers (no NaN values)*/
			void add(const std::vector < site_seq >& toAdd);
			
			/*! \brief Get the number of samples in the 
			 \link MCTMicrosatLocus MCTMicrosatLocuss\endlink 
			 in the collection.	*/
			size_t nsam() const;
			
			/*! \brief Get the population structure.*/
			const boost::shared_ptr< const hudson_ms::PopulationStructure > 
				populationStructure() const;
			
			
			/*! \brief Get the number of 
			 \link MCTMicrosatLocus MCTMicrosatLocuss\endlink
			 in the collection.	*/
			size_t nloci() const;
			
			/*! \brief Get the collection of the number of 
			sites in each locus .	*/
			std::vector < size_t > nsites() const;
			
			/*! \brief Get the total of the number of 
			sites in each locus .	*/
			size_t totalNsites() const;
			
			/*! \brief Return whether this contains no loci data.	*/
			bool empty() const;
			
			/*! \brief Get the 
			 \link SummaryStatistic SummaryStatistic\endlink
			 for this.
			 
			 \note Returns a SummaryStatistic made up of the average
			 repeats*/
			virtual boost::shared_ptr < const SummaryStatistic  > 
									getSummaryStatistic() const;
			
			/*! \brief Get the average repeats over all sites for this.
			 
			 One average is calculated for each individual in the sample,
			 taking the average over all sites in all loci.
			 
			 \note If there is no loci data (empty (), or nloci() == 0)
			 then this will return a vector of nsam nan values. 
			 
			\return a SummaryStatistic holding the average repeats, one
			 for each individual in the sample, in the same
			 order as the individuals are in the sample.*/
			virtual boost::shared_ptr < const SummaryStatistic  >
				getAverageRepeats() const;
			
			/*! \brief Get the variance-covariance matrix for the N
			sites in the collection.
			
			Each individual is an 'observation', a vector of N values; 
			the variance-covariance matrix is nsam x nsam, where 
			nsam is the number of individuals we have data for for 
			each site.
			 
			The SummaryStatistic returned 
			holds (nsam+1)*nsam/2 elements
			representing the top diagonal half of the 
			variance - covariance matrix between individuals in
			in the sample. 
			
			\f$ ( q_{1,1}, q_{1,2}, ... q_{1,nsam},
			q_{2,2}, x_{2,3}, .. q_{2,nsam},
			... 
			q_{nsam-1,nsam-1},q_{nsam-1,xnsam},
			q_{nsam,nsam} ) \f$
		  
			Where the element \f$ q_{ij} \f$ would have been 
			in row \f$ i \f$ and column \f$ j \f$ of the 
			variance-covariance
			matrix and 
			
			\f[
			q_{ij} = \frac{1}{N}\sum_{k=1}^N(x_{ik}-\bar{x}_i)(x_{jk}-\bar{x}_j)
			\f]
			
			where N is the total number of sites over all the loci
			and \f$ \bar{x}_i \f$ is the average repeat for individula 
			\f$ i \f$ in the sample 
			
			\note The variance-covariance matrix is calculated using the MLE 
			for the variance (covariance), i.e. division by N not N-1
			where N is the number of sites covered by the loci.
			\note If there is no loci data (empty(), or nloci() == 0)
			then this will return a vector of (nsam+1)nsam/2) 0 values.
			If there are no individuals the returned vector will be empty.  
			\return A SummaryStatistic representing the top diagonal half of the 
			variance - covariance matrix between individuals in
			in the sample.*/
			virtual boost::shared_ptr < const SummaryStatistic  >
				getCovariances() const;	
			
			/** @name Population comparison statistics.
			  
			 These functions return values for statistics 
			 used to compare and evaluate 'distances' between the
			 subpopulations with the population represented
			 by this.  
			  
			 Where within-subpopulation statistics are returned in an
			 ordered collection, each collection will be ordered so 
			 that if there are \f$ P > 0 \f$ subpopulations in this,
			  indexed 
			 \f$ j = 0, 1, 2, \ldots, P-1 \f$, and \f$ stat_{(j,j)} \f$
			 is the \f$ stat \f$ is the 
			 appropriate statistic calculated for the
			  subpopulation indexed \f$ j \f$ , then 
			  the returned collection will contain
			  \f$ stat_{(0,0)},stat_{(1,1)}, \ldots, stat_{(P-1,P-1)} \f$.
			
			
			 Where pairwise comparison statistics are returned in an
			 ordered collection, each collection will be ordered so 
			 that if there are \f$ P > 0 \f$ subpopulations in this,
			  indexed 
			 \f$ j = 0, 1, 2, \ldots, P-1 \f$, and \f$ stat_{(j,k)} \f$
			 , \f$ j \neq k \f$ is the \f$ stat \f$ is the 
			 appropriate statistic calculated between
			  subpopulations indexed \f$ j \f$ and \f$ k \f$, 
			  then the returned collection will contain
			  \f$ stat_{(0,0)},stat_{(0,1)},\ldots, stat_{(0,P-1)},stat_{(1,1)},stat_{(1,2)}, \ldots, stat_{(1,P-1)}, \ldots, stat_{(P-2,P-2)},stat_{(P-2,P-1)},stat_{(P-1,P-1)} \f$.
			.  If there is just one subpopulation in the population 
			represented by this, any collection of pairwise comparison
			statistics will be contain just one element. 
			
			Where between-subpopulation comparison statistics are returned in an
			 ordered collection, each collection will be ordered so 
			 that if there are \f$ P > 0 \f$ subpopulations in this,
			  indexed 
			 \f$ j = 0, 1, 2, \ldots, P-1 \f$, and \f$ stat_{(j,k)} \f$
			 , \f$ j \neq k \f$ is the \f$ stat \f$ is the 
			 appropriate statistic calculated between
			  subpopulations indexed \f$ j \f$ and \f$ k \f$, 
			  then the returned collection will contain
			  \f$ stat_{(0,1)},\ldots, stat_{(0,P-1)},stat_{(1,2)}, \ldots, stat_{(1,P-1)}, \ldots, stat_{(P-2,P-1)} \f$.
			.  If there is just one subpopulation in the population 
			represented by this, any collection of between-population comparison
			statistics will be empty. 
			
			Note that no exceptions or errors will be reported if
			between-subpopulation comparison statistics are requested
			when one or more subpopulations 
			contains no samples but that in this case any
			between-subpopulation comparison statistic
			involving one of the empty subpopulations should be treated
			with caution.  (No errors are reported to allow scenarios 
			that include one of more subpopulations from which no
			samples are taken but where past demographic events, including
			migration and population joins/splits link those subpopulations
			with the subpopulations from which samples have been taken.)
			 */
			//@{
			
			/*! \brief Get the average average repeats for each sub-population.
			  
			 Returned SummaryStatistic contains one entry for every 
			 sub-population, and each entry is the average over all sites 
			 in all loci of the average repeat size over all individuals in 
			 in that sub-population.
			 
			 If the number of sites is 0 the entry for each subpopulation
			 will be 0.0.  If any subpopulation has no samples, the average
			 repeat size for that subpopulation is 0.0.
			 
			 \return a SummaryStatistic holding the population average
			 repeats, one
			 for each population in the sample, in the same
			 order as the subpopulations are in the population.
			 \pre There is at least one sub-population containing at 
			 least one individual.*/
			virtual boost::shared_ptr < const SummaryStatistic  > 
							averageAverageRepeatSizeWithinPop() const;
			
			/*! \brief Get Goldstein's Delta Mu distance between
			  each pair of sub-populations i, j (i < j).\
			  
			 Each entry is the average over all sites 
			 in all loci of the site-by-site delta-mu distances 
			 for that pair of sub-populations.
			 
			 For each site, the delta-mu distance between a pair of 
			 subpopulations is the square of the difference in the average
			 repeat size at that site between the two subpopulations, ie
			 
			 \f$ \left( \mu_A - \mu_B \right)^2 \f$
			 
			 \note Goldstein's Delta Mu Squared distance is related to 
			 Slatkin's Rst distance.  
			 Delta Mu Distance = ASD - within population variances
			  
			 The within-population variances removed here are the 
			 biased variance estimates, as in Goldstein's formula. See
			 http://hpgl.stanford.edu/projects/microsat/dgboxc.html
			 
			 The method should return the same values as if they were
			 calculated by Arlequin 3.5.
			 (see documentation for Arlequin 3.5, section 8.2.4.5, 
			 except in
			 cases where a subpopulation has sample size 0 when Arlequin
			 will abort with an error.
			 
			 \return a SummaryStatistic holding the average delta mu distances, one
			 for each pair of sub-populations.
			 \pre There are at least two sub-populations with at least
			 one individual in each and at least one
			 site containing data.*/
			virtual boost::shared_ptr < const SummaryStatistic  > 
						distanceGoldsteinDeltaMuSqBetweenPop() const;
			
			/*! \brief Get Slatkin's \f$ Rst \f$ distance between
			  each pair of sub-populations.
			  
			 For each site, 
			 * 
			 \f$ Rst = \frac{(\bar{S} - S_w)}{\bar{S}} \f$
			 * 
			 where $\f$ \bar{S} \f$ is the average squared difference
			 in repeat size between pairs of genes taken from the genes 
			 of all subpopulations and $\f$ S_w \f$ is the average
			 squared differnce in repeat size between pairs of genes
			 within populations.  See Slatkin, M., (1995), "A measure
			 of population subdivision based on microsatellite allele
			 frequencies", Genetics, 139, pp. 457-462 and Michalakis, Y. 
			 and Excoffier, L. (1996), "A generic estimation of population
			 subdivision using distances between alleles with specail
			 reference to microsatellite loci", Genetics, 142, pp:1061-1064.  
			 
			 When \f$ Rst \f$ is calculated over multiple sites
			 the variance components are averaged over loci before 
			 \f$ Rst \f$ is calculated (rather than averaging the
			 site-by-site \f$ Rst \f$'s).  See Goodman, S.J., (1997),
			 "Rst Calc: a collection of computer programs for calculating
			 estimates of genetic differentiation from microsatellite
			 data and determining their significance", Molecular Ecology, 
			 6, pp: 881-885, 
			 * 
			 
			 The result should be the same as when using Arlequin 3.5
			 to calculate 
			 between-population " \f$ Fst \f$ " values using the sum
			 of squared size difference as the distance metric 
			 * (see documentation for Arlequin 3.5, section 8.1.2.7.2,
			  except in
			 cases where a subpopulation has sample size 0 when Arlequin
			 will abort with an error. 
			 			 
			 \return a SummaryStatistic holding the between-population
			 \f$ Rst \f$ distances, one
			 for each pair of sub-populations.*/
			virtual boost::shared_ptr < const SummaryStatistic  >
						distanceSlatkinRstBetweenPop() const;
			
			
								
			/*! \brief Get the pairwise average sum of squared
			 differences between repeat sizes between
			  each pair of sub-populations.
			  
			 Each entry is the average over possible pairs (within
			 or between subpopulations, as required) of the sum of the 
			 squared differences in repeat sizes between the pairs.
			 * 
			 The result should be the same as the population average
			 pairwise distances calculated using Arlequin 3.5 when the 
			 distance metric is the sum of squared size differences 
			 (see documentation for Arlequin 3.5, section 8.1.2.7.2),
			 except in cases where a subpopulation has sample size 0 
			 when Arlequin will abort with an error.
			 * 
			 \return a SummaryStatistic holding the average sum
			 of squared differences, one
			 for each pair of sub-populations.*/
			
			virtual boost::shared_ptr < const SummaryStatistic  > 
						pairwiseAverageSumSqDiffs() const;
			
			/*! \brief Get Slatkin's \f$ Rst \f$ distance over
			 the set of subpopulations.
			  
			   
			 For each site, 
			 * 
			 \f$ Rst = \frac{(\bar{S} - S_w)}{\bar{S}} \f$
			 * 
			 where $\f$ \bar{S} \f$ is the average squared difference
			 in repeat size between pairs of genes taken from the genes 
			 of all subpopulations and $\f$ S_w \f$ is the average
			 squared differnce in repeat size between pairs of genes
			 within populations.  See Slatkin, M., (1995), "A measure
			 of population subdivision based on microsatellite allele
			 frequencies", Genetics, 139, pp. 457-462 and Michalakis, Y. 
			 and Excoffier, L. (1996), "A generic estimation of population
			 subdivision using distances between alleles with specail
			 reference to microsatellite loci", Genetics, 142, pp:1061-1064.  
			 
			 When \f$ Rst \f$ is calculated over multiple sites
			 the variance components are averaged over loci before 
			 \f$ Rst \f$ is calculated (rather than averaging the
			 site-by-site \f$ Rst \f$'s).  See Goodman, S.J., (1997),
			 "Rst Calc: a collection of computer programs for calculating
			 estimates of genetic differentiation from microsatellite
			 data and determining their significance", Molecular Ecology, 
			 6, pp: 881-885, 
			 * 
			 
			 The result should be the same as when using Arlequin 3.5
			 to calculate 
			 between-population " \f$ Fst \f$ " values using the sum
			 of squared size difference as the distance metric 
			 * (see documentation for Arlequin 3.5, section 8.1.2.7.2,
			  except in
			 cases where a subpopulation has sample size 0 when Arlequin
			 will abort with an error. 
			 			 
			 \return a SummaryStatistic holding the overall
			 \f$ Rst \f$ value over all subpopulations.*/
			virtual boost::shared_ptr < const SummaryStatistic  > 
						rst() const;
	
			
			//@}
			
			/*! \brief Print the repeats information to an output
			stream, separating the subpopulation data. 
			* 
			Each line is information over all loci for one individual in
			the sample.	*/
			std::ostream& printSubpopData(std::ostream& out) const;
			
			/*! \brief Get a string representation of this.*/
			virtual std::string toString() const;
					
		protected:
		
					
			MultiLociMicrosat();
			
			explicit MultiLociMicrosat(
			 						const MultiLociMicrosat& other);
			
			MultiLociMicrosat& operator=(MultiLociMicrosat tmp);
			
			// typedefs for calculating allele sharing measures
			typedef std::map < int, size_t > freqmap;
			typedef std::pair < int, size_t > freqmap_pair;
			typedef freqmap::iterator freqmap_it;
			typedef std::pair < freqmap_it, bool > freqmap_check;
		
			// typedefs for differences in repeat size (diff or not diff)
			typedef std::map < size_t, std::vector < site_seq > >  popDataMap;
			typedef popDataMap::iterator popDataMapIt;
			typedef popDataMap::const_iterator popDataMapConstIt;

					
			/* Get repeats data for individuals in a 
			 given sub-population only.
			 
			 Returns a collection of collections of ints.  Each internal 
			 collection is the repeats, site by site, for an individual 
			 in the given sub-population. The outer collection holds 
			 these for all individuals in the given sub-population.*/			
			std::vector < site_seq > 
					getSubPopLoci(
						const std::vector< site_seq >& data,
						size_t popIndex) const;
			
			/* Get averages over all sites for a given data
			 collection.
			 
			 Return an ordered collection of averages, the average at
			 the ith position being the average of the contents of
			 the ith inner collection in the input.
			 If in is empty, out is empty.  If an inner vector of in
			 is empty, the corresponding entry in out is 0.0 */
			static std::vector < double >& averageOverVecs(
					std::vector < double >& out,
					const std::vector < std::vector < double > >& in);

			static std::vector < double >& averageOverVecs(
					std::vector < double >& out,
					const std::vector < std::vector < int > >& in);

			/* Get ordered collection of average repeat sizes for each site 
			 for each population.
			  
			 In returned object, inner collections are average repeats 
			 over all individuals for each site for a population, i.e.
			 each inner collection is totalNsites long.  The ith postion in
			 jth inner collection is the average over all the individuals
			 in sub-population j of repeat sizes for site i.
			 
			 If any subpopulation has nsam == 0 the inner vector
			 for that subpop will consist of totalNsites 0's.
			  
			 There are as many 
			 inner collections as there are sub-populations and the 
			 overall container is the outer collection.  */
			std::vector < std::vector < double > > 
					averageRepeatSizeWithinPopBySite() const;
			
			/* delta mu distances between different pops by site.
			one outer vector per population pair i,j with i < j
			inner vectors are all  nsites long and each entry 
			is the squared difference of average repeat size
			for corresponding site in subpops i and j.
			
			If nsites == 0, each entry will be 0.0.
			If any subpop has nsam == 0, the entry for any pair including
			that subpop will be the squared average 
			repeat size for the other subpop in the pair
			
			This the 'biased' version estimator, 
			ie variances removed from ASD are biased var estimates
			see http://hpgl.stanford.edu/projects/microsat/dgboxc.html */
			std::vector < std::vector < double > > 
					distanceGoldsteinDeltaMuSqBetweenPopsBySite() const;
			
			
			/* a vector of average repeat size for each individual
			ie nsam entries in the vector, each the average repeat size
			for that individual.
			Returns an empty vector if there are no individuals.  
			Returns a vector of 0's if nsites == 0. */
			virtual std::vector < double> getAverageRepeatsVec() const;
			
			/* a vector of variance-covariance of repeat size within-between
			individuals for each individual
			ie the top half of the variance-covariance matrix including
			variances on diagonals (total length nsam*(nsam+1)/2).
			Returns an empty vector if there are no individuals.
			Returns a vector of nsam*(nsam+1)/2 0's if nsites == 0. */
			virtual std::vector < double > getCovariancesVec() const;
			
			/* a vector of npops entries, each entry being the 
			 average repeat size over all individuals over all sites in
			 all loci for the corresponding subpopulation.
			 * 
			 If nsites == 0, the vector is npop 0's.
			 If any subpopulation has nsam = 0 the corresponding 
			 average average repeat size for that subpopulation will be 0.*/
			virtual std::vector < double > 
					averageAverageRepeatSizeWithinPopVec() const;
			
			
			/* a vector of npops*(npops-1)/2 entries, each entry 
			 corresponding to a pair of sub-populations i,j, i < j
			 with value = to the average over all sites of the per-site 
			 delta mus squared for that pair of sub-populations
			 * 
			 If nsites == 0, the vector is all 0's.
			 If any subpopulation has nsam = 0 any pair including that
			 subpop will have deltamusquared = average of the squared
			 per site average repeats for the other subpopulation in the 
			 pair.*/
			virtual std::vector < double > 
					distanceGoldsteinDeltaMuSqBetweenPopVec() const;
			
			
			virtual double getRst() const;
			
			// one entry for each pair of pops i, j, j > i
			virtual std::vector < double > getBetweenPopRst() const;

			
			// return true if cached differences are considered to be up to date
			bool checkCachedDifferences() const;
			
			void recalcCachedDifferences() const;
			
			
			MultiLociMicrosat::popDataMap& fillPopDataMap(
						MultiLociMicrosat::popDataMap& mp) const;
			
			static size_t countDifferences(const site_seq& base, 
						const std::vector < site_seq >& against );
			
			void countAllDiffsBetweenPops(
					const MultiLociMicrosat::popDataMap& mp) const;
							
			static size_t countDiffsBetweenPops(
					const MultiLociMicrosat::popDataMap& mp,
					size_t pop1, size_t pop2);
			
			void countAllDiffsWithinPops(
					const MultiLociMicrosat::popDataMap& mp) const;
							
			static size_t countDiffsWithinPop(
					const MultiLociMicrosat::popDataMap& mp,
					size_t pop);
			
			// sum squared difference over all sites for base and against
			static int sumSqDifferencesTotalOverSites(
						const site_seq& base, 
						const std::vector < site_seq >& against );
			
			// one vector of the total over sites of the sum sq differences per site for every pair of pops
			std::vector < int > sumSqAllDiffsBetweenPops(
				const MultiLociMicrosat::popDataMap& mp) const;
			
			// one figure for sum of squared differences in repeat size
			// summed over all sites.
			// All 0 if < 2 individuals in pop
			int sumSqDiffsBetweenPopsTotalOverSites(
							const MultiLociMicrosat::popDataMap& mp,
							size_t pop1, size_t pop2) const;

			// one vector of the total over sites of the sum sq differences per site for every pop
			std::vector < int > sumSqAllDiffsWithinPops(
			const MultiLociMicrosat::popDataMap& mp) const;
			
			// one figure for sum of squared differences in repeat size
			// summed over all sites.
			// All 0 if < 2 individuals in pop
			int sumSqDiffsWithinPopTotalOverSites(
							const MultiLociMicrosat::popDataMap& mp,
							size_t pop) const;
			
			double calcPairwiseRst(
						const std::vector < int >& nsams,
						const std::vector < double >& wp_diff,
						const std::vector < int >& sumSqAllDiffsWithinPops,
						const std::vector < int >& sumSqAllDiffsBetweenPops,
						size_t pop1, size_t pop2) const;
			
			std::vector < double > getBetweenPopAverageSumSqDiffs() const;
				
			std::vector < double > getBetweenPopAverageSumSqDiffs(
				const MultiLociMicrosat::popDataMap& mp) const;
			
			std::vector < double > getCorrectedBetweenPopAverageSumSqDiffs() const;
			
			std::vector < double > getWithinPopAverageSumSqDiffs() const;
			
			std::vector < double > getWithinPopAverageSumSqDiffs(
					const MultiLociMicrosat::popDataMap& mp) const;
			
			std::vector < double > getPairwiseAverageSumSqDiffs() const;			
			
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
			
			/*! \brief The collection of microsat repeats data held.
			
			Each inner vector represents an individual, containing
			repeats data for that individual for each site in order.
			*/
			std::vector< site_seq > container;
			
			//cache the number of loci last time we did population comparison calcs
			mutable size_t cachedLociInPopCompCalcs;
			
			//cache counts of differences between populations
			mutable std::vector < size_t > cachedCountAllDiffsBetweenPops;
			
			//cache counts of difference within populations
			mutable std::vector < size_t > cachedCountAllDiffsWithinPops;
			
			#if(0)
			// not sure that this is correct for microsatellites
			
			/*! \brief Get Nei's \f$ D_a \f$ distance between
			  each pair of sub-populations.
			  
			 Returned SummaryStatistic contains one value for every 
			 pair of sub-populations (i.e. m(m-1) entries for
			 m > 1 sub-populations).
			 
			 If there are m sub-populations, the pairs are, in order
			 sub-pop 1 x sub-pop 2, .. , sub-pop 1 x sub-pop m-1,
			 ..., sub-pop m-1 x sub-pop m.
			  
			 Each entry is the average over all sites 
			 in all loci of the site-by site \f$ D_a \f$ distances 
			 for that pair of sub-populations.
			 
			 \return a SummaryStatistic holding the average \f$ D_a \f$ distances, one
			 for each pair of sub-populations.
			 \pre There are at least two sub-populations containing at 
			 least one individual each and at least one
			 site containing data.*/
			 virtual boost::shared_ptr < const SummaryStatistic  >
						distanceNeiDaBetweenPop() const;
			
			virtual std::vector < double > 
					distanceNeiDaBetweenPopVec() const;
			
			
			/* components for calculating allele sharing measures
			inner vectors are allele sharing measures for each site 
			for a pair of populations i,j i > j
			outer vector holds these over all pairs of populations.
			If either subpop has no samples then the allele sharing
			vector for that pair is nsites 0's.  If nsites == 0 then
			the inner vectors of the returned vector are all empty.
			If npop == 1, returned vector will be empty */
			std::vector < std::vector < double > > 
						alleleSharingBetweenPairedPopBySite() const;
			
			// min and max repeat sizes over all samples and sites
			// returns 0,0 if there are no pops or no sites
			std::pair <int, int > getMinMax() const;
			
			/* make a collection of frequency maps from a population
			If there are no sites then the vector is empty.
			If the subpop has no individual then the returned vector
			contains nsites empty maps*/
			std::vector < freqmap > makeFreqMap( 
			const std::vector < site_seq >& popsites ) const;
			
			
			#endif			
			
	};
	
	// does not need to be friend
	std::ostream& operator<< (std::ostream &out, const MultiLociMicrosat &m);

} // end namespace mct

#endif 


