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
\brief GeneticPopStructureAnalysable declarations.
*/

#ifndef _INC_GENETIC_POP_STRUCT_ANALYSABLE_H
#define _INC_GENETIC_POP_STRUCT_ANALYSABLE_H


#include "data_object.hpp"

#include <vector>


namespace mct {
	
	
	/** @brief An abstract class representing 
	a type representing genetic data samples and capable of producing 
	general statistics for population 
	structure analysis:
	
	<ul>
	<li>An overall \f$ Fst \f$ statistic.</li>
	<li>Between population and pairwise \f$ Fst \f$ statistics for 
	each possible 2 out of all subpopulations represented.</li>
	<li>Population average pairwise differences (ie between and within
	populations) and correct average pairwise differences.</li>
	</ul>
	
	The underlying genetic sample data represented in concrete 
	subtypes of the %GeneticPopStructureAnalysable type may be of many 
	varieties (nucleotide 
	sequences, polymorphic nucleotide sites, single nucleotide 
	polymorphisms, microsatellite repeats etc): the only requirement
	is that these concrete types should be capable of producing the 
	required statistics for population structure analysis.
	*/
	
	class GeneticPopStructureAnalysable : public aabc::DataObject {
		public:
		
			
			virtual ~GeneticPopStructureAnalysable();
			
			
			/** @name Population comparison statistics.
			  
			 These functions return values for statistics 
			 used to compare and evaluate 'distances' between the
			 subpopulations with the population represented
			 by this.  
			  
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
			* 
			Note also that because the overall \f$ Fst \f$ value 
			is calculated using only \b non-empty subpopulations, the
			presence of empty sub-populations should not distort this
			statistic.  
			 */
			
			//@{
			
			/*! \brief Get the value of the \f$ Fst \f$ statistic for this.
			 * 
			The \f$ Fst \f$ is only calculated over subpopulations with 
			<b> sample sizes > 0 </b>.  This avoids distortion from the degrees
			of freedom using the number of subpopulations P 
			in the calculations:  this method will take P to be 
			the number of subpopulations with sample size > 0.
			 
			 If there are < 2 subpopulations with sample size > 0 in this,
			 the returned value will be 0.0.
			 
			 \f$ Fst \f$ is calculated over all loci and all subpopulations 
			 (with sample size > 0) represented in this.
			 
			 The method should return the same value as if it were
			 calculated by Arlequin 3.5 using the pairwise difference 
			 distance metric (see documentation for Arlequin 3.5, section 8.2.1.1 
			 Haplotypic data, one group of populations), except in
			 cases where a subpopulation has sample size 0 when Arlequin
			 will abort with an error. */
			virtual boost::shared_ptr < const SummaryStatistic  > fst() const;
			
			/*! \brief Get a collection of the pairwise subpopulation
			 comparison \f$ Fst \f$ values for the subpopulations represented
			 in this.
			 
			The \f$ Fst \f$ value for a pair of the same subpopulation
			(ie the diagonals of a pairwise Fst matrix) is 0.0.
			
			Otherwise, if the sample size is 0 for any subpopulation in a pair,
			of different subpopulations, then the \f$ Fst \f$ value
			 for that pair will be 0.0.
			  
			Each \f$ Fst \f$ for pairs of different subpopulations
			with both sample sizes > 0
			is calculated over all loci represented in this.
			 
			 The method should return the same values as if they were
			 calculated by Arlequin 3.5. using the pairwise difference
			 distance metric
			 (see documentation for Arlequin 3.5, section 8.2.1.1 
			 Haplotypic data, one group of populations), except in
			 cases where a subpopulation has sample size 0 when Arlequin
			 will abort with an error.*/
			virtual boost::shared_ptr < const SummaryStatistic  > 
												pairwiseFst() const;
			
			/*! \brief Get a collection of the between subpopulation
			 comparison \f$ Fst \f$ values for the subpopulations represented
			 in this.
			 
			The between subpopulation collection does not include comparison
			of a subpopulation with itself.
			
			If the sample size is 0 for any subpopulation in a pair,
			of different subpopulations, then the \f$ Fst \f$ value
			 for that pair will be 0.0.
			  
			Each \f$ Fst \f$ for pairs of different subpopulations
			with both sample sizes > 0
			is calculated over all loci represented in this.
			 
			 The method should return the same values as if they were
			 calculated by Arlequin 3.5.
			 (see documentation for Arlequin 3.5, section 8.2.1.1 
			 Haplotypic data, one group of populations), except in
			 cases where a subpopulation has sample size 0 when Arlequin
			 will abort with an error.*/
			virtual boost::shared_ptr < const SummaryStatistic  > 
												betweenPopFst() const;
			
			/*! \brief Get a collection of the between-subpopulation 
			 average differences
			 (Nei's raw \f$ D \f$ statistic) for subpopulations 
			 represented in this.
			 
			Each average difference is calculated over all loci represented in this.
			 
			 The method should return the same values as if they were
			 calculated by Arlequin 3.5.
			 (see documentation for Arlequin 3.5, section 8.2.4.4 .*/
			virtual boost::shared_ptr < const SummaryStatistic  > 
									betweenPopAverageDifferences() const;
			
			/*! \brief Get a collection of the corrected between-subpopulation
			  average differences (Nei's \f$ D_A \f$ statistic) for
			  subpopulations represented in this.
			 
			Each correct average difference is calculated over all loci represented in this.
			 
			 The method should return the same values as if they were
			 calculated by Arlequin 3.5.
			 (see documentation for Arlequin 3.5, section 8.2.4.4 .*/
			virtual boost::shared_ptr < const SummaryStatistic  > 
									correctedBetweenPopAverageDifferences() const;
			
			
			/*! \brief Get a collection of the within-subpopulation average 
			 differences for the subpopulations represented in this.
			 
			Each average difference is calculated over all loci represented in this.
			 
			 The method should return the same values as if they were
			 calculated by Arlequin 3.5.
			 (see documentation for Arlequin 3.5, section 8.2.4.4 .*/
			virtual boost::shared_ptr < const SummaryStatistic  > 
									withinPopAverageDifferences() const;
			
			/*! \brief Get a collection of the pairwise-subpopulation average 
			 differences for the subpopulations represented in this.
			 
			Each pairwise average difference is calculated over all loci 
			represented in this.
			 
			 The method should return the same values as if they were
			 calculated by Arlequin 3.5.
			 (see documentation for Arlequin 3.5, section 8.2.4.4 .*/
			virtual boost::shared_ptr < const SummaryStatistic  > 
									pairwiseAverageDifferences() const;
			
			//@}
									
		protected:
		
			virtual std::vector < int > getPopnsams() const = 0;
			
			virtual size_t getNpops() const = 0 ;
			
			virtual std::vector < size_t > getCountAllDiffsBetweenPops() const = 0;
			
			virtual std::vector < size_t > getCountAllDiffsWithinPops() const = 0;
		
			virtual double getFst() const;
			
			virtual std::vector < double > getPairwiseFst() const;
			
			virtual std::vector < double > getBetweenPopFst() const;
			
			virtual std::vector < double > getBetweenPopAverageDifferences() const;
			
			virtual std::vector < double > getCorrectedBetweenPopAverageDifferences() const;
			
			virtual std::vector < double > getWithinPopAverageDifferences() const;
			
			virtual std::vector < double > getPairwiseAverageDifferences() const;
			
			virtual double calcPairwiseFst(
						const std::vector < int >& nsams,
						const std::vector < double >& wp_diff,
						size_t pop1, size_t pop2) const;
			
			static double protectedDivision(double num, double denom);
			
						
		private:
		
			
			
			
	};
	
	
} // end namespace mct

#endif 

