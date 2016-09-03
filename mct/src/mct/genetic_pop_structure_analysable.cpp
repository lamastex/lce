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
\brief GeneticPopStructureAnalysable definitions.
*/

#include "genetic_pop_structure_analysable.hpp"


#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/lambda/lambda.hpp>

#include <algorithm> // for transform
#include <functional> // for binary and unary functions
#include <numeric> // for accumulate

//#define MYDEBUGPOPSTATS
//#define MYDEBUGPOPSTATS1
#if defined (MYDEBUGPOPSTATS) || defined (MYDEBUGPOPSTATS1)
	#include <iterator> // for ostream_iterator
#endif

#include <iostream>

using namespace mct;



GeneticPopStructureAnalysable::~GeneticPopStructureAnalysable()
{
	//std::cout << "Deleting GeneticPopStructureAnalysable" << std::endl;
		
}

boost::shared_ptr < const SummaryStatistic  > 
	GeneticPopStructureAnalysable::fst() const
{
	boost::shared_ptr < SummaryStatistic  > result
		(new SummaryStatistic( std::vector < double > (1, getFst()) ));
	return result;
}

boost::shared_ptr < const SummaryStatistic  > 
	GeneticPopStructureAnalysable::pairwiseFst() const
{
	boost::shared_ptr < SummaryStatistic  > result
		(new SummaryStatistic( getPairwiseFst() ));
	return result;
}

boost::shared_ptr < const SummaryStatistic  > 
	GeneticPopStructureAnalysable::betweenPopFst() const
{
	boost::shared_ptr < SummaryStatistic  > result
		(new SummaryStatistic( getBetweenPopFst() ));
	return result;
}

boost::shared_ptr < const SummaryStatistic  > 
	GeneticPopStructureAnalysable::betweenPopAverageDifferences() const
{
	boost::shared_ptr < SummaryStatistic  > result
		(new SummaryStatistic( getBetweenPopAverageDifferences() ));
	return result;
}

boost::shared_ptr < const SummaryStatistic  > 
	GeneticPopStructureAnalysable::correctedBetweenPopAverageDifferences() const
{
	boost::shared_ptr < SummaryStatistic  > result
		(new SummaryStatistic( getCorrectedBetweenPopAverageDifferences() ));
	return result;
}

boost::shared_ptr < const SummaryStatistic  > 
	GeneticPopStructureAnalysable::withinPopAverageDifferences() const
{
	boost::shared_ptr < SummaryStatistic  > result
		(new SummaryStatistic( getWithinPopAverageDifferences() ));
	return result;
}

boost::shared_ptr < const SummaryStatistic  > 
	GeneticPopStructureAnalysable::pairwiseAverageDifferences() const
{
	boost::shared_ptr < SummaryStatistic  > result
		(new SummaryStatistic( getPairwiseAverageDifferences() ));
	return result;
}							

double GeneticPopStructureAnalysable::getFst() const
{
	#ifdef MYDEBUGPOPSTATS
		std::cout << "\nIn getFst()" <<std::endl;
	#endif
	
	double fst = 0.0;
	
	std::vector < size_t > countAllDiffsWithinPops 
				= getCountAllDiffsWithinPops(); // pure virtual method
	
	std::vector < int > nsams = getPopnsams(); // pure virtual method
	
	std::vector < int > nsams_adj;
	for (std::vector < int >::iterator it = nsams.begin();
		it < nsams.end(); ++it) {
		if ((*it) > 0) nsams_adj.push_back(*it);
	}
	size_t npops = nsams_adj.size(); // how many pops have nsam > 0
		
	if ( npops > 1) {
		
		std::vector < double > wp_diff;
		wp_diff.assign(countAllDiffsWithinPops.begin(),
									countAllDiffsWithinPops.end());
		
		boost::lambda::placeholder1_type a1;
		boost::lambda::placeholder2_type n1;
		
		std::transform(wp_diff.begin(), wp_diff.end(), 
						nsams.begin(), wp_diff.begin(),
						boost::lambda::bind(GeneticPopStructureAnalysable::protectedDivision, a1, n1)); 
		
		#ifdef MYDEBUGPOPSTATS
			std::cout << "\nwp_diffs are " <<std::endl;
			{
				std::ostream_iterator<double> out_it (std::cout,"\t");
				copy ( wp_diff.begin(), wp_diff.end(), out_it );
				std::cout << std::endl;
			} 
		#endif
		
		size_t totnsam = 0;
		totnsam = std::accumulate(nsams_adj.begin(), nsams_adj.end(), totnsam);
		
		//work out SSD(WP) and the df N-P for this
		
		double SSD_wp = std::accumulate(wp_diff.begin(), wp_diff.end(), 0.0);
			
		//hence get MSD(WP) == SSD(WP)/N-P which is sigmasq(wp)
		double MSD_wp = protectedDivision(SSD_wp, (totnsam - npops));
		double sigmasq_wp = MSD_wp;
		
		#ifdef MYDEBUGPOPSTATS
			std::cout << "\nSSD_wp = " << SSD_wp << " totnsam = " << totnsam << " npops = " << npops << std::endl;
			std::cout << "MSD_wp = SSD_wp/(totnsam - npops) = " << MSD_wp <<std::endl;
			std::cout << "sigmasq_wp = " << sigmasq_wp <<std::endl;
		#endif
		
		//work out SSD(total) using accumulations of pop comparisons
		
		std::vector < size_t > countAllDiffsBetweenPops 
					= getCountAllDiffsBetweenPops(); // pure virtual method
		
		double SSD_total = 0.0;
		SSD_total = std::accumulate(countAllDiffsBetweenPops.begin(),
									countAllDiffsBetweenPops.end(),
									SSD_total);
		SSD_total = std::accumulate(countAllDiffsWithinPops.begin(),
									countAllDiffsWithinPops.end(),
									SSD_total);
		#ifdef MYDEBUGPOPSTATS
			std::cout << "\nSSD_total before div = " 
			<< SSD_total << " totnsam = " << totnsam << std::endl;
		#endif
		
		SSD_total /= totnsam;
		
		#ifdef MYDEBUGPOPSTATS
			std::cout << "\nSSD_total / totnsam = " << SSD_total << std::endl;
		#endif
		
		//hence work out SSD(AP) = SSD(total) - SSD(WP) and df P-1
		double SSD_ap = SSD_total - SSD_wp;
		
		#ifdef MYDEBUGPOPSTATS
			std::cout << "\nSSD_ap = SSD_total - SSD_wp = " << SSD_ap << std::endl;
			std::cout << "npops - 1 = " << (npops-1) << std::endl;
		#endif
		
		//hence get MSD(AP) = SSD(AP)/P-1
		double MSD_ap = SSD_ap/(npops - 1);
		
		#ifdef MYDEBUGPOPSTATS
			std::cout << "\nMSD_ap = SSD_ap/(npops - 1) = " << MSD_ap << std::endl;
		#endif
		
		// work out n = (N - sum_over_pops(nsam^2/N))/(P-1)
		int sum = 0;	
		boost::lambda::placeholder1_type p1;
		sum = std::inner_product ( nsams_adj.begin(), nsams_adj.end(), 
									nsams_adj.begin(), sum);
		
		#ifdef MYDEBUGPOPSTATS
			std::cout << "\nsum of squared subpop sizes  = " << sum << std::endl;
		#endif
		
		double n = (totnsam - sum/static_cast<double>(totnsam))/(static_cast<double>(npops-1)); 		
		
		#ifdef MYDEBUGPOPSTATS
			std::cout << "n = (totnsam - sum/totnsam)/(npops-1) = " << n << std::endl;
		#endif
		
		// hence work out sigmasq(ap) = (MSD(AP) - sigmasq(wp))/n
		double sigmasq_ap = protectedDivision(MSD_ap - sigmasq_wp, n);
		
		#ifdef MYDEBUGPOPSTATS
			std::cout << "\nsigmasq_ap = (MSD_ap - sigmasq_wp)/n = " << sigmasq_ap << std::endl;
		#endif
		
		// hence work out sigmasq(total) = sigmasq(wp) + sigmasq(ap)
		double sigmasq_total = sigmasq_wp + sigmasq_ap;
		
		#ifdef MYDEBUGPOPSTATS
			std::cout << "sigmasq_total = sigmasq_wp + sigmasq_ap = " << sigmasq_total << std::endl;
		#endif
		
		// hence work out Fst = sigmasq(ap)/sigmasq(total)
		fst = protectedDivision(sigmasq_ap, sigmasq_total); // update cache
		#ifdef MYDEBUGPOPSTATS
			std::cout << "fst = sigmasq_ap/sigmasq_total = " << fst << std::endl;
		#endif
	}
	return fst;
}

			
std::vector < double > GeneticPopStructureAnalysable::getPairwiseFst() const
{
	
	std::vector < double > betweenFsts = getBetweenPopFst();
	
	size_t npops = getNpops();
	
	std::vector < double > pairwise;
	
	//index into between pops is pop1*(npops-1) + pop2-1 - ((pop1+1)*pop1)/2
	if (npops > 0) {
		pairwise.reserve(npops + betweenFsts.size());
		
		size_t pop1 = 0;
		while (pop1 < npops) {
			for (size_t pop2 = pop1; pop2 < npops; ++pop2) {
				
				if (pop1 == pop2) pairwise.push_back(0.0);
				else pairwise.push_back(
					betweenFsts[pop1*(npops-1) + pop2-1 - ((pop1+1)*pop1)/2]);
			}
			pop1++;
		}
	}
	
	#ifdef MYDEBUGPOPSTATS
		std::cout << "\nEnd of getPairwiseFst, Fsts are " <<std::endl;
		{
			std::ostream_iterator<double> out_it (std::cout,"\t");
			copy ( pairwise.begin(), pairwise.end(), out_it );
			std::cout << std::endl;
		} 
	#endif
	
	return pairwise;
}

std::vector < double > GeneticPopStructureAnalysable::getBetweenPopFst() const
{
	std::vector < int > nsams = getPopnsams(); // pure virtual method
	
	size_t npops = nsams.size();
	
	std::vector < double > fsts;
	
	if ( npops > 1) {
	
		fsts.reserve(npops*(npops+1)/2);
		
		std::vector < size_t > countAllDiffsWithinPops 
					= getCountAllDiffsWithinPops();
		
		std::vector < double > wp_diff;
		wp_diff.assign(countAllDiffsWithinPops.begin(),
									countAllDiffsWithinPops.end());
		
		boost::lambda::placeholder1_type a1;
		boost::lambda::placeholder2_type n1;
		
		std::transform(wp_diff.begin(), wp_diff.end(), 
			nsams.begin(), wp_diff.begin(),
			boost::lambda::bind(GeneticPopStructureAnalysable::protectedDivision, a1, n1)); 
		
		size_t pop1 = 0;
		while (pop1 < npops-1) {
			
			for (size_t pop2 = pop1+1; pop2 < npops; ++pop2) {
				fsts.push_back( calcPairwiseFst(nsams, wp_diff,
								pop1, pop2) );
			}
			pop1++;
		}
		
		#ifdef MYDEBUGPOPSTATS
			std::cout << "\nEnd of getBetweenPopFst, Fsts are " <<std::endl;
			{
				std::ostream_iterator<double> out_it (std::cout,"\t");
				copy ( fsts.begin(), fsts.end(), out_it );
				std::cout << std::endl;
			} 
		#endif
	}
	return fsts;
}
			
std::vector < double > GeneticPopStructureAnalysable::getBetweenPopAverageDifferences() const
{
	#ifdef MYDEBUGPOPSTATS
		std::cout << "\nIn getBetweenPopAverageDifferences " <<std::endl;
	#endif
	
	std::vector < double > avs;
	std::vector < int > nsams = getPopnsams(); // pure virtual method
	size_t npops = nsams.size();
		
	if (npops > 1) {
		
		std::vector < size_t > countAllDiffsBetweenPops 
									= getCountAllDiffsBetweenPops();
		
		//subpopsizes
		avs.assign ( countAllDiffsBetweenPops.begin(), 
									countAllDiffsBetweenPops.end() );
		size_t pop1 = 0;
		std::vector < double >::iterator it = avs.begin();
		while (pop1 + 1 < npops) {
			for (size_t pop2 = pop1 + 1; pop2 < npops; ++pop2) {
				
				#ifdef MYDEBUGPOPSTATS
					std::cout << "pop1 = " << pop1 << " pop2 = " 
						<< pop2 << " total diffs = " << (*it)
							<< " nsams[pop1] = " << nsams[pop1] 
							<< " nsams[pop2] = " << nsams[pop2] << std::endl; 
				#endif
				
				if ((*it) > 0) (*it) /= ( nsams[pop1] * nsams[pop2] );
				
				#ifdef MYDEBUGPOPSTATS
					std::cout << "totaldiff / ( nsams[pop1] * nsams[pop2] ) = " << (*it) << std::endl; 
				#endif
				
				++it;
			}
			++pop1;
			
			//no checking on length avs vs length nsams
		}
	}
	#ifdef MYDEBUGPOPSTATS
		std::cout << "\nEnd of getBetweenPopAverageDifferences, average diffs are " <<std::endl;
		{
			std::ostream_iterator<double> out_it (std::cout,"\t");
			copy ( avs.begin(), avs.end(), out_it );
			std::cout << std::endl;
		} 
	#endif
	return avs;
}
			
std::vector < double > GeneticPopStructureAnalysable::getCorrectedBetweenPopAverageDifferences() const
{
	#ifdef MYDEBUGPOPSTATS
		std::cout << "\nIn getCorrectedBetweenPopAverageDifferences " <<std::endl;
	#endif
	
	std::vector < double > avs = getBetweenPopAverageDifferences();
	
	if (!avs.empty()) {
		
		std::vector < double > withinPopAvs = getWithinPopAverageDifferences();
		size_t npops = withinPopAvs.size();
	
		size_t pop1 = 0;
		std::vector < double >::iterator it = avs.begin();
		while (pop1 + 1 < npops) {
			for (size_t pop2 = pop1 + 1; pop2 < npops; ++pop2) {
				
				#ifdef MYDEBUGPOPSTATS
					std::cout << "pop1 = " << pop1 << " pop2 = " 
						<< pop2 << " av diffs = " << (*it) << std::endl;
					std::cout << " withinPopAvs[pop1] = " 
						<< withinPopAvs[pop1] << " withinPopAvs[pop2] = " 
						<< withinPopAvs[pop2] << std::endl; 
				#endif
				
				
				(*it) -= 0.5*( withinPopAvs[pop1] + withinPopAvs[pop2] );
				
				#ifdef MYDEBUGPOPSTATS
					std::cout << "av diff- 0.5*( withinPopAvs[pop1] + withinPopAvs[pop2] ) = " << (*it) << std::endl; 
				#endif
				
				++it;
			}
			++pop1;
			
			//no checking on length avs vs length withinPopAvs
		}
	}
	
	#ifdef MYDEBUGPOPSTATS
		std::cout << "\nEnd of getCorrectedBetweenPopAverageDifferences, average diffs are " <<std::endl;
		{
			std::ostream_iterator<double> out_it (std::cout,"\t");
			copy ( avs.begin(), avs.end(), out_it );
			std::cout << std::endl;
		} 
	#endif
	return avs;
	
}

// divisor is nsams[pop]*(nsams[pop]-1)/2			
std::vector < double > GeneticPopStructureAnalysable::getWithinPopAverageDifferences() const
{
	#ifdef MYDEBUGPOPSTATS
		std::cout << "\nIn getWithinPopAverageDifferences " <<std::endl;
	#endif
	
	//subpopsizes
	std::vector < int > nsams = getPopnsams();
	
	std::vector < size_t > countAllDiffsWithinPops = getCountAllDiffsWithinPops();
	
	std::vector < double > avs ( countAllDiffsWithinPops.begin(), 
								countAllDiffsWithinPops.end() );
	
	boost::lambda::placeholder1_type a1;
	boost::lambda::placeholder2_type n1;
	
	std::transform(avs.begin(), avs.end(), 
		nsams.begin(), avs.begin(),
		boost::lambda::bind( GeneticPopStructureAnalysable::protectedDivision, a1, (n1*(n1-1)/2.0) ) ); 
		
	
	#ifdef MYDEBUGPOPSTATS
		std::cout << "\nEnd of getWithinPopAverageDifferences, average diffs are " <<std::endl;
		{
			std::ostream_iterator<double> out_it (std::cout,"\t");
			copy ( avs.begin(), avs.end(), out_it );
			std::cout << std::endl;
		} 
	#endif
	
	return avs;
}

std::vector < double > GeneticPopStructureAnalysable::getPairwiseAverageDifferences() const			
{
	std::vector < double > withinPops = getWithinPopAverageDifferences();
	std::vector < double > betweenPops = getBetweenPopAverageDifferences();
	
	size_t npops = withinPops.size();
	std::vector < double > pairwise;
	
	//index into between pops is pop1*(npops-1) + pop2-1 - ((pop1+1)*pop1)/2
	if (npops) {
		pairwise.reserve(npops + betweenPops.size());
		
		size_t pop1 = 0;
		while (pop1 < npops) {
			for (size_t pop2 = pop1; pop2 < npops; ++pop2) {
				
				if (pop1 == pop2) pairwise.push_back(withinPops[pop1]);
				else pairwise.push_back(
					betweenPops[pop1*(npops-1) + pop2-1 - ((pop1+1)*pop1)/2]);
			}
			pop1++;
		}
	}
	return pairwise;
	
}

double GeneticPopStructureAnalysable::calcPairwiseFst(
						const std::vector < int >& nsams,
						const std::vector < double >& wp_diff,
						size_t pop1, size_t pop2) const
{
	double Fst = 0.0; 
	//std::vector < double > wp_diff = getWithinPopAverageDifferences();
	size_t nsam1 = nsams[pop1];
	size_t nsam2 = nsams[pop2];
	size_t totnsam = nsam1 + nsam2;
	
	#ifdef MYDEBUGPOPSTATS
		std::cout << "\nIn calcPairwiseFst, pop1 = " << pop1 << " pop2 = " << pop2 << std::endl;	
		std::cout << "nsam1 = " << nsam1 << " nsam2 = " << nsam2 << " totnsam = " << totnsam << std::endl;		 
	#endif
		
	
	if (nsam1 > 0 && nsam2 > 0) {
		
		size_t npops = 2; // we are always comparing 2 pops
		
		//work out SSD(WP) and the df N-P for this
		
		double SSD_wp = wp_diff[pop1] + wp_diff[pop2];;
			
		//hence get MSD(WP) == SSD(WP)/N-2 (P=2) which is sigmasq(wp)
		double MSD_wp = protectedDivision(SSD_wp, (totnsam - npops));
		double sigmasq_wp = MSD_wp;
		
		//work out SSD(total) using accumulations of pop comparisons
		
		double SSD_total = 0.0;
		
		std::vector < size_t > countAllDiffsBetweenPops 
								= getCountAllDiffsBetweenPops();
		std::vector < size_t > countAllDiffsWithinPops 
								= getCountAllDiffsWithinPops();
		
		SSD_total = countAllDiffsBetweenPops[pop1*(wp_diff.size()-1) + pop2-1 - ((pop1+1)*pop1)/2];
		
		#ifdef MYDEBUGPOPSTATS
			std::cout << "Index into countAllDiffsBetweenPops = [pop1*(wp_diff.size()-1) + pop2-1 - ((pop1+1)*pop1)/2] = " 
				<< (pop1*(wp_diff.size()-1) + pop2-1 - ((pop1+1)*pop1)/2) << std::endl;	
			std::cout << "countAllDiffsBetweenPops[] = " << SSD_total << std::endl;		 
		#endif
		
		SSD_total += countAllDiffsWithinPops[pop1]
						+ countAllDiffsWithinPops[pop2];
		
		#ifdef MYDEBUGPOPSTATS
			std::cout << "countAllDiffsWithinPops[pop1] + countAllDiffsWithinPops[pop2] = " 
				<< (countAllDiffsWithinPops[pop1] + countAllDiffsWithinPops[pop2]) << std::endl;	
			std::cout << "SSD_total = " << SSD_total << std::endl;		 
		#endif
		
		SSD_total /= totnsam;
		
		#ifdef MYDEBUGPOPSTATS
			std::cout << "SSD_total / totnsam = " << SSD_total << std::endl;		 
		#endif
		
		//hence work out SSD(AP) = SSD(total) - SSD(WP) and df P-1
		double SSD_ap = SSD_total - SSD_wp;
		
		//hence get MSD(AP) = SSD(AP)/P-1 where P-1 = 2-1 = 1;
		double MSD_ap = SSD_ap;
		
		// work out n = (N - sum_over_pops(nsam^2/N))/(P-1) remembering P-1 = 1
		double n = totnsam - (nsam1*nsam1 + nsam2*nsam2)/static_cast<double>(totnsam); 		
		
		// hence work out sigmasq(ap) = (MSD(AP) - sigmasq(wp))/n
		double sigmasq_ap = protectedDivision(MSD_ap - sigmasq_wp, n);
		
		// hence work out sigmasq(total) = sigmasq(wp) + sigmasq(ap)
		double sigmasq_total = sigmasq_wp + sigmasq_ap;
		
		// hence work out Fst = sigmasq(ap)/sigmasq(total)
		Fst = protectedDivision(sigmasq_ap, sigmasq_total); // 
	}
	
	#ifdef MYDEBUGPOPSTATS
		std::cout << "Fst = " << Fst << std::endl;		 
	#endif
	return Fst;
}

double GeneticPopStructureAnalysable::protectedDivision(double num, double denom)
{
	if (num != 0.0) num /= denom;
	return num; 
}





