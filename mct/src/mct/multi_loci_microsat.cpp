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
\brief MultiLociMicrosat definitions.
*/

#include "multi_loci_microsat.hpp"
#include "poly_data_sequence_maker.hpp"

#include <boost/lambda/bind.hpp>
#include <boost/lambda/lambda.hpp>

#include <gsl/gsl_math.h> // GSL_NAN and gsl_isnan()

#include <sstream>  // to be able to manipulate strings as streams
#include <algorithm> // for transform
#include <numeric> // for accumulate
#include <cmath> // squareroot
#include <iterator> // for ostream_iterator
#include <stdexcept>
#include <cassert>

//#define MYDEBUG
//#define POPDEBUG
//#define MYDEBUGPOPSTATS
#if defined (MYDEBUG) || defined (POPDEBUG) || defined (MYDEBUGPOPSTATS)
	#include <iostream>
	
#endif

using namespace mct;

MultiLociMicrosat::MultiLociMicrosat(
			const boost::shared_ptr< const hudson_ms::PopulationStructure > _p)
		: 	_pop(_p),
			container(std::vector < site_seq > (nsam())) 
{
	if ( !(_pop) ) {
		throw std::invalid_argument(
			"MultiLociMicrosat::MultiLociMicrosat(...) :\nPopulation pointer null" );
	}
	if (_pop->empty() ) {
		throw std::runtime_error(
			"MultiLociMicrosat::MultiLociMicrosat(...) :\nPopulation empty" );
	}
}


MultiLociMicrosat::~MultiLociMicrosat()
{
		
}

		
boost::shared_ptr <MCTMicrosatLocus > MultiLociMicrosat::at(
										const size_t index) const
{
	// return a MCTMicrosatLocus for a particular locus
	// make the vector of site_seqs and use this to make the MCTMicrosatLocus
	if (index > nloci()) {
		throw std::invalid_argument(
			"MultiLociMicrosat::at(const size_t) :\nindex >= number of loci");
	}
	
	// find the start and stop indexes
	size_t start = 0;
	size_t stop = 0;
	if ( index > 0 ) {
		for (size_t i = 0; i < index; ++i) {
			start += _nsites.at(i); // at will throw exception if out of range
		}
	}
	stop = start + _nsites.at(index); // i now indexes number of sites for locus at index
	
	std::vector < site_seq > seqs(nsam());
	
	// go through each of the nsam and take out these sites
	std::vector < site_seq >::iterator it = seqs.begin();
	
	for (std::vector < site_seq >::const_iterator cit = container.begin();
		cit < container.end(), it < seqs.end();
		++cit, ++it) {
			
		for (size_t i = start; i < stop; ++i) {
			it->push_back(cit->at(i));
		}
			
	}
	
	boost::shared_ptr< MCTMicrosatLocus > result(new MCTMicrosatLocus(seqs));
	return result;
}


/*
void MultiLociMicrosat::add(
				const boost::shared_ptr <MCTMicrosatLocus >& toAdd)
{
	// internally data is kept as a vector of site_seq
	if (toAdd->size() != _nsam) {
		throw MCTexception("nsam incompatible\n");
	}
	_nsites.push_back(toAdd->nsites());
	container.push_back(toAdd);
	
}
*/
void MultiLociMicrosat::add(const std::vector < site_seq >& toAdd) 
{
	// copy of container and _nsites
			
	std::vector < site_seq > temp = container;
	std::vector < size_t > tnsites = _nsites;
	
	if (toAdd.size() != nsam()) {
		throw std::invalid_argument(
			"MultiLociMicrosat::add(const std::vector < site_seq >&) :\nnsam incompatible");
	}
	if (!toAdd.empty() ) {
		size_t this_nsites = (toAdd.at(0)).size();
		
		if (this_nsites == 0) {
			throw std::invalid_argument(
				"MultiLociMicrosat::add(const std::vector < site_seq >&) :\nnsites = 0");
		}
		
		std::vector < site_seq >::iterator it = temp.begin();
		for (std::vector < site_seq >::const_iterator cit = toAdd.begin();
			cit < toAdd.end(), it < temp.end();
			++cit, ++it) {
			
			if (cit->size() == this_nsites) {
				it->insert( it->end(), cit->begin(), cit->end() );
				
			}
			else {
				throw std::invalid_argument(
					"MultiLociMicrosat::add(const std::vector < site_seq >&) :\nnsites uneven");
			}
		}
		tnsites.push_back(this_nsites);
		
		/* swap container and temp
		   this ensures that container is unchanged if an 
		   exception is thrown during transfer of data */
		 std::swap ( container, temp); 
		 std::swap ( _nsites, tnsites);
		 
		 #ifdef MYDEBUG
			std::cout << "Data for individuals held is " << std::endl;
			
			for (std::vector < site_seq >::iterator it = container.begin();
			it < container.end();
			++it) {
				std::cout << "[\t";
				copy(it->begin(), it->end(), 
				std::ostream_iterator< size_t >(std::cout, "\t"));
				std::cout << "]"<< std::endl;
			}
		#endif
		  
	}
}


bool MultiLociMicrosat::empty() const
{
	return _nsites.empty();
}

size_t MultiLociMicrosat::nsam() const
{
	return _pop->totalNsam();
}

const boost::shared_ptr< const hudson_ms::PopulationStructure > 
				MultiLociMicrosat::populationStructure() const
{
	return _pop;
}

std::vector < size_t > MultiLociMicrosat::nsites() const
{
	return _nsites;
}

size_t MultiLociMicrosat::totalNsites() const
{
	return std::accumulate(_nsites.begin(),
									_nsites.end(), 0);
}


size_t MultiLociMicrosat::nloci() const
{
	
	return _nsites.size();
}


boost::shared_ptr < const SummaryStatistic  > 
						MultiLociMicrosat::getSummaryStatistic() const
{
	return getAverageRepeats();
		
}


boost::shared_ptr < const SummaryStatistic  > 
				MultiLociMicrosat::getAverageRepeats() const
{
	boost::shared_ptr < SummaryStatistic  > result
		(new SummaryStatistic(getAverageRepeatsVec() ));
	return result;
}


/* returns summary stat 
   that is essentially a vector of (nsam+1)*nsam/2 elements representing the 
	top diagonal half of the variance covariance matrix where
	each individual is treated as an observation
	( (x[1][1], x[1][2], ... x[1][nsam]),
	(x[2][2], x[2][3], .. x[2][nsam]),
	.. 
	(x[nsam-1][nsam-1],x[nsam-1][xnsam]),
	(x[nsam][nsam]) ) */

boost::shared_ptr < const SummaryStatistic  >
				MultiLociMicrosat::getCovariances() const
{
	boost::shared_ptr < SummaryStatistic  > result
		(new SummaryStatistic(getCovariancesVec() ));
	return result;
}
	
boost::shared_ptr < const SummaryStatistic  >
		MultiLociMicrosat::averageAverageRepeatSizeWithinPop() const
{
	boost::shared_ptr < SummaryStatistic  > result
		(new SummaryStatistic(averageAverageRepeatSizeWithinPopVec() ));
	return result;
}



boost::shared_ptr < const SummaryStatistic  > 
	MultiLociMicrosat::distanceGoldsteinDeltaMuSqBetweenPop() const
{
	boost::shared_ptr < SummaryStatistic  > result
		(new SummaryStatistic( distanceGoldsteinDeltaMuSqBetweenPopVec() ));
	return result;
}

boost::shared_ptr < const SummaryStatistic  > 
		MultiLociMicrosat::distanceSlatkinRstBetweenPop() const
{
	boost::shared_ptr < SummaryStatistic  > result
		(new SummaryStatistic(getBetweenPopRst() ));
	return result;
}


boost::shared_ptr < const SummaryStatistic  > 
	MultiLociMicrosat::pairwiseAverageSumSqDiffs() const
{
	boost::shared_ptr < SummaryStatistic  > result
		(new SummaryStatistic( getPairwiseAverageSumSqDiffs() ));
	return result;
}							


boost::shared_ptr < const SummaryStatistic  > 
	MultiLociMicrosat::rst() const
{
	boost::shared_ptr < SummaryStatistic  > result
		(new SummaryStatistic( std::vector < double > (1, getRst()) ));
	return result;
}



std::ostream& MultiLociMicrosat::printSubpopData(std::ostream& out) const
{
	
	size_t npops = _pop->nSubPops();
	
	int counter = 0;
	for (size_t i = 0; i < npops; ++i) {
		
		std::vector < site_seq > this_data = getSubPopLoci(container, i);
		
		for (std::vector < site_seq >::iterator it = this_data.begin();
				it < this_data.end();
				++it) {
		
			out << "\t" << counter << "\t1\t"; // freq 1 for each
			
			std::ostream_iterator<int> out_it (out," ");
			copy ( it->begin(), it->end(), out_it );
			
			out << "\n"; 
			counter++;
		}
		
		out << std::endl; // new line after each pop
	}
	return out; // return by reference
}

std::string MultiLociMicrosat::toString() const
{
	std::ostringstream stm;
	
	stm << "lo.";
	for (size_t i = 0; i < _nsites.size(); ++i) {
		stm << '\t' << i;
		for (size_t j = 0; j < _nsites[i]-1; ++j) {
			stm << '\t';
		}
	}
	
	// for each nsam print his stuff, after a \t
	for (std::vector< std::vector< int  > >::const_iterator it 
									= container.begin();
			it < container.end(); ++it) {
		stm << "\n\t";	
		copy(it->begin(), it->end(), 
		std::ostream_iterator<int>(stm, "\t"));
				
	}
	stm << std::endl;
	
	return stm.str();
}


std::vector < site_seq > 
		MultiLociMicrosat::getSubPopLoci(
						const std::vector< site_seq >& data,
						size_t popIndex) const
{
	// population structure knows about pops
	std::vector < int > nsams = _pop->makePopnsams();
	size_t totnsam = _pop->totalNsam();
	
	assert (data.size() == totnsam);
	
	std::vector < site_seq >::const_iterator first_it 
							= data.begin();
	std::vector < site_seq >::const_iterator last_it 
							= data.end();
	advance( last_it, -(totnsam) );

	for (size_t i = 0; i < popIndex; ++i) {
		advance( first_it, nsams[i] );
		advance( last_it, nsams[i] );
	}
	size_t thisNsam = nsams[popIndex];
	advance( last_it, thisNsam );
	
	return std::vector < site_seq >(first_it, last_it);
}



std::vector < double >& MultiLociMicrosat::averageOverVecs(
		std::vector < double >& out,
		const std::vector < std::vector < double > >& in)
{
	if (!in.empty() ) {
		#ifdef POPDEBUG
			for (std::vector< std::vector< double > >::const_iterator it = in.begin();
						it < in.end(); ++it) {
							copy(it->begin(), it->end(), 
							std::ostream_iterator<double>(std::cout, "\t"));
							std::cout << std::endl;
						}
		
		#endif
		
		
		size_t len = in.at(0).size();
		if (len > 0 ) {
			
		
			for (std::vector < std::vector < double > >::const_iterator it = in.begin();
				it < in.end();
				++it ) {
				
				out.push_back( accumulate(it->begin(), it->end(), 0.0)/len );
								
			}
		}
		else out.push_back(0.0);
	}
	
	return out;
}

std::vector < double >& MultiLociMicrosat::averageOverVecs(
		std::vector < double >& out,
		const std::vector < std::vector < int > >& in)
{
	if (!in.empty() ) {
		#ifdef POPDEBUG
			for (std::vector< std::vector< double > >::const_iterator it = in.begin();
						it < in.end(); ++it) {
							copy(it->begin(), it->end(), 
							std::ostream_iterator<double>(std::cout, "\t"));
							std::cout << std::endl;
						}
		
		#endif
		
		
		size_t len = in.at(0).size();
		if (len > 0 ) {
			
		
			for (std::vector < std::vector < int > >::const_iterator it = in.begin();
				it < in.end();
				++it ) {
				
				out.push_back( accumulate(it->begin(), it->end(), 0.0)/len );
								
			}
		}
		else out.push_back(0.0);
	}
	
	return out;
}

// needed by Goldstein distance
std::vector < std::vector < double > > 
				MultiLociMicrosat::averageRepeatSizeWithinPopBySite() const
{
	std::vector < std::vector < double > > res;
	
	// population structure knows about pops
	std::vector < int > nsams = _pop->makePopnsams();
	//size_t totnsam = _pop->totalNsam();
	size_t npop = _pop->nSubPops();
	
	size_t ns = totalNsites();
	
	for (size_t i = 0; i < npop; ++i) {
	
		using namespace boost::lambda;
		
		std::vector < site_seq > thisPop = getSubPopLoci(container, i);
	
		if (ns == 0) res.push_back( std::vector < double > ());
		else {
			res.push_back( std::vector < double >( ns, 0.0) );
			
			size_t thisNsam = thisPop.size();
					
			if (thisNsam > 0) {
			
				for (std::vector < site_seq >::iterator it = thisPop.begin();
					it < thisPop.end();
					++it) {
						
						transform( res.back().begin(), res.back().end(), it->begin(), 
															res.back().begin(),
															(_1 + _2) );
					
				}
			
				transform( res.back().begin(), res.back().end(),  
												res.back().begin(),
												(_1/thisNsam ) );
			} // end if thisNsam > 0
			// if thisNsam == 0 the entry will be a vector of nsites 0's
		}
	}
	
	return res;
}


std::vector < std::vector < double > > 
		MultiLociMicrosat::distanceGoldsteinDeltaMuSqBetweenPopsBySite() const
{
	std::vector < std::vector < double > > res;
	
	#ifdef MYDEBUG
		std::cout << "In distanceGoldsteinDeltaMuSqBetweenPopsBySite()" << std::endl;
	#endif
	
	std::vector < std::vector < double > > avRepeats 
						= averageRepeatSizeWithinPopBySite();
	
	#ifdef MYDEBUG
		std::cout << "size of avRepeats is " << avRepeats.size() << std::endl;
	#endif
	
	
	for (size_t p1 = 0; p1 < avRepeats.size() - 1; ++p1 ) {
		for (size_t p2 = p1+1; p2 < avRepeats.size(); ++p2 ) {
			
			#ifdef MYDEBUG
				std::cout << "p1 = " << p1 << " p2 = " << p2 << std::endl;
			#endif
			res.push_back( avRepeats[p1] );
			std::vector < double > pop2 = avRepeats[p2];
			
			using namespace boost::lambda;
			transform( res.back().begin(), res.back().end(), 
						pop2.begin(), 
						res.back().begin(),
						(_1 - _2)*(_1 - _2) );
		}
	}
	
	#ifdef MYDEBUG
		std::cout << "at end size of res is " << res.size() << std::endl;
	#endif
	// if npop == 1, res will be empty
	return res;
}





std::vector < double> 
				MultiLociMicrosat::getAverageRepeatsVec() const
{
	std::vector < double> averages;
	
	size_t total_nsites = totalNsites();
	
	if (total_nsites > 0) {
		for (std::vector < site_seq >::const_iterator it =
											container.begin();
			it < container.end();
			++it) {
			
			averages.push_back ( std::accumulate(it->begin(),
									it->end(), 0.0)/total_nsites );
		}
	}
	else if (nsam() > 0) {
		averages = std::vector < double > (container.size(), 0.0);
	}

	return averages;
}

std::vector < double > MultiLociMicrosat::getCovariancesVec() const
{
	// each inner vector in the container
	// is a vector of repeats for one individual (one row)
	// there are nsam of these, for nsam observations (nsam rows)
	// each observation consists of N values, the sites covered
	
	std::vector < double > result;

	if (!container.empty()) {
	
		size_t cols = totalNsites(); // N
		
		#ifdef DEBUG
			size_t check_cols = container.empty() ? 0 : 
									(container.at(0)).size(); // N
			assert(cols == check_cols);
		#endif
		 
		size_t rows = nsam(); // nsam
		
		#ifdef DEBUG
			assert(rows == container.size());
		#endif
		
		std::vector < std::vector < double > > work(rows);
		
		if (cols > 0) {
			for (size_t i = 0; i < rows; ++i) {
				
				std:: vector < int > a = container.at(i);
				std:: vector < double > aa(
										cols, std::accumulate(a.begin(),
										a.end(), 0.0)/cols);
				// subtract values from means from values
				std::transform(aa.begin(), aa.end(), a.begin(),
								aa.begin(), std::minus< double >());
				
				work[i] = aa;
			}
				
			for (size_t i = 0; i < rows; ++i) {	
				for (size_t j = i; j < rows; ++j) {
					
					double v = inner_product((work[i]).begin(), (work[i]).end(),
					(work[j]).begin(), 0.0)/(cols); 
					
					result.push_back(v);
					
				}
			}
		}
		else if (rows > 0) result = std::vector < double > ((rows * (rows+1))/2, 0.0);
	}
	return result;
}

std::vector < double > MultiLociMicrosat::averageAverageRepeatSizeWithinPopVec() const
{
	std::vector < double > res;
		
	std::vector < std::vector < double > > avBySitesByPop 
									= averageRepeatSizeWithinPopBySite();
	// an inner vector could have only 0's in if nsam_pop = 0, or could be empty if nsites == 0
	
	res = averageOverVecs(res, avBySitesByPop);
	
	return res;
}


std::vector < double > 
	MultiLociMicrosat::distanceGoldsteinDeltaMuSqBetweenPopVec() const
{
	std::vector < double > res;
		
	std::vector < std::vector < double > > goldsteinBySite 
					= distanceGoldsteinDeltaMuSqBetweenPopsBySite();
	
	res = averageOverVecs(res, goldsteinBySite);
	
	return res;
}






double MultiLociMicrosat::getRst() const
{
	#ifdef MYDEBUGPOPSTATS
		std::cout << "\nIn getRst()" <<std::endl;
	#endif
	
	std::vector < int > nsams = getPopnsams();
	
	std::vector < int > nsams_adj;
	for (std::vector < int >::iterator it = nsams.begin();
		it < nsams.end(); ++it) {
		if ((*it) > 0) nsams_adj.push_back(*it);
	}
	size_t npops = nsams_adj.size(); // how many pops have nsam > 0
	double rst = 0.0;
	
	if ( npops > 1) {
		
		// do the map
		MultiLociMicrosat::popDataMap mp;
		fillPopDataMap(mp);
	
		std::vector < int > sumSqWithin
					= sumSqAllDiffsWithinPops(mp);
		
		std::vector < double > wp_diff;
		wp_diff.assign(sumSqWithin.begin(),
									sumSqWithin.end());
		
		
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
		
		std::vector < int > sumSqBetween
					= sumSqAllDiffsBetweenPops(mp);
		
		double SSD_total = 0.0;
		SSD_total = std::accumulate(sumSqBetween.begin(),
									sumSqBetween.end(),
									SSD_total);
		SSD_total = std::accumulate(sumSqWithin.begin(),
									sumSqWithin.end(),
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
		
		// hence work out Rst = sigmasq(ap)/sigmasq(total)
		rst = protectedDivision(sigmasq_ap, sigmasq_total); // update cache
		#ifdef MYDEBUGPOPSTATS
			std::cout << "rst = sigmasq_ap/sigmasq_total = " << rst << std::endl;
		#endif
	}
	return rst;
}


// one entry for each pair of pops i, j, j > i
std::vector < double > MultiLociMicrosat::getBetweenPopRst() const
{
	std::vector < int > nsams = getPopnsams();
	
	size_t npops = nsams.size();
	
	std::vector < double > rsts;
	
	if ( npops > 1) {
		
		// do the map
		MultiLociMicrosat::popDataMap mp;
		fillPopDataMap(mp);
	
		rsts.reserve(npops*(npops+1)/2);
		
		std::vector < int > sumSqWithin
					= sumSqAllDiffsWithinPops(mp);
		
		std::vector < int > sumSqBetween
					= sumSqAllDiffsBetweenPops(mp);
		
		std::vector < double > wp_diff;
		wp_diff.assign(sumSqWithin.begin(),
									sumSqWithin.end());
		
		boost::lambda::placeholder1_type a1;
		boost::lambda::placeholder2_type n1;
		
		std::transform(wp_diff.begin(), wp_diff.end(), 
			nsams.begin(), wp_diff.begin(),
			boost::lambda::bind(GeneticPopStructureAnalysable::protectedDivision, a1, n1)); 
		
		size_t pop1 = 0;
		while (pop1 < npops-1) {
			
			for (size_t pop2 = pop1+1; pop2 < npops; ++pop2) {
				rsts.push_back( calcPairwiseRst(nsams, wp_diff,
								sumSqWithin, sumSqBetween,
								pop1, pop2) );
			}
			pop1++;
		}
		
		#ifdef MYDEBUGPOPSTATS
			std::cout << "\nEnd of getBetweenPopRst, Rsts are " <<std::endl;
			{
				std::ostream_iterator<double> out_it (std::cout,"\t");
				copy ( rsts.begin(), rsts.end(), out_it );
				std::cout << std::endl;
			} 
		#endif
	}
	return rsts;
}


/* return true if cached differences are up to date*/
bool MultiLociMicrosat::checkCachedDifferences() const
{
	return (cachedLociInPopCompCalcs == _nsites.size());
	
}

/* Recalculate the the cached differences*/
void MultiLociMicrosat::recalcCachedDifferences() const // stats updated are mutable
{
	#ifdef MYDEBUGPOPSTATS
		std::cout << "\nIn recalcCachedDifferences" <<std::endl;
	#endif
	
	// recalculate and cache
	MultiLociMicrosat::popDataMap mp;
	fillPopDataMap(mp);
	
	// make sure we have up to date stats
	countAllDiffsWithinPops(mp);
	countAllDiffsBetweenPops(mp);
		
	cachedLociInPopCompCalcs = _nsites.size();

}


MultiLociMicrosat::popDataMap& MultiLociMicrosat::fillPopDataMap(
						MultiLociMicrosat::popDataMap& mp) const
{
	#ifdef MYDEBUGPOPSTATS
		std::cout << "In fillPopDataMap" << std::endl;
	#endif
	
	//get just the polymorphic sites in the data
	
	std::vector< site_seq > polyData = container;
	PolyDataSequenceMaker::reduceToPolyData(polyData);
	
	size_t totnsam = _pop->totalNsam();
	
	assert(polyData.size() == totnsam);
	
	if (totnsam > 0) {
	
		// population structure knows about pops
		std::vector < int > nsams = _pop->makePopnsams();
		size_t npops = nsams.size();
			
		size_t nPolySites = polyData.front().size();

		// divvie up by population
		
		#ifdef MYDEBUGPOPSTATS1
			std::cout << "npops = " << npops << std::endl;
				 
		#endif
		
		if (nPolySites == 0) {
			popDataMapIt mit = mp.begin();
			site_seq tmp; // empty vector
			
			for (size_t i = 0; i < npops; ++i) {
				
				// max efficiency inserting
				mit = mp.insert (mit, 
							std::pair< size_t, std::vector < site_seq > >
								(i, std::vector< site_seq >( nsams[i], tmp ) ) );
			}
		}

		else {
			
			std::vector < site_seq >::const_iterator first_it 
									= polyData.begin();
			std::vector < site_seq >::const_iterator last_it 
									= polyData.end();
			advance( last_it, -(totnsam-nsams[0]) ); // we know there is at least one pop

			popDataMapIt mit = mp.begin();
			
			for (size_t i = 0; i < npops; ++i) {
				
				#ifdef MYDEBUGPOPSTATS
					std::cout << "pop index " << i << " nsams[i] = " << nsams[i] << std::endl;
				#endif
				
				// after first pop, set the iterator around that ones strings
				if (i > 0 ) { 
					advance( first_it, nsams[i-1] );
					advance( last_it, nsams[i] );
				}
				
				// add entry to map for contents of for pop index i
				// max efficiency inserting
				
				mit = mp.insert (mit, 
							std::pair< size_t, std::vector < site_seq > >
								(i, std::vector <  site_seq >(first_it, last_it) ) );  
				#ifdef MYDEBUGPOPSTATS
				{
					std::cout << "map is now" << std::endl;
					for (MultiLociMicrosat::popDataMapIt debug_it = mp.begin();
						debug_it != mp.end(); ++debug_it) {
						
						for (std::vector < site_seq >:: iterator 
										dit = debug_it->second.begin();
										dit < debug_it->second.end();
										++dit ) {
						
							std::ostream_iterator< int > out_it (std::cout," ");
							copy ( dit->begin(), dit->end(), out_it );
							std::cout << std::endl;
						}
					}
				}	
				#endif
				
			}
		}
	}
	
	return mp; // return by reference
}

size_t MultiLociMicrosat::countDifferences(const site_seq& base, 
						const std::vector < site_seq >& against )
{
	#ifdef MYDEBUGPOPSTATS
	{
		std::cout << "In countDifferences, base = " << std::endl;
		{
			std::ostream_iterator<int> out_it (std::cout," ");
			copy ( base.begin(), base.end(), out_it );
			std::cout << std::endl;
		}
		
		std::cout << " and against = " << std::endl;
		for (std::vector < site_seq >::const_iterator it = against.begin();
				it < against.end();
				++it) {
			std::ostream_iterator<int> out_it (std::cout," ");
			copy ( it->begin(), it->end(), out_it );
			std::cout << std::endl;
		}
	}
	#endif 
	size_t diffCount = 0;
	
	if (!against.empty() ) {
		
		size_t nsites = base.size();
		
		for (std::vector < site_seq >::const_iterator it = against.begin();
				it < against.end();
				++it) {
				
				site_seq tmp = *it;
				
				using namespace boost::lambda;
				
				// transform tmp to tmp - base
				std::transform(tmp.begin(), tmp.end(), 
							base.begin(), tmp.begin(), 
							(_1 -_2) );
				
				diffCount += (nsites - std::count_if( tmp.begin(), tmp.end(), (_1 == 0) ));						
				
		}
		
	}	
	return diffCount;
}

void MultiLociMicrosat::countAllDiffsBetweenPops(
				const MultiLociMicrosat::popDataMap& mp) const
{
	#ifdef MYDEBUGPOPSTATS
		std::cout << "\nIn countAllDiffsBetweenPops" <<std::endl;
	#endif
	
	size_t npops = mp.size();
	
	std::vector < size_t > countDiffs;
	
	if (npops > 1) {
		countDiffs.reserve(npops*(npops+1)/2);
		
		size_t pop1 = 0;
		while (pop1 + 1 < npops) {
			for (size_t pop2 = pop1 + 1; pop2 < npops; ++pop2) {
				countDiffs.push_back( countDiffsBetweenPops(mp, pop1, pop2) );
			}
			pop1++;
		}
	}
	// update cached stats
	countDiffs.swap(cachedCountAllDiffsBetweenPops);
	
	#ifdef MYDEBUGPOPSTATS
		std::cout << "\nEnd of countAllDiffsBetweenPops, diffs are " <<std::endl;
		{
			std::ostream_iterator<size_t> out_it (std::cout,"\t");
			copy ( cachedCountAllDiffsBetweenPops.begin(), cachedCountAllDiffsBetweenPops.end(), out_it );
			std::cout << std::endl;
		} 
	#endif
	
}

size_t MultiLociMicrosat::countDiffsBetweenPops(
							const MultiLociMicrosat::popDataMap& mp,
							size_t pop1, size_t pop2)
{
	#ifdef MYDEBUGPOPSTATS
		std::cout << "\nIn countDiffsBetweenPops, pop1 = " << pop1 << " pop2 = " << pop2 <<std::endl;
		
	#endif
	//iterator to pair for first pop
	popDataMapConstIt cit = mp.find(pop1);
	//get the pop 2 site_seqs out once and work with the local copy
	std::vector < site_seq > seq2 = mp.find(pop2)->second;
	
	size_t countDiff = 0;
	for (std::vector < site_seq >::const_iterator it = cit->second.begin();
			it < cit->second.end();
			++it) {
		countDiff += countDifferences(*it, seq2);
	}
	
	#ifdef MYDEBUGPOPSTATS
		std::cout << "countDiff = " << countDiff <<std::endl;		 
	#endif
	
	return countDiff;
}


void MultiLociMicrosat::countAllDiffsWithinPops(
			const MultiLociMicrosat::popDataMap& mp) const
{
	#ifdef MYDEBUGPOPSTATS
		std::cout << "\nIn countAllDiffsWithinPops" <<std::endl;
	#endif
	
	size_t npops = mp.size();
	
	std::vector < size_t > countDiffs;
	
	if (npops > 0) {
		countDiffs.reserve(npops);
		
		for (size_t i = 0; i < npops; ++i) {
			countDiffs.push_back( countDiffsWithinPop(mp, i) );
		}
	}
	
	// update cached stats
	countDiffs.swap(cachedCountAllDiffsWithinPops);
	#ifdef MYDEBUGPOPSTATS
		std::cout << "\nEnd of countAllDiffsWithinPops, diffs are " <<std::endl;
		{
			std::ostream_iterator<size_t> out_it (std::cout,"\t");
			copy ( cachedCountAllDiffsWithinPops.begin(), cachedCountAllDiffsWithinPops.end(), out_it );
			std::cout << std::endl;
		} 
	#endif
}

size_t MultiLociMicrosat::countDiffsWithinPop(
							const MultiLociMicrosat::popDataMap& mp,
							size_t pop)
{
	#ifdef MYDEBUGPOPSTATS
		std::cout << "\nIn countDiffsWithinPop, pop = " << pop <<std::endl;
	#endif
	//get the data out once and work with the local copy
	std::vector < site_seq > popdata = mp.find(pop)->second;
	
	size_t countDiff = 0;
	size_t sz = popdata.size();
	
	while (sz > 1) {
		site_seq lastone= popdata.back(); // last one 
		// remove that last one
		popdata.resize( sz-1 );
		
		sz--;
		
		countDiff += countDifferences(lastone, popdata);
	}
	
	#ifdef MYDEBUGPOPSTATS
		
		std::cout << "countDiff = " << countDiff <<std::endl;		 
	#endif
		
	return countDiff;
}

// distance method square of difference of repeat size, site by site, eg for Slatkin


// sum squared difference over all sites for base and against
int MultiLociMicrosat::sumSqDifferencesTotalOverSites(
						const site_seq& base, 
						const std::vector < site_seq >& against )
{
	#ifdef MYDEBUGPOPSTATS
	{
		std::cout << "In sumSqDifferences, base = " << std::endl;
		{
			std::ostream_iterator<int> out_it (std::cout," ");
			copy ( base.begin(), base.end(), out_it );
			std::cout << std::endl;
		}
		
		std::cout << " and against = " << std::endl;
		for (std::vector < site_seq >::const_iterator it = against.begin();
				it < against.end();
				++it) {
			std::ostream_iterator<int> out_it (std::cout," ");
			copy ( it->begin(), it->end(), out_it );
			std::cout << std::endl;
		}
	}
	#endif 
	
	int sumSqDiff = 0;
	
	if (!against.empty() ) {
			
		for (std::vector < site_seq >::const_iterator it = against.begin();
				it < against.end();
				++it) {
				
				site_seq tmp = *it;
				
				using namespace boost::lambda;
				
				// transform tmp to tmp - base
				std::transform(tmp.begin(), tmp.end(), 
							base.begin(), tmp.begin(), 
							(_1 -_2)*(_1 -_2) );
				
				// accumulate tmp onto sumSqDiff
				sumSqDiff = std::accumulate(tmp.begin(), tmp.end(), 
							sumSqDiff );
				
		}
		
	}	
	return sumSqDiff;
}

// one vector of the total over sites of the sum sq differences per site for every pair of pops
std::vector < int > MultiLociMicrosat::sumSqAllDiffsBetweenPops(
				const MultiLociMicrosat::popDataMap& mp) const
{
	#ifdef MYDEBUGPOPSTATS
		std::cout << "\nIn sumSqAllDiffsBetweenPops" <<std::endl;
	#endif
	
	size_t npops = mp.size();
	
	std::vector < int > sumSqDiffs;
	
	if (npops > 1) {
		sumSqDiffs.reserve(npops*(npops+1)/2);
		
		size_t pop1 = 0;
		while (pop1 + 1 < npops) {
			for (size_t pop2 = pop1 + 1; pop2 < npops; ++pop2) {
				sumSqDiffs.push_back( sumSqDiffsBetweenPopsTotalOverSites(mp, pop1, pop2) );
			}
			pop1++;
		}
	}
	
	
	#ifdef MYDEBUGPOPSTATS
		
		std::cout << "total sum sq diffs over sites are " <<std::endl;
		{
			std::ostream_iterator< int > out_it (std::cout,"\t");
			copy ( sumSqDiffs.begin(), sumSqDiffs.end(), out_it );
			std::cout << std::endl;

		} 
	#endif
	
	return sumSqDiffs;
	
}

// one figure for sum of squared differences in repeat size
// summed over all sites.
// All 0 if no individuals in pop1 or pop2
int MultiLociMicrosat::sumSqDiffsBetweenPopsTotalOverSites(
				const MultiLociMicrosat::popDataMap& mp,
				size_t pop1, size_t pop2) const
{
	#ifdef MYDEBUGPOPSTATS
		std::cout << "\nIn sumSqDiffsBetweenPopsTotalOverSite, pop1 = " 
			<< pop1 << " pop2 = " << pop2 <<std::endl;
		
	#endif
	//iterator to pair for first pop
	popDataMapConstIt cit = mp.find(pop1);
	//get the pop 2 site_seqs out once and work with the local copy
	std::vector < site_seq > seq2 = mp.find(pop2)->second;
	
	int sumSqDiffs = 0.0;
	
	size_t ns = totalNsites();
	
	if ( ns > 0) {
		
		for (std::vector < site_seq >::const_iterator it = cit->second.begin();
				it < cit->second.end();
				++it) {
			sumSqDiffs += sumSqDifferencesTotalOverSites(*it, seq2);
		}
	}	
	
	#ifdef MYDEBUGPOPSTATS
		std::cout << "sumSqDiffs = " << sumSqDiffs << std::endl;		 
	#endif
	
	return sumSqDiffs;
}



// one vector of the sum over sites of the sum sq differences per site for every pop
std::vector < int > MultiLociMicrosat::sumSqAllDiffsWithinPops(
			const MultiLociMicrosat::popDataMap& mp) const
{
	#ifdef MYDEBUGPOPSTATS
		std::cout << "\nIn sumSqAllDiffsWithinPops" <<std::endl;
	#endif
	
	size_t npops = mp.size();
	
	std::vector < int > sumSqDiffs;
	
	if (npops > 0) {
		sumSqDiffs.reserve(npops);
		
		for (size_t i = 0; i < npops; ++i) {
			sumSqDiffs.push_back( sumSqDiffsWithinPopTotalOverSites(mp, i) );
		}
	}
	
	#ifdef MYDEBUGPOPSTATS
		
		std::cout << "total sum sq diffs over sites are " <<std::endl;
		{
			std::ostream_iterator< int > out_it (std::cout,"\t");
			copy ( sumSqDiffs.begin(), sumSqDiffs.end(), out_it );
			std::cout << std::endl;

		} 
	#endif
	
	return sumSqDiffs;
}

// one figure for sum of squared differences in repeat size
// summed over all sites.
// 0 if < 2 individuals in pop or no sites
int MultiLociMicrosat::sumSqDiffsWithinPopTotalOverSites(
				const MultiLociMicrosat::popDataMap& mp,
				size_t pop) const
{
	#ifdef MYDEBUGPOPSTATS
		std::cout << "\nIn sumSqDiffsWithinPopTotalOverSites, pop = " << pop <<std::endl;
	#endif
	//get the data out once and work with the local copy
	std::vector < site_seq > popdata = mp.find(pop)->second;
	
	int sumSqDiffs = 0;
	
	size_t ns = totalNsites();
	
	if ( ns > 0) {
		
		size_t sz = popdata.size();
	
		while (sz > 1) {
			site_seq lastone= popdata.back(); // last one 
			// remove that last one
			popdata.resize( sz-1 );
			
			sz--;
			
			sumSqDiffs += sumSqDifferencesTotalOverSites(lastone, popdata);
		}
	}
	
	#ifdef MYDEBUGPOPSTATS
		std::cout << "sumSqDiffs = " << sumSqDiffs << std::endl;		 
	#endif
		
	return sumSqDiffs;
}

std::vector < double > MultiLociMicrosat::getBetweenPopAverageSumSqDiffs() const
{
	MultiLociMicrosat::popDataMap mp;
	fillPopDataMap(mp);
	
	return getBetweenPopAverageSumSqDiffs(mp);
}

std::vector < double > MultiLociMicrosat::getBetweenPopAverageSumSqDiffs(
				const MultiLociMicrosat::popDataMap& mp) const
{
	#ifdef MYDEBUGPOPSTATS
		std::cout << "\nIn getBetweenPopAverageSumSqDiffs " <<std::endl;
	#endif
	
	std::vector < double > avs;
	std::vector < int > nsams = getPopnsams();
	size_t npops = nsams.size();
		
	if (npops > 1) {
		
		std::vector < int > sumSqBetweenPops 
									= sumSqAllDiffsBetweenPops(mp);
		
		//subpopsizes
		avs.assign ( sumSqBetweenPops.begin(), 
									sumSqBetweenPops.end() );
		size_t pop1 = 0;
		std::vector < double >::iterator it = avs.begin();
		while (pop1 + 1 < npops) {
			for (size_t pop2 = pop1 + 1; pop2 < npops; ++pop2) {
				
				#ifdef MYDEBUGPOPSTATS
					std::cout << "pop1 = " << pop1 << " pop2 = " 
						<< pop2 << " total sum sq diffs = " << (*it)
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
		std::cout << "\nEnd of getBetweenPopAverageSumSqDiffs, average diffs are " <<std::endl;
		{
			std::ostream_iterator<double> out_it (std::cout,"\t");
			copy ( avs.begin(), avs.end(), out_it );
			std::cout << std::endl;
		} 
	#endif
	return avs;
}

			
std::vector < double > MultiLociMicrosat::getCorrectedBetweenPopAverageSumSqDiffs() const
{
	#ifdef MYDEBUGPOPSTATS
		std::cout << "\nIn getCorrectedBetweenPopAverageSumSqDiffs " <<std::endl;
	#endif
	
	MultiLociMicrosat::popDataMap mp;
	fillPopDataMap(mp);
	
	std::vector < double > avs = getBetweenPopAverageSumSqDiffs(mp);
	
	if (!avs.empty()) {
		
		std::vector < double > withinPopAvs = getWithinPopAverageSumSqDiffs(mp);
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
		std::cout << "\nEnd of getCorrectedBetweenPopAverageSumSqDiffs, average diffs are " <<std::endl;
		{
			std::ostream_iterator<double> out_it (std::cout,"\t");
			copy ( avs.begin(), avs.end(), out_it );
			std::cout << std::endl;
		} 
	#endif
	return avs;
	
}

std::vector < double > MultiLociMicrosat::getWithinPopAverageSumSqDiffs() const
{
	MultiLociMicrosat::popDataMap mp;
	fillPopDataMap(mp);
	
	return getWithinPopAverageSumSqDiffs(mp);
}

// divisor is nsams[pop]*(nsams[pop]-1)/2			
std::vector < double > MultiLociMicrosat::getWithinPopAverageSumSqDiffs(
					const MultiLociMicrosat::popDataMap& mp) const
{
	#ifdef MYDEBUGPOPSTATS
		std::cout << "\nIn getWithinPopAverageSumSqDiff " <<std::endl;
	#endif
	
	//subpopsizes
	std::vector < int > nsams = getPopnsams();
	
	std::vector < int > sumSqDiffsWithinPops = sumSqAllDiffsWithinPops(mp);
	
	std::vector < double > avs ( sumSqDiffsWithinPops.begin(), 
								sumSqDiffsWithinPops.end() );
	
	boost::lambda::placeholder1_type a1;
	boost::lambda::placeholder2_type n1;
	
	std::transform(avs.begin(), avs.end(), 
		nsams.begin(), avs.begin(),
		boost::lambda::bind( GeneticPopStructureAnalysable::protectedDivision, a1, (n1*(n1-1)/2.0) ) ); 
		
	
	#ifdef MYDEBUGPOPSTATS
		std::cout << "\nEnd of getWithinPopAverageSumSqDiff, average diffs are " <<std::endl;
		{
			std::ostream_iterator<double> out_it (std::cout,"\t");
			copy ( avs.begin(), avs.end(), out_it );
			std::cout << std::endl;
		} 
	#endif
	
	return avs;
}

std::vector < double > MultiLociMicrosat::getPairwiseAverageSumSqDiffs() const			
{
	
	MultiLociMicrosat::popDataMap mp;
	fillPopDataMap(mp);
	
	std::vector < double > withinPops = getWithinPopAverageSumSqDiffs(mp);
	std::vector < double > betweenPops = getBetweenPopAverageSumSqDiffs(mp);
	
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


double MultiLociMicrosat::calcPairwiseRst(
						const std::vector < int >& nsams,
						const std::vector < double >& wp_diff,
						const std::vector < int >& sumSqWithin,
						const std::vector < int >& sumSqBetween,
						size_t pop1, size_t pop2) const
{
	double Rst = 0.0; 
	size_t nsam1 = nsams[pop1];
	size_t nsam2 = nsams[pop2];
	size_t totnsam = nsam1 + nsam2;
	
	#ifdef MYDEBUGPOPSTATS
		std::cout << "\nIn calcPairwiseRst, pop1 = " << pop1 << " pop2 = " << pop2 << std::endl;	
		std::cout << "nsam1 = " << nsam1 << " nsam2 = " << nsam2 << " totnsam = " << totnsam << std::endl;		 
	#endif
		
	
	if (nsam1 > 0 && nsam2 > 0) {
		
		size_t npops = 2; // we are always comparing 2 pops
		
		//work out SSD(WP) and the df N-P for this
		
		#ifdef MYDEBUGPOPSTATS
			std::cout << "wp_diff[pop1] = ";
			std::cout << wp_diff[pop1] << std::endl;
			std::cout << "wp_diff[pop2] = "; 
			std::cout << wp_diff[pop2] << std::endl;
		#endif
		
		double SSD_wp = wp_diff[pop1] + wp_diff[pop2];
			
		//hence get MSD(WP) == SSD(WP)/N-2 (P=2) which is sigmasq(wp)
		double MSD_wp = protectedDivision(SSD_wp, (totnsam - npops));
		double sigmasq_wp = MSD_wp;
		
		//work out SSD(total) using accumulations of pop comparisons
		
		double SSD_total = 0.0;
		
		
		SSD_total = sumSqBetween[pop1*(wp_diff.size()-1) + pop2-1 - ((pop1+1)*pop1)/2];
		
		#ifdef MYDEBUGPOPSTATS
			std::cout << "Index into sumSqBetween = [pop1*(wp_diff.size()-1) + pop2-1 - ((pop1+1)*pop1)/2] = " 
				<< (pop1*(wp_diff.size()-1) + pop2-1 - ((pop1+1)*pop1)/2) << std::endl;	
			std::cout << "sumSqBetween[] = " << SSD_total << std::endl;		 
		#endif
		
		SSD_total += sumSqWithin[pop1]
						+ sumSqWithin[pop2];
		
		#ifdef MYDEBUGPOPSTATS
			std::cout << "sumSqAllWithin[pop1] + sumSqWithin[pop2] = " 
				<< (sumSqWithin[pop1] + sumSqWithin[pop2]) << std::endl;	
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
		Rst = protectedDivision(sigmasq_ap, sigmasq_total); // 
	}
	
	#ifdef MYDEBUGPOPSTATS
		std::cout << "Rst = " << Rst << std::endl;		 
	#endif
	return Rst;
}




std::vector < int > MultiLociMicrosat::getPopnsams() const
{
	return _pop->makePopnsams();
}

size_t MultiLociMicrosat::getNpops() const
{
	return _pop->nSubPops();
}

std::vector < size_t > MultiLociMicrosat::getCountAllDiffsBetweenPops() const
{
	if ( !checkCachedDifferences() ) recalcCachedDifferences();
	return cachedCountAllDiffsBetweenPops;
}

std::vector < size_t > MultiLociMicrosat::getCountAllDiffsWithinPops() const
{
	if ( !checkCachedDifferences() ) recalcCachedDifferences();
	return cachedCountAllDiffsWithinPops;
}



// non member functions
std::ostream& mct::operator<< (std::ostream &out, const MultiLociMicrosat &m)
{
	out << m.toString();
	return out;
}



#if(0)
// not sure that this is correct for microsatellites

boost::shared_ptr < const SummaryStatistic  > 
	MultiLociMicrosat::distanceNeiDaBetweenPop() const
{
	boost::shared_ptr < SummaryStatistic  > result
		(new SummaryStatistic(distanceNeiDaBetweenPopVec() ));
	return result;
}


std::vector < double > 
	MultiLociMicrosat::distanceNeiDaBetweenPopVec() const
{
	std::vector < double > res;
	
	std::vector < std::vector < double > > alleleSharingBetweenPop 
				= alleleSharingBetweenPairedPopBySite();
	
	res = averageOverVecs(res, alleleSharingBetweenPop);
	
	using namespace boost::lambda;
	// 1- average summed sqrt product of freqs
	transform(res.begin(), res.end(), res.begin(), 1 - _1);
	
	return res;
}

std::vector < std::vector < double > > 
		MultiLociMicrosat::alleleSharingBetweenPairedPopBySite() const
{
	std::vector < std::vector < double > > res;
	
	size_t ns = totalNsites();
	
	// population structure knows about pops
	size_t npop = _pop->nSubPops();
	
	assert(npop > 0); // see constructor
	
	std::pair<int, int> min_max = getMinMax();
	
	for (size_t p1 = 0; p1 < npop - 1; ++p1) {
		
		std::vector < site_seq > thisPop1 
										= getSubPopLoci(container, p1);
		size_t thisNsam1 = thisPop1.size();
		// could be 0
		
		std::vector < freqmap > fmapvec1 = makeFreqMap(thisPop1);
				
		for (size_t p2 = p1+1; p2 < npop; ++p2) {
			
			std::vector < site_seq > thisPop2 
										= getSubPopLoci(container, p2);
			size_t thisNsam2 = thisPop2.size();
			
			if (ns > 0) {
			
				res.push_back( std::vector < double >( ns, 0.0 ) );
				
				if ( (thisNsam1 > 0) && (thisNsam2 > 0)) {

					std::vector < freqmap > fmapvec2 = makeFreqMap(thisPop2);
					
					double divisor = static_cast < double > (thisNsam1 * thisNsam2);
					
					
					
					// site index s
					for (size_t s = 0; s < ns; ++s) {
						
						double sumsqrprod = 0.0;
						
						// for each possible repeat size
						// if it is in both pops, take product of freq
						
						for (int rep = min_max.first; rep <= min_max.second ; ++rep ) {
							
							freqmap_it it1 = fmapvec1[s].find(rep);
							freqmap_it it2 = fmapvec2[s].find(rep);
							
							if ( (it1 != fmapvec1[s].end()) && (it2 != fmapvec2[s].end()) ) {
								
								sumsqrprod += std::sqrt( 
									static_cast < double >((it1->second) * (it2->second)) );
								
							}
							
						} // end loop over possible rep values
						
						res.back()[s] = sumsqrprod/std::sqrt( divisor) ;
						
					} // end locus loop
				} // end if both subpops have samples
				// if at least one of the subpops has no samples the allele sharing values are all 0
			} // end if ns > 0
			else res.push_back( std::vector < double >() ); // empty vector
		} // end pop2 loop
	} // end pop1 loop
	
	// if npops == 1, res will be empty
	return res;
}

std::pair <int, int > MultiLociMicrosat::getMinMax() const
{
	int min = 0; 
	int max = 0;
	
	if ( !empty() && !(container.at(0).empty())) {
		
		int min = -1; // dummy value
		int max = 0;
		
		//size_t ns = totalNsites();
		
		// we are going to have to rearrange elements
		for (std::vector< site_seq >::const_iterator it = container.begin();
			it < container.end();
			++it ) {
				site_seq tmp = *it;
				
				// use nth element to get max in range
				std::nth_element (tmp.begin(), tmp.end(), tmp.end());
				if ( tmp.back() > max ) max = tmp.back();
				if ( min < 0) min = max; // ie only on first run
				// use nth element to get min in range
				std::nth_element (tmp.begin(), tmp.begin(), tmp.end());
				if ( tmp.front() < min ) min = tmp.front();
		}
	}
	return std::pair < int, int > (min, max);
}
	
		
std::vector < MultiLociMicrosat::freqmap > MultiLociMicrosat::makeFreqMap( 
			const std::vector < site_seq >& popsites ) const
{
	// for each site, 
	size_t ns = totalNsites();
	size_t nsam = popsites.size();
	
	std::vector < freqmap > fmap;
	if (ns > 0) {
		// stuff the vector with empty maps
		fmap = std::vector < freqmap >(ns);
		
		for ( size_t s = 0; s < ns; ++s ) {
			
			for (size_t i = 0; i < nsam; ++i) { // individual
				
				freqmap_check check = (fmap[s]).insert(
							freqmap_pair( popsites[i][s],1 ) );
				if (!check.second) { // value was already there
					(check.first)->second ++;
				}
				
			}
			// if nsam == 0, each map will have no entries in it
		}
		
		#ifdef POPDEBUG
		{
			std::cout << "--- The map is ---" << std:: endl;
			for (std::vector < freqmap >::iterator it = fmap.begin(); 
							it < fmap.end();
							++it ) {
				std::cout << "locus:" << std:: endl;
				for (freqmap_it mit = it->begin(); mit != it->end(); ++mit ) {
					
					std::cout << "( " << mit->first << " -> " << mit->second << " )";
				}
				std::cout << std:: endl;
			}
			std::cout << "-------------------" << std:: endl;
			
		}
		#endif
	}
	// if ns == 0 fmap will be empty
	return fmap;
}



#endif
