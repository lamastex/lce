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
\brief MultiLociPolyTable definitions.
*/

#include "multi_loci_polytable.hpp"

#include <Sequence/PolySNP.hpp>

#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/lambda/lambda.hpp>

#include <gsl/gsl_math.h> //gsl_isnan()

#include <sstream>  // to be able to manipulate strings as streams
#include <algorithm> // for transform
#include <functional> // for binary and unary functions
#include <numeric> // for accumulate
#include <stdexcept>

//#define MYDEBUGPOPSTATS
//#define MYDEBUGPOPSTATS1
#if defined (MYDEBUGPOPSTATS) || defined (MYDEBUGPOPSTATS1)
	#include <iterator> // for ostream_iterator
#endif

#include <iostream>

using namespace mct;

MultiLociPolyTable::MultiLociPolyTable(
			const boost::shared_ptr< const hudson_ms::PopulationStructure > _p)
			: 	_pop(_p), cachedLociInPopCompCalcs(0)
{
	if ( !(_pop) ) {
		throw std::invalid_argument(
			"MultiLociPolyTable::MultiLociPolyTable(...) :\nPopulation pointer null" );
	}
	if (_pop->empty() ) {
		throw std::runtime_error(
			"MultiLociPolyTable::MultiLociPolyTable(...) :\nPopulation empty" );
	}
}



MultiLociPolyTable::~MultiLociPolyTable()
{
	//std::cout << "Deleting MultiLociPolyTable" << std::endl;
		
}

			
const boost::shared_ptr <Sequence::PolyTable >& MultiLociPolyTable::at(
										const size_t index) const
{
	return container.at(index);
}

size_t MultiLociPolyTable::nsitesAt(const size_t index) const
{
	return _nsites.at(index);
}

void MultiLociPolyTable::add(
				const boost::shared_ptr <Sequence::PolyTable >& toAdd,
				size_t _sites)
{
	if (toAdd->size() == nsam() || toAdd->numsites() == 0) {
		container.push_back(toAdd);
		_nsites.push_back(_sites);
	}
	
	else {
		
		std::ostringstream stm;
		toAdd->print(stm);
		std::string msg = ":\nThe problem PolyTable is at " + stm.str();
	
		throw std::invalid_argument(
			"MultiLociPolyTable::at(const size_t index)" + msg );
	}
}

bool MultiLociPolyTable::empty() const
{
	return container.empty();
}

const boost::shared_ptr< const hudson_ms::PopulationStructure > 
					MultiLociPolyTable::populationStructure() const
{
	return _pop;
}
								
size_t MultiLociPolyTable::nsam() const
{
	return _pop->totalNsam();
}

std::vector < size_t > MultiLociPolyTable::nsites() const
{
	return _nsites;
}

size_t MultiLociPolyTable::nloci() const
{
	return container.size();
}


boost::shared_ptr < const SummaryStatistic  > 
						MultiLociPolyTable::getSummaryStatistic() const
{
	double ss[] = {getMeanSegSitesValue(), getMeanHeterozygosityValue()}; 
	
	boost::shared_ptr < SummaryStatistic  > result
		(new SummaryStatistic(ss, sizeof(ss) / sizeof(double) ));
	return result;
		
}


boost::shared_ptr < const SummaryStatistic  > 
						MultiLociPolyTable::getAllAvailableStatistic() const
{
	double ss[] = {getMeanSegSitesValue(), getMeanHeterozygosityValue(), 
										getMeanTajimasDValue()}; 
	boost::shared_ptr < SummaryStatistic  > result
		(new SummaryStatistic(ss, sizeof(ss) / sizeof(double) ));
	return result;
		
}

boost::shared_ptr < const SummaryStatistic  > 
						MultiLociPolyTable::getMeanTajimasD() const
{
	boost::shared_ptr < SummaryStatistic  > result
		(new SummaryStatistic( std::vector < double > (1, getMeanTajimasDValue() )) );
	return result;
}

boost::shared_ptr < const SummaryStatistic  > 
						MultiLociPolyTable::getMeanHeterozygosity() const
{
	boost::shared_ptr < SummaryStatistic  > result
		(new SummaryStatistic( std::vector < double > (1, getMeanHeterozygosityValue() )) );
	return result;
}

boost::shared_ptr < const SummaryStatistic  > 
						MultiLociPolyTable::getMeanSegSites() const
{
	boost::shared_ptr < SummaryStatistic  > result
		(new SummaryStatistic( std::vector < double > (1, getMeanSegSitesValue() )) );
	return result;
}

double MultiLociPolyTable::getMeanTajimasDValue() const
{
	std::vector < double > Ds = MultiLociPolyTable::getTajimasD();
		
		return accumulate(Ds.begin(), Ds.end(), 0.0)/nloci();
}

double MultiLociPolyTable::getMeanHeterozygosityValue() const
{
		std::vector < double > Hs = MultiLociPolyTable::getHeterozygosity();
		return accumulate(Hs.begin(), Hs.end(), 0.0)/nloci();
}
			
double MultiLociPolyTable::getMeanSegSitesValue() const
{
	std::vector < size_t > Ss = MultiLociPolyTable::getSegSites();
	return (accumulate(Ss.begin(), Ss.end(), 0.0))/nloci();
}
		
			
std::string MultiLociPolyTable::toString() const
{
	std::ostringstream stm;
	
	stm << "{ ";
	size_t i = 0;
	for (std::vector < 
				boost::shared_ptr <Sequence::PolyTable > 
				>::const_iterator it = container.begin();
		it < container.end();
		++it, ++i) {
		stm << "(locus " << i << ")\n";
		(*it)->print(stm);
		stm << ": nsites " << _nsites.at(i) << "\n";
	}
	stm << "}\n";
	
	return stm.str();
}


std::vector < double > MultiLociPolyTable::getTajimasD() const
{
	std::vector<double> result(nloci());
	
	std::transform(container.begin(), container.end(), result.begin(),
	getD );
	return result;
}



std::vector < size_t > MultiLociPolyTable::getSegSites() const
{
	std::vector<size_t> result(nloci());

	std::transform(container.begin(), container.end(), result.begin(),
	boost::bind(&Sequence::PolyTable::numsites, _1) );
	//getS );
	
	return result;
}


std::vector < double > MultiLociPolyTable::getHeterozygosity() const
{
	std::vector<double> result(nloci());

	std::transform(container.begin(), container.end(), result.begin(),
	getThetaPi );
	return result;
}

std::ostream& MultiLociPolyTable::printSubpopData(std::ostream& out) const
{
	// popstrings
	MultiLociPolyTable::popStringsMap mp;
	fillPopStringsMap(mp);
	
	#ifdef MYDEBUGPOPSTATS
		std::cout << "Have map" << std::endl;
	#endif
	
	int counter = 0;
	for (MultiLociPolyTable::popStringsMapIt mit = mp.begin();
			mit != mp.end();
			++mit) {
		
		for (std::vector < std::string >::iterator it = mit->second.begin();
				it < mit->second.end();
				++it) {
		
			out << "\t" << counter << "\t1\t" << (*it) << "\n"; // freq 1 for each
			counter++;
		}
		
		out << std::endl; // new line after each pop
	}
	return out; // return by reference
}



// protected

double MultiLociPolyTable::getD (const boost::shared_ptr < Sequence::PolyTable >& ptr)
{
	Sequence::PolySNP polysnp(&(*ptr));
	
	double d = polysnp.TajimasD();
	if (gsl_isnan(d)) d = 0.0; // change nans to 0.0
	return d;
}


double MultiLociPolyTable::getThetaPi (const boost::shared_ptr < Sequence::PolyTable >& ptr)
{
	Sequence::PolySNP polysnp(&(*ptr));
	
	return polysnp.ThetaPi();
}

// implement from base class
std::vector < int > MultiLociPolyTable::getPopnsams() const
{
	return _pop->makePopnsams();
}

size_t MultiLociPolyTable::getNpops() const
{
	return _pop->nSubPops();
}

std::vector < size_t > MultiLociPolyTable::getCountAllDiffsBetweenPops() const
{
	if ( !checkCachedDifferences() ) recalcCachedDifferences();
	return cachedCountAllDiffsBetweenPops;
}

std::vector < size_t > MultiLociPolyTable::getCountAllDiffsWithinPops() const
{
	if ( !checkCachedDifferences() ) recalcCachedDifferences();
	return cachedCountAllDiffsWithinPops;
}


/* return true if cached differences are up to date*/
bool MultiLociPolyTable::checkCachedDifferences() const
{
	return (cachedLociInPopCompCalcs == _nsites.size());
	
}

/* Recalculate the the cached differences*/
void MultiLociPolyTable::recalcCachedDifferences() const // stats updated are mutable
{
	#ifdef MYDEBUGPOPSTATS
		std::cout << "\nIn recalcCachedDifferences" <<std::endl;
	#endif
	
	// recalculate and cache
	MultiLociPolyTable::popStringsMap mp;
	fillPopStringsMap(mp);
	
	// make sure we have up to date stats
	countAllDiffsWithinPops(mp);
	countAllDiffsBetweenPops(mp);
		
	cachedLociInPopCompCalcs = _nsites.size();

}

MultiLociPolyTable::popStringsMap& MultiLociPolyTable::fillPopStringsMap(
						MultiLociPolyTable::popStringsMap& mp) const
{
	#ifdef MYDEBUGPOPSTATS1
		std::cout << "In fillPopStringsMap" << std::endl;
	#endif
	// go through all the loci, get their strings data, bolt together;
	size_t totnsam = _pop->totalNsam();
	
	if (totnsam > 0) {
	
		std::vector < std::string > strs(totnsam);
		
		for (std::vector< boost::shared_ptr <Sequence::PolyTable > >::const_iterator 
			it = container.begin();
			it < container.end();
			++it) {
			
			#ifdef MYDEBUGPOPSTATS1
				std::cout << "Starting on a locus" << std::endl;
			#endif
			
			std::vector < std::string > this_strs= (*it)->GetData();
			
			// add the new strings if any onto the old ones
			if (!this_strs.empty()) {
			
				transform( strs.begin(), strs.end(), this_strs.begin(), 
							strs.begin(),
							std::plus < std::string >() );
			}
			
			#ifdef MYDEBUGPOPSTATS1
			{
				std::cout << "Added a locus, strings are now" << std::endl;
				std::ostream_iterator<std::string> out_it (std::cout,"\n");
				copy ( strs.begin(), strs.end(), out_it );
 
			}	
			#endif
							
		}
		
		// population structure knows about pops
		std::vector < int > nsams = _pop->makePopnsams();
		size_t npops = nsams.size();
			
		if (strs.front().empty()) {
			popStringsMapIt mit = mp.begin();
			std::vector < std::string > tmp; // empty vector
			
			for (size_t i = 0; i < npops; ++i) {
				
				// max efficiency inserting
				mit = mp.insert (mit, 
							std::pair< size_t, std::vector < std::string > >
								(i, tmp ) );
			}
		}
		else {

			// divvie up by population
			
			#ifdef MYDEBUGPOPSTATS1
				std::cout << "npops = " << npops << std::endl;
					 
			#endif
			
			std::vector < std::string >::const_iterator first_it 
									= strs.begin();
			std::vector < std::string >::const_iterator last_it 
									= strs.end();
			advance( last_it, -(totnsam-nsams[0]) ); // we know there is at least one pop

			popStringsMapIt mit = mp.begin();
			
			for (size_t i = 0; i < npops; ++i) {
				
				#ifdef MYDEBUGPOPSTATS1
					std::cout << "pop index " << i << " nsams[i] = " << nsams[i] << std::endl;
				#endif
				
				// after first pop, set the iterator around that ones strings
				if (i > 0 ) { 
					advance( first_it, nsams[i-1] );
					advance( last_it, nsams[i] );
				}
				
				// add entry to map for strings for pop index i
				// max efficiency inserting
				mit = mp.insert (mit, 
							std::pair< size_t, std::vector < std::string > >
								(i, std::vector < std::string >(first_it, last_it) ) );  
				#ifdef MYDEBUGPOPSTATS1
				{
					std::cout << "map is now" << std::endl;
					for (MultiLociPolyTable::popStringsMapIt debug_it = mp.begin();
						debug_it != mp.end(); ++debug_it) {
					
						std::ostream_iterator<std::string> out_it (std::cout,"\n");
						copy ( debug_it->second.begin(), debug_it->second.end(), out_it );
						std::cout << std::endl;
					}
				}	
				#endif
				
			}
		}
	}
	
	return mp; // return by reference
}


size_t MultiLociPolyTable::countDifferences(const std::string& str, 
						const std::vector < std::string >& strs )
{
	#ifdef MYDEBUGPOPSTATS2
	{
		std::cout << "In countDifferences, str = " << str << " and strs = " << std::endl;
		std::ostream_iterator<std::string> out_it (std::cout,"\n");
		copy ( strs.begin(), strs.end(), out_it );
		std::cout << std::endl;
	}
	#endif 
	size_t diffCount = 0;
	
	if (!strs.empty() ) {
		std::vector < std::string > tmp;
		tmp.reserve(strs.size());
	
		// not sure if this is worth it or not
		for (std::vector < std::string  >::const_iterator it = strs.begin();
			it < strs.end();
			++it) {
				if ( str != (*it) ) tmp.push_back(*it);
		}
		
		
		/* could go from the end of str
		 * and erase one from the end of every string in strs every time it's checked
		 * and recheck if rest is == str in every loop? */
		
		// assume the strings are of compatible length
		for (size_t i = 0; i < str.size(); ++i) {
			char c = str[i];
			
			#ifdef MYDEBUGPOPSTATS2
				std::cout << "i = " << i << " and c = " << c << std::endl;
			#endif 
			
			// char by char comparison
			for (std::vector < std::string  >::iterator it = tmp.begin();
				it < tmp.end();
				++it) {
					
					#ifdef MYDEBUGPOPSTATS2
						std::cout << "comp string = " << (*it) << " it->at(i) = ";
						std::cout << it->at(i) << std::endl;
					#endif 
					if ( c != it->at(i) ) diffCount++;
			}
		}
	}	
	return diffCount;
}

void MultiLociPolyTable::countAllDiffsBetweenPops(
				const MultiLociPolyTable::popStringsMap& mp) const
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

size_t MultiLociPolyTable::countDiffsBetweenPops(
							const MultiLociPolyTable::popStringsMap& mp,
							size_t pop1, size_t pop2)
{
	#ifdef MYDEBUGPOPSTATS
		std::cout << "\nIn countDiffsBetweenPops, pop1 = " << pop1 << " pop2 = " << pop2 <<std::endl;
		
	#endif
	
	//get the strings out once and work with the local copies
	popStringsMapConstIt cit = mp.find(pop1);
	std::vector < std::string > strs2 = mp.find(pop2)->second;
	
	size_t countDiff = 0;
	for (std::vector < std::string >::const_iterator it = cit->second.begin();
			it < cit->second.end();
			++it) {
		countDiff += countDifferences(*it, strs2);
	}
	
	#ifdef MYDEBUGPOPSTATS
		std::cout << "countDiff = " << countDiff <<std::endl;		 
	#endif
	
	return countDiff;
}

void MultiLociPolyTable::countAllDiffsWithinPops(
			const MultiLociPolyTable::popStringsMap& mp) const
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

size_t MultiLociPolyTable::countDiffsWithinPop(
							const MultiLociPolyTable::popStringsMap& mp,
							size_t pop)
{
	#ifdef MYDEBUGPOPSTATS
		std::cout << "\nIn countDiffsWithinPop, pop = " << pop <<std::endl;
	#endif
	//get the strings out once and work with the local copy
	std::vector < std::string > strs = mp.find(pop)->second;
	
	size_t countDiff = 0;
	size_t sz = strs.size();
	
	while (sz > 1) {
		std::string str = strs.back(); // last one 
		// remove that last one
		strs.resize( sz-1 );
		
		sz--;
		
		countDiff += countDifferences(str, strs);
	}
	
	#ifdef MYDEBUGPOPSTATS
		
		std::cout << "countDiff = " << countDiff <<std::endl;		 
	#endif
		
	return countDiff;
}



// non member functions
std::ostream& mct::operator<< (std::ostream &out, const MultiLociPolyTable &m)
{
	out << m.toString();
	return out;
}

/*
 * MultiLociPolyTable::MultiLociPolyTable(const MultiLociPolyTable& other)
		: 	_nsam(other._nsam), _nsites(other._nsites),
			container(other.container) {}


MultiLociPolyTable& MultiLociPolyTable::operator=(MultiLociPolyTable tmp)
{
	std::swap( _nsam, tmp._nsam );
	std::swap( _nsites, tmp._nsites );
    std::swap( container, tmp.container );

	return *this;
}
* */

