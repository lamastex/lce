/* mscplusplus - ms for c++.

   Copyright (C) 2011, 2012 Jenny Harlow
   
   This code was adapted for C++ from ms (Richard R. Hudson).
 
   ms - a program for generating samples under neutral models.
  
   Hudson, R. R. (2002).  Generating samples under a Wright-Fisher neutral 
   model of genetic variation.  Bioinformatics 18: 337-338.
  
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



#include "segsites.hpp"

#include <cassert>
//#include <numeric>
//#include <algorithm>
#include <iterator>
#include <stdexcept>
#include <iomanip>
#include <sstream>


//#define MYDEBUG

//#define MYDEBUG1

//#define MYDEBUG2

//for reserving space for seg positions
#define MAXSITES 100

using namespace hudson_ms;


// constructor
HudsonMSSegSites::HudsonMSSegSites(int _nsam, 
									bool _showProbss)
	: 	showProbss(_showProbss), probss(0.0)
{
	
	if (_nsam < 1) {
		throw std::invalid_argument(
			"HudsonMSSegSites::HudsonMSSegSites(...) : _nsam < 1");
	}
	nsam = static_cast<size_t> (_nsam);
		
	std::vector < double > tmp1;
	tmp1.reserve(MAXSITES);
	tmp1.swap(positions);
	// positions is empty
	
	std::vector< int > tmp2;
	tmp2.reserve(MAXSITES);
	std::vector < std::vector < int > > tmp3(nsam, tmp2);
	tmp3.swap(segsites);
	// segsites contains nsam empty vectors
}

HudsonMSSegSites::~HudsonMSSegSites() {}

void HudsonMSSegSites::setProbss(double _probss)
{
	probss = _probss;
}

void HudsonMSSegSites::addPositions(const std::vector < double >& newPos,
			const std::vector < std::vector < int > >& newSS)
{
	if (!newPos.empty()) {
		if (newSS.size() != nsam) {
			throw std::invalid_argument(
				"HudsonMSSegSites::addPositions(...) : newSS.size() != nsam");
		}
		
		size_t nNewPos = newPos.size();
		for (std::vector < std::vector < int > >::const_iterator it = newSS.begin();
			it < newSS.end();
			++it) {
			
			if (it->size() != nNewPos) {
				throw std::invalid_argument(
				"HudsonMSSegSites::addPositions(...) : (element in newSS).size != newPos.size() ");
			}
		}
				
		positions.insert(positions.end(), newPos.begin(), newPos.end());
		
		for (size_t i = 0; i < nsam; ++i) {
			
			segsites[i].insert(segsites[i].end(), newSS[i].begin(), newSS[i].end());
			
		}
		
		assert( segsites.size() == nsam
				&& positions.size() == segsites.front().size()
				&& positions.size() == segsites.back().size() );
	}
}

std::ostream& HudsonMSSegSites::outputSegSites(std::ostream& os, 
									int precPos, int precProb) const
{
	os << toString(precPos, precProb);
	
	return os;
	
}

std::string HudsonMSSegSites::toString(int precPos, int precProb) const
{
	std::ostringstream oss1;
	std::ostringstream oss2;
	
	if (showProbss) { // segsitesin > 0 && theta > 0
		//oss1 << std::scientific << std::setprecision(precProb);
		if (precProb < 0) oss1.unsetf(std::ios::floatfield);
		else oss1 << std::scientific << std::setprecision(precProb);
		oss1 << "prob: " << probss << "\n";
	}
	oss1 << "segsites: " << positions.size() << std::endl;
	
	if (positions.empty()) {
		oss1 << " " << std::endl;
	}
	else { //theta > 0
		oss1 << std::fixed << std::setprecision(precPos);
		oss1 << "positions: ";
		std::ostream_iterator<double> pos_it (oss1, " ");
		copy ( positions.begin(), positions.end(), pos_it );
	}
	std::string result = oss1.str();
	result = result.substr(0, result.size()-1); // take off last " "
	
	oss2 << result << std::endl;
	
	if (!positions.empty()) {
	
		// display the ss vertically side by side
		//string_index = 0 - nsam
		//size_t pos_index = 0 - segsites.size() = positions.size();
		
		
		
		for (std::vector < std::vector < int > >::const_iterator nit = segsites.begin();
				nit < segsites.end();
				++nit) {
			
			std::ostream_iterator<int> ss_it (oss2, "");
			copy ( nit->begin(), nit->end(), ss_it );
			oss2 << std::endl;
		}
	}
	return oss2.str(); 
}
