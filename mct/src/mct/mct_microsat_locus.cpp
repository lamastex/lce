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
\brief MCTMicrosatLocus definitions.
*/
#include "mct_microsat_locus.hpp"

#include <gsl/gsl_math.h> // GSL_NAN

#include <algorithm> // swap
#include <numeric> // accumulate
#include <iterator> 
#include <sstream>  // to be able to manipulate strings as streams
#include <stdexcept>


using namespace mct;


MCTMicrosatLocus::MCTMicrosatLocus() : _nsites(0) {}

MCTMicrosatLocus::MCTMicrosatLocus( const std::vector < site_seq >& seqs )
{
	for (std::vector < site_seq >::const_iterator it = seqs.begin();
		it < seqs.end();
		++it) {
			add(*it);
			// will throw an exception if any of the seqs are incompatible
			
		}
	
}

#if (0)
MCTMicrosatLocus::MCTMicrosatLocus(
				const MCTMicrosatLocus& other)
		: container(other.container) {};


MCTMicrosatLocus& MCTMicrosatLocus::operator=(MCTMicrosatLocus tmp)
{
	std::swap( container, tmp.container );

	return *this;
}
#endif

MCTMicrosatLocus::~MCTMicrosatLocus()
{
		
}
			
const site_seq& MCTMicrosatLocus::at(const size_t index) const
{
	return container.at(index);
}

void MCTMicrosatLocus::add(
			const site_seq& seq)
{
	if (!empty() && seq.size() != _nsites) { 
		throw std::invalid_argument(
				"MCTMicrosatLocus::add(const site_seq&)");
	}
	if (empty()) {
		_nsites = seq.size();
	}
	container.push_back(seq);
	
}


size_t MCTMicrosatLocus::size() const
{
	return container.size();
}

bool MCTMicrosatLocus::empty() const
{
	return container.empty();
}

size_t MCTMicrosatLocus::nsites() const
{
	return _nsites;
}

/*
std::vector < double > MCTMicrosatLocus::getAverageRepeats() const
{
	try {
		std::vector < double > result;
		result.reserve(size());
		
		if ( _nsites > 0 ) {
			for (MCTMicrosatLocus::const_iterator it = container.begin();
				it < container.end();
				++it) {
				double init = 0.0; // type double to avoid integer div
				std::accumulate( it->begin(), it->end(), init);
				
				result.push_back(init/_nsites);
			}	
		}
		else {
				result = std::vector< double > (size(), GSL_NAN);
		} 
		return result;

	}
	catch (std::exception& e) {
		throw MCTexception(
				"Error in MCTMicrosatLocus::getAverageRepeats :\n"
						+ std::string(e.what()) );
		
	}	
		
}
*/

MCTMicrosatLocus::const_iterator MCTMicrosatLocus::begin() const
{
	return container.begin();
}

MCTMicrosatLocus::const_iterator MCTMicrosatLocus::end() const
{
	return container.end();
}

std::string MCTMicrosatLocus::toString() const
{
	std::ostringstream stm;
	
	for (std::vector < site_seq >::const_iterator it
			= container.begin();
			it < container.end();
			++it) {

		stm << "(\t";
		
		std::ostream_iterator<int> out_it (stm, "\t");
		copy ( it->begin(), it->end(), out_it );
	  
		stm << ")\n";
	}
			
	return stm.str();

}


std::string MCTMicrosatLocus::plainString() const
{
	std::ostringstream stm;
	
	for (std::vector < site_seq >::const_iterator it
			= container.begin();
			it < container.end();
			++it) {

		std::ostream_iterator<int> out_it (stm, "\t");
		copy ( it->begin(), it->end(), out_it );
	  
		stm << "\n";
	}
			
	return stm.str();

}


// non member functions
std::ostream& mct::operator<< (std::ostream &out, const MCTMicrosatLocus &m)
{
		out << m.toString();
		return out;
}



