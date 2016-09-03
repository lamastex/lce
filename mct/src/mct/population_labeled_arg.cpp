/* MCT - Markov Chains on Trees.

   Copyright (C) 2009 Brendan Bycroft <brb44@student.canterbury.ac.nz>

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
\brief Definitions for a class holding ARG data population labels.
*/

#include "population_labeled_arg.hpp"
#include "libseq.hpp" // for toString(hist)


#include <sstream>
#include <iterator>
#include <algorithm> // for copy


using namespace mct;


PopulationLabeledARG::PopulationLabeledARG(
			boost::shared_ptr < Sequence::arg > _h,
			const std::vector < std::string >& _pl)
		: hist(_h), populationLabels(_pl) {}

std::vector< std::string > 
				PopulationLabeledARG::getPopulationLabels() const
{
	
	return populationLabels;
}

std::vector< size_t > 
				PopulationLabeledARG::getPopulationNumbers() const
{
	
	std::vector < size_t > result;
	
	if (!populationLabels.empty()) {
		
		std::vector < std::string >::const_iterator it = populationLabels.begin();
		std::string currentLabel = (*it);
		size_t current = 0;
		do {
			if ( (*it) != currentLabel ) {
				currentLabel = (*it);
				current++;}
			result.push_back(current);
			++it;
			
		} while (it < populationLabels.end() );
	}
	assert (result.size() == populationLabels.size() );
	
	return result;
}

boost::shared_ptr < const Sequence::arg > 
							PopulationLabeledARG::getARG() const
{
	return hist;
}

size_t PopulationLabeledARG::getMarginalTreeCount() const
{
	return hist->size();
}

size_t PopulationLabeledARG::getNsam() const
{
	return populationLabels.size();
}

std::string PopulationLabeledARG::toString() const
{
	std::ostringstream stm;
	stm << LibSeq::argNewickString(*hist);
	stm << "Population labels:\t";
	copy( populationLabels.begin(), populationLabels.end(), 
						std::ostream_iterator< std::string >(stm, "\t") );
	stm << std::endl;
	
	return stm.str();
}

PopulationLabeledARG::arg_iterator PopulationLabeledARG::trees_begin() const
{
	return hist->begin();
}


PopulationLabeledARG::arg_iterator PopulationLabeledARG::trees_end() const
{
	return hist->end();
}

// Get outstream representation of a PopulationLabeledARG 
std::ostream& mct::operator<<(std::ostream & s, const PopulationLabeledARG &pa)
{
    s << pa.toString();
    return s;
}

