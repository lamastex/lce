/* MCT - Markov Chains on Trees.

   Copyright (C) 2012 Jenny Harlow

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
\brief MultiLociModelSetSet definitions.
*/

#include "multi_loci_modelset_set.hpp"


#include <sstream>  // to be able to manipulate strings as streams


using namespace mct;

MultiLociModelSetSet::MultiLociModelSetSet() {}

MultiLociModelSetSet::MultiLociModelSetSet(
					const boost::shared_ptr < ModelSet >& ms,
					size_t nloci) 
	: container (std::vector< boost::shared_ptr < ModelSet > >(nloci, ms))
{
	
}

MultiLociModelSetSet::~MultiLociModelSetSet()
{
		
}
			
const boost::shared_ptr < ModelSet >& MultiLociModelSetSet::at(
										const size_t index) const
{
	return container.at(index);
	
}

void MultiLociModelSetSet::add(
			const boost::shared_ptr < ModelSet >& ms)
{
	
	container.push_back(ms);
}




size_t MultiLociModelSetSet::size() const
{
	return container.size();
}

bool MultiLociModelSetSet::empty() const
{
	return container.empty();
}



std::string MultiLociModelSetSet::toString() const
{
	std::ostringstream stm;
	
	for (std::vector < 
				boost::shared_ptr < ModelSet > 
				>::const_iterator it = container.begin();
		it < container.end();
		++it) {
		stm << (*it)->toString() << "\n";
	}
	
	return stm.str();

}

void MultiLociModelSetSet::resetThetas(
			const std::vector < std::vector < double > >& new_thetas)
{
	if( new_thetas.size()!= size() ) {
		throw std::invalid_argument(
		"MultiLociModelSetSet::resetThetas(const std::vector < std::vector < double > >&)");
	}
	std::vector < std::vector < double > >::const_iterator cit = new_thetas.begin();
	for (std::vector< boost::shared_ptr < ModelSet > >::iterator it = container.begin();
			it < container.end();
			++it, ++cit) {
				(*it)->resetThetas(*cit);
	}
}

void MultiLociModelSetSet::resetThetas(double new_theta)
{
	for (std::vector< boost::shared_ptr < ModelSet > >::iterator it = container.begin();
			it < container.end();
			++it) {
				(*it)->resetThetas(new_theta);
	}
}


// non member functions
std::ostream& mct::operator<< (std::ostream &out, const MultiLociModelSetSet &m)
{
	out << m.toString();
	return out;
}

