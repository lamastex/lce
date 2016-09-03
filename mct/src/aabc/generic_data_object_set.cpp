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
\brief GenericDataObjectSet definitions.
*/

#include "generic_data_object_set.hpp"
#include "summary_statistic_set.hpp"

#include <boost/function.hpp>
#include <boost/bind.hpp>

#include <sstream>  // to be able to manipulate strings as streams


using namespace aabc;

GenericDataObjectSet::GenericDataObjectSet() {}

GenericDataObjectSet::~GenericDataObjectSet()
{}
			
const boost::shared_ptr < const DataObject >& GenericDataObjectSet::at(
										const size_t index) const
{
	return container.at(index);
}

void GenericDataObjectSet::add(
			const boost::shared_ptr < const DataObject >& doptr)
{
	
	container.push_back(doptr);
}

void GenericDataObjectSet::add(
			const boost::shared_ptr < const GenericDataObjectSet >& sptr)
{
	std::vector< boost::shared_ptr < const DataObject > > otherContainer
		= sptr->_getContainer();
	container.insert( container.end(),
					otherContainer.begin(), otherContainer.end() );
}

size_t GenericDataObjectSet::size() const
{
	return container.size();
}

bool GenericDataObjectSet::empty() const
{
	return container.empty();
}

boost::shared_ptr < mct::SummaryStatisticSet > 
				GenericDataObjectSet::getSummaryStatisticSet() const
{
	boost::shared_ptr < mct::SummaryStatisticSet >
			ss_set( new mct::SummaryStatisticSet( size() ) );
	
	boost::function<boost::shared_ptr < const mct::SummaryStatistic > 
			(boost::shared_ptr < const DataObject >)>
				f1 = boost::bind(&DataObject::getSummaryStatistic, _1);
	
	boost::function< void (boost::shared_ptr < const mct::SummaryStatistic >) >
			// need cast to specify exactly which add to use - args are not automatically used to resolve function address
			f2 = boost::bind(static_cast<void (mct::SummaryStatisticSet::*)( const boost::shared_ptr < const mct::SummaryStatistic >& )>(&mct::SummaryStatisticSet::add),
						ss_set, _1);
	
	std::for_each( container.begin(), container.end(),
		boost::bind(f2, boost::bind(f1, _1) ) );
	
	return ss_set;	
}

boost::shared_ptr < const mct::SummaryStatistic >
					GenericDataObjectSet::getSummaryStatistic() const
{
	return getSummaryStatisticSet()->getSummaryStatisticMeans(); 
}

std::string GenericDataObjectSet::toString() const
{
	std::ostringstream stm;
	
	stm << "[ ";
	size_t i = 0;
	for (std::vector < 
				boost::shared_ptr < const DataObject > 
				>::const_iterator it = container.begin();
		it < container.end();
		++it, ++i) {
		stm << "(sequences " << i << ")\n" << std::endl;
		stm << (*it)->toString() << "\n";
	}
	stm << "]\n";
	
	return stm.str();

}


// protected

const std::vector< boost::shared_ptr < const DataObject > >&
				GenericDataObjectSet::_getContainer() const
{
	return container;
}



// non member functions
std::ostream& aabc::operator<< (std::ostream &out, const GenericDataObjectSet &m)
{
	out << m.toString();
	return out;
}

