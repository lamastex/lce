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
\brief MultiLociPolyTableSet definitions.
*/

#include "multi_loci_polytable_set.hpp"
#include "summary_statistic_set.hpp"



#include <boost/function.hpp>
#include <boost/bind.hpp>

#include <sstream>  // to be able to manipulate strings as streams


using namespace mct;

MultiLociPolyTableSet::MultiLociPolyTableSet() {}

MultiLociPolyTableSet::~MultiLociPolyTableSet()
{
		
}
			
const boost::shared_ptr < const MultiLociPolyTable >& MultiLociPolyTableSet::at(
										const size_t index) const
{
	return container.at(index);
	
}

void MultiLociPolyTableSet::add(
			const boost::shared_ptr < const MultiLociPolyTable >& mlptr)
{
	
	container.push_back(mlptr);
}

void MultiLociPolyTableSet::add(
			const boost::shared_ptr < const MultiLociPolyTableSet >& sptr)
{
	std::vector< boost::shared_ptr < const MultiLociPolyTable > > otherContainer
		= sptr->_getContainer();
	container.insert( container.end(),
					otherContainer.begin(), otherContainer.end() );
	
}


size_t MultiLociPolyTableSet::size() const
{
	return container.size();
}

bool MultiLociPolyTableSet::empty() const
{
	return container.empty();
}

boost::shared_ptr < SummaryStatisticSet > 
				MultiLociPolyTableSet::getSummaryStatisticSet() const
{
	boost::shared_ptr < SummaryStatisticSet >
			ss_set( new SummaryStatisticSet( size() ) );
	
	boost::function<boost::shared_ptr < const SummaryStatistic > 
			(boost::shared_ptr < const MultiLociPolyTable >)>
				f1 = boost::bind(&MultiLociPolyTable::getSummaryStatistic, _1);
	
	boost::function< void (boost::shared_ptr < const SummaryStatistic >) >
			// need cast to specify exactly which add to use - args are not automatically used to resolve function address
			f2 = boost::bind(static_cast<void (SummaryStatisticSet::*)( const boost::shared_ptr < const SummaryStatistic >& )>(&SummaryStatisticSet::add),
						ss_set, _1);
	
	std::for_each( container.begin(), container.end(),
		boost::bind(f2, boost::bind(f1, _1) ) );
	
	return ss_set;	
}


boost::shared_ptr < SummaryStatisticSet > 
				MultiLociPolyTableSet::getAllAvailableStatisticsSet() const
{
	boost::shared_ptr < SummaryStatisticSet >
			ss_set( new SummaryStatisticSet( size() ) );
	
	boost::function<boost::shared_ptr < const SummaryStatistic > 
			(boost::shared_ptr < const MultiLociPolyTable >)>
				f1 = boost::bind(&MultiLociPolyTable::getAllAvailableStatistic, _1);
	
	boost::function< void (boost::shared_ptr < const SummaryStatistic >) >
			// need cast to specify exactly which add to use - args are not automatically used to resolve function address
			f2 = boost::bind(static_cast<void (SummaryStatisticSet::*)( const boost::shared_ptr < const SummaryStatistic >& )>(&SummaryStatisticSet::add),
						ss_set, _1);
	
	std::for_each( container.begin(), container.end(),
		boost::bind(f2, boost::bind(f1, _1) ) );
	
	return ss_set;	
}


boost::shared_ptr < const SummaryStatistic >
					MultiLociPolyTableSet::getSummaryStatistic() const
{
	return getSummaryStatisticSet()->getSummaryStatisticMeans(); 
	
}

boost::shared_ptr < const SummaryStatistic >
		MultiLociPolyTableSet::getTajimasDSummaryStatistic() const
{
	return getTajimasDMeans()->getSummaryStatisticMeans(); 
	
}

boost::shared_ptr < const SummaryStatistic >
		MultiLociPolyTableSet::getHeterozygositySummaryStatistic() const
{
	return getHeterozygosityMeans()->getSummaryStatisticMeans(); 
	
}

boost::shared_ptr < const SummaryStatistic >
		MultiLociPolyTableSet::getSegSitesSummaryStatistic() const
{
	return getSegSitesMeans()->getSummaryStatisticMeans(); 
	
}


std::string MultiLociPolyTableSet::toString() const
{
	std::ostringstream stm;
	
	stm << "[ ";
	size_t i = 0;
	for (std::vector < 
				boost::shared_ptr < const MultiLociPolyTable > 
				>::const_iterator it = container.begin();
		it < container.end();
		++it, ++i) {
		stm << "(sequences " << i << ")\n" << std::endl;
		stm << (*it)->toString() << "\n";
	}
	stm << "]\n";
	
	return stm.str();

}


boost::shared_ptr < SummaryStatisticSet > 
					MultiLociPolyTableSet::getTajimasDMeans() const
{
	boost::shared_ptr < SummaryStatisticSet >
			ss_set( new SummaryStatisticSet( size() ) );
	
	boost::function<boost::shared_ptr < const SummaryStatistic > 
			(boost::shared_ptr < const MultiLociPolyTable >)>
				f1 = boost::bind(&MultiLociPolyTable::getMeanTajimasD, _1);
	
	boost::function< void (boost::shared_ptr < const SummaryStatistic >) >
			// need cast to specify exactly which add to use - args are not automatically used to resolve function address
			f2 = boost::bind(static_cast<void (SummaryStatisticSet::*)( const boost::shared_ptr < const SummaryStatistic >& )>(&SummaryStatisticSet::add),
						ss_set, _1);
	
	std::for_each( container.begin(), container.end(),
		boost::bind(f2, boost::bind(f1, _1) ) );
	
	return ss_set;	
}


boost::shared_ptr < SummaryStatisticSet > 
					MultiLociPolyTableSet::getSegSitesMeans() const
{
	boost::shared_ptr < SummaryStatisticSet >
			ss_set( new SummaryStatisticSet( size() ) );
	
	boost::function<boost::shared_ptr < const SummaryStatistic > 
			(boost::shared_ptr < const MultiLociPolyTable >)>
				f1 = boost::bind(&MultiLociPolyTable::getMeanSegSites, _1);
	
	boost::function< void (boost::shared_ptr < const SummaryStatistic >) >
			// need cast to specify exactly which add to use - args are not automatically used to resolve function address
			f2 = boost::bind(static_cast<void (SummaryStatisticSet::*)( const boost::shared_ptr < const SummaryStatistic >& )>(&SummaryStatisticSet::add),
						ss_set, _1);
	
	std::for_each( container.begin(), container.end(),
		boost::bind(f2, boost::bind(f1, _1) ) );
	
	return ss_set;	
}

boost::shared_ptr < SummaryStatisticSet > 
					MultiLociPolyTableSet::getHeterozygosityMeans() const
{
	boost::shared_ptr < SummaryStatisticSet >
			ss_set( new SummaryStatisticSet( size() ) );
	
	boost::function<boost::shared_ptr < const SummaryStatistic > 
			(boost::shared_ptr < const MultiLociPolyTable >)>
				f1 = boost::bind(&MultiLociPolyTable::getMeanHeterozygosity, _1);
	
	boost::function< void (boost::shared_ptr < const SummaryStatistic >) >
			// need cast to specify exactly which add to use - args are not automatically used to resolve function address
			f2 = boost::bind(static_cast<void (SummaryStatisticSet::*)( const boost::shared_ptr < const SummaryStatistic >& )>(&SummaryStatisticSet::add),
						ss_set, _1);
	
	std::for_each( container.begin(), container.end(),
		boost::bind(f2, boost::bind(f1, _1) ) );
	
	return ss_set;	
}


// protected

std::vector < std::vector < double > > MultiLociPolyTableSet::getTajimasD() const
{
	std::vector< std::vector < double > > result;
	
	for (std::vector< 
		boost::shared_ptr < const MultiLociPolyTable > 
		>::const_iterator it = container.begin();
		it < container.end();
		++it)
	{
		result.push_back( (*it)->getTajimasD() );
	}
	return result;
}


std::vector < std::vector < size_t > > MultiLociPolyTableSet::getSegSites() const
{
	std::vector<std::vector < size_t > > result;

	for (std::vector< 
		boost::shared_ptr < const MultiLociPolyTable > 
		>::const_iterator it = container.begin();
		it < container.end();
		++it)
	{
		result.push_back( (*it)->getSegSites() );
	}
	return result;
}


std::vector < std::vector < double > > 
						MultiLociPolyTableSet::getHeterozygosity() const
{
	std::vector<std::vector < double > > result;
	
	for (std::vector< 
		boost::shared_ptr < const MultiLociPolyTable > 
		>::const_iterator it = container.begin();
		it < container.end();
		++it)
	{
		result.push_back( (*it)->getHeterozygosity() );
	}
	return result;
}


const std::vector< boost::shared_ptr < const MultiLociPolyTable > >&
				MultiLociPolyTableSet::_getContainer() const
{
	return container;
}



// non member functions
std::ostream& mct::operator<< (std::ostream &out, const MultiLociPolyTableSet &m)
{
	out << m.toString();
	return out;
}

