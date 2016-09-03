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
\brief MultiLociMicrosatSet definitions.
*/

#include "multi_loci_microsat_set.hpp"
#include "summary_statistic_set.hpp"

#include <boost/function.hpp>
#include <boost/bind.hpp>

#include <sstream>  // to be able to manipulate strings as streams


using namespace mct;

MultiLociMicrosatSet::MultiLociMicrosatSet() {}

MultiLociMicrosatSet::~MultiLociMicrosatSet()
{
		
}
			
const boost::shared_ptr < const MultiLociMicrosat >& 
		MultiLociMicrosatSet::at(const size_t index) const
{
	return container.at(index);
}

void MultiLociMicrosatSet::add(
			const boost::shared_ptr < const MultiLociMicrosat >& mlptr)
{
	container.push_back(mlptr);
}

void MultiLociMicrosatSet::add(
			const boost::shared_ptr < const MultiLociMicrosatSet >& sptr)
{
	std::vector< boost::shared_ptr < const MultiLociMicrosat > > otherContainer
		= sptr->_getContainer();
	container.insert( container.end(),
					otherContainer.begin(), otherContainer.end() );
	
}


size_t MultiLociMicrosatSet::size() const
{
	return container.size();
}

bool MultiLociMicrosatSet::empty() const
{
	return container.empty();
}

boost::shared_ptr < SummaryStatisticSet > 
				MultiLociMicrosatSet::getSummaryStatisticSet() const
{
	boost::shared_ptr < SummaryStatisticSet >
			ss_set( new SummaryStatisticSet( size() ) );
	
	boost::function<boost::shared_ptr < const SummaryStatistic > 
			(boost::shared_ptr < const MultiLociMicrosat >)>
				f1 = boost::bind(&MultiLociMicrosat::getSummaryStatistic, _1);
	
	boost::function< void (boost::shared_ptr < const SummaryStatistic >) >
			// need cast to specify exactly which add to use - args are not automatically used to resolve function address
			f2 = boost::bind(static_cast<void (SummaryStatisticSet::*)( const boost::shared_ptr < const SummaryStatistic >& )>(&SummaryStatisticSet::add),
						ss_set, _1);
	
	std::for_each( container.begin(), container.end(),
		boost::bind(f2, boost::bind(f1, _1) ) );
	
	return ss_set;	
}


boost::shared_ptr < const SummaryStatistic >
					MultiLociMicrosatSet::getSummaryStatistic() const
{
	return getSummaryStatisticSet()->getSummaryStatisticMeans(); 
	
}

std::string MultiLociMicrosatSet::toString() const
{
	std::ostringstream stm;
	
	size_t i = 0;
	for (std::vector < 
				boost::shared_ptr < const MultiLociMicrosat > 
				>::const_iterator it = container.begin();
		it < container.end();
		++it, ++i) {
		stm << "(sequences " << i << ")\n" << std::endl;
		stm << (*it)->toString() << "\n";
	}
	stm << "\n";
	
	return stm.str();

}

boost::shared_ptr < SummaryStatisticSet >
				MultiLociMicrosatSet::getAverageRepeats() const
{
	boost::shared_ptr < SummaryStatisticSet >
			ss_set( new SummaryStatisticSet( size() ) );
	
	boost::function<boost::shared_ptr < const SummaryStatistic > 
			(boost::shared_ptr < const MultiLociMicrosat >)>
				f1 = boost::bind(&MultiLociMicrosat::getAverageRepeats, _1);
	
	boost::function< void (boost::shared_ptr < const SummaryStatistic >) >
			// need cast to specify exactly which add to use - args are not automatically used to resolve function address
			f2 = boost::bind(static_cast<void (SummaryStatisticSet::*)( const boost::shared_ptr < const SummaryStatistic >& )>(&SummaryStatisticSet::add),
						ss_set, _1);
	
	std::for_each( container.begin(), container.end(),
		boost::bind(f2, boost::bind(f1, _1) ) );
	
	return ss_set;	
}


boost::shared_ptr < SummaryStatisticSet > 
			MultiLociMicrosatSet::getCovariances() const
{
	boost::shared_ptr < SummaryStatisticSet >
			ss_set( new SummaryStatisticSet( size() ) );
	
	boost::function<boost::shared_ptr < const SummaryStatistic > 
			(boost::shared_ptr < const MultiLociMicrosat >)>
				f1 = boost::bind(&MultiLociMicrosat::getCovariances, _1);
	
	boost::function< void (boost::shared_ptr < const SummaryStatistic >) >
			// need cast to specify exactly which add to use - args are not automatically used to resolve function address
			f2 = boost::bind(static_cast<void (SummaryStatisticSet::*)( const boost::shared_ptr < const SummaryStatistic >& )>(&SummaryStatisticSet::add),
						ss_set, _1);
	
	std::for_each( container.begin(), container.end(),
		boost::bind(f2, boost::bind(f1, _1) ) );
	
	return ss_set;	
}


boost::shared_ptr < SummaryStatisticSet > 
			MultiLociMicrosatSet::getDistanceDeltaMuSqBetweenPop() const
{
	boost::shared_ptr < SummaryStatisticSet >
			ss_set( new SummaryStatisticSet( size() ) );
	
	boost::function<boost::shared_ptr < const SummaryStatistic > 
			(boost::shared_ptr < const MultiLociMicrosat >)>
				f1 = boost::bind(&MultiLociMicrosat::distanceGoldsteinDeltaMuSqBetweenPop, _1);
	
	boost::function< void (boost::shared_ptr < const SummaryStatistic >) >
			// need cast to specify exactly which add to use - args are not automatically used to resolve function address
			f2 = boost::bind(static_cast<void (SummaryStatisticSet::*)( const boost::shared_ptr < const SummaryStatistic >& )>(&SummaryStatisticSet::add),
						ss_set, _1);
	
	std::for_each( container.begin(), container.end(),
		boost::bind(f2, boost::bind(f1, _1) ) );
	
	return ss_set;	
}




boost::shared_ptr < SummaryStatisticSet > 
		MultiLociMicrosatSet::getDistanceSlatkinRstBetweenPop() const
{
	boost::shared_ptr < SummaryStatisticSet >
			ss_set( new SummaryStatisticSet( size() ) );
	
	boost::function<boost::shared_ptr < const SummaryStatistic > 
			(boost::shared_ptr < const MultiLociMicrosat >)>
				f1 = boost::bind(&MultiLociMicrosat::distanceSlatkinRstBetweenPop, _1);
	
	boost::function< void (boost::shared_ptr < const SummaryStatistic >) >
			// need cast to specify exactly which add to use - args are not automatically used to resolve function address
			f2 = boost::bind(static_cast<void (SummaryStatisticSet::*)( const boost::shared_ptr < const SummaryStatistic >& )>(&SummaryStatisticSet::add),
						ss_set, _1);
	
	std::for_each( container.begin(), container.end(),
		boost::bind(f2, boost::bind(f1, _1) ) );
	
	return ss_set;	
}


boost::shared_ptr < SummaryStatisticSet > 
			MultiLociMicrosatSet::getPairwiseAverageSumSqDifferences() const
{
	boost::shared_ptr < SummaryStatisticSet >
			ss_set( new SummaryStatisticSet( size() ) );
	
	boost::function<boost::shared_ptr < const SummaryStatistic > 
			(boost::shared_ptr < const MultiLociMicrosat >)>
				f1 = boost::bind(&MultiLociMicrosat::pairwiseAverageSumSqDiffs, _1);
	
	boost::function< void (boost::shared_ptr < const SummaryStatistic >) >
			// need cast to specify exactly which add to use - args are not automatically used to resolve function address
			f2 = boost::bind(static_cast<void (SummaryStatisticSet::*)( const boost::shared_ptr < const SummaryStatistic >& )>(&SummaryStatisticSet::add),
						ss_set, _1);
	
	std::for_each( container.begin(), container.end(),
		boost::bind(f2, boost::bind(f1, _1) ) );
	
	return ss_set;	
}

boost::shared_ptr < SummaryStatisticSet > 
			MultiLociMicrosatSet::getRst() const
{
	boost::shared_ptr < SummaryStatisticSet >
			ss_set( new SummaryStatisticSet( size() ) );
	
	boost::function<boost::shared_ptr < const SummaryStatistic > 
			(boost::shared_ptr < const MultiLociMicrosat >)>
				f1 = boost::bind(&MultiLociMicrosat::rst, _1);
	
	boost::function< void (boost::shared_ptr < const SummaryStatistic >) >
			// need cast to specify exactly which add to use - args are not automatically used to resolve function address
			f2 = boost::bind(static_cast<void (SummaryStatisticSet::*)( const boost::shared_ptr < const SummaryStatistic >& )>(&SummaryStatisticSet::add),
						ss_set, _1);
	
	std::for_each( container.begin(), container.end(),
		boost::bind(f2, boost::bind(f1, _1) ) );
	
	return ss_set;	
}

boost::shared_ptr < SummaryStatisticSet > 
			MultiLociMicrosatSet::getFst() const
{
	boost::shared_ptr < SummaryStatisticSet >
			ss_set( new SummaryStatisticSet( size() ) );
	
	boost::function<boost::shared_ptr < const SummaryStatistic > 
			(boost::shared_ptr < const MultiLociMicrosat >)>
				f1 = boost::bind(&MultiLociMicrosat::fst, _1);
	
	boost::function< void (boost::shared_ptr < const SummaryStatistic >) >
			// need cast to specify exactly which add to use - args are not automatically used to resolve function address
			f2 = boost::bind(static_cast<void (SummaryStatisticSet::*)( const boost::shared_ptr < const SummaryStatistic >& )>(&SummaryStatisticSet::add),
						ss_set, _1);
	
	std::for_each( container.begin(), container.end(),
		boost::bind(f2, boost::bind(f1, _1) ) );
	
	return ss_set;	
}

boost::shared_ptr < SummaryStatisticSet > 
			MultiLociMicrosatSet::getBetweenPopFst() const
{
	boost::shared_ptr < SummaryStatisticSet >
			ss_set( new SummaryStatisticSet( size() ) );
	
	boost::function<boost::shared_ptr < const SummaryStatistic > 
			(boost::shared_ptr < const MultiLociMicrosat >)>
				f1 = boost::bind(&MultiLociMicrosat::betweenPopFst, _1);
	
	boost::function< void (boost::shared_ptr < const SummaryStatistic >) >
			// need cast to specify exactly which add to use - args are not automatically used to resolve function address
			f2 = boost::bind(static_cast<void (SummaryStatisticSet::*)( const boost::shared_ptr < const SummaryStatistic >& )>(&SummaryStatisticSet::add),
						ss_set, _1);
	
	std::for_each( container.begin(), container.end(),
		boost::bind(f2, boost::bind(f1, _1) ) );
	
	return ss_set;	
}

// protected

const std::vector< boost::shared_ptr < const MultiLociMicrosat > >&
				MultiLociMicrosatSet::_getContainer() const
{
	return container;
}



// non member functions
std::ostream& mct::operator<< (std::ostream &out, const MultiLociMicrosatSet &m)
{
	out << m.toString();
	return out;
}

