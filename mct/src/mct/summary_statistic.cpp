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
\brief SummaryStatistic definitions.
*/

#include "summary_statistic.hpp"


#include <boost/lambda/bind.hpp>
#include <boost/lambda/lambda.hpp>

#include <cmath> //sqrt

#include <sstream>  // to be able to manipulate strings as streams
#include <iterator> // for ostream_iterator
#include <algorithm> // transform, for_each
#include <numeric> // accumulate
#include <stdexcept>
#include <iomanip>


using namespace mct;

SummaryStatistic::SummaryStatistic(){}

	
SummaryStatistic::SummaryStatistic(const std::vector<double>& v)
	: values(v) {}

SummaryStatistic::SummaryStatistic(const double a[],
									const size_t n)
{
	if (NULL == a) {
		throw std::invalid_argument(
			"SummaryStatistic::SummaryStatistic(const double [], const size_t)");
	}
	values.insert (values.begin(), a, a + n);
}

SummaryStatistic::SummaryStatistic(const SummaryStatistic& other)
		: values(other.values)
		{}


SummaryStatistic& SummaryStatistic::operator=(SummaryStatistic tmp)
{
	 std::swap( values, tmp.values );

	return *this;
}

SummaryStatistic::~SummaryStatistic()
{}


const std::vector<double>& SummaryStatistic::getValues() const
{
	
	return values;	
}

const double& SummaryStatistic::at(const size_t index) const
{
	return values.at(index);
}
			
			
size_t SummaryStatistic::size() const
{
	return values.size();
}	

void SummaryStatistic::pushBack(double d)
{
	values.push_back(d);
}
	
bool SummaryStatistic::empty() const
{
	return values.empty();
}			

std::string SummaryStatistic::toString(int prec) const
{
	std::ostringstream stm;
	
	stm << "(" << plainString(prec) << ")";
	
	return stm.str();
}


std::string SummaryStatistic::plainString(int prec) const
{
	
	std::ostringstream stm;
	
	stm << std::fixed << std::setprecision(prec);
		
	std::ostream_iterator<double> out_it (stm, "\t");
	copy ( values.begin(), values.end(), out_it );
  
	std::string s = stm.str();


	size_t pos = s.find_last_of('\t');
	if (pos < std::string::npos) s.erase(pos, 1);
	

	return s;
}


double SummaryStatistic::getEuclideanDistance(
		const SummaryStatistic& other)
		const
{	
	std::vector < double > vals1 = getValues();
	std::vector < double > vals2 = other.getValues();
			
	//  sum of ((values-ref)^2)/normals
	
	//using namespace boost::lambda;
	
	boost::lambda::placeholder1_type v1;
	boost::lambda::placeholder2_type v2;
	
	//size_t n = (stats.front())->size();	// expected size
	transform ( vals1.begin(), vals1.end(), vals2.begin(), 
						vals1.begin(), (v1 - v2)*(v1 - v2) );
	
	return std::sqrt(std::accumulate(vals1.begin(), vals1.end(), 0.0));
	
}



// protected


// non member functions
std::ostream& mct::operator<< (std::ostream &out, const SummaryStatistic &ss)
{
	out << ss.toString();
	return out;
}

