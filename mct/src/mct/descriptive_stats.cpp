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
\brief DescriptiveStats definitions.
*/

#include "descriptive_stats.hpp"

#include <gsl/gsl_statistics_double.h>

#include <cmath>
#include <limits> 
#include <algorithm>
#include <stdexcept>
#include <iomanip>
#include <sstream>

using namespace mct;

DescriptiveStats::DescriptiveStats(const std::vector < double >& vec,
									double c)
									: criticalPercentileRange(c)
{
	try {
		if ( vec.empty() )  {
			throw std::invalid_argument("vec.empty()");
		}
		if (c <= 0.0 || c > 1.0) {
			throw std::invalid_argument("c");
		}
		
		double mean, sd, med, min, max, UQ, LQ, _Lcrit, _Ucrit;
		size_t n;
		boost::tie(mean, sd) = calcMeanAndSD(vec);
		boost::tie(min, _Lcrit, LQ, med, UQ, _Ucrit, max, n) = 
						calcRangeStats(vec);
		stats = boost::make_tuple(mean, sd, 
				min, _Lcrit, LQ, med, UQ, _Ucrit, max, n);
	}
	catch (std::invalid_argument& ia) {
	
		throw std::invalid_argument(
			std::string("DescriptiveStats::DescriptiveStats(...) : ") + ia.what());
	}
}

DescriptiveStats::DescriptiveStats(
				const std::vector < std::vector < double > >& vecvec,
				int index, 
				double c)
{
	try {
		if ( vecvec.empty()) {
			throw std::invalid_argument("vecvec.empty()");
		}
		
		size_t n_els = vecvec.size();
		
		if ( index < 0 || static_cast<size_t>(index) >= n_els) {
			throw std::invalid_argument("index");
		}
		
		DescriptiveStats tmp( vecvec[index], c );
		
		tmp.swap(*this);
		
	}
	catch (std::invalid_argument& ia) {
		throw std::invalid_argument(
			std::string(
			"DescriptiveStats::DescriptiveStats(...) :\n")
			+ ia.what() );
	}		
	
}

DescriptiveStats::DescriptiveStats(
				const SummaryStatisticSet& ss_set,
				int index, 
				double c)
{
	try {
		
		DescriptiveStats tmp(*(ss_set.getAllValues()), index, c);
		
		tmp.swap(*this);
		
	}
	catch (std::invalid_argument& ia) {
		throw std::invalid_argument(
			std::string(
			"DescriptiveStats::DescriptiveStats(...) :\n")
			+ ia.what() );
	}		
}

DescriptiveStats::DescriptiveStats(
				const boost::shared_ptr < const SummaryStatisticSet > ss_ptr,
				int index, 
				double c)
{
	try {
		
		if (!ss_ptr) {
			throw std::invalid_argument("ss_ptr empty");
		}
		
		DescriptiveStats tmp(*(ss_ptr->getAllValues()), index, c);
		
		tmp.swap(*this);
		
	}
	catch (std::invalid_argument& ia) {
		throw std::invalid_argument(
			std::string(
			"DescriptiveStats::DescriptiveStats(...) :\n")
			+ ia.what() );
	}		
}


DescriptiveStats::~DescriptiveStats() {}
				
double DescriptiveStats::mean() const
{
	return stats.get<0>();
}

double DescriptiveStats::sampleSD() const
{
	return stats.get<1>();
}

double DescriptiveStats::min() const
{
	return stats.get<2>();
}

double DescriptiveStats::lowerCriticalValue() const
{
	return stats.get<3>();
}

double DescriptiveStats::LQ() const
{
	return stats.get<4>();
}
double DescriptiveStats::median() const
{
	return stats.get<5>();
}

double DescriptiveStats::UQ() const
{
	return stats.get<6>();
}

double DescriptiveStats::upperCriticalValue() const
{
	return stats.get<7>();
}

double DescriptiveStats::max() const
{
	return stats.get<8>();
}

size_t DescriptiveStats::n() const
{
	return stats.get<9>();
}

DescriptiveStats::DStats DescriptiveStats::all() const
{
	return stats;
}

double DescriptiveStats::getCriticalPercentileRange() const
{
	return criticalPercentileRange;
}

std::string DescriptiveStats::toString(int prec) const
{
	
	std::ostringstream stm;
	
	stm << std::fixed << std::setprecision(prec) << all();
		
	return stm.str();
}

//protected


DescriptiveStats::expStats DescriptiveStats::calcMeanAndSD(
							const std::vector < double >& vec) 
{
	size_t n = vec.size();
	
	double sum_sqs = 0.0;
	double sum = 0.0;
	for (std::vector < double >::const_iterator it = vec.begin();
			it < vec.end();
			++it)
	{
		sum_sqs += (*it)*(*it);
		sum += (*it);
	}
	
	double mean = sum/n;
	double var = 0.0;
	
	if (n > 1) {
		var = (sum_sqs - sum * sum/n)/(n-1);
		
	}
	return boost::make_tuple( mean, std::sqrt(var) );
	
}

DescriptiveStats::rngStats DescriptiveStats::calcRangeStats(
							const std::vector < double >& vec) 
{
	double percentileLQ = (1.0-criticalPercentileRange)/2.0;
	double percentileUQ = 1.0-percentileLQ;
	
	
	std::vector < double > sortedvec = vec;
	std::sort(sortedvec.begin(), sortedvec.end());
	double min = sortedvec.front();
	double _Lcrit = getPercentile(sortedvec, percentileLQ);
	double LQ = getPercentile(sortedvec, 0.25);
	double med = getPercentile(sortedvec, 0.5);
	double UQ = getPercentile(sortedvec, 0.75);
	double _Ucrit = getPercentile(sortedvec, percentileUQ);
	double max = sortedvec.back();
	size_t n = sortedvec.size();
	
	return boost::make_tuple(min, _Lcrit, LQ,
								med, UQ, _Ucrit, max, n);
	

}


double DescriptiveStats::getPercentile(
							const std::vector < double >& sortedvec,
							double percentile) 
{
	return ( gsl_stats_quantile_from_sorted_data(&sortedvec[0],
						1, // stride
						sortedvec.size(),
						percentile) );
}

//private

DescriptiveStats::DescriptiveStats(const DescriptiveStats& other)
	: criticalPercentileRange(other.criticalPercentileRange), 
		stats (other.stats)
{}
			 
DescriptiveStats& DescriptiveStats::operator=(DescriptiveStats rhs)
{
	rhs.swap(*this);
	return *this;
}

void DescriptiveStats::swap (DescriptiveStats& ds) //throw()
{
	std::swap(criticalPercentileRange, ds.criticalPercentileRange); 
	std::swap(stats, ds.stats);

}
