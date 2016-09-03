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
\brief DataObjectSet definitions.
*/

#include "data_object_set.hpp"

using namespace aabc;


	
boost::shared_ptr < const mct::SummaryStatistic > 
					DataObjectSet::getSummaryStatisticMeans() const
{
	boost::shared_ptr < mct::SummaryStatisticSet > ss_set
									= getSummaryStatisticSet();
	
	return ss_set->getSummaryStatisticMeans();	
}
						
boost::shared_ptr < const mct::SummaryStatistic > 
					DataObjectSet::getSummaryStatisticVars() const
{
	boost::shared_ptr < mct::SummaryStatisticSet > ss_set
									= getSummaryStatisticSet();
	
	return ss_set->getSummaryStatisticVars();	
}

/*			
double DataObjectSet::getDistance(
					const boost::shared_ptr < DataObject >& ref)
					const
{	
	//std::vector < double > Ds = getTajimasDMeans();
	std::vector < double > Hs = getHeterozygosityMeans();
	std::vector < double > Ss = getSegSitesMeans();
	
	//double refD = ref->getMeanTajimasD();
	double refH = ref->getMeanHeterozygosity();
	double refS = ref->getMeanSegSites();
	
	//double distanceD = normalisedDistance(refD, Ds);
	double distanceH = normalisedDistance(refH, Hs);
	double distanceS = normalisedDistance(refS, Ss);
	double distance = sqrt(+ distanceH*distanceH
							+ distanceS*distanceS);
	
	#ifdef DEBUG
		
		
		
	#endif
	
					
	return distance;
}
*/


// protected


	
/*		
double DataObjectSet::normalisedDistance(const double& ref, 
						const std::vector < double >& values)
{
	size_t n = values.size();
	
	double sum_sqs = 0.0;
	double sum = 0.0;
	for (std::vector < double >::const_iterator it = values.begin();
			it < values.end();
			++it)
	{
		sum_sqs += (*it)*(*it);
		sum += (*it);
	}
	
	double mean = sum/n;
	double var = 0.0;
	double result = (ref > mean ? ref - mean: mean - ref);
	
	if (n > 1) {
		var = (sum_sqs - sum * sum/n)/(n-1);
		
	}
	if (var > 0.0) {
		result /= sqrt(var);
	}
	
	#ifdef DEBUG
		if (isnan(result) || isinf(result)) {
						
			std::cout << "Ref is " << ref << " and values are " << std::endl;
			std::ostream_iterator<double> out_it (std::cout, "\t");
			copy ( values.begin(), values.end(), out_it );
			std::cout << std::endl;
		}
	#endif
	
	return result;	
}

*/

// non member functions
std::ostream& aabc::operator<< (std::ostream &out, const DataObjectSet &d)
{
	out << d.toString();
	return out;
}


