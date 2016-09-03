/* MCT - Markov Chains on Trees.

   Copyright (C) 2011 Jenny Harlow, Kevin Thornton
   Copyright (C) 2008-2010 Kevin Thornton	

   This program is based on the transformations code in ABCreg,
   Thornton, K.R. (2009) Automating approximate Bayesian computation
   by local linear regression.  BMC Genetics 10: 35.	

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
\brief TanTransformer definitions.
*/

#include "reg_transformer_tan.hpp"


#include <boost/function.hpp>
#include <boost/bind.hpp>

#include <algorithm>

#include <cmath>
#include <gsl/gsl_math.h>


using namespace reg;

TanTransformer::TanTransformer(const std::vector < double >& _mins,
				const std::vector < double >& _maxs)
				 : Transformer("tangent"),
				 mins(_mins),
				 maxs(_maxs)
				 {}

TanTransformer::TanTransformer(const TanTransformer& other)
: Transformer(other.type), mins(other.mins), maxs(other.maxs) 
	{}
	
TanTransformer::~TanTransformer() {}

std::vector< std::vector <double> >& TanTransformer::data_transform(
				std::vector< std::vector <double> >& prior)
{
	size_t index = 0;	// to index into mins and maxes
	for (std::vector< std::vector <double> > ::iterator it = prior.begin();
			it < prior.end();
			++it, index++) { // move the index as well
		double this_min = mins.at(index);
		double this_max = maxs.at(index);
		
		boost::function < double (double) > tangent 
			= boost::bind(&tan_transform, _1, this_min, this_max);
			
		transform( it->begin(), it->end(), it->begin(), tangent );
	}
	
	return prior;
}

								
std::vector< std::vector <double> >&  TanTransformer::data_untransform(
				std::vector< std::vector <double> >& posterior)
{
	size_t index = 0;	// to index into mins and maxes
	
	for (std::vector< std::vector <double> > ::iterator it = posterior.begin();
			it < posterior.end();
			++it, index++) {
		double this_min = mins.at(index);
		double this_max = maxs.at(index);
				
		boost::function < double (double) > untangent 
			= boost::bind(&tan_untransform, _1, this_min, this_max);
		
		transform( it->begin(), it->end(), it->begin(), untangent );
	}
	
	return posterior;
}

									
double TanTransformer::tan_transform( const double x,
										const double minval,
										const double maxval)
//Apply the tangent transformation of Hamilton et al. 2005 PNAS 7476
{
  return -std::log(1./(std::tan( ((x-minval)/(maxval-minval))*(M_PI/2.) )));
}

double TanTransformer::tan_untransform( const double y,
										const double minval,
										const double maxval)
//undo the tangent transformation
{
  return minval + (2./M_PI)*(maxval-minval)*std::atan(std::exp(y));
}
