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
\brief ParameterPriorNormal definitions.
*/

#include "parameter_prior_normal.hpp"


#include <boost/bind.hpp>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

#include <functional> // for plus and minus


#include <sstream>  // to be able to manipulate strings as streams

using namespace aabc;

ParameterPriorNormal::ParameterPriorNormal(
					double _m, double _s,
					const boost::shared_ptr < mct::PRNGen >& _r_ptr)
					: m(_m), s(_s),
					  cdfFunc (boost::function<double (double)>
							( boost::bind(gsl_cdf_gaussian_P, 
							boost::bind(std::minus<double>(), _1, _m), _s) ) 
							),
					  pdfFunc (boost::function<double (double)>
							( boost::bind(gsl_ran_gaussian_pdf, 
							boost::bind(std::minus<double>(), _1, _m), _s) )
							),
					  r_rand_ptr(_r_ptr )
					  {}
	
/*	
ParameterPriorNormal::ParameterPriorNormal(
									const ParameterPriorNormal& other)
										: m(other.m), s(other.s),
										cdfFunc(other.cdfFunc), 
										pdfFunc(other.pdfFunc), 
										randFunc(other.randFunc) {}	
*/			

ParameterPriorNormal::~ParameterPriorNormal() {}
			
double ParameterPriorNormal::cdf(const double x) const
{
	return cdfFunc(x);
	
}

double ParameterPriorNormal::pdf(const double x) const
{
	return pdfFunc(x);
	
}
			
double ParameterPriorNormal::rand() const
{
	return r_rand_ptr->gaussian_double(m, s);
	
}


std::string ParameterPriorNormal::toString() const
{
	std::ostringstream stm;
	
	stm << "Normal prior: mean = " << m << ", sigma = " << s;
	return stm.str();
}
