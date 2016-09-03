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
\brief ParameterPriorUniform definitions.
*/

#include "parameter_prior_uniform.hpp"


#include <boost/bind.hpp>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>


#include <sstream>  // to be able to manipulate strings as streams

using namespace aabc;

ParameterPriorUniform::ParameterPriorUniform(
					double _a, double _b,
					const boost::shared_ptr < mct::PRNGen >& _r_ptr)
					: 	a(_a), b(_b),
						cdfFunc (boost::function<double (double)>
							( boost::bind(gsl_cdf_flat_P, _1, _a, _b) )
							),
						pdfFunc (boost::function<double (double)>
							( boost::bind(gsl_ran_flat_pdf, _1, _a, _b) )
							),
						r_rand_ptr(_r_ptr )
						{}
/*	
ParameterPriorUniform::ParameterPriorUniform(
								const ParameterPriorUniform& other)
										: a(other.a), b(other.b),
										cdfFunc(other.cdfFunc), 
										pdfFunc(other.pdfFunc), 
										randFunc(other.randFunc) {}	
*/			

ParameterPriorUniform::~ParameterPriorUniform() {}
			
double ParameterPriorUniform::cdf(const double x) const
{
	return cdfFunc(x);
	
}

double ParameterPriorUniform::pdf(const double x) const
{
	return pdfFunc(x);
	
}
			
double ParameterPriorUniform::rand() const
{
	return r_rand_ptr->uniform_double(a ,b);
	
}

std::string ParameterPriorUniform::toString() const
{
	std::ostringstream stm;
	
	stm << "Uniform prior: U(" << a << ", " << b << ")";
	return stm.str();
}
