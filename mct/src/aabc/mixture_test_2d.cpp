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
\brief MixtureTest2D definitions.
*/



#include "mixture_test_2d.hpp"


#include <sstream>  // to be able to manipulate strings as streams
#include <iterator> // for ostream_iterator

#include <cmath>

using namespace aabc;

MixtureTest2D::MixtureTest2D(double v0, double v1)
{
	v = NULL;
	try {
		v = gsl_vector_calloc(2);
		gsl_vector_set(v, 0, v0);
		gsl_vector_set(v, 1, v1);
	}
	catch (...) {
		try {
			if (v != NULL) gsl_vector_free(v);
		}
		catch (...) {}
		throw;
	}
}

MixtureTest2D::~MixtureTest2D()
{
	
	try {
		if (NULL != v) gsl_vector_free(v);
	}
	catch (...) {}
	
}

boost::shared_ptr < const mct::SummaryStatistic  > 
							MixtureTest2D::getSummaryStatistic() const
{
	boost::shared_ptr < mct::SummaryStatistic  > s(new mct::SummaryStatistic());
	s->pushBack(gsl_vector_get(v, 0));
	s->pushBack(gsl_vector_get(v, 1));
	return s;
}
			
std::string MixtureTest2D::toString() const
{
		std::ostringstream stm;
		
		stm << "(" << gsl_vector_get(v, 0) << "\t" << gsl_vector_get(v, 1) << ")";
			
		return stm.str();
}

/*
const gsl_vector* MixtureTest2D::getValue() const
{
	
	return v;
}
*/

