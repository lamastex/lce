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
\brief MixtureTest definitions.
*/

#include "mixture_test.hpp"

#include <sstream>  // to be able to manipulate strings as streams
#include <iterator> // for ostream_iterator


using namespace aabc;

MixtureTest::MixtureTest(double v)
		: value(v) {}

MixtureTest::~MixtureTest()
{
	
}

		
boost::shared_ptr < const mct::SummaryStatistic  > 
							MixtureTest::getSummaryStatistic() const
{
	boost::shared_ptr < mct::SummaryStatistic  > s(new mct::SummaryStatistic());
	s->pushBack(value);
	return s;
}
			
std::string MixtureTest::toString() const
{
	std::ostringstream stm;
	
	stm << value;
		
	return stm.str();
}
/*
double MixtureTest::getValue() const
{
	return value;
}
*/

