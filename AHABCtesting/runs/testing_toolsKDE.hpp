/* MCT - Markov Chains on Trees.


   Copyright (C) 2011, 2012 Jennifer Harlow

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
\brief Tools for doing various KDE routines

*/

#ifndef _INC_TESTINGTOOLS_KDE_H
#define _INC_TESTINGTOOLS_KDE_H


#include "piecewise_constant_function.hpp"
#include "cxsc.hpp"

#include <boost/shared_ptr.hpp>

#include <vector>
#include <string>



// get an approximation to a kde
boost::shared_ptr < subpavings::PiecewiseConstantFunction > getKDEapprox(
				const std::vector < std::vector < double > >& container,
				const cxsc::ivector& box,
				const std::string& prefix,
				const std::string& spec,
				cxsc::real maxMeasure,
				const std::string& logFilename);

// just do a kde and output
void doKDE(
				const std::vector < std::vector < double > >& container,
				const cxsc::ivector& box,
				const std::string& prefix,
				const std::string& spec,
				cxsc::real maxMeasure,
				size_t id,
				const std::string& logFilename);			


	
#endif
