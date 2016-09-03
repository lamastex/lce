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
\brief Tools for doing various mcmc routines

*/

#ifndef _INC_TESTINGTOOLS_GR_H
#define _INC_TESTINGTOOLS_GR_H



#include "piecewise_constant_function.hpp"

#include "automcmc_factory.hpp"


#include <boost/shared_ptr.hpp>


#include <string>
#include <vector>

	boost::shared_ptr < subpavings::PiecewiseConstantFunction > getMCMCGRaverage(
				const std::vector < std::vector < double > >& container,
				const cxsc::ivector& box,
				int np,
				int ns,
				const std::string& prefix,
				const std::string& spec,
				ahabc::AutoMCMCFactory::AUTOTYPE autotype,
				const std::string& logFilename,
				const std::string& memLogFilename);
	
	boost::shared_ptr < subpavings::PiecewiseConstantFunction > getMCMCGRaverage(
				const std::vector < std::vector < double > >& container,
				const cxsc::ivector& box,
				real minVol,
				int ns,
				const std::string& prefix,
				const std::string& spec,
				ahabc::AutoMCMCFactory::AUTOTYPE autotype,
				const std::string& logFilename,
				const std::string& memLogFilename);
				
	boost::shared_ptr < subpavings::PiecewiseConstantFunction > getMCMCGRaverage(
				const std::vector < std::vector < double > >& container,
				const cxsc::ivector& box,
				int ns,
				const std::string& prefix,
				const std::string& spec,
				ahabc::AutoMCMCFactory::AUTOTYPE autotype,
				const std::string& logFilename,
				const std::string& memLogFilename);
				
	boost::shared_ptr < subpavings::PiecewiseConstantFunction > getMCMCGRaverage(
				const std::vector < std::vector < double > >& container,
				const cxsc::ivector& box,
				real minVol,
				int np,
				int ns,
				const std::string& prefix,
				const std::string& spec,
				ahabc::AutoMCMCFactory::AUTOTYPE autotype,
				const std::string& logFilename,
				const std::string& memLogFilename);				

	
	
#endif
