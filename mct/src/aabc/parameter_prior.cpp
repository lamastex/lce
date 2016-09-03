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
\brief ParameterPrior definitions.
*/

#include "parameter_prior.hpp"

#include <iostream> // for testing output

using namespace aabc;



void ParameterPrior::test() const
{
	for (size_t i = 0; i < 10; ++i) {

		double randResult = rand(); 

		// call the bound functions 
		double pdfResult = pdf(randResult); 
		double cdfResult = cdf(randResult); 

		std::cout << "rand() = " << randResult << "\t";
		std::cout << "pdf(" << randResult << ") = " << pdfResult << "\t";
		std::cout << "cdf(" << randResult << ") = " << cdfResult << "\n";
		
	}
	std::cout << std::endl;
}
