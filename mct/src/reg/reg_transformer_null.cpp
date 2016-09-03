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
\brief NullTransformer definitions.
*/

#include "reg_transformer_null.hpp"


using namespace reg;

NullTransformer::NullTransformer() : Transformer("null") {}

NullTransformer::NullTransformer(const NullTransformer& other)
: Transformer(other.type) {}

	
NullTransformer::~NullTransformer() {}

std::vector< std::vector <double> >& NullTransformer::data_transform(
				std::vector< std::vector <double> >& prior)
{
		return prior;
}
									
std::vector< std::vector <double> >&  NullTransformer::data_untransform(
				std::vector< std::vector <double> >& posterior)
{
		return posterior;
}
	
