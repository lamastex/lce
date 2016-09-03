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
\brief Transformer declarations.
*/
   
#ifndef __REG_TRANSFORMER_HPP__
#define __REG_TRANSFORMER_HPP__

#include <vector>
#include <string>

namespace reg {
	
	/*! \brief An abstract class to provide a strategy for 
	 bulk transformations 
	 and untransformations on blocks of data.
	 
	The data to be transformed is presented as vectors of vectors.*/
	class Transformer {
		
		public:
		
			explicit Transformer(const std::string& ts);
			
			explicit Transformer(const Transformer& other);
		
			virtual ~Transformer();
	
			/*! \brief Transform the given block of data.
			 
			 \note This method is not \c const to allow
			 concrete Tansformer classes to use it to
			 set data members.
			 
			\param prior A reference to the data to be transformed.
			\return \a prior by reference.*/	
			virtual std::vector< std::vector <double> >& data_transform( 
				std::vector< std::vector <double> >& prior) = 0;
			
			/*! \brief Apply the inverse of the transformation 
			 process to untransform the given block of data.
			 
			 \note This method is not \c const to allow
			 concrete Tansformer classes to use it to
			 set (or unset) data members.
			 			 
			\param posterior A reference to the data to be untransformed.
			\return \a posterior by reference.*/	
			virtual std::vector< std::vector <double> >&  data_untransform(
				std::vector< std::vector <double> >& posterior) = 0;
				
			virtual std::string getType() const;
			
		protected:
		
			explicit Transformer();
			
			std::string type;
	};
}
#endif
