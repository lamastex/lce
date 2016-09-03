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
\brief NullTransformer declarations.
*/
   
#ifndef __REG_TRANSFORMER_NULL_HPP__
#define __REG_TRANSFORMER_NULL_HPP__

#include "reg_transformer.hpp"



namespace reg {
	
	/*! \brief A class to perform a null transformation.
	 
	 This class can be used when no transformation is required
	 
	 See the interface Transformer.*/
	class NullTransformer : public Transformer {
		
		public:
		
			/*! \brief No-args constructor.*/
			explicit NullTransformer();
		
			/*! \brief Copy constructor.*/
			explicit NullTransformer(const NullTransformer& other);
		
			virtual ~NullTransformer();
	
			/*! \brief No transform of the given block of data.
			 
			\param prior A reference to the data.
			\return \a prior by reference.
			\post The data is unchanged.*/	
			std::vector< std::vector <double> >& data_transform( 
				std::vector< std::vector <double> >& prior);
			
			/*! \brief No transform of the given block of data.
			 
			\param posterior A reference to the data.
			\return \a prior by reference.
			\post The data is unchanged.*/	
			std::vector< std::vector <double> >&  data_untransform(
				std::vector< std::vector <double> >& posterior);
				
		protected:
			NullTransformer& operator=(NullTransformer tmp);
				
		private :
			
	};
}
#endif
