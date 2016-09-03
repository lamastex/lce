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
\brief LogTransformer declarations.
*/
   
#ifndef __REG_TRANSFORMER_LOG_HPP__
#define __REG_TRANSFORMER_LOG_HPP__

#include "reg_transformer.hpp"



namespace reg {
	
	/*! \brief A class to perform bulk natural log transformations 
	 and untransformations on blocks of data.
	 
	 See the interface Transformer.*/
	class LogTransformer : public Transformer {
		
		public:
		
			/*! \brief No-args constructor.*/
			explicit LogTransformer();
		
			/*! \brief Copy constructor.*/
			explicit LogTransformer(const LogTransformer& other);
		
			virtual ~LogTransformer();
	
			/*! \brief Transform the given block of data.
			 
			 The transformation applied is the natural log of the data.
			 
			\param prior A reference to the data to be transformed.
			\return \a prior by reference.
			\pre The data must all be log-able.*/	
			std::vector< std::vector <double> >& data_transform( 
				std::vector< std::vector <double> >& prior);
			
			/*! \brief Apply the inverse of the transformation 
			 process to untransform the given block of data.
			 
			\param posterior A reference to the data to be untransformed.
			\return \a posterior by reference.
			\pre The data must all be exponentiable.*/		
			std::vector< std::vector <double> >&  data_untransform(
				std::vector< std::vector <double> >& posterior);
				
		protected:
			
			LogTransformer& operator=(LogTransformer tmp);
			
				
		private :
			
	};
}
#endif
