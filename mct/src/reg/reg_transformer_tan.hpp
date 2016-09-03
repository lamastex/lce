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
\brief TanTransformer declarations.
*/
   
#ifndef __REG_TRANSFORMER_TAN_HPP__
#define __REG_TRANSFORMER_TAN_HPP__

#include "reg_transformer.hpp"


namespace reg {
	
	/*! \brief A class to perform bulk transformations
	 and untransformations on blocks of data
	 using the tangent transformation.
	 
	 The tangent transformation is the transformation proposed by 
	 Hamilton, G., Stoneking, M., Excoffier., L. (2005), "Molecular
	 analysis reveals tighter social regulation of immigration in
	 patrilocal populations than in matrilocal populations."
	 PNAS, 102, pp. 746-7480.
	 
	 Note that Thornton uses mins and maxs from data referred to as 'prior', ie
	 the parameter values to be transformed before being regressed.
	 This is not what Hamilton et al use - they use (or say they do) the 
	 actual min and max bounds on the prior as the distribution
	 from which possible parameter values are generated.  
	 
	 This transformation guarantees that the posterior values
	 remain within the bounds of the prior.
	 
	 See the interface Transformer.*/
	class TanTransformer : public Transformer {
		
		public:
		
			/*! \brief Parameterised constructor.
			 
			 If all priors had bounds it would be easier to send in
			 the actual parameter prior set
			 object here - would be easier for user - but this 
			 creates a dependency of the transformer on the parameter
			 prior, and - most importantly - not all priors do
			 have bounds.
			 */
			explicit TanTransformer(const std::vector < double >& _mins,
				const std::vector < double >& _maxs);
			
			/*! \brief Copy constructor.*/
			explicit TanTransformer(const TanTransformer& other);
		
			virtual ~TanTransformer();
	
			/*! \brief Transform the given block of data.
			 
			 The transformation applied is the tangent transformation.
			 
			\param prior A reference to the data to be transformed.
			\return \a prior by reference.
			\pre The data must all be tangent-transformable.*/	
			std::vector< std::vector <double> >& data_transform( 
				std::vector< std::vector <double> >& prior);
			
			/*! \brief Apply the inverse of the transformation 
			 process to untransform the given block of data.
			 
			 \note the mins and maxs used for the un-transform assume are
			 those for the last use of this transformer to perform a 
			 tangent transformation; i.e. it is assumed that this
			 transformer was last used to perform the transformation
			 of the data.  If this assumption is violated, unsuitable
			 mins and maxes may be used in the un-transformation and
			 the results will therefore be incorrect. 
			 
			\param posterior A reference to the data to be untransformed.
			\return \a posterior by reference.
			\pre The data must all be inverse-tangent-transformable.*/
			std::vector< std::vector <double> >&  data_untransform(
				std::vector< std::vector <double> >& posterior);
				
				
		protected:
		
			/*! \brief No-args constructor.*/
			explicit TanTransformer();
			
		
			TanTransformer& operator=(TanTransformer tmp);
				
				
			static double tan_transform( const double x,
										const double minval,
										const double maxval);

			static double tan_untransform( 	const double y,
													const double minval,
													const double maxval);
			
			const std::vector < double > mins;
			
			const std::vector < double > maxs;
			
		private:
		
			
		};
}
#endif
