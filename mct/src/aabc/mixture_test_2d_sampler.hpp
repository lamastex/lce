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
\brief MixtureTest2DSampler declarations.
*/

#ifndef _MIXTURE_TEST_2D_SAMPLER_H
#define _MIXTURE_TEST_2D_SAMPLER_H

#include "generic_data_object_sampler.hpp"
#include "mixture_test_2d.hpp"

#include "prng.hpp"


namespace aabc {
	
	/** @brief A class capable of 
	 generating a sample of MixtureTest2D tests.
	 
	 */

	class MixtureTest2DSampler : public GenericDataObjectSampler {
		public:
		
			/*! \brief Constructor. 
			 
			\param _r_ptr A pointer to a random number generator.
			 */
			explicit MixtureTest2DSampler(
						const boost::shared_ptr < mct::PRNGen >& _r_ptr);
		
		
			virtual ~MixtureTest2DSampler();
			
			
		protected:
		
			virtual boost::shared_ptr < GenericDataObjectSet >& 
				_sample(
					boost::shared_ptr < GenericDataObjectSet >& 
							gsetptr,
					const std::size_t n_reps, 
					const boost::shared_ptr< const ParameterParticle>& p) const;

		
		
			virtual boost::shared_ptr < MixtureTest2D > genNew(
					const boost::shared_ptr < const ParameterParticle >& p)
					const;
		
			virtual void initialisePRNGs(
							const boost::shared_ptr < mct::PRNGen >& _r_ptr);
		
		private:
		
			boost::shared_ptr < mct::PRNGen > rand;
		
	};

} // end namespace aabc

#endif /* _SAMPLER_TEST_2D_H */

