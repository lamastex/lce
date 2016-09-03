/* MCT - Markov Chains on Trees.

   Copyright (C) 2009 Brendan Bycroft <brb44@student.canterbury.ac.nz>

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
\brief Declarations for ARGmsConverter class to convert ms args to ours.
*/

#ifndef __ARG_MS_CONVERTER_NEW_HPP__
#define __ARG_MS_CONVERTER_NEW_HPP__

//#include "arg_parser.hpp"
#include <mscplusplus/arg.hpp>

#include <mscplusplus/population_structure.hpp>

#include <mscplusplus/event.hpp>

#include <Sequence/Coalescent/SimTypes.hpp>

#include <boost/smart_ptr.hpp>

#include <vector>

/* note: At present ARGmsConverter is a static class.  Each use 
	is a separate communication with the code (hudson_ms_arg_builder...)
	used to generate the ARG.  We may be able 
	to alter this so that ARGmsConverter is not static but holds
	the HudsonMSArgFactory with most of its properties
	already set up, and then use this with additional parameters in each 
	call.
	*/
	


namespace mct {
	
	/*! \brief Static class to convert from ms representation of args to
	 the libsequence args.
	 */
	
	class ARGmsConverter {
		
		public:
		
			
			/*! \brief Use given arguments to pass to ms-style 
			 routine to make an ms-style ARG and then convert
			 this to the libsequence-style ARG.
			 * 
			 \param nsites The total number of sites this ARG represents.
			 \param rho crossover parameter for a locus.
			 @todo specify rho better - relationship with nsites - ms docs.
			 \param popnsams An ordered collection of sample sizes
			 for the sub-populations for which the ARG is to 
			 be generated.
			 \param popsizeratios An ordered collection of population size
			 ratios for the sub-populations for which the ARG is to 
			 be generated.
			 \param popgrowthrates An ordered collection of growth rates
			 for the sub-populations for which the ARG is to 
			 be generated.
			 \param migration_matrix A vector representation of the 
			 migration matrix for the sub-populations for which the ARG 
			 is to be generated.  Diagonals are ignored. 
			 \param eList A sequence of demographic events to happen to
			 the sub-populations.
			 \param seed A seed to be used for random number generation
			 for simulating the ARG.
			 \return A pointer to the new Sequence::arg.*/
			static boost::shared_ptr < Sequence::arg > makeArg(
				size_t nsites,
				double rho,
				const boost::shared_ptr< const hudson_ms::PopulationStructure > pop,
				int seed);
		
		protected:
		
		private:
			ARGmsConverter();
	};

} //end namespace mct

#endif
