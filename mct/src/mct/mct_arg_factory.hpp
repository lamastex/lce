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
\brief ArgFactory declarations.
*/

#ifndef _INC_MCT_ARG_FACTORY_H
#define _INC_MCT_ARG_FACTORY_H

#include <mscplusplus/population_structure.hpp>

#include "population_labeled_arg.hpp"

#include <boost/smart_ptr.hpp>

#include "prng.hpp"

#include <string>

namespace mct {

	/*! \brief Class to create 
	\link PopulationLabeledARG PopulationLabeledARGs\endlink.
	
	The ARGFactory acts as the single common intermediary between many
	types wanting to make 
	PopulationLabeledARGs, like the samplers, and the many methods
	which may be used to make ancestral recombination graphs (ARGs)
	in different ways:
	<ul> 
	<li>Using the %ARGFactory, a type wanting to
	make a %PopulationLabeledARG only has to know about the interface
	of the %ARGFactory, not the details of the different ways 
	in which an ARG can be be created.</li>
	<li>The %ARGFactory is the common way in which information needed
	to make an ARG can be extracted from a population structure for
	any other type needing that information to make the ARG.  i.e.
	only the %ARGFactory needs to know about the methods a 
	hudson_ms::PopulationStructure provides for accessing this information,
	and the %ARGFactory can if necessary manipulate this into the format
	needed by the types that use the information directly.</li>
	<li>Only the %ARGFactory needs to know how to construct a 
	%PopulationLabeledARG given some kind of representation of an ARG 
	(created in one of 
	many possible ways) and a population structure.</li>
	
	At present the %ARGFactory can construct a PopulationLabeledARG
	using the mscplusplus library, \ref libsequence, or using a newick-encoded
	description of the required ARG in a file.
	
	\note ARGFactory is a static class, i.e. a ARGFactory object 
	cannot be instantiated and the class only contains static members.
	 */
	class ARGFactory {
		
		
		public:
			
			/*! \brief Create an ARG using Libsequence.
			 * 
			This method can only create an ARG for a population structure
			consisting of a single population.
			
			\param  pop The population structure of the sample for the ARG.
			\param  nsites The total number of sites covered by the ARG.
			\param  growth The growth rate per locus 
					(i.e. over all nsites).
			\param  rho The crossover parameter for a locus 
					(i.e. over all nsites).
			\param r_ptr A pointer to a random number generator.
			\return A pointer to a PopulationLabeledARG.
			\pre The population structure pointed to by \a pop must have
			only one sub-population.*/
			static boost::shared_ptr < PopulationLabeledARG > 
				makeHistFromLibseq(
				const boost::shared_ptr< hudson_ms::PopulationStructure > pop,
				const size_t nsites,
				const double growth,
				const double rho,
				const boost::shared_ptr < const PRNGen >& r_ptr);

			/*! \brief Create an ARG using MS.
			
			\note The population growth rates are assumed to be set
			in the population pointed to by \a pop.
			
			\param  pop The population structure of the sample for the ARG.
			\param  nsites The total number of sites covered by the ARG.
			\param  rho The crossover parameter for locus 
					(i.e. over all nsites).
			\param r_ptr A pointer to a random number generator.
			\return A pointer to a PopulationLabeledARG. 	*/
			static boost::shared_ptr < PopulationLabeledARG > 
				makeHistFromMS(
				const boost::shared_ptr< const hudson_ms::PopulationStructure > pop,
				const size_t nsites,
				const double rho,
				const boost::shared_ptr < const PRNGen >& r_ptr);

			/*! \brief Create an ARG from a given string.
			
			\param  newicktreestring A string giving an ARG in 
			newick format.
			\param  pop The population structure of the sample for the ARG.
			\param  nsites The total number of sites covered by the ARG.
			\return A pointer to a PopulationLabeledARG. 	*/
			static boost::shared_ptr < PopulationLabeledARG > 
				makeHistFromNewickString(
				const std::string& newicktreestring,
				const boost::shared_ptr< hudson_ms::PopulationStructure > pop,
				const size_t nsites);
			
		protected:
		
		/*
		static std::string makeMSTreesString(
				const boost::shared_ptr< hudson_ms::PopulationStructure > pop,
				const double growth,
				const double rho,
				const size_t nsites,
				int seed);
		*/
			
		private:
		
		ARGFactory(); // can't instantiate class

	};

} // end namespace mct

#endif /* _INC_ARG_FACTORY_H */

