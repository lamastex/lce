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
\brief LibSeq declarations.
*/

#ifndef _INC_LIBSEQ_H
#define _INC_LIBSEQ_H

#include "prng.hpp"

#include <Sequence/Coalescent/SimTypes.hpp>
#include <Sequence/SimData.hpp>
#include <Sequence/PolyTable.hpp>


namespace mct {
	
	/** @brief A class of statics for useful routines using
	 Libsequence. 
	 */
	 

	class LibSeq {
		public:
		
			/*! \brief Create and return a Sequence::arg.
			
			\param nsam The number of samples in each ARG.
			\param nsites 	The number of sites at which mutations occur. 
			\param r_ptr A pointer to a random number generator.
			\param g 	The rate of exponential change in effective size.
			\param t_begin 	The time in the past (in units of 
			4Ne generations) at which population size change begins 
			(i.e., ends, moving forward in time).
			\param t_end 	The time in the past (in units of 
			4Ne generations) at which populations size change 
			ends (begins forward in time).
			\param rho 	The population recombination rate 4N0r. 
			\param size_at_end 	At time t_end in the past, 
			the population size is set to size_at_end. 
			\note If size_at_end = 1, the population is set to the same size 
				that is was at t=0 (i.e. the beginning of the simulation).  
				If size_at_and < 0, the population size is not adjusted at 
				t_end. In other words, it is left at whatever it grew or shrank to. 
			\return A Sequence::arg.
			*/
			static Sequence::arg initialise_hist(
				const size_t nsam,
				const int nsites,
				const boost::shared_ptr < const mct::PRNGen >& r_ptr,
				const double g,
				const double t_begin,
				const double t_end,
				const double rho,
				const double size_at_end);
			
			/*! \brief Create and return a Sequence::SimData.
			
			\param nsites 	The length of the region begin simulated.
			\param r_ptr 	A pointer to a random number generator.
			\param hist 	The ARG.
			\param theta 	The coalescent-scaled mutation rate per
			per locus (ie over all nsites).
			\return A Sequence::SimData.
			*/
			static Sequence::SimData infinite_sites_sim_data(
				const int	nsites,
				const boost::shared_ptr < const mct::PRNGen >& r_ptr,
				const Sequence::arg& hist,
				const double& theta);
			
			static void startConversion (Sequence::PolyTable& poly,
				bool haveOutgroup, unsigned outgroup, bool strictInfSites,
				size_t this_nsites);
			
			/*! \brief String representation of a 
			 Sequence::node in a marginal tree.*/
			static std::string nodeToString(const Sequence::node& snode);

			/*! \brief String representation of a 
			 Sequence::marginal tree.*/
			static std::string margToString(const Sequence::marginal& marg);

			/*! \brief String representation of a 
			 Sequence::arg.*/
			static std::string argToString(const Sequence::arg& ARG);
		
			/*! \brief String representation of a 
			 Sequence::arg in newick format.
			 
			 \param ARG The ARG to get the string representation of.
			 \param prec The precision to use for doubles in string.
			 \return A string representation of \a ARG.*/
			static std::string argNewickString(
							const Sequence::arg& ARG,
							const size_t prec = 3);
			
			
			/*! \brief String representation of an 
			 Sequence::PolyTable.*/
			static std::string polytableToString(
									const Sequence::PolyTable& simdata);			
			
		private:
		
			LibSeq();
					
			explicit LibSeq(const LibSeq& other);
			
			LibSeq& operator=(LibSeq tmp);
			
			static std::ostream & makeNewick(Sequence::arg::const_iterator m,
						std::ostream & o);
			
			static std::ostream & parens( 
							Sequence::marginal::const_iterator mi,
							const std::vector < int >& left,
							const std::vector < int >& right,
							const int & noden,
							std::ostream & o);
			
			
	};
	
	// non_member function
	bool seq_node_comp_height(const Sequence::node& n1, const Sequence::node& n2);
	
	
} // end namespace mct

#endif /* _INC_LIBSEQ_H */

