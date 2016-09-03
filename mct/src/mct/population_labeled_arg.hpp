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
\brief Declarations for PopulationLabeledARG class holding ARG
and population labels for each individual in the ARG.
*/

#ifndef _INC_POPULATION_ARG_H
#define _INC_POPULATION_ARG_H

#include <boost/smart_ptr.hpp>

#include <Sequence/Coalescent/SimTypes.hpp>

#include <string>
#include <vector>

/* */

namespace mct {
	
	/*! \brief Class holding an ARG and labels for the individual
	 in the ARG.
	 * 
	 An ARG is an Ancestral Recombination Graph for a sample 
	of nsam individuals over nsite sites.  An ARG consists of one
	or more trees showing an ancestry for the present day individuals.
	An ARG can have more than one tree if recombination events can
	can occur between sites.  

	The %PopulationLabeledARG also contains information about 
	the populations that the individuals in the ARG come from.
	
	Objects using the %PopulationLabeledARG to generate sequences
	can use the labels to label the sequences. 
	* 
	\internal 
	The PopulationLabeledARG  could have held a pointer to the whole
	 population structure rather than the collection of labels.  This 
	 might actually be better.  I did it this way on the basis that the 
	 individual samples only need to be labeled by sub-population, not 
	 to know about the entire population structure.  In the present 
	 set up (June 2011) we don't actually use the individual sequence
	 labels directly at all anyway because sequences just go into 
	 the multi-loci objects (after compression to SNP arrays in the case
	 of nucleotide sequences) and the multi-loci objects do have a pointer
	 to the entire population structure.  However, using labeled sequences 
	 does mean that for different uses we can still access the population
	 labels.  We may decide tht the whole labelling thing is unnecessary 
	 but I'll leave it in for the moment. 
	 */
	class PopulationLabeledARG {
		
		public:
		
		typedef Sequence::arg::iterator arg_iterator;
		
		/*! \brief Constructor taking a Libsequence::ARG
		and an ordered collection of population labels for
		each individual in the ARG.
		
		\param _h A pointer to a Sequence::arg ARG.
		\param _pl An ordered collection of population labels, one for 
		each individual in the total sample, showing which 
		sub-population each individual comes from.
		\pre The population labels should be appropriate for 
		the ARG:  no checks are done on this in the constructor
		but subsequent use of the %PopulationLabeledARG will 
		probably result in rubbish, exceptions or segmentation faults 
		if this pre-condition is violated.
		*/
		PopulationLabeledARG(
			boost::shared_ptr < Sequence::arg > _h,
			const std::vector < std::string >& _pl);
		
		virtual ~PopulationLabeledARG() {}
		
		/*! \brief Return an ordered collection of the population
		labels for the individuals in the population whose genealogical
		history is contained in the arg.
		
		if subpop i has label \f$ label_i \f$ and \f$ n_i \f$ individuals, \f$ i = 0,\,\ldots,\,m-1 \f$
		for \f$ m \f$ subpopulations, then this method returns an
		ordered collection of labels
		\f$ label_0^{0},\,\ldots,\,label_0^{[n_0-1]},\,\, label_1^{0},\,\ldots,\,label_1^{[n_1-1]},\,\,\dots,\,\,label_{m-1}^{0},\,\ldots,\,label_{m-1}^{[n_{m-1}-1]} \f$
		i.e. \f$ n_0 \f$ label_0 then \f$ n_1 \f$ label_1, then \f$ n_2 \f$ label_2 , etc, and finally \f$ n_{m-1} \f$ label_{m-1}.
		
		\return a collection of population labels for the individuals in the arg. 
		*/
		std::vector< std::string > getPopulationLabels() const;
		
		/*! \brief Return an ordered collection of the population
		numbers for the individuals in the population whose genealogical
		history is contained in the arg.
		
		if subpop i has \f$ n_i \f$ individuals, \f$ i = 0,\,\ldots,\,m-1 \f$
		for \f$ m \f$ subpopulations, then this method returns an
		ordered collection of numbers
		\f$ 0^{0},\,\ldots,\,0^{[n_0-1]},\,\, 1^{0},\,\ldots,\,1^{[n_1-1]},\,\,\dots,\,\,(m-1)^{0},\,\ldots,\,(m-1)^{[n_{m-1}-1]} \f$
		i.e. \f$ n_0 \f$ 0' then \f$ n_1 \f$ 1's, then \f$ n_2 \f$ 2's , etc, and finally \f$ n_{m-1} \f$ (m-1)'s.
		
		\return a collection of population numbers for the individuals in the arg. 
		*/
		std::vector< size_t > getPopulationNumbers() const;
		
		/*! \brief Get a copy of pointer to the ARG. */
		boost::shared_ptr < const Sequence::arg > getARG() const;
		
		/*! \brief Get the number of individuals in the sample
		covered by the ARG (or, equivalently, covered by the
		 population labels). */
		size_t getNsam() const;
		
		/*! \brief Get number of marginal trees in the ARG. */
		size_t getMarginalTreeCount() const;
		
		/*! \brief Get a string representation of this. */
		std::string toString() const;
		
		/*! \brief Return a const iterator to the first 
		marginal tree in the ARG.*/
		arg_iterator trees_begin() const;
		
		/*! \brief Return a const iterator to the
		'just past the end' marginal tree in the ARG.
		
		\note The idiom for iterating through the margingal
		trees of a %PopulationLabeledARG  pla is 
		for (PopulationLabeledARG::const_arg_iterator it = pla.begin();
					it < pla.end();
					++it) {...}.
		*/
		arg_iterator trees_end() const;
		
		
		
		protected:
		
		PopulationLabeledARG();
		
		PopulationLabeledARG(const PopulationLabeledARG& other);
		
		PopulationLabeledARG& operator=(PopulationLabeledARG tmp);
		
		const boost::shared_ptr < Sequence::arg > hist;
		
		const std::vector < std::string > populationLabels;
		
		private:
		
		
	};
	
	/* Get outstream representation of a PopulationLabeledARG. */
	std::ostream& operator<<(std::ostream & s, const PopulationLabeledARG &pa);


} // end namespace mct

#endif
