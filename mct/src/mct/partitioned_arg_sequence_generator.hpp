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
\brief Declarations for PartionedARG class holding ARG and site 
meta-data.
*/

#ifndef _INC_PARTITIONED__SEQ_GEN_ARG_H
#define _INC_PARTITIONED__SEQ_GEN_ARG_H


#include "partition.hpp"
#include "modelset.hpp"
#include "mct_sample_sequence.hpp"
#include "population_labeled_arg.hpp"

#include <boost/smart_ptr.hpp>

#include <Sequence/Coalescent/SimTypes.hpp>


/*
A graphical description of how the site data is layed out.

- = site
p = partition
t = tree
| = partition boundary
s = sample
n = node

A sample size of 3, with 4 different trees (partitions) over 25 sites {1}
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
--------|--------|---|------ s
--------|--------|---|------ s
--------|--------|---|------ s
    p       p      p     p
    t       t      t     t

The 5 nodes in the first partition {2}
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
-------- s      <-- the first three are the same as the sample in the
-------- s          first partition
-------- s
-------- n      <-- the last two are internal nodes of the tree. i.e. at
-------- n          coalescent events

======================================
* The Partition struct holds the data in {2} and also data about where the
  partition starts and finishes.

* The markov models will run through the trees, and generate the values for 
  all the `-'s in {2}.
  
* Once all of the trees have been filled with data, the sample nodes (e.g. the
  first 3 in {2}, will have their dna strand concat'ed with the others to
  result in {1}, with one strand per sample.

* Sites, `-', are encoded by an integer, which represents an arbitrary value in
  the state-space of a model. e.g. 0=`A', 1=`C', 2=`G', 3=`T'. This mapping is
  defined by the model used at that site.

* Each site (i.e. a column of `-'s in {1}) has a model that defines how the
  site in a child node is generated given the site in its parent node and the
  time between the nodes.

* The models can be placed arbitrarily, and are independent of the partition
  boundaries.

*/

namespace mct {
	
	/*! \brief Class holding an ARG and site meta-data for a section of DNA
	(locus) using this with a ModelSet to generate a sequence for a sample.

	The evolution of a section of DNA (locus), including possible recombination
	events, is described by an ARG.  
	* 
	A PartitionedARGSequenceGenerator contains the ARG and 
	meta-data related to the total DNA section.  Recombination events 
	divide the overall section of DNA into a
	PartionSet of \link Partition mct Partitions\endlink.  The PartitionSet
	describes the \link Partition mct Partitions\endlink of the 
	overall DNA section.  Each Partition relates to one marginal tree in the 
	total ARG for the entire DNA section.

	A PartitionedARGSequenceGenerator is used to generate a sample of sequence data for
	the present-day individuals in the ARG given a set of mct models 
	(ModelSet) describing the model to be used for each site in the locus. 

	 */

	class PartitionedARGSequenceGenerator {
		
		public:
		
		/*! \brief Constructor taking a number of sites
		 and a population structure and a pointer to an ARG.
		
		\param  nsites The total number of sites covered by the ARG.
		\param _lh A pointer to a  PopulationLabeledARG.
		\pre _lh is assumed to be a non-empty pointer.
		*/
		PartitionedARGSequenceGenerator(
			int nsites,
			const boost::shared_ptr < const PopulationLabeledARG > _lh);
		
		
		virtual ~PartitionedARGSequenceGenerator() {}
		
		/*! \brief Create a sample of DNA data sequences.
		
		Creates a MCTSampleSequence with one sequence for each 
		present-day
		individual in an ARG.
		
		\param  models The ModelSet of \link Model Models\endlink used
		to generate the sample sequences.
		\return A pointer to the new MCTSampleSequence.
		*/							
		boost::shared_ptr <MCTSampleSequence> createSampleSeqs(
				const boost::shared_ptr < ModelSet >& models);

		/*! \brief Return a representation of this PartitionedARGSequenceGenerator in 
		Newick format. */
		std::string newick_string() const;
		
		/*! \brief Get the sample size covered by the ARG represented by 
		this PartitionedARGSequenceGenerator. */
		int getNsam() const;
		
		/*! \brief Get the number of sites covered by the ARG 
		 represented by this PartitionedARGSequenceGenerator. */
		int getNsites() const;
		
		/*! \brief Get a string representation of this. */
		std::string toString() const;
		
		protected:
		
		PartitionedARGSequenceGenerator();
		
		PartitionedARGSequenceGenerator(const PartitionedARGSequenceGenerator& other);
		
		PartitionedARGSequenceGenerator& operator=(PartitionedARGSequenceGenerator tmp);
		
		/*! \brief An ordered collection of Partitions. */
		typedef std::vector<Partition> PartitionSet;
		
		void evolve_arg(
				const boost::shared_ptr < ModelSet >& models);
		
		void initialise_partitions(int nsites);
		
		
		int n_sites;
		
		const boost::shared_ptr < const PopulationLabeledARG > labeled_hist;
		
		PartitionSet p_set;
		
		
		private:
		
		
	};
	
	/* Get outstream representation of a PartitionSet. */
	std::ostream& operator<<(std::ostream & s, const PartitionedARGSequenceGenerator &p);


} // end namespace mct

#endif
