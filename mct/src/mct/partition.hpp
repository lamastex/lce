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

/*! \file     partition.hpp

\brief Declarations for a Partition structure holding data for one 
section of DNA in which no recombination takes place.
*/

#ifndef _INC_PARTITION_H
#define _INC_PARTITION_H

#include "site_seq.hpp"
#include "modelset.hpp"

#include <Sequence/Coalescent/SimTypes.hpp>


#include <boost/smart_ptr.hpp>

#include <iostream>
#include <vector>


namespace mct {
	
	//forward declaration
	class ModelSet;


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

	* Sites, `-', are encoded by an integer, which represents an arbitrary value in
	  the state-space of a model. e.g. 0=`A', 1=`C', 2=`G', 3=`T'. This mapping is
	  defined by the model used at that site.

	* Each site (i.e. a column of `-'s in {1}) has a model that defines how the
	  site in a child node is generated given the site in its parent node and the
	  time between the nodes.
	*/

	/*!
	\brief A structure storing meta-data for one marginal tree in an ARG.

	Each partition relates to a section of DNA in which no recombination
	takes place.  Each partition has a one-to-one relationship
	with a marginal tree in an ARG.

	The meta-data includes the position of the contiguous sub-section of sites
	within the entire DNA sequence to which the Partition relates.
	*/
	/*
	Site numbering is from 0 - nsites-1 where nsites is the total number of
		sites in the entire section of DNA.*/
	struct Partition {

		/*! \brief The marginal tree associated with this partition.
		
		(see \ref libsequence_marginal "Sequence::marginal"). */
	 	Sequence::marginal * tree; 
									
		
		/* The numerical value of the starting site
		for the partition within the entire section of DNA.*/
		int site_start;   /* \note site_start = tree->beg */
		
		int site_end;     /* \note site_end = next_tree->beg */
		
		int num_sites;    /* \note num_sites = site_end - site_start */
		
		/*! \brief An ordered collection of 
		 \link mct::site_seq site_seqs sequences\endlink,
		each sequence 
		corresponding to a node in the \link mct::Partition::tree tree\endlink 
		for this Partition and each sequence covering as many sites as
		this Partition relates to.
		
		\note Ordering is from leaves to root, so that the site_seqs
		relating to the leaf nodes (sample individuals) are at the start
		of the collection. */
		/* 
		The following should hold for seqs: \n
		 - length(seqs) = number of nodes in marginal tree \n
		 - length(seqs[i]) = num_sites
		*/
		std::vector<mct::site_seq> seqs;
		
		/*! \brief Return a string representation of this Partition. */
		std::string toString() const;
		
		/*! \brief Return a string representation of the 
		\link mct::Partition::tree tree\endlink for 
		this Partition in Newick format. */
		std::string newick_string() const;

		/*! \brief Evolve the DNA sequences \link mct::Partition::seqs seqs\endlink
		for the nodes in the \link mct::Partition::tree tree\endlink using
		the \link mct::Model Models\endlink in the ModelSet \a models which apply to the sites
		covered by this Partition.
		
		\param models The total ModelSet covering the whole section of DNA
		of which this Partition relates to a part.
		\post The DNA sequences \link mct::Partition::seqs seqs\endlink have
		been filled in using the \link mct::Model Models\endlink in the ModelSet \a models which apply 
		to the sites covered by this Partition.
		*/
		void evolve_tree(
				const boost::shared_ptr < ModelSet >& models);
		
		/*! \brief Initialise the sequences \link mct::Partition::seqs seqs\endlink
		for the nodes in the \link mct::Partition::tree tree\endlink.
		
		Sets the size of \link mct::Partition::seqs seqs\endlink
		and the length of each site_seq within it. 
		*/
		void initialise_sequences();
	};


	

	/* Get outstream representation of a Partition. */
	std::ostream& operator<<(std::ostream & s, const Partition &p);

	/* Get outstream representation of a list of site_seq. */
	std::ostream & operator<<(std::ostream & s, 
						const std::vector<mct::site_seq> &p);

} // end namespace mct

#endif /* _INC_PARTITION_H */
