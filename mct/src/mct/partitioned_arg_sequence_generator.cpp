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
\brief Definitions for a class holding ARG data and additional site
* information.
*/

#include "partitioned_arg_sequence_generator.hpp"

#include "mct_locus_map.hpp"
#include "ind_site_sequence.hpp"

#include <sstream>
#include <stdexcept>

//#define MYDEBUG
#if defined (MYDEBUG)
	#include <iostream>
	
#endif

using namespace mct;

// turn the arg into partitions
PartitionedARGSequenceGenerator::PartitionedARGSequenceGenerator(
			int nsites,
			const boost::shared_ptr < const PopulationLabeledARG > _lh) 
		: n_sites(nsites), labeled_hist(_lh)
{
		if (!_lh) { // _lh is null
			throw std::invalid_argument(
			"PartitionedARGSequenceGenerator::PartitionedARGSequenceGenerator(int, boost::shared_ptr < PopulationLabeledARG >)");
		}
		
		initialise_partitions(nsites);
}

#if(0)
PartitionedARGSequenceGenerator::PartitionedARGSequenceGenerator(const PartitionedARGSequenceGenerator& other)
	: n_sites(other.n_sites), _pop(other._pop), 
		hist(other.hist), p_set(other.p_set)
	{};

PartitionedARGSequenceGenerator& PartitionedARGSequenceGenerator::operator=(PartitionedARGSequenceGenerator tmp)
{
	std::swap( n_sites, tmp.n_sites );
	std::swap( _pop, tmp._pop );
    std::swap( hist, tmp.hist );
	std::swap( p_set, tmp.p_set );

	return *this;
	
}
#endif

/*
 * Collate sequence data from Partitions into a per-sample sequence.
 *
 * This calls evolve_arg to make sure all the sites have been evolved through.
 */
boost::shared_ptr <MCTSampleSequence> PartitionedARGSequenceGenerator::createSampleSeqs(
			const boost::shared_ptr < ModelSet >& models)
{
	
	if ( getNsites() != static_cast<int>(models->getNsites()) ) {
		throw std::invalid_argument(
		std::string("PartitionedARGSequenceGenerator::createSampleSeqs(const boost::shared_ptr < ModelSet >&)")
		+ std::string(":\nSites covered by models != sites in ARG") );
	}
	
	evolve_arg(models);
		
	PartitionSet::iterator i_part;
	
	size_t nsam = getNsam();
	std::vector<IndividualSiteSequence> sample_seq(nsam);
	
	std::vector < size_t > popNumbers 
							= labeled_hist->getPopulationNumbers();

	// each site sequence is nsites long - start by making sure we
	// have all this space to avoid reallocations
	// and then for each site sequence, makes sure we clear it
	// thus we will have enough space in each but no data
	
	for (size_t a=0; a < nsam; ++a) {
		
		// do the population labels
		sample_seq[a].pop = popNumbers[a];
					
		sample_seq[a].seq = site_seq(n_sites);
		(sample_seq[a].seq).clear();
		
		//mct::site_seq temp_seq;
		for (i_part = p_set.begin(); i_part != p_set.end(); ++i_part)
		{
			/* concat i_part->seqs[a] onto the current sequence */
			((sample_seq[a]).seq).insert(((sample_seq[a]).seq).end(),
							i_part->seqs[a].begin(), i_part->seqs[a].end());
		}
		//sample_seq.push_back(temp_seq);
	}

	// construct the sample sequence supplying the locus map from the models
	return boost::shared_ptr < MCTSampleSequence> (
				new MCTSampleSequence(sample_seq, 
				models->getLocusMap() ) );
}



std::string PartitionedARGSequenceGenerator::newick_string() const
{
	PartitionSet::const_iterator it;
	std::ostringstream stm;
	for (it = p_set.begin(); it < p_set.end(); ++it) {
		stm << it->newick_string();
	}
	return stm.str();
}

int PartitionedARGSequenceGenerator::getNsam() const
{
	return labeled_hist->getNsam();
}

int PartitionedARGSequenceGenerator::getNsites() const
{
	return n_sites;
}

std::string PartitionedARGSequenceGenerator::toString() const
{
	PartitionSet::const_iterator it;
	std::ostringstream stm;
	stm << "PartitionedARGSequenceGenerator: nsam = " << getNsam() << "\n";
	stm << "Partitions [\n";
	for (it = p_set.begin(); it < p_set.end(); ++it) {
		stm << *it << "\n";
	}
	stm << "]" << std::endl;
	return stm.str();
}


/*
 * Makes repeated calls to partition::evolve_tree().
 */
void PartitionedARGSequenceGenerator::evolve_arg(
				const boost::shared_ptr < ModelSet >& models)
{
	
	PartitionSet::iterator i_part;

    /* Evolve each tree in the arg */
    for (i_part = p_set.begin(); i_part < p_set.end(); ++i_part)
    {
        i_part->evolve_tree(models);

		#ifdef MYDEBUG
			std::cout << "This sequenced partition is " << i_part->toString() << std::endl;
		#endif
    
    }
}

void PartitionedARGSequenceGenerator::initialise_partitions(int nsites)
{									
	Sequence::arg::iterator itera, itera_next;

	//efficient clear and size
	PartitionSet((labeled_hist->getARG())->size()).swap(p_set);

	//p_set.clear();
	//p_set.resize(hist.size());
	
	PartitionSet::iterator i_part = p_set.begin();

    for (itera = labeled_hist->trees_begin(); itera != labeled_hist->trees_end(); ++itera, ++i_part)
    {
        //need a better way of dereferencing the (n+1)th value
        itera_next = itera;
        ++itera_next;

        if (itera_next == labeled_hist->trees_end())
            i_part->site_end = nsites;
        else
            i_part->site_end = itera_next->beg;

        i_part->site_start = itera->beg;
        i_part->num_sites = i_part->site_end - i_part->site_start;

        i_part->tree = &*itera;

        /* clear and allocate correct space to sequences
		This is now done by the partitions themselves */
        //i_part->initialise_sequences();
		
    }
}


// Get outstream representation of a PartitionedARGSequenceGenerator 
std::ostream& mct::operator<<(std::ostream & s, const PartitionedARGSequenceGenerator &p)
{
    s << p.toString();
    return s;
}

