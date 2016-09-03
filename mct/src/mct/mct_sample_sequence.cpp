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
\brief MCTSampleSequence definitions.
*/

#include "mct_sample_sequence.hpp"
#include "mct_model_data_types.hpp"
#include <Sequence/Alignment.hpp>

#include <boost/bind.hpp>

#include <algorithm> // for transform
#include <sstream>  // to be able to manipulate strings as streams
#include <iterator> // for ostream_iterator
#include <stdexcept>

using namespace mct;

//declare static data member
std::string MCTSampleSequence::nucleotides = "ACGT";

MCTSampleSequence::MCTSampleSequence(
					const std::vector < IndividualSiteSequence >& seqs,
					const boost::shared_ptr <MCTLocusMap>& lm)
		: locusMap(lm) 
{
	size_t ns = seqs.size();
	sample.resize( ns );
	popLabels.resize( ns );
	
	std::transform(seqs.begin(), seqs.end(), sample.begin(),
	boost::bind(&IndividualSiteSequence::getSeq, _1) );
	
	std::transform(seqs.begin(), seqs.end(), popLabels.begin(),
	boost::bind(&IndividualSiteSequence::getPop, _1) );
	
}

MCTSampleSequence::MCTSampleSequence(
					const std::vector < size_t >& pls,
					const std::vector < site_seq >& ss,
					const boost::shared_ptr < MCTLocusMap >& lm)
		: popLabels(pls), sample(ss), locusMap(lm) 
{
	if (sample.size() != popLabels.size() ) {
		throw std::invalid_argument(
		std::string("MCTSampleSequence::MCTSampleSequence(")
		+ std::string("const std::vector < size_t >&, ")
		+ std::string("const std::vector < site_seq >&, ")
		+ std::string("const boost::shared_ptr < MCTLocusMap >&)")
		+ std::string(":\nIncompatible number of sequences and labels") );
	}
}


MCTSampleSequence::~MCTSampleSequence()
{
	//std::cout << "Deleting MCTSampleSequence" << std::endl;
		
}
			

IndividualSiteSequence MCTSampleSequence::getIndividualSequenceAt(
										const size_t seq_index) const
{
	return IndividualSiteSequence( popLabels.at(seq_index),
									sample.at(seq_index) );
}

std::string MCTSampleSequence::toString() const
{
	std::ostringstream stm;
	
	std::vector < size_t >::const_iterator pit = popLabels.begin();
	for (MCTSampleSequence::const_iterator it
		= begin();
		it < end();
		++it, ++pit) {
	
		stm << *(pit) << "( ";
		
		std::ostream_iterator<int> out_it (stm, "\t");
		copy ( it->begin(), it->end(), out_it );
	  
		stm << ")\n";
	}
	
    return stm.str();
}

size_t MCTSampleSequence::nsam() const
{
	return sample.size();
}

std::string MCTSampleSequence::stringPhylipFormat() const
{
	std::vector < std::string > strs;
	
	strs = codeSeqs(strs);
			
	std::string s;
	for (std::vector < std::string >::iterator it = strs.begin();
			it < strs.end();
			++it) {
				
		s.append(*it);
		s += "\n";
	}

	return s;
}


boost::shared_ptr < Sequence::PolySites > 
			MCTSampleSequence::polyNucTableFormatPtr() const
{
	std::vector < std::string > strs;
	strs = codeSeqs(strs);
	
	if ( !Sequence::Alignment::validForPolyAnalysis
				(strs.begin(), strs.end()) ) {
			throw std::domain_error(
			"MCTSampleSequence::polyNucTableFormatPtr() :\nnot a valid aligned pure DNA sequence");
	}
	
	boost::shared_ptr < Sequence::PolySites > 
			polyPtr( new Sequence::PolySites(strs) ); 
			
	return polyPtr;
}

Sequence::PolySites MCTSampleSequence::polyNucTableFormat() const
{
	return *(polyNucTableFormatPtr());
}


std::vector < site_seq > MCTSampleSequence::microsatRepeatsFormat() const
{
	// check the locus map only has microsat data in it
	if ( (locusMap->size() > 1) 
						|| 
		( !locusMap->empty() && ((locusMap->at(0))->type != MICROSAT)) ) {
		throw std::domain_error(
			"MCTSampleSequence::microsatSiteSeqFormat() :\nSequence is not pure MICROSAT data");
	}
	
	std::vector < site_seq > m_format(sample.size());
	std::vector < site_seq >::iterator it = m_format.begin();
	
	for (MCTSampleSequence::const_iterator sit = begin();
			sit < end();
			++sit) {
		it->reserve(sit->size());
		it->insert(it->end(), sit->begin(), sit->end());
		++it;
				
	}
	
	return m_format;
}


std::vector < site_seq > MCTSampleSequence::siteSeqFormat() const
{
		return sample;
}

//protected

std::vector < std::string >& MCTSampleSequence::codeSeqs(
		std::vector < std::string >& strs) const
{
	if ( !sample.empty() ) {
		/*reset the vector of strings to be "" for each sample
		 but so that each "" has enough space reserved for 
		 a large whole sequence (remember we don't know that
		 integer codes will translate to single-characters though. */
		std::string tmp = "";
		tmp.reserve(2 * (sample.front()).size());
		
		std::vector < std::string >( nsam(),tmp).swap(strs);
		
		//get each block from the locus map
		// get the start, the end, and how to treat it
		
		size_t blocks = locusMap->size();
		
		for (size_t i = 0; i < blocks; ++i) {
			
			DATATYPE thisType = (locusMap->at(i))->type;
					
			if (thisType == DNA ) {
				strs = dnaCodeToString( strs, 
							(locusMap->at(i))->from, 
							(locusMap->at(i))->to );
			}
			if (thisType == MICROSAT ) {
				strs = microsatCodeToString( strs, 
							(locusMap->at(i))->from, 
							(locusMap->at(i))->to );
			}
		}
		for (std::vector < std::string >::iterator it = strs.begin();
				it < strs.end();
				++it) {
			// shrink-to-fit each string
			std::string(*it).swap(*it);
		}
	}
	
	return strs;
}


std::vector < std::string >& MCTSampleSequence::dnaCodeToString(
				std::vector < std::string>& strs, 
				const size_t from, 
				const size_t to) const
{
	MCTSampleSequence::const_iterator sam_it;
	std::vector < std::string>::iterator strs_it;
	
	for (size_t i = from; i <= to; ++i) {
		for (sam_it = begin(), strs_it = strs.begin();
				sam_it < end(), strs_it < strs.end();
				++sam_it, ++strs_it) {

			strs_it->append( 1, nucleotides.at( sam_it->at(i) ) );
		} 
					
	}
	
	return strs; // return by reference
}

std::vector < std::string >& MCTSampleSequence::rawToString(
				std::vector < std::string>& strs) const
{
	MCTSampleSequence::const_iterator sam_it;
	
	for (sam_it = begin(); sam_it < end(); ++sam_it) {
	
		std::ostringstream stm;
		std::ostream_iterator<int> out_it (stm, "");
		copy ( sam_it->begin(), sam_it->end(), out_it );
		strs.push_back(stm.str());
							
	}
	
	return strs; // return by reference
}

std::vector < std::string >& MCTSampleSequence::microsatCodeToString(
				std::vector < std::string>& strs, 
				const size_t from, 
				const size_t to) const
{
	MCTSampleSequence::const_iterator sam_it;
	std::vector < std::string>::iterator strs_it;
	
	for (sam_it = begin(), strs_it = strs.begin();
				sam_it < end(), strs_it < strs.end();
				++sam_it, ++strs_it) {
					
		std::ostringstream stm;
		stm << "."; // start with a dot
			
		for (size_t i = from; i <= to; ++i) {
			
			stm << ( sam_it->at(i) ) << ".";
		}
		strs_it->append(stm.str()); 
					
	}
	
	return strs; // return by reference
}

MCTSampleSequence::const_iterator MCTSampleSequence::begin() const
{
	return sample.begin();
}
			
MCTSampleSequence::const_iterator MCTSampleSequence::end() const
{
	return sample.end();
}


// non member functions
std::ostream& mct::operator<< (std::ostream &out, 
								const MCTSampleSequence &seq)
{
	out << seq.toString();
	return out;
}
