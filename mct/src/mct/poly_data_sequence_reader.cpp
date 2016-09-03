/* MCT - Markov Chains on Trees.

   Copyright (C) 2011, 2012 Jenny Harlow

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
\brief PolyDataSequenceReader definitions.
*/


#include "poly_data_sequence_reader.hpp"

#include "libseq.hpp"


#include <mscplusplus/population_structure.hpp>



#include <Sequence/PolySites.hpp>

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <cassert>

using namespace std;
using namespace mct;


//#define MYDEBUG

#if defined (MYDEBUG)
	#include <iostream>
	
#endif

boost::shared_ptr < mct::MultiLociPolyTable > 
	PolyDataSequenceReader::readMultiLociPolytableFromArlequinDNAfileSingleLocus(
			const std::string& s)
{
	try {
		std::vector < std::map < int, int > > idFreqMapVec;
		std::map < int, std::string > idsAndNucStrings;
			
			
		bool success = readArlequinDNAFile(s,
						idFreqMapVec,
						idsAndNucStrings) ;
		
		if (!success) throw std::runtime_error(
			"Unsuccessful read");
		
		boost::shared_ptr < MultiLociPolyTable > mlpt 
			= makeMultiLociPolytable(idFreqMapVec, idsAndNucStrings);
		
		return mlpt;
	}
	catch (std::runtime_error& re) {
		std::string msg("PolyDataSequenceReader::");
		msg += "readMultiLociPolytableFromArlequinDNAfile(const std;:string) :\n";
		msg += re.what();
		throw;
	}
}


boost::shared_ptr < mct::MultiLociPolyTable > 
			PolyDataSequenceReader::makeMultiLociPolytable(
			const std::vector < std::map < int, int > >& idFreqMapVec,
			const std::map < int, std::string >& idsAndNucStrings)
{
			
	std::vector < std::string > strs;
	// pops is number of entries in the vec
	// each pop has number of individuals calculated from map in vec
	/* for each individual, need to get a string of their nuc data*/
	
	size_t n_pops = idFreqMapVec.size();
	std::vector < size_t > individuals(n_pops, 0);
	
	for (size_t i = 0; i < n_pops; ++i) {
		for (std::map < int, int >::const_iterator it = idFreqMapVec[i].begin();
				it != idFreqMapVec[i].end(); ++it) {
			
			int reps = it->second;
			individuals[i] += reps;
			
			int id = it->first;
			std::map < int, std::string >::const_iterator nit = idsAndNucStrings.find(id);
			if (nit != idsAndNucStrings.end()) {
				string str = nit->second;
				while (reps > 0) {
					strs.push_back(str);
					--reps;
				}
			}
			else throw std::runtime_error("Cannot find id entry");
		}
		
	}
	
	
	if ( !Sequence::Alignment::validForPolyAnalysis
		(strs.begin(), strs.end()) ) {
			throw std::runtime_error(
			"Not a valid aligned pure DNA sequence");
	}

	boost::shared_ptr < Sequence::PolySites > 
			polyPtr( new Sequence::PolySites(strs) ); 
			
	size_t n_sites = strs.back().size();
	
	boost::shared_ptr< hudson_ms::PopulationStructure > pop 
			= makePop(individuals);
	
	boost::shared_ptr < mct::MultiLociPolyTable > mlptr(new mct::MultiLociPolyTable(pop));
	mlptr->add(polyPtr, n_sites);
	
	return mlptr;
	
}
			



bool PolyDataSequenceReader::readArlequinDNAFile(const std::string& s,
		std::vector < std::map < int, int > >& idFreqMapVec,
		std::map < int, std::string >& idsAndNucStrings) 
{
	bool readOkay = true;
	ifstream dataFile(s.c_str());
	
	string line;
	
	if (dataFile.is_open())
	{
		if ( !checkDataTypeDNA(dataFile) ) {
			cerr << "readArlequinDNAFile(...) : DataType not DNA" << endl;
			readOkay = false;
		}
		while (readOkay && dataFile.good()) {
			readOkay = findAndReadArlequinDNASampleData(dataFile,
							idFreqMapVec,
							idsAndNucStrings);
			
		}
		if (readOkay && (idFreqMapVec.empty() || idsAndNucStrings.empty()) ) {
						cerr << "No population sample data found" << endl;
						readOkay = false;
		}
	}
	else { // dataFile not open
	// todo should have this as our own subpavings io exception 
		std::cerr << "Error in "
			<< "_readRvectorsFromTxt."
			<< "Unable to open file " << s << std::endl;
		readOkay = false;
	}
	return readOkay;
}

bool PolyDataSequenceReader::findAndReadArlequinDNASampleData(ifstream& dataFile,
		std::vector < std::map < int, int > >& idFreqMapVec,
		std::map < int, std::string >& idsAndNucStrings) 
{
	string line;
	string findStartPhrase("SampleData");
	string findEndPhrase("}");
	bool foundStartPhrase = false;
	bool foundEndPhrase = false;
	bool readOkay = true;
	
	while ( !foundStartPhrase && (dataFile.good())) {
		getline (dataFile,line);
		
		if ( line.find(findStartPhrase) != std::string::npos ) 
										foundStartPhrase = true;
		
	}
	if (foundStartPhrase) {
		
		idFreqMapVec.push_back(std::map < int, int >());
			
		while ( !foundEndPhrase && (dataFile.good() && readOkay) ) {
			getline (dataFile,line);
			if ( line.find(findEndPhrase) != std::string::npos ) 
											foundEndPhrase = true;
			else if ( line.find(findStartPhrase) != std::string::npos ) {
							cerr << "Found another start before this end" << endl;
							readOkay = false;
			}
			else { // process the line
				
				readOkay = processArlequinDNALine(line,
										idFreqMapVec,
										idsAndNucStrings);
			}
		}
		if (readOkay && !foundEndPhrase) {
			cerr << "findAndReadArlequinDNASampleData(...) : cannot find end }" << endl;
			readOkay = false;
		}
		if (!readOkay) {
			std::vector < std::map < int, int > >().swap(idFreqMapVec);
			std::map < int, std::string >().swap(idsAndNucStrings);
		}
	}
	return readOkay;
}

bool PolyDataSequenceReader::processArlequinDNALine(std::string& line,
				std::vector < std::map < int, int > >& idFreqMapVec,
				std::map < int, std::string >& idsAndNucStrings)
{
	string nucCodesUC("ACGT");
	string nucCodesLC("acgt");
	string nucCodes = nucCodesUC + nucCodesLC;
	string whitespace(" \t");
	string allowable = nucCodes+whitespace;
	bool readOkay = true;
	
	int id;
	int freq;
	string nucs;
	istringstream sin(line);
	if (line.find_first_of(nucCodes) != std::string::npos) {
		sin >> id >> freq >> nucs;
		
	}
	else {
		// read the ID and look up if in map
		sin >> id >> freq;
	}
					
	// if ID not in map, add id and nuc to map
	if (!idsAndNucStrings.count(id)) {
		if (nucs.find_first_not_of(allowable) != std::string::npos) {
			readOkay = false;
		}
		if (readOkay) {
			size_t len = nucs.size();
			size_t pos = nucs.find_first_of(whitespace);
			while (pos != std::string::npos) {
				nucs.erase(pos, 1);
				--len;
				if (pos < len) pos = nucs.find_first_of(whitespace, pos);
				else pos = std::string::npos;
			}
		}
		if (nucs.empty()) {
			cerr << "Expected to find nucleotide string" << endl;
			readOkay = false;		
		}
		if (readOkay) {
			
			size_t posLC = nucs.find_first_of(nucCodesLC);
				
			while (posLC != std::string::npos) {
				nucs[posLC] = nucCodesUC [nucCodesLC.find(nucs[posLC])];
				posLC = nucs.find_first_of(nucCodesLC, posLC);
				
			}
			idsAndNucStrings[id] = nucs;
			
		}
		else cerr << "Error with line " << line << endl;
	}
	else if (!nucs.empty()) {
		if (idsAndNucStrings[id] != nucs) {
			readOkay = false;
			cerr << "Incompatible nucleotide string in line " << line << endl;
		}
	}
	
	// add info to thisIdFreqMap
	if (readOkay) {
		std::pair < std::map < int, int >::iterator, bool >insPair;
		insPair = idFreqMapVec.back().insert( 
			std::pair < int, int >(id, freq) );
		if (!insPair.second) insPair.first->second += freq;
	}
	return readOkay;
}




bool PolyDataSequenceReader::checkDataTypeDNA(ifstream& dataFile)
{
	string line;
	string findPhrase("DataType");
	bool foundPhrase = false;
	string findPhraseDNA("DNA");
	bool isDNA = false;
	
	while ( !foundPhrase && (dataFile.good() ) ) {
		getline (dataFile,line);
		if ( line.find(findPhrase) != std::string::npos ) 
										foundPhrase = true;
	}
	if ( foundPhrase &&
		(line.find(findPhraseDNA) != std::string::npos) )
										isDNA = true;
	
	return isDNA;
}


boost::shared_ptr< hudson_ms::PopulationStructure > 
		PolyDataSequenceReader::makePop(std::vector < size_t > nsams)
{
	boost::shared_ptr< hudson_ms::PopulationStructure > pop ( new hudson_ms::PopulationStructure );
	
	for (size_t i = 0; i < nsams.size(); ++i) {
		// add a sub-pop
		pop->add(nsams[i]);
	}
	
	return pop;	
}
