/* MCT - Markov Chains on Trees.

   Copyright (C) 2012 Jenny Harlow

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
\brief PolyDataSequenceReader declarations.
*/

#ifndef _INC_POLY_SITE_SEQUENCE_READER_H
#define _INC_POLY_SITE_SEQUENCE_READER_H

#include "multi_loci_polytable.hpp"

#include <boost/smart_ptr.hpp>

#include <map>
#include <vector>
#include <fstream>

namespace mct {
	
	
	/** @brief A static class for 
	 reading in sequence data from a file, such as an Arlequin 
	 project file.
	 * 
	
	Various methods to read in sequence information from a file and
	turn it into one of the mct sequence objects containing
	only polymorphic site information. 
	
	A 'polymorphic site' means any site position where there
	is a difference between two or more members of the sample
	in the element at that position. 
	
	\internal
	Only method to make polytable type inplemented at present. I'll add 
	more as I need them. 

	*/
	 
	

	class PolyDataSequenceReader {
		public:
			
			/*! \brief Read in nucleotide sequence information 
			from a file and
			turn it into a MultiLociPolyTable containing
			just one locus. 
			
			All DNA informaion in the file is read into a single locus
			in the MultiLociPolyTable.
			
			A 'polymorphic site' means any site position where there
			is a difference between two or more members of the sample
			in the element at that position. 
			
			\param s name of file to read in.
			\return A pointer to a MultiLociPolyTable containing
			just one locus, with population structure and sample data
			read in from the file named \a s.
			\pre The file named \a s is in Arlequin format for haplotypic
			DNA data.*/
			static boost::shared_ptr < mct::MultiLociPolyTable > 
				readMultiLociPolytableFromArlequinDNAfileSingleLocus(
						const std::string& s);
			
					
		private:
		
					
			static boost::shared_ptr< hudson_ms::PopulationStructure > 
					makePop(std::vector < size_t > nsams);

			static bool readArlequinDNAFile(const std::string& s,
					std::vector < std::map < int, int > >& idFreqMapVec,
					std::map < int, std::string >& idsAndNucStrings) ;
					
			static bool findAndReadArlequinDNASampleData(
					std::ifstream& dataFile,
					std::vector < std::map < int, int > >& idFreqMapVec,
					std::map < int, std::string >& idsAndNucStrings);
					
			static bool checkDataTypeDNA(std::ifstream& dataFile);

			static bool processArlequinDNALine(std::string& line,
							std::vector < std::map < int, int > >& idFreqMapVec,
							std::map < int, std::string >& idsAndNucStrings);

			static boost::shared_ptr < mct::MultiLociPolyTable > makeMultiLociPolytable(
						const std::vector < std::map < int, int > >& idFreqMapVec,
						const std::map < int, std::string >& idsAndNucStrings);

						
	};
	
	
} // end namespace mct

#endif 

