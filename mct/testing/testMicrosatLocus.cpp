/* MCT - Markov Chains on Trees.


   Copyright (C) 2011, 2012 Jennifer Harow

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
\brief Testing of MicrosatLocus class.

*/

#include "config.h"

#include "mct_microsat_locus.hpp"
#include "site_seq.hpp"

#include <iostream>
#include <iterator>


using namespace std;
using namespace mct;


void testMCTMicrosatLocus();


int main()
{
	testMCTMicrosatLocus();
	
	
	return 0;
	
}

	
void testMCTMicrosatLocus()
{
	try {
		
		// make an empty table
		mct::MCTMicrosatLocus mtab;
		cout << "mtab has size " << mtab.size() << " and nsites " << mtab.nsites() << std::endl;
		cout << "mtab empty =  " << mtab.empty() << std::endl;
				
		try {
			cout << "toString is" << endl;
			cout << mtab.toString() << endl;
		}
		catch (exception& ee) {
			std::cout << "That resulted in exception\n" 
				+ std::string(ee.what()) << std::endl;
		}
		
		std::cout << "Try adding a site_seq: adding" << endl; 
		try {
			int myints[] = {30}; 
			mct::site_seq seq (myints, myints + sizeof(myints) / sizeof(int) );
			copy(seq.begin(), seq.end(), ostream_iterator<int>(std::cout, "\t"));
			std::cout << std::endl;
			
			mtab.add(seq);
			cout << "mtab has size " << mtab.size() << " and nsites " << mtab.nsites() << std::endl;
			cout << "mtab empty =  " << mtab.empty() << std::endl;
			cout << "toString is" << endl;
			cout << mtab.toString() << endl;
		
		}
		catch (exception& ee) {
			std::cout << "That resulted in exception\n" 
				+ std::string(ee.what()) << std::endl;
		}		
		
		std::cout << "Try adding a site_seq: adding" << endl; 
		try {
			int myints[] = {30, 25}; 
			mct::site_seq seq (myints, myints + sizeof(myints) / sizeof(int) );
			copy(seq.begin(), seq.end(), ostream_iterator<int>(std::cout, "\t"));
			std::cout << std::endl;
			
			mtab.add(seq);
			cout << "mtab has size " << mtab.size() << " and nsites " << mtab.nsites() << std::endl;
			cout << "toString is" << endl;
			cout << mtab.toString() << endl;
		
		}
		catch (exception& ee) {
			std::cout << "That resulted in exception\n" 
				+ std::string(ee.what()) << std::endl;
		}		
		
		std::cout << "Try adding a site_seq: adding" << endl; 
		try {
			int myints[] = {25}; 
			mct::site_seq seq (myints, myints + sizeof(myints) / sizeof(int) );
			copy(seq.begin(), seq.end(), ostream_iterator<int>(std::cout, "\t"));
			std::cout << std::endl;
			
			mtab.add(seq);
			cout << "mtab has size " << mtab.size() << " and nsites " << mtab.nsites() << std::endl;
			cout << "toString is" << endl;
			cout << mtab.toString() << endl;
		
		}
		catch (exception& ee) {
			std::cout << "That resulted in exception\n" 
				+ std::string(ee.what()) << std::endl;
		}		
			
		
		std::cout << "Try getting element at " << endl; 
		try {
			size_t i = 0;
			cout << i << endl;
			mct::site_seq seq = mtab.at(i);
			
				cout <<  "This element is (\t";
				copy(seq.begin(), seq.end(), ostream_iterator<int>(std::cout, "\t"));
			std::cout << ")\n" << std::endl;
		
		}
		catch (exception& ee) {
			std::cout << "That resulted in exception\n" 
				+ std::string(ee.what()) << std::endl;
		}	
		std::cout << "Try getting element at " << endl; 
		try {
			size_t i = 1;
			cout << i << endl;
			mct::site_seq seq = mtab.at(i);
			
				cout << "This element is (\t";
				copy(seq.begin(), seq.end(), ostream_iterator<int>(std::cout, "\t"));
			std::cout << ")\n" << std::endl;
		
		}
		catch (exception& ee) {
			std::cout << "That resulted in exception\n" 
				+ std::string(ee.what()) << std::endl;
		}	
		std::cout << "Try getting element at " << endl; 
		try {
			size_t i = 2;
			cout << i << endl;
			mct::site_seq seq = mtab.at(i);
			
				cout << "This element is (\t";
				copy(seq.begin(), seq.end(), ostream_iterator<int>(std::cout, "\t"));
			std::cout << ")\n" << std::endl;
		
		}
		catch (exception& ee) {
			std::cout << "That resulted in exception\n" 
				+ std::string(ee.what()) << std::endl;
		}	
		
		std::cout << "Try iterating through" << endl; 
		try {
			
			for (mct::MCTMicrosatLocus::const_iterator it = mtab.begin();
			it < mtab.end();
			++it) {
				cout << "(\t";
				copy(it->begin(), it->end(), ostream_iterator<int>(std::cout, "\t"));
			std::cout << ")\n" << std::endl;
		
			}
		}
		catch (exception& ee) {
			std::cout << "That resulted in exception\n" 
				+ std::string(ee.what()) << std::endl;
		}
		std::cout << "\n\nEnd of testing\n" << std::endl;
	
	}
	catch (std::exception& e) {
	
		std::cout << "std::exception::\n" + std::string(e.what()) << std::endl;
	}
	catch (...) {
	
		std::cout << "Unknown error\n" << std::endl;
	}
	
	
}
