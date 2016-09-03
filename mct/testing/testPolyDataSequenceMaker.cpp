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
\brief Small test of PolyDataSequenceMaker.

*/

#include "config.h"

#include "poly_data_sequence_maker.hpp"

#include "site_seq.hpp"


#include <iostream>
#include <vector>
#include <iterator>



using namespace std;
using namespace mct;


void polydataSequenceMakerTesting();

void testSub(const vector < site_seq >& data);
	

int main()
{
	polydataSequenceMakerTesting();
	
	
	return 0;
	
}

	

void polydataSequenceMakerTesting()
{
	try {
		int myPid = getpid();
		std::cout << "testing PolyDataSequenceMaker" << std::endl;
		std::cout << "This process id is " << myPid << std::endl;
		
		std::vector <site_seq> data;
		
		{
			std::cout << "\nTry to reduce empty vector" << std::endl;
			testSub(data);
		}
		
		
		{
			size_t nsam = 1;
			site_seq tmp;
			data = std::vector < site_seq > (nsam, tmp);
			std::cout << "\nTry to reduce vector of " << data.size() << " empty site_seq" << std::endl;
			testSub(data);
		}
		
		{
			size_t nsam = 3;
			site_seq tmp;
			data = std::vector < site_seq > (nsam, tmp);
			std::cout << "\nTry to reduce vector of " << data.size() << " empty site_seq" << std::endl;
			testSub(data);
		}
		
		{
			size_t nsam = 1;
			site_seq tmp;
			tmp.push_back(1);
			data = std::vector < site_seq > (nsam, tmp);
			
			std::cout << "\nTry to reduce vector of " << data.size() 
					<< " site_seq"
					<< " with size " << data.front().size() << std::endl;
			testSub(data);
		}
		
		{
			data.back().push_back(1);
			data.back().push_back(2);
			
			std::cout << "\nTry to reduce vector of " << data.size()
					<< " site_seq"
					<< " with size " << data.front().size() << std::endl;
			testSub(data);
		}
		
		{
			data.push_back(data.front());
			
			std::cout << "\nTry to reduce vector of " << data.size()
					<< " site_seq"
					<< " with size " << data.front().size() << std::endl;
			testSub(data);
		}
		
		{
			data.push_back(data.front());
			
			std::cout << "\nTry to reduce vector of " << data.size()
					<< " site_seq"
					<< " with size " << data.front().size() << std::endl;
			testSub(data);
		}
		{
			data.back().back() = 1 + data.front().back();
			
			std::cout << "\nTry to reduce vector of " << data.size()
					<< " site_seq"
					<< " with size " << data.front().size() << std::endl;
			testSub(data);
		}
		{
			data.push_back(data.back());
			data.back().front() = 1 + data.front().front();
			
			std::cout << "\nTry to reduce vector of " << data.size()
					<< " site_seq"
					<< " with size " << data.front().size() << std::endl;
			testSub(data);
		}
		{
			data.push_back(data.front());
			for (size_t i = 1; i < data.size()-1 ; ++i) {
				data[i] = data.front();
			}
			for (size_t i = 0; i < data.front().size(); ++i) {
				data[(i% data.size()) + 1][i] = data.front()[i]+1;
			}
			
			std::cout << "\nTry to reduce vector of " << data.size()
					<< " site_seq"
					<< " with size " << data.front().size() << std::endl;
			testSub(data);
		}
		
	}
	catch (std::exception& e) {
	
		std::cout << "std::exception::\n" + std::string(e.what()) << std::endl;
	}
	catch (...) {
	
		std::cout << "Unknown error\n" << std::endl;
	}
}

void testSub(const vector < site_seq >& data)
{
	std::cout << "Before reduction, data = " << std::endl;
	if (!data.empty()) {
		for (std::vector < site_seq >::const_iterator it = data.begin();
				it < data.end();
				++it) {
			std::ostream_iterator<int> out_it (std::cout," ");
			copy ( it->begin(), it->end(), out_it );
			std::cout << std::endl;
		}
		if (data.front().empty()) std::cout << "(inner sequences are empty)" << std::endl;
	}
	else std::cout << "(data is empty)" << std::endl;
	std::cout << std::endl;
	
	std::vector <site_seq> dataCopy = data;
	PolyDataSequenceMaker::reduceToPolyData(dataCopy);
	
	std::cout << "\nAfter reduction, data = " << std::endl;
	if (!data.empty()) {
		for (std::vector < site_seq >::const_iterator it = dataCopy.begin();
				it < dataCopy.end();
				++it) {
			std::ostream_iterator<int> out_it (std::cout," ");
			copy ( it->begin(), it->end(), out_it );
			std::cout << std::endl;
		}
		if (data.front().empty()) std::cout << "(inner sequences are empty)" << std::endl;
	}
	else std::cout << "(data is empty)" << std::endl;
}
