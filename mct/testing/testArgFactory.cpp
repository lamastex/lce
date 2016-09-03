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
\brief Test our mct::ArgFactory.

The arg factory produces args using (indirectly) my mscplusplus library.

*/

#include "config.h"

#include "mct_arg_factory.hpp"
#include "population_labeled_arg.hpp"

#include "prng.hpp"
#include "prng_gsl.hpp"

#include "make_test_pops.hpp"

#include <mscplusplus/population_structure.hpp>

#include <boost/shared_ptr.hpp>

#include <iostream>


using namespace std;
using namespace mct;

void argSub(boost::shared_ptr< hudson_ms::PopulationStructure > pop,
	double g, double rho, size_t nsites, int seed);

void testARGFactory();

int main()
{
	testARGFactory();
	
	
	return 0;
	
}

	
void testARGFactory()

{ 
	try {
		try {
		
			cout << "Try making one-sub-pop population" << endl;
			boost::shared_ptr< hudson_ms::PopulationStructure > pop  = makePopSingle();
			
			double g = 5.0;
			double rho = 0.0; // 4No r where r is cross proby over between sites in nsites 
			size_t nsites = 2500;
			int seed = 1245;
			
			argSub(pop, g, rho, nsites, seed);	
		}
		catch (exception& ee) {
			std::cout << "That resulted in exception\n" 
				+ std::string(ee.what()) << std::endl;
				throw;
		}

		try {	
			cout << "Try making a more complex population" << endl;
			size_t Nzero = 1000000;
			boost::shared_ptr< hudson_ms::PopulationStructure > pop  = makePopIslandComplex(Nzero);
			
			double g = 5.0;
			double rho = 0.0; // 4No r where r is cross proby over between sites in nsites 
			size_t nsites = 2500;
			int seed = 1245;
			
			argSub(pop, g, rho, nsites, seed);
		}
		catch (exception& ee) {
			std::cout << "That resulted in exception\n" 
				+ std::string(ee.what()) << std::endl;
				throw;
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

void argSub(boost::shared_ptr< hudson_ms::PopulationStructure > pop,
	double g, double rho, size_t nsites, int seed)
{
	
	cout << "Try using this to make an ARG " << endl;
	cout << "growth " << g << ", rho " << rho << ", nsites " << nsites << " and seed " << seed << endl;
	
	pop->setAllGrowthRates(g);
	// use it to make an ARG
	
	boost::shared_ptr < PRNGen > rr ( new PRNGenGSL(seed) );

	boost::shared_ptr < const PopulationLabeledARG > labeled_hist =
		mct::ARGFactory::makeHistFromMS(pop, nsites, rho, rr);
	
	
	std::string labeled_arg = labeled_hist->toString();
	cout << "The labeled tree is\n" << labeled_arg << endl; 
		
	
}
