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
\brief Test our mct::ArgConverter.

args to convert come from my mscplusplus library.


*/

#include "config.h"

#include "arg_ms_converter.hpp"
#include "libseq.hpp"

#include "prng.hpp"
#include "prng_gsl.hpp"

#include "make_test_pops.hpp"

#include <mscplusplus/population_structure.hpp>
#include <mscplusplus/recomb_params.hpp>

#include <Sequence/Coalescent/SimTypes.hpp>

#include <boost/shared_ptr.hpp>

#include <iostream>


using namespace std;
using namespace mct;

void testArgConverterSub(
	const hudson_ms::HudsonMSRecombParams& recombParams,
	boost::shared_ptr< hudson_ms::PopulationStructure > pop);

void testARGconverter1();
void testARGconverter2();

int main()
{
	testARGconverter1();
	testARGconverter2();
	
	
	return 0;
	
}

/*
 * You can compile arg_ms_converter.cpp with MYDEBUG_OUT defined to and compile library before 
 * running this test for extra output*/


void testARGconverter1() {
	try {
		
		cout << "\n\nTry making one-sub-pop population" << endl;
		boost::shared_ptr< hudson_ms::PopulationStructure > pop  = makePopSingle();
					
		double g = 5.0;
		double rho = 5.0; // 4No r where r is cross proby over between sites in nsites 
		size_t nsites = 2500;
		
		pop->setAllGrowthRates( g );
	
		hudson_ms::HudsonMSRecombParams recombParams(rho, nsites);
			
		testArgConverterSub(recombParams, pop);
	}
	catch (exception& ee) {
		std::cout << "That resulted in exception\n" 
			+ std::string(ee.what()) << std::endl;
	}

	try {
		
		cout << "\n\nTry making a more complex population" << endl;
		size_t Nzero = 10000;
		boost::shared_ptr< hudson_ms::PopulationStructure > pop  = makePopIslandComplex(Nzero);
		
		double g = 5.0;
		double rho = 5.0; // 4No r where r is cross proby over between sites in nsites 
		size_t nsites = 2500;
		
		pop->setAllGrowthRates( g );
			
		hudson_ms::HudsonMSRecombParams recombParams(rho, nsites);
			
		testArgConverterSub(recombParams, pop);

	}
	catch (exception& ee) {
		std::cout << "That resulted in exception\n" 
			+ std::string(ee.what()) << std::endl;
	}	
	
}

void testARGconverter2()

{ 
	try {
		
		try {
			
			cout << "\n\nTest -eN:  Try making instantaneous size change" << endl;
			boost::shared_ptr< hudson_ms::PopulationStructure > pop  = makePopInstantSizeChange();
		
			double rho = 0.0; // 4No r where r is cross proby over between sites in nsites 
			size_t nsites = 500;
			
			hudson_ms::HudsonMSRecombParams recombParams(rho, nsites);
			
			testArgConverterSub(recombParams, pop);
			
		}
		catch (exception& ee) {
			std::cout << "That resulted in exception\n" 
				+ std::string(ee.what()) << std::endl;
		}
		
		
		try {
			
			cout << "\n\nTry making outgroup sequence" << endl;
			boost::shared_ptr< hudson_ms::PopulationStructure > pop  = makePopOutgroupSequence();
		
			double rho = 0.0; // 4No r where r is cross proby over between sites in nsites 
			size_t nsites = 500;
			
			hudson_ms::HudsonMSRecombParams recombParams(rho, nsites);
			
			testArgConverterSub(recombParams, pop);

		}
		catch (exception& ee) {
			std::cout << "That resulted in exception\n" 
				+ std::string(ee.what()) << std::endl;
		} 

		try {
			
			cout << "\n\nTry making outgroup sequence with recombination" << endl;
			boost::shared_ptr< hudson_ms::PopulationStructure > pop  = makePopOutgroupSequence();
		
			int nZero = 1.0E4;			
			size_t nsites = 500;//double g = 5.0;
			double rho_per_base_pair = 1.0E-8;
			double rho =  (nsites-1)*rho_per_base_pair*4*nZero; // probability of cross over in locus per generation
			
			hudson_ms::HudsonMSRecombParams recombParams(rho, nsites);
			
			testArgConverterSub(recombParams, pop);
			
		}
		catch (exception& ee) {
			std::cout << "That resulted in exception\n" 
				+ std::string(ee.what()) << std::endl;
		} 
		
		
		try {
			
			cout << "\n\nTry making instantaneous size change followed by expontial growth" << endl;
			boost::shared_ptr< hudson_ms::PopulationStructure > pop  = makePopInstantSizeChangeAndExpGrowth();
		
			double rho = 0.0; // 4No r where r is cross proby over between sites in nsites 
			size_t nsites = 500;
			
			hudson_ms::HudsonMSRecombParams recombParams(rho, nsites);
			
			testArgConverterSub(recombParams, pop);
			
		}
		catch (exception& ee) {
			std::cout << "That resulted in exception\n" 
				+ std::string(ee.what()) << std::endl;
		}
		
		
		
		try {
			
			cout << "\n\nTry making two species with population size differences" << endl;
			boost::shared_ptr< hudson_ms::PopulationStructure > pop 
				= makePopTwoSpeciesWithPopSizeDiff();
		
			double rho = 0.0; // 4No r where r is cross proby over between sites in nsites 
			size_t nsites = 7000;
			
			hudson_ms::HudsonMSRecombParams recombParams(rho, nsites);
			
			testArgConverterSub(recombParams, pop);
			
		}
		catch (exception& ee) {
			std::cout << "That resulted in exception\n" 
				+ std::string(ee.what()) << std::endl;
		} 
		
		
		try {
			// this gives slightly different results to classic MS because I 
			// have made all demographic events at the same time be processed together
			cout << "\n\nTry making stepping stone model with recent barrier:" << endl;
			cout << "this gives slightly different results to classic MS because I" << endl; 
			cout << "have made all demographic events at the same time be processed together" << endl;
			
			boost::shared_ptr< hudson_ms::PopulationStructure > pop 
				= makePopSteppingStoneRecentBarrier1();
		
			double rho = 0.0; // 4No r where r is cross proby over between sites in nsites 
			size_t nsites = 7000;
			
			hudson_ms::HudsonMSRecombParams recombParams(rho, nsites);
			
			testArgConverterSub(recombParams, pop);
			
		}
		catch (exception& ee) {
			std::cout << "That resulted in exception\n" 
				+ std::string(ee.what()) << std::endl;
		}
		
		
		try {
			
			cout << "\n\nTry making population with backwards split (forward admixture)" << endl;
			boost::shared_ptr< hudson_ms::PopulationStructure > pop 
				= makePopBackwardsSplit();
		
			double rho = 0.0; // 4No r where r is cross proby over between sites in nsites 
			size_t nsites = 7000;
			
			hudson_ms::HudsonMSRecombParams recombParams(rho, nsites);
			
			testArgConverterSub(recombParams, pop);
		
		}
		catch (exception& ee) {
			std::cout << "That resulted in exception\n" 
				+ std::string(ee.what()) << std::endl;
		}
		//boost::shared_ptr< hudson_ms::PopulationStructure > makePopBackwardsSplit()
		
		std::cout << "\n\nEnd of testing\n" << std::endl;
	
	}
	catch (std::exception& e) {
	
		std::cout << "std::exception::\n" + std::string(e.what()) << std::endl;
	}
	catch (...) {
	
		std::cout << "Unknown error\n" << std::endl;
	}
	
}


void testArgConverterSub(
	const hudson_ms::HudsonMSRecombParams& recombParams,
	boost::shared_ptr< hudson_ms::PopulationStructure > pop)
{
	
	cout << "Try using this to make an ARG " << endl;
	
	size_t initseed = 1234;

	int seed = 0;
	{	boost::shared_ptr < PRNGen > rr ( new PRNGenGSL(initseed) );

		size_t max_seed = 999999998;
		seed = 1 + rr->uniform_int( max_seed );
		cout << "seed is " << seed << endl;
	}
	
	boost::shared_ptr < PRNGen > rr ( new PRNGenGSL(initseed) );

	cout << recombParams.toString() << endl;
	cout << "seed " << seed << endl;
	
	std::string cmd = pop->makeMScommand(recombParams, seed);
	
	std::cout << "The ms command line would be" << std::endl;
	std::cout << cmd << std::endl;
	
	cout << "\nabout to use arg converter" << endl;
			
	boost::shared_ptr < Sequence::arg >  arg_ptr = ARGmsConverter::makeArg(
							recombParams.getNsites(),
							recombParams.getRho(),
							pop,
							seed);

	cout << "\nthe arg in newick format is" << endl;
	cout << LibSeq::argNewickString(*arg_ptr) << endl;
}
