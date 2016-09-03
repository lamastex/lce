
// definitions for testing emulation of ms args code

#include "config.h"
#include "test_args.hpp"
#include "make_test_pops.hpp"

#include <mscplusplus/arg_factory.hpp>
#include <mscplusplus/arg.hpp>
#include <mscplusplus/prng_gsl_glibc2.hpp>

#include <gsl/gsl_rng.h>


#include <iostream>

using namespace std;

void testArgsSub(boost::shared_ptr< hudson_ms::PopulationStructure > pop,
		const hudson_ms::HudsonMSRecombParams& recombParams);

int main()
{
		testNewArgs1();
		testNewArgs2();
		//testMSArgPrinting();
		
		return 0;
	
}




void testNewArgs1()

{ 
	try {
		
		try {
			
			cout << "\n\nTry making one-sub-pop population" << endl;
			boost::shared_ptr< hudson_ms::PopulationStructure > pop  = makePopSingle();
			
			hudson_ms::HudsonMSRecombParams recombParams;
		
			testArgsSub(pop, recombParams);
		
		}
		catch (exception& ee) {
			std::cout << "That resulted in exception\n" 
				+ std::string(ee.what()) << std::endl;
		}
		try {
			
			cout << "\n\nTry making one-sub-pop population" << endl;
			boost::shared_ptr< hudson_ms::PopulationStructure > pop  = makePopSingle();
			
			double g = 5.0;
			//double rho = 20.0; // 4No r where r is cross proby over between sites in nsites 
			double rho = 5.0; // 4No r where r is cross proby over between sites in nsites 
			size_t nsites = 2500;
			
			pop->setAllGrowthRates( g );
		
			hudson_ms::HudsonMSRecombParams recombParams(rho, nsites);
		
			testArgsSub(pop, recombParams);
		
		}
		catch (exception& ee) {
			std::cout << "That resulted in exception\n" 
				+ std::string(ee.what()) << std::endl;
		}
		
		try {
			
			cout << "\n\nTry making one-sub-pop population with gene conversion" << endl;
			boost::shared_ptr< hudson_ms::PopulationStructure > pop  = makePopSingle();
		
			double g = 5.0;
			//double rho = 20.0; // 4No r where r is cross proby over between sites in nsites 
			double rho = 5.0; // 4No r where r is cross proby over between sites in nsites 
			size_t nsites = 2500;
			
			//rho = r*(nsites-1)*4*Nzero
			#if(0)
			int nZero = 1.0E4;			
			double r = rho/(4*nZero*(nsites-1));
			double gene_con = 3.0E-8;
			double q = 3.0E-3;
			double conversion_parameter = gene_con/r;
			double conversion_track_len = 1/q;
			#endif
			
			double conversion_parameter = 0.6;
			double conversion_track_len = 300.0;
			
			pop->setAllGrowthRates( g );
		
			hudson_ms::HudsonMSRecombParams recombParams(rho, nsites, 
							conversion_parameter, conversion_track_len);
			
			testArgsSub(pop, recombParams);
		
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
			//double rho = 20.0; // 4No r where r is cross proby over between sites in nsites 
			double rho = 5.0; // 4No r where r is cross proby over between sites in nsites 
			size_t nsites = 2500;
			
			pop->setAllGrowthRates( g );
			
			hudson_ms::HudsonMSRecombParams recombParams(rho, nsites);
			
			testArgsSub(pop, recombParams);
		}
		catch (exception& ee) {
			std::cout << "That resulted in exception\n" 
				+ std::string(ee.what()) << std::endl;
		}
		
	}
	catch (std::exception& e) {
	
		std::cout << "Error:\n" + std::string(e.what()) << std::endl;
	}
	
	
}

void testNewArgs2()

{ 
	try {
		
		try {
			
			cout << "\n\nTest -eN:  Try making instantaneous size change" << endl;
			boost::shared_ptr< hudson_ms::PopulationStructure > pop  = makePopInstantSizeChange();
		
			double rho = 0.0; // 4No r where r is cross proby over between sites in nsites 
			size_t nsites = 500;
			
			hudson_ms::HudsonMSRecombParams recombParams(rho, nsites);
			
			testArgsSub(pop, recombParams);
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
			
			testArgsSub(pop, recombParams);
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
			
			testArgsSub(pop, recombParams);
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
			
			testArgsSub(pop, recombParams);
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
			
			testArgsSub(pop, recombParams);
		}
		catch (exception& ee) {
			std::cout << "That resulted in exception\n" 
				+ std::string(ee.what()) << std::endl;
		} 
		
		
		try {
			// this gives slightly different results to classic MS because I 
			// have made all demographic events at the same time be processed together
			cout << "\n\nTry making stepping stone model with recent barrier" << endl;
			cout << "this gives slightly different results to classic MS because I" << endl; 
			cout << "have made all demographic events at the same time be processed together" << endl;
			
			boost::shared_ptr< hudson_ms::PopulationStructure > pop 
				= makePopSteppingStoneRecentBarrier1();
		
			double rho = 0.0; // 4No r where r is cross proby over between sites in nsites 
			size_t nsites = 7000;
			
			hudson_ms::HudsonMSRecombParams recombParams(rho, nsites);
			
			testArgsSub(pop, recombParams);
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
			
			testArgsSub(pop, recombParams);
		}
		catch (exception& ee) {
			std::cout << "That resulted in exception\n" 
				+ std::string(ee.what()) << std::endl;
		}
		//boost::shared_ptr< hudson_ms::PopulationStructure > makePopBackwardsSplit()
	
	}
	catch (std::exception& e) {
	
		std::cout << "Error:\n" + std::string(e.what()) << std::endl;
	}
	
}


void testMSArgPrinting()

{ 
	try {
		
		cout << "\n\nTest arg output" << endl;
		cout << "\n\nTry making one-sub-pop population with recombination" << endl;
		
		boost::shared_ptr< hudson_ms::PopulationStructure > pop  = makeTestPopForSS1();
		
		
		try {
			
			size_t seed = 205642444;

			
			boost::shared_ptr < hudson_ms::HudsonMSPRNGen > 
				rr ( new hudson_ms::HudsonMSPRNGenGLS_glibc2(seed) );
			
						
			//double g = 5.0;
			double rho = 0.2; // 4No r where r is cross proby over between sites in nsites 
			size_t nsites = 2501;
			
			//pop->setAllGrowthRates( g );
		
			hudson_ms::HudsonMSRecombParams recombParams(rho, nsites);
			
			hudson_ms::HudsonMSArgFactory factory(recombParams,
													pop,
													rr);
			
			
			boost::shared_ptr < hudson_ms::HudsonMSArg > arg_ptr = factory.makeArg();
			
			std::string cmd = pop->makeMScommand(recombParams, seed);
			std::cout << "The ms command line would be" << std::endl;
			std::cout << cmd << std::endl;
			
			cout << "\n\narg_ptr->outputArg(std::cout);\n" << endl;
			arg_ptr->outputArg(std::cout);
			cout << "\n\narg_ptr->outputArg(std::cout, 5, true);\n" << endl;
			arg_ptr->outputArg(std::cout, 5, true);
			cout << "\n\narg_ptr->outputArg(std::cout, 5, true, 3);\n" << endl;
			arg_ptr->outputArg(std::cout, 5, true, 3);
			cout << "\n\narg_ptr->outputArg(std::cout, true);\n" << endl;
			arg_ptr->outputArg(std::cout, true);
			cout << "\n\narg_ptr->outputArg(std::cout, true, 3);\n" << endl;
			arg_ptr->outputArg(std::cout, true, 3);
			cout << "\n\narg_ptr->outputArg(std::cout, 5);\n" << endl;
			arg_ptr->outputArg(std::cout, 5);
			cout << "\n\narg_ptr->outputArg(std::cout, 2, 5);\n" << endl;
			arg_ptr->outputArg(std::cout, 2, 5);
			cout << "\n\narg_ptr->outputArg(std::cout, 2, 5, 3);\n" << endl;
			arg_ptr->outputArg(std::cout, 2, 5, 3);
			
		}
		catch (exception& ee) {
			std::cout << "That resulted in exception\n" 
				+ std::string(ee.what()) << std::endl;
		}
		
	}
	catch (std::exception& e) {
	
		std::cout << "Error:\n" + std::string(e.what()) << std::endl;
	}
	
	
}

void testArgsSub(boost::shared_ptr< hudson_ms::PopulationStructure > pop,
		const hudson_ms::HudsonMSRecombParams& recombParams)
{

	cout << "Try using this to make an ARG " << endl;

	size_t seed = 205642444;
	cout << "seed is " << seed << endl;

	boost::shared_ptr < hudson_ms::HudsonMSPRNGen > 
		rr ( new hudson_ms::HudsonMSPRNGenGLS_glibc2(seed) );
	
				
	cout << recombParams.toString() << endl;
	
	std::string cmd = pop->makeMScommand(recombParams, seed);
	
	std::cout << "The revised ms command line would be" << std::endl;
	std::cout << cmd << std::endl;
		
	hudson_ms::HudsonMSArgFactory factory(recombParams,
											pop,
											rr);
	
	boost::shared_ptr < hudson_ms::HudsonMSArg > arg = factory.makeArg();
	
	arg->outputArg(std::cout);
	
}
