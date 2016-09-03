
// definitions for testing emulation of ms segsites code

#include "config.h"
#include "test_segsites.hpp"
#include "make_test_pops.hpp"

#include <mscplusplus/segsites_factory.hpp>
#include <mscplusplus/segsites.hpp>
#include <mscplusplus/arg_factory.hpp>
#include <mscplusplus/arg.hpp>
#include <mscplusplus/prng_gsl_glibc2.hpp>


#include <gsl/gsl_rng.h>



#include <iostream>

using namespace std;



int main()
{
		testMSSegsites1();
		testMSSegsites2();
		
		return 0;
	
}


void testSegSitesSub(boost::shared_ptr< hudson_ms::PopulationStructure > pop,
	const hudson_ms::HudsonMSRecombParams& recombParams,
	const hudson_ms::HudsonMSSegSiteParams& ssParams)
{
	size_t seed = 205642444;
	cout << "seed is " << seed << endl;

	
	boost::shared_ptr < hudson_ms::HudsonMSPRNGen > 
		rr ( new hudson_ms::HudsonMSPRNGenGLS_glibc2(seed) );
	
				
	cout << recombParams.toString() << endl;
	
				
	hudson_ms::HudsonMSArgFactory factory(recombParams,
											pop,
											rr);
	
	
	boost::shared_ptr < hudson_ms::HudsonMSArg > arg_ptr = factory.makeArg();
	
	hudson_ms::HudsonMSSegSitesFactory ssfactory(rr);
	
	cout << ssParams.toString() << endl;
	
	bool showTime = true;
	
	std::string cmd = pop->makeMScommand(recombParams, ssParams, seed, showTime);
	std::cout << "The ms command line would be" << std::endl;
	std::cout << cmd << std::endl;
	
	boost::shared_ptr < hudson_ms::HudsonMSSegSites > ss_ptr 
			= ssfactory.makeSegSites(arg_ptr, ssParams);
	
	int mfreq = ssParams.getMfreq();
		
	if (mfreq < 2) {
		arg_ptr->outputArg(std::cout, true);
	}
	else {
		arg_ptr->outputArg(std::cout, mfreq, 3); // default precisions
	}
	
	
	ss_ptr->outputSegSites(std::cout);
	
}


void testMSSegsitesFramework(boost::shared_ptr< hudson_ms::PopulationStructure > pop)

{ 
	try {
		
		// theta only
		try {
			
			cout << "\n\nTest theta only\n" << endl;
			
			//double g = 5.0;
			double rho = 0.2; // 4No r where r is cross proby over between sites in nsites 
			size_t nsites = 2501;
			
			//pop->setAllGrowthRates( g );
		
			hudson_ms::HudsonMSRecombParams recombParams(rho, nsites);
			
			double theta = 2.0; //theta per locus, ie 4N0 * mu_per_site * nsites 
								// where N0 = effective pop size and nsites = sites in locus
			
			hudson_ms::HudsonMSSegSiteParams ssParams(theta);
			
			testSegSitesSub(pop, recombParams, ssParams);
		}
		catch (exception& ee) {
			std::cout << "That resulted in exception\n" 
				+ std::string(ee.what()) << std::endl;
		}
		
		// theta only, no recombination
		try {
			
			cout << "\n\nTest theta only, no recombination\n" << endl;
			
			//double g = 5.0;
			//double rho = 0.2; // 4No r where r is cross proby over between sites in nsites 
			//size_t nsites = 2501;
			
			//pop->setAllGrowthRates( g );
		
			hudson_ms::HudsonMSRecombParams recombParams;
			
			double theta = 2.0;
			
			hudson_ms::HudsonMSSegSiteParams ssParams(theta);
			
			testSegSitesSub(pop, recombParams, ssParams);
		}
		catch (exception& ee) {
			std::cout << "That resulted in exception\n" 
				+ std::string(ee.what()) << std::endl;
		}
		
		// theta and segsites in
		try {
			
			cout << "\n\nTest theta and segsitesin only\n" << endl;
			
			//double g = 5.0;
			double rho = 0.2; // 4No r where r is cross proby over between sites in nsites 
			size_t nsites = 2501;
			
			//pop->setAllGrowthRates( g );
		
			hudson_ms::HudsonMSRecombParams recombParams(rho, nsites);
			double theta = 2.0;
			int segsitesin = 5;
			int mfreq = 1;
			
			hudson_ms::HudsonMSSegSiteParams ssParams(theta, segsitesin, mfreq);
			
			testSegSitesSub(pop, recombParams, ssParams);
		}
		catch (exception& ee) {
			std::cout << "That resulted in exception\n" 
				+ std::string(ee.what()) << std::endl;
		}
		//segsites in only
		try {
			
			cout << "\n\nTest segsitesin only\n" << endl;
			
			//double g = 5.0;
			double rho = 0.2; // 4No r where r is cross proby over between sites in nsites 
			size_t nsites = 2501;
			
			//pop->setAllGrowthRates( g );
		
			hudson_ms::HudsonMSRecombParams recombParams(rho, nsites);
			double theta = 0.0;
			int segsitesin = 5;
			int mfreq = 1;
			
			hudson_ms::HudsonMSSegSiteParams ssParams(theta, segsitesin, mfreq);
			
			testSegSitesSub(pop, recombParams, ssParams);
		}
		catch (exception& ee) {
			std::cout << "That resulted in exception\n" 
				+ std::string(ee.what()) << std::endl;
		}
		
		// theta and mfreq only
		try {
			
			cout << "\n\nTest theta and mfreq only\n" << endl;
			
			//double g = 5.0;
			double rho = 0.2; // 4No r where r is cross proby over between sites in nsites 
			size_t nsites = 2501;
			
			//pop->setAllGrowthRates( g );
		
			hudson_ms::HudsonMSRecombParams recombParams(rho, nsites);
			
			double theta = 2.0;
			int mfreq = 3;
			
			hudson_ms::HudsonMSSegSiteParams ssParams(theta, mfreq);
			
			testSegSitesSub(pop, recombParams, ssParams);
		}
		catch (exception& ee) {
			std::cout << "That resulted in exception\n" 
				+ std::string(ee.what()) << std::endl;
		}
		
		
		// theta and mfreq only, with no recombination but with growth
		try {
			
			cout << "\n\nTest theta and mfreq only, no recombination but with growth\n" << endl;
			
			double g = 5.0;
			
			pop->setAllGrowthRates( g );
		
			hudson_ms::HudsonMSRecombParams recombParams;
			
			double theta = 2.0;
			int mfreq = 3;
			
			hudson_ms::HudsonMSSegSiteParams ssParams(theta, mfreq);
			
			testSegSitesSub(pop, recombParams, ssParams);
		}
		catch (exception& ee) {
			std::cout << "That resulted in exception\n" 
				+ std::string(ee.what()) << std::endl;
		}
		
		// theta and mfreq only, with gene cross-over
		try {
			
			cout << "\n\nTest theta and mfreq only, with gene conversion\n" << endl;
			
			//double g = 5.0;
			double rho = 0.2; // 4No r where r is cross proby over between sites in nsites 
			size_t nsites = 2501;
			
			/* nZero = 10000;
			 * r = rho/(4*nZero*(nsites-1))
			 * gene_con = 3.0E-8
			 * q = 3.0E-3
			 * double conversion_parameter = gene_con/r
			 * double conversion_track_len = 1/q */
			double conversion_parameter = 15;
			double conversion_track_len = 300.0;
			
			
			hudson_ms::HudsonMSRecombParams recombParams(rho, nsites, 
							conversion_parameter, conversion_track_len);
			
		
			double theta = 2.0;
			int mfreq = 3;
			
			hudson_ms::HudsonMSSegSiteParams ssParams(theta, mfreq);
			
			testSegSitesSub(pop, recombParams, ssParams);
		}
		catch (exception& ee) {
			std::cout << "That resulted in exception\n" 
				+ std::string(ee.what()) << std::endl;
		}
		
		// segsitesin and mfreq only
		try {
			
			cout << "\n\nTest segsitesin and mfreq only\n" << endl;
			
			//double g = 5.0;
			double rho = 0.2; // 4No r where r is cross proby over between sites in nsites 
			size_t nsites = 2501;
			
			//pop->setAllGrowthRates( g );
		
			hudson_ms::HudsonMSRecombParams recombParams(rho, nsites);
			
			int segsitesin = 5;
			int mfreq = 2;
			
			hudson_ms::HudsonMSSegSiteParams ssParams(segsitesin, mfreq);
			
			testSegSitesSub(pop, recombParams, ssParams);
		}
		catch (exception& ee) {
			std::cout << "That resulted in exception\n" 
				+ std::string(ee.what()) << std::endl;
		}
		
		// theta and segsitesin and mfreq
		try {
			
			cout << "\n\nTest theta and segsitesin and mfreq\n" << endl;
			
			//double g = 5.0;
			double rho = 0.2; // 4No r where r is cross proby over between sites in nsites 
			size_t nsites = 2501;
			
			//pop->setAllGrowthRates( g );
		
			hudson_ms::HudsonMSRecombParams recombParams(rho, nsites);
			
			double theta = 2.0;
			int segsitesin = 5;
			int mfreq = 2;
			
			hudson_ms::HudsonMSSegSiteParams ssParams(theta, segsitesin, mfreq);
			
			testSegSitesSub(pop, recombParams, ssParams);
		}
		catch (exception& ee) {
			std::cout << "That resulted in exception\n" 
				+ std::string(ee.what()) << std::endl;
		}
		
		// theta and segsitesin and mfreq, with growth and recombination
		try {
			
			cout << "\n\nTest theta and segsitesin and mfreq, with growth but no recombination\n" << endl;
			
			double g = 5.0;
			//double rho = 0.2; // 4No r where r is cross proby over between sites in nsites 
			//size_t nsites = 2501;
			
			pop->setAllGrowthRates( g );
		
			hudson_ms::HudsonMSRecombParams recombParams;
			
			double theta = 2.0;
			int segsitesin = 5;
			int mfreq = 2;
			
			hudson_ms::HudsonMSSegSiteParams ssParams(theta, segsitesin, mfreq);
			
			testSegSitesSub(pop, recombParams, ssParams);
		}
		catch (exception& ee) {
			std::cout << "That resulted in exception\n" 
				+ std::string(ee.what()) << std::endl;
		}
		
		// theta and segsitesin and mfreq, with gene conversion as well
		try {
			
			cout << "\n\nTest theta and segsitesin and mfreq,  with gene conversion as well\n" << endl;
			
			double rho = 0.2; // 4No r where r is cross proby over between sites in nsites 
			size_t nsites = 2501;
			
			/* nZero = 10000
			 * r = rho/(4*nZero*(nsites-1))
			 * gene_con = 3.0E-8
			 * q = 3.0E-3
			 * conversion_parameter = gene_con/r
			 * conversion_track_len = 1/q */
			double conversion_parameter = 15;
			double conversion_track_len = 300.0;
			
			
			hudson_ms::HudsonMSRecombParams recombParams(rho, nsites, 
							conversion_parameter, conversion_track_len);
			
			double theta = 2.0;
			int segsitesin = 5;
			int mfreq = 2;
			
			hudson_ms::HudsonMSSegSiteParams ssParams(theta, segsitesin, mfreq);
			
			testSegSitesSub(pop, recombParams, ssParams);
		}
		catch (exception& ee) {
			std::cout << "That resulted in exception\n" 
				+ std::string(ee.what()) << std::endl;
		}
			
		
	}
	catch (...) {
	
		std::cout << "Non-standard exception error:\n" << std::endl;
		throw;
	}
	
	
}


void testMSSegsites1()

{ 
	try {
		cout << "\n\n** Try making one-sub-pop population with recombination **" << endl;
		
		boost::shared_ptr< hudson_ms::PopulationStructure > pop  = makeTestPopForSS1();
		
		testMSSegsitesFramework(pop);
	}
	catch (exception& ee) {
	
		std::cout << "std::exception error:\n" << ee.what() << std::endl;
	}
	catch (...) {
	
		std::cout << "Non-standard exception error:\n" << std::endl;
	}
}

void testMSSegsites2()

{ 
	try {
		cout << "\n\nTry making outgroup sequence" << endl;
		
		boost::shared_ptr< hudson_ms::PopulationStructure > pop  = makePopOutgroupSequence();
		
		testMSSegsitesFramework(pop);
	}
	catch (exception& ee) {
	
		std::cout << "std::exception error:\n" << ee.what() << std::endl;
	}
	catch (...) {
	
		std::cout << "Non-standard exception error:\n" << std::endl;
	}
}


