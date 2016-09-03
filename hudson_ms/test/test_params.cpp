
// definitions for testing recombination and segsites parameter objects

#include "config.h"
#include "test_params.hpp"
#include "make_test_pops.hpp"

#include <mscplusplus.hpp>


#include <boost/smart_ptr.hpp>


#include <gsl/gsl_rng.h>


#include <iostream>
#include <stdexcept>

using namespace std;


int main()
{
		testMSParams();
		
		return 0;
	
}


void testMSParams()

{
	//testing recombParams
	
	try {
		
		cout << "\n\nTest hudson_ms::HudsonMSRecombParams" << endl;
		
		bool ok = false;
		try {
			
			double rho = -0.1;
			int nsites = 2;
			
			std::cout << "\nrho = " << rho << " nsites = " << nsites << std::endl;
			ok = false;
			
			hudson_ms::HudsonMSRecombParams recombParams(rho, nsites);
			
			if (!ok) throw std::logic_error("Should not be able to do that");
		
		}
		catch (std::invalid_argument& ia) {
		
			std::cout << "invalid_argument:\n" + std::string(ia.what()) << std::endl;
			if (ok) throw;
		}
		
		try {
			
			double rho = 0.1;
			int nsites = 1;
			
			std::cout << "\nrho = " << rho << " nsites = " << nsites << std::endl;
			ok = false;
			
			hudson_ms::HudsonMSRecombParams recombParams(rho, nsites);
			
			if (!ok) throw std::logic_error("Should not be able to do that");
		
		}
		catch (std::invalid_argument& ia) {
		
			std::cout << "invalid_argument:\n" + std::string(ia.what()) << std::endl;
			if (ok) throw;
		}
		try {
			
			double rho = 0.0;
			int nsites = 1;
			
			std::cout << "\nrho = " << rho << " nsites = " << nsites << std::endl;
			ok = true;
			
			hudson_ms::HudsonMSRecombParams recombParams(rho, nsites);
			
			if (!ok) throw std::logic_error("Should not be able to do that");
		
		}
		catch (std::invalid_argument& ia) {
		
			std::cout << "invalid_argument:\n" + std::string(ia.what()) << std::endl;
			if (ok) throw;
		}
		try {
			
			double rho = 0.1;
			int nsites = 2;
			
			std::cout << "\nrho = " << rho << " nsites = " << nsites << std::endl;
			ok = true;
						
			hudson_ms::HudsonMSRecombParams recombParams(rho, nsites);
			
			if (!ok) throw std::logic_error("Should not be able to do that");
		
		}
		catch (std::invalid_argument& ia) {
		
			std::cout << "invalid_argument:\n" + std::string(ia.what()) << std::endl;
			if (ok) throw;
		}
		
		try {
			
			double rho = 0.0;
			int nsites = 2;
			double cp = -0.1;
			double ctl = 0.0;
			
			std::cout << "\nrho = " << rho << " nsites = " << nsites << " cp  " << cp << " ctl = " << ctl << std::endl;
			ok = false;
			
			hudson_ms::HudsonMSRecombParams recombParams(rho, nsites, cp, ctl);
			
			if (!ok) throw std::logic_error("Should not be able to do that");
		
		}
		catch (std::invalid_argument& ia) {
		
			std::cout << "invalid_argument:\n" + std::string(ia.what()) << std::endl;
			if (ok) throw;
		}
		try {
			
			double rho = 0.0;
			int nsites = 2;
			double cp = 0.0;
			double ctl = -0.1;
			
			std::cout << "\nrho = " << rho << " nsites = " << nsites << " cp  " << cp << " ctl = " << ctl << std::endl;
			ok = false;
			
			hudson_ms::HudsonMSRecombParams recombParams(rho, nsites, cp, ctl);
			
			if (!ok) throw std::logic_error("Should not be able to do that");
		
		}
		catch (std::invalid_argument& ia) {
		
			std::cout << "invalid_argument:\n" + std::string(ia.what()) << std::endl;
			if (ok) throw;
		}
		try {
			
			double rho = 0.0;
			int nsites = 2;
			double cp = 0.1;
			double ctl = 0.0;
			
			std::cout << "\nrho = " << rho << " nsites = " << nsites << " cp  " << cp << " ctl = " << ctl << std::endl;
			ok = false;
			
			hudson_ms::HudsonMSRecombParams recombParams(rho, nsites, cp, ctl);
			
			if (!ok) throw std::logic_error("Should not be able to do that");
		
		}
		catch (std::invalid_argument& ia) {
		
			std::cout << "invalid_argument:\n" + std::string(ia.what()) << std::endl;
			if (ok) throw;
		}
		try {
			
			double rho = 0.0;
			int nsites = 1;
			double cp = 0.1;
			double ctl = 0.1;
			
			std::cout << "\nrho = " << rho << " nsites = " << nsites << " cp  " << cp << " ctl = " << ctl << std::endl;
			ok = false;
			
			hudson_ms::HudsonMSRecombParams recombParams(rho, nsites, cp, ctl);
			
			if (!ok) throw std::logic_error("Should not be able to do that");
		
		}
		catch (std::invalid_argument& ia) {
		
			std::cout << "invalid_argument:\n" + std::string(ia.what()) << std::endl;
			if (ok) throw;
		}
		try {
			
			double rho = 0.0;
			int nsites = 2;
			double cp = 0.0;
			double ctl = 0.0;
			
			std::cout << "\nrho = " << rho << " nsites = " << nsites << " cp  " << cp << " ctl = " << ctl << std::endl;
			ok = true;
			
			hudson_ms::HudsonMSRecombParams recombParams(rho, nsites, cp, ctl);
			
			if (!ok) throw std::logic_error("Should not be able to do that");
		
		}
		catch (std::invalid_argument& ia) {
		
			std::cout << "invalid_argument:\n" + std::string(ia.what()) << std::endl;
			if (ok) throw;
		}
		try {
			
			double rho = 0.0;
			int nsites = 2;
			double cp = 0.1;
			double ctl = 0.1;
			
			std::cout << "\nrho = " << rho << " nsites = " << nsites << " cp  " << cp << " ctl = " << ctl << std::endl;
			ok = true;
			
			hudson_ms::HudsonMSRecombParams recombParams(rho, nsites, cp, ctl);
			
			if (!ok) throw std::logic_error("Should not be able to do that");
		
		}
		catch (std::invalid_argument& ia) {
		
			std::cout << "invalid_argument:\n" + std::string(ia.what()) << std::endl;
			if (ok) throw;
		}
		
		
		cout << "\n\nTest hudson_ms::HudsonMSSegSiteParams" << endl;
		
		try {
			
			double theta = -0.1;
			int segsitesin = 0;
			int mfreq = 1;
			
			std::cout << "\ntheta = " << theta << "segnsitesin = " << segsitesin << " mfreq " << mfreq << std::endl;
			ok = false;			
			
			hudson_ms::HudsonMSSegSiteParams ssParams(theta, segsitesin, mfreq);
			
			if (!ok) throw std::logic_error("Should not be able to do that");
			
		}
		catch (std::invalid_argument& ia) {
		
			std::cout << "invalid_argument:\n" + std::string(ia.what()) << std::endl;
			if (ok) throw;
		}
		
		try {
			
			double theta = 0.0;
			int segsitesin = -1;
			int mfreq = 1;
			
			std::cout << "\ntheta = " << theta << "segnsitesin = " << segsitesin << " mfreq " << mfreq << std::endl;
			ok = false;			
			
			hudson_ms::HudsonMSSegSiteParams ssParams(theta, segsitesin, mfreq);
			
			if (!ok) throw std::logic_error("Should not be able to do that");
			
		}
		catch (std::invalid_argument& ia) {
		
			std::cout << "invalid_argument:\n" + std::string(ia.what()) << std::endl;
			if (ok) throw;
		}
		
		try {
			
			double theta = 0.0;
			int segsitesin = 1;
			int mfreq = 0;
			
			std::cout << "\ntheta = " << theta << "segnsitesin = " << segsitesin << " mfreq " << mfreq << std::endl;
			ok = false;			
			
			hudson_ms::HudsonMSSegSiteParams ssParams(theta, segsitesin, mfreq);
			
			if (!ok) throw std::logic_error("Should not be able to do that");
			
		}
		catch (std::invalid_argument& ia) {
		
			std::cout << "invalid_argument:\n" + std::string(ia.what()) << std::endl;
			if (ok) throw;
		}
		
		cout << "\n\nTest hudson_ms::HudsonMSArgFactory" << endl;
						
		try {
			cout << "\n\nEmpty pop pointer" << endl;
			ok = false;
			
			double rho = 0.0;
			int nsites = 1;
			
			int seed = 1234;
			
			hudson_ms::HudsonMSRecombParams recombParams(rho, nsites);
			boost::shared_ptr < hudson_ms::PopulationStructure > pop;
			
			boost::shared_ptr < hudson_ms::HudsonMSPRNGen > 
				rr ( new hudson_ms::HudsonMSPRNGenGLS_glibc2(seed) );
			
			hudson_ms::HudsonMSArgFactory factory(recombParams,
											pop,
											rr);
			
			if (!ok) throw std::logic_error("Should not be able to do that");
	
		}
		catch (std::invalid_argument& ia) {
		
			std::cout << "invalid_argument:\n" + std::string(ia.what()) << std::endl;
			if (ok) throw;
		}
		
		try {
			
			cout << "\n\nDefault constructed pop pointer: npop = 0, nsam = 0" << endl;
			ok = false;
			
			double rho = 0.0;
			int nsites = 1;
			
			int seed = 1234;
			
			hudson_ms::HudsonMSRecombParams recombParams(rho, nsites);
			boost::shared_ptr < hudson_ms::PopulationStructure > pop(new hudson_ms::PopulationStructure());
			
			boost::shared_ptr < hudson_ms::HudsonMSPRNGen > 
				rr ( new hudson_ms::HudsonMSPRNGenGLS_glibc2(seed) );
			
			hudson_ms::HudsonMSArgFactory factory(recombParams,
											pop,
											rr);
			
			if (!ok) throw std::logic_error("Should not be able to do that");
	
		}
		catch (std::invalid_argument& ia) {
		
			std::cout << "invalid_argument:\n" + std::string(ia.what()) << std::endl;
			if (ok) throw;
		}
		
		try {
			
			cout << "\n\nPop pointer with npop = 1, nsam = 0" << endl;
			ok = false;
			
			double rho = 0.0;
			int nsites = 1;
			
			int seed = 1234;
			
			hudson_ms::HudsonMSRecombParams recombParams(rho, nsites);
			boost::shared_ptr < hudson_ms::PopulationStructure > pop(new hudson_ms::PopulationStructure());
			
			size_t ns = 0;
			
			pop->add(ns);
			
			boost::shared_ptr < hudson_ms::HudsonMSPRNGen > 
				rr ( new hudson_ms::HudsonMSPRNGenGLS_glibc2(seed) );
			
			hudson_ms::HudsonMSArgFactory factory(recombParams,
											pop,
											rr);
			
			if (!ok) throw std::logic_error("Should not be able to do that");
	
		}
		catch (std::invalid_argument& ia) {
		
			std::cout << "invalid_argument:\n" + std::string(ia.what()) << std::endl;
			if (ok) throw;
		}
		
		
		try {
			
			cout << "\n\nEmpty prng pointer" << endl;
			ok = false;
			
			double rho = 0.0;
			int nsites = 1;
			
			hudson_ms::HudsonMSRecombParams recombParams(rho, nsites);
			boost::shared_ptr < hudson_ms::PopulationStructure > pop = makePopSingle();
			
			boost::shared_ptr < hudson_ms::HudsonMSPRNGen > 
				rr;
						
			hudson_ms::HudsonMSArgFactory factory(recombParams,
											pop,
											rr);
			
			if (!ok) throw std::logic_error("Should not be able to do that");
	
		}
		catch (std::invalid_argument& ia) {
		
			std::cout << "invalid_argument:\n" + std::string(ia.what()) << std::endl;
			if (ok) throw;
		}
		
		cout << "\n\nTest hudson_ms::HudsonMSArg" << endl;
		
		try {
			
			cout << "\n\nnsam < 1 in arg constructor" << endl;
			ok = false;
			
			double theta = 0.0;
			int segsitesin = 0;
			int mfreq = 1;
			
			hudson_ms::HudsonMSSegSiteParams ssParams(theta, segsitesin, mfreq);
			
			int seed = 1234;
			
			boost::shared_ptr < hudson_ms::PopulationStructure > pop = makePopSingle();
			
			boost::shared_ptr < hudson_ms::HudsonMSPRNGen > 
				rr( new hudson_ms::HudsonMSPRNGenGLS_glibc2(seed) );
			
						
			hudson_ms::HudsonMSSegSitesFactory factory(rr);
			
			// cannot test nsegs in arg = 0 because arg constructor won't allow that
			boost::shared_ptr < hudson_ms::HudsonMSArg > arg_ptr(new hudson_ms::HudsonMSArg(0));
			
			if (!ok) throw std::logic_error("Should not be able to do that");
			
		}
		catch (std::invalid_argument& ia) {
		
			std::cout << "invalid_argument:\n" + std::string(ia.what()) << std::endl;
			if (ok) throw;
		}
		
		cout << "\n\nTest hudson_ms::HudsonMSSegSitesFactory" << endl;
		
		try {
			
			cout << "\n\nEmpty prng pointer" << endl;
			ok = false;
			
			boost::shared_ptr < hudson_ms::HudsonMSPRNGen > 
				rr;
			
						
			hudson_ms::HudsonMSSegSitesFactory factory(rr);
	
			if (!ok) throw std::logic_error("Should not be able to do that");
	
		}
		catch (std::invalid_argument& ia) {
		
			std::cout << "invalid_argument:\n" + std::string(ia.what()) << std::endl;
			if (ok) throw;
		}
		
		try {
			cout << "\n\nEmpty arg pointer" << endl;
			ok = false;
			
			double theta = 0.0;
			int segsitesin = 0;
			int mfreq = 1;
			
			hudson_ms::HudsonMSSegSiteParams ssParams(theta, segsitesin, mfreq);
			
			int seed = 1234;
			
			boost::shared_ptr < hudson_ms::PopulationStructure > pop = makePopSingle();
			
			boost::shared_ptr < hudson_ms::HudsonMSPRNGen > 
				rr( new hudson_ms::HudsonMSPRNGenGLS_glibc2(seed) );
			
			hudson_ms::HudsonMSSegSitesFactory factory(rr);
			
			boost::shared_ptr < hudson_ms::HudsonMSArg > arg_ptr;
			
			boost::shared_ptr < hudson_ms::HudsonMSSegSites > ss_ptr
					= factory.makeSegSites(arg_ptr, ssParams);
	
			if (!ok) throw std::logic_error("Should not be able to do that");
	
		}
		catch (std::invalid_argument& ia) {
		
			std::cout << "invalid_argument:\n" + std::string(ia.what()) << std::endl;
			if (ok) throw;
		}
		
		try {
			
			cout << "\n\nnsites < 2" << endl;
			ok = false;
			
			double theta = 0.0;
			int segsitesin = 0;
			int mfreq = 1;
			
			hudson_ms::HudsonMSSegSiteParams ssParams(theta, segsitesin, mfreq);
			
			int seed = 1234;
			
			boost::shared_ptr < hudson_ms::PopulationStructure > pop = makePopSingle();
			
			boost::shared_ptr < hudson_ms::HudsonMSPRNGen > 
				rr( new hudson_ms::HudsonMSPRNGenGLS_glibc2(seed) );
			
			hudson_ms::HudsonMSSegSitesFactory factory(rr);
			
			// cannot test nsegs in arg = 0 because arg constructor won't allow that
			boost::shared_ptr < hudson_ms::HudsonMSArg > arg_ptr(new hudson_ms::HudsonMSArg(1));
			
			
			boost::shared_ptr < hudson_ms::HudsonMSSegSites > ss_ptr
					= factory.makeSegSites(arg_ptr, ssParams);
	
			if (!ok) throw std::logic_error("Should not be able to do that");
			
		}
		catch (std::invalid_argument& ia) {
		
			std::cout << "invalid_argument:\n" + std::string(ia.what()) << std::endl;
			if (ok) throw;
		}
	
		try {
			
			cout << "\n\nmreq > nsam/2" << endl;
			ok = false;
			
			double rho = 0.0;
			int nsites = 1;
			
			hudson_ms::HudsonMSRecombParams recombParams(rho, nsites);
			double theta = 0.0;
			int segsitesin = 0;
			int mfreq = 3;
			
			hudson_ms::HudsonMSSegSiteParams ssParams(theta, segsitesin, mfreq);
			
			int seed = 1234;
			int nZero = 10000;
			boost::shared_ptr < hudson_ms::PopulationStructure > pop = makePopIslandSmall(nZero);
			
			boost::shared_ptr < hudson_ms::HudsonMSPRNGen > 
				rr( new hudson_ms::HudsonMSPRNGenGLS_glibc2(seed) );
			
						
			hudson_ms::HudsonMSArgFactory argFactory(recombParams,
											pop,
											rr);
	
			hudson_ms::HudsonMSSegSitesFactory ssFactory(rr);
			
			boost::shared_ptr < hudson_ms::HudsonMSArg > arg_ptr = argFactory.makeArg();
			
			boost::shared_ptr < hudson_ms::HudsonMSSegSites > ss_ptr
					= ssFactory.makeSegSites(arg_ptr, ssParams);
			
			if (!ok) throw std::logic_error("Should not be able to do that");
	
		}
		catch (std::invalid_argument& ia) {
		
			std::cout << "invalid_argument:\n" + std::string(ia.what()) << std::endl;
			if (ok) throw;
		}
	
		cout << "\n\nend of test \n\n" << endl;
	}
	catch (std::exception& e) {
	
		std::cout << "std::exception:\n" + std::string(e.what()) << std::endl;
	}
	catch (...) {
	
		std::cout << "non-standard exception" << std::endl;
	}
	
	
	
}

