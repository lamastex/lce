
/*
 * definitions for making populations for runs
*/


#include "config.h"

#include "make_run_pops.hpp"

#include <cmath>
#include <stdexcept>

using namespace std;




boost::shared_ptr< hudson_ms::PopulationStructure > makePopSingle()
{
	try {
	
		boost::shared_ptr< hudson_ms::PopulationStructure > pop ( new hudson_ms::PopulationStructure );
		
		{
			// add a sub-pop
			size_t subNsamSize = 5;
			pop->add(subNsamSize);
		}
			
		return pop;	
	}
	catch (std::exception& e) {
	
		std::cout << "Error in makePopIsland:\n" + std::string(e.what()) << std::endl;
		throw;
	}
}

boost::shared_ptr< hudson_ms::PopulationStructure > makePopIslandSmall(size_t Nzero)
{
	try {
	
		boost::shared_ptr< hudson_ms::PopulationStructure > pop ( new hudson_ms::PopulationStructure );
		
		{
			// add a sub-pop
			size_t subNsamSize = 2;
			pop->add(subNsamSize);
		}
		{
			// add a sub-pop
			size_t subNsamSize = 2;
			pop->add(subNsamSize);
		}
		
		size_t i = 0;
		size_t j = 1;
		double _m = 1.0/(4.0*1000000); // fudge so that we get 4nZero*_m = 1.0 when Nzero = 1000000
		pop->setMigration(i,j,_m*4*Nzero);
		
		std::string migCheckResult = pop->doMigCheck();
		
		if (!migCheckResult.empty()) {
			throw std::logic_error( migCheckResult);
		}
		
		return pop;	
	}
	catch (std::exception& e) {
	
		std::cout << "Error in makePopIslandSmall:\n" + std::string(e.what()) << std::endl;
		throw;
	}
}


boost::shared_ptr< hudson_ms::PopulationStructure > makePopIslandBasic(size_t Nzero)
{
	try {
	
		boost::shared_ptr< hudson_ms::PopulationStructure > pop ( new hudson_ms::PopulationStructure );
		
		{
			// add a sub-pop
			size_t subNsamSize = 3;
			pop->add(subNsamSize);
		}
		{
			// add a sub-pop
			size_t subNsamSize = 2;
			pop->add(subNsamSize);
		}
		
		size_t i = 0;
		size_t j = 1;
		double _m = 5.0/(4.0*1000000); // fudge so that we get 4nZero*_m = 5.0 when Nzero = 1000000
		pop->setMigration(i,j,_m*4*Nzero);
		
		std::string migCheckResult = pop->doMigCheck();
		
		if (!migCheckResult.empty()) {
			throw std::logic_error( migCheckResult);
		}
		
		return pop;	
	}
	catch (std::exception& e) {
	
		std::cout << "Error in makePopIslandBasic:\n" + std::string(e.what()) << std::endl;
		throw;
	}
}

boost::shared_ptr< hudson_ms::PopulationStructure > makePopIslandComplex(size_t Nzero)
{ 
	try {
		
	
		// set up an empty population structure
		boost::shared_ptr< hudson_ms::PopulationStructure > pop ( new hudson_ms::PopulationStructure );
				
		// add a sub-pop
		{
			size_t subNsamSize = 3;
			pop->add(subNsamSize);
		}
		
		{
			// add another sub-pop
			size_t subNsamSize = 5;
			double _sr = 0.5;
			double _gr = 0.8;
			pop->add(subNsamSize, _sr, _gr);
		}
			
		// set migration matrix elements
		{
			size_t i = 0;
			size_t j = 1;
			double _m = 5.0/(4.0*1000000); // fudge so that we get 4nZero*_m = 5.0 when Nzero = 1000000
			pop->setMigration(i,j,_m*4*Nzero);
		}
		
		// set migration matrix elements
		{
			size_t i = 1;
			size_t j = 0;
			double _m = 6.0/(4.0*1000000); // fudge so that we get 4nZero*_m = 6.0 when Nzero = 1000000
			pop->setMigration(i,j,_m*4*Nzero);
		}
		
		std::string migCheckResult = pop->doMigCheck();
		
		if (!migCheckResult.empty()) {
			throw std::logic_error( migCheckResult);
		}
		
		return pop;	
	}
	catch (std::exception& e) {
	
		std::cout << "Error in makePopIslandComplex:\n" + std::string(e.what()) << std::endl;
		throw;
	}
	
	
}

//alpha is mig rate into pop 0 from pop 1, beta is mig rate into pop 1 from pop 0
boost::shared_ptr< hudson_ms::PopulationStructure > makePopRaazRuns(
				size_t subNsamSizeEach,
				double alpha, double beta)
{
	try {
	
		boost::shared_ptr< hudson_ms::PopulationStructure > pop ( new hudson_ms::PopulationStructure );
		
		{
			// add a sub-pop
			size_t subNsamSize = subNsamSizeEach;
			pop->add(subNsamSize);
		}
		{
			// add a sub-pop
			size_t subNsamSize = subNsamSizeEach;
			pop->add(subNsamSize);
		}
		
		{
			size_t i = 0;
			size_t j = 1;
			double _m = alpha;
			pop->setMigration(i,j,_m);
		}
		
		{
			size_t i = 1;
			size_t j = 0;
			double _m = beta;
			pop->setMigration(i,j,_m);
		}
		
		std::string migCheckResult = pop->doMigCheck();
		
		if (!migCheckResult.empty()) {
			throw std::logic_error( migCheckResult);
		}
		
		return pop;	
	}
	catch (std::exception& e) {
	
		std::cout << "Error in makePopRaazRuns:\n" + std::string(e.what()) << std::endl;
		throw;
	}
}


//mct only
// 3 subpops and 3 samples in each
boost::shared_ptr< hudson_ms::PopulationStructure > makePopIslandBalanced(size_t Nzero)
{
	try {
	
		boost::shared_ptr< hudson_ms::PopulationStructure > pop ( new hudson_ms::PopulationStructure );
		
		{
			// add a sub-pop
			size_t subNsamSize = 3;
			pop->add(subNsamSize);
		}
		{
			// add a sub-pop
			size_t subNsamSize = 3;
			pop->add(subNsamSize);
		}
		
		{
			// add a sub-pop
			size_t subNsamSize = 3;
			pop->add(subNsamSize);
		}
		
		for (size_t i = 0; i < 3; ++i) {
			for (size_t j = 0; j < 3; ++j) {
				double _m = 1.0/(4.0*1000000); // fudge so that we get 4nZero*_m = 1.0 when Nzero = 1000000
				if ( i != j )pop->setMigration(i,j,_m*4*Nzero);
			}
		}
		
		std::string migCheckResult = pop->doMigCheck();
		
		if (!migCheckResult.empty()) {
			throw std::logic_error( migCheckResult);
		}
		
		return pop;	
	}
	catch (std::exception& e) {
	
		std::cout << "Error in makePopIslandBalanced:\n" + std::string(e.what()) << std::endl;
		throw;
	}
}

// 3 subpops and 3 samples in each but communication only between 1 and 2, 1 and 3
boost::shared_ptr< hudson_ms::PopulationStructure > makePopIslandBalLessComm(size_t Nzero)
{
	try {
	
		boost::shared_ptr< hudson_ms::PopulationStructure > pop ( new hudson_ms::PopulationStructure );
		
		{
			// add a sub-pop
			size_t subNsamSize = 3;
			pop->add(subNsamSize);
		}
		{
			// add a sub-pop
			size_t subNsamSize = 3;
			pop->add(subNsamSize);
		}
		
		{
			// add a sub-pop
			size_t subNsamSize = 3;
			pop->add(subNsamSize);
		}
		
		{
			size_t i = 0;
			size_t j = 1;
			double _m = 1.0/(4.0*1000000); // fudge so that we get 4nZero*_m = 1.0 when Nzero = 1000000
			pop->setMigration(i,j,_m*4*Nzero);

		}
		{
			size_t i = 0;
			size_t j = 2;
			double _m = 1.0/(4.0*1000000); // fudge so that we get 4nZero*_m = 1.0 when Nzero = 1000000
			pop->setMigration(i,j,_m*4*Nzero);

		}
		{
			size_t i = 1;
			size_t j = 0;
			double _m = 1.0/(4.0*1000000); // fudge so that we get 4nZero*_m = 1.0 when Nzero = 1000000
			pop->setMigration(i,j,_m*4*Nzero);
		}
		
		std::string migCheckResult = pop->doMigCheck();
		
		if (!migCheckResult.empty()) {
			throw std::logic_error( migCheckResult);
		}
		
		return pop;	
	}
	catch (std::exception& e) {
	
		std::cout << "Error in makePopIslandBalLessComm:\n" + std::string(e.what()) << std::endl;
		throw;
	}
}


// a test population for microsat popn stats
boost::shared_ptr< hudson_ms::PopulationStructure > makePopIslandMSstats(size_t Nzero)
{
	try {
	
		boost::shared_ptr< hudson_ms::PopulationStructure > pop ( new hudson_ms::PopulationStructure );
		
		{
			// add a sub-pop
			size_t subNsamSize = 6;
			pop->add(subNsamSize);
		}
		{
			// add a sub-pop
			size_t subNsamSize = 4;
			pop->add(subNsamSize);
		}
		
		{
			// add a sub-pop
			size_t subNsamSize = 5;
			pop->add(subNsamSize);
		}
		{
			size_t i = 0;
			size_t j = 1;
			double _m = 5.0/(4.0*1000000); // fudge so that we get 4nZero*_m = 5.0 when Nzero = 1000000
			pop->setMigration(i,j,_m*4*Nzero);
		}
		{
			size_t i = 1;
			size_t j = 0;
			double _m = 5.0/(4.0*1000000); // fudge so that we get 4nZero*_m = 5.0 when Nzero = 1000000
			pop->setMigration(i,j,_m*4*Nzero);
		}
		{
			size_t i = 0;
			size_t j = 2;
			double _m = 5.0/(4.0*1000000); // fudge so that we get 4nZero*_m = 5.0 when Nzero = 1000000
			pop->setMigration(i,j,_m*4*Nzero);
		}
		{
			size_t i = 2;
			size_t j = 0;
			double _m = 5.0/(4.0*1000000); // fudge so that we get 4nZero*_m = 5.0 when Nzero = 1000000
			pop->setMigration(i,j,_m*4*Nzero);
		}
		{
			size_t i = 1;
			size_t j = 2;
			double _m = 5.0/(4.0*1000000); // fudge so that we get 4nZero*_m = 5.0 when Nzero = 1000000
			pop->setMigration(i,j,_m*4*Nzero);
		}
		{
			size_t i = 2;
			size_t j = 1;
			double _m = 5.0/(4.0*1000000); // fudge so that we get 4nZero*_m = 5.0 when Nzero = 1000000
			pop->setMigration(i,j,_m*4*Nzero);
		}

		std::string migCheckResult = pop->doMigCheck();
		
		if (!migCheckResult.empty()) {
			throw std::logic_error( migCheckResult);
		}

		return pop;	

	}
	catch (std::exception& e) {
	
		std::cout << "Error in makePopIslandMSstats:\n" + std::string(e.what()) << std::endl;
		throw;
	}
}
