
/*
 * definitions for making test populations for testing
*/


#include "config.h"

#include "make_test_pops.hpp"

#include <cmath>
#include <stdexcept>

using namespace std;



boost::shared_ptr< hudson_ms::PopulationStructure > makeTestPop1()
{ 
	try {
		
		// set up an empty population structure
		cout << "Setting up a population structure using PopulationStructure(_sn)" << endl;
		size_t _sn = 3;
		cout << "_sn =" << _sn << endl;
		boost::shared_ptr< hudson_ms::PopulationStructure > pop ( new hudson_ms::PopulationStructure(_sn) );
		cout << "Pop number of sub-pops is " << pop->nSubPops() << endl;
		
		try {
			size_t i = 0;
			cout << "trying to access properties of subpop at " << i << endl;
			size_t sns = pop->getSubNsamAt(i);
			cout << "subpop nsam is " << sns << endl;
			double ssr = pop->getSizeRatioAt(i);
			cout << "subpop size ratio is " << ssr << endl;
			double sgr = pop->getGrowthRateAt(i);
			cout << "subpop growth rate is " << sgr << endl;
			
		}
		catch (invalid_argument& ia) {
			std::cout << "That resulted in an invalid argument exception\n" 
				+ std::string(ia.what()) << std::endl;
		}
		
		
		// print the migration matrix
		try {
			// print the migration matrix
		cout << "Pop migration matrix is " << endl;
		cout << pop->migrationMatrixToString() << endl;
		
			
		}
		catch (invalid_argument& ia) {
			std::cout << "That resulted in an invalid argument exception\n" 
				+ std::string(ia.what()) << std::endl;
		}
		
		
		return pop;
	}		
	catch (std::exception& e) {
	
		std::cout << "std::exception in makeTestPop1:\n" + std::string(e.what()) << std::endl;
		throw;
	}
	catch (...) {
	
		std::cout << "Unknown exception type in makeTestPop1:\n" << std::endl;
		throw;
	}
	
	
}		

boost::shared_ptr< hudson_ms::PopulationStructure > makeTestPop2()
{ 
	try {
		
	
		// set up an empty population structure
		cout << "\n\nSetting up empty population structure" << endl;
		boost::shared_ptr< hudson_ms::PopulationStructure > pop ( new hudson_ms::PopulationStructure );
		cout << "Pop number of sub-pops is " << pop->nSubPops() << endl;
		
		
		// print the migration matrix
		cout << "Pop migration matrix is " << endl;
		cout << pop->migrationMatrixToString() << endl;
		
		try {
			pop->getSubNsamAt(1);
			
		}
		catch (invalid_argument& ia) {
			std::cout << "That resulted in an invalid argument exception\n" 
				+ std::string(ia.what()) << std::endl;
		}
		
		// add a sub-pop
		size_t subNsamSize = 3;
		std::cout << "Adding subpopulation with " 
				<< subNsamSize << " samples" << std::endl;
		pop->add(subNsamSize);
		
		cout << "Pop number of sub-pops is " << pop->nSubPops() << endl;
		cout << "Pop number total samples is " << pop->totalNsam() << endl;
		
		try {
			size_t i = 0;
			cout << "trying to access properties of subpop at " << i << endl;
			size_t sns = pop->getSubNsamAt(i);
			cout << "subpop nsam is " << sns << endl;
			double ssr = pop->getSizeRatioAt(i);
			cout << "subpop size ratio is " << ssr << endl;
			double sgr = pop->getGrowthRateAt(i);
			cout << "subpop growth rate is " << sgr << endl;
			
		}
		catch (invalid_argument& ia) {
			std::cout << "That resulted in an invalid argument exception\n" 
				+ std::string(ia.what()) << std::endl;
		}
		
		// print the migration matrix
		cout << "Pop migration matrix is " << endl;
		cout << pop->migrationMatrixToString() << endl;
		
		
		// add another sub-pop
		subNsamSize = 5;
		double _sr = 0.5;
		double _gr = 0.8;
		std::cout << "Adding subpopulation with " 
				<< subNsamSize << " samples, "  
				<< _sr << " size ratio, " 
				<< _gr << " growth rate ratio" << std::endl;
		pop->add(subNsamSize, _sr, _gr);
		
		cout << "Pop number of sub-pops is " << pop->nSubPops() << endl;
		cout << "Pop number total samples is " << pop->totalNsam() << endl;
		
		try {
			size_t i = 1;
			cout << "trying to access properties of subpop at " << i << endl;
			size_t sns = pop->getSubNsamAt(i);
			cout << "subpop nsam is " << sns << endl;
			double ssr = pop->getSizeRatioAt(i);
			cout << "subpop size ratio is " << ssr << endl;
			double sgr = pop->getGrowthRateAt(i);
			cout << "subpop growth rate is " << sgr << endl;
			
		}
		catch (invalid_argument& ia) {
			std::cout << "That resulted in an invalid argument exception\n" 
				+ std::string(ia.what()) << std::endl;
		}
		try {
			size_t i = 2;
			cout << "trying to access properties of subpop at " << i << endl;
			size_t sns = pop->getSubNsamAt(i);
			cout << "subpop nsam is " << sns << endl;
			double ssr = pop->getSizeRatioAt(i);
			cout << "subpop size ratio is " << ssr << endl;
			double sgr = pop->getGrowthRateAt(i);
			cout << "subpop growth rate is " << sgr << endl;
			
		}
		catch (invalid_argument& ia) {
			std::cout << "That resulted in an invalid argument exception\n" 
				+ std::string(ia.what()) << std::endl;
		}
		
		// print the migration matrix
		cout << "Pop migration matrix is " << endl;
		cout << pop->migrationMatrixToString() << endl;
		
		// set migration matrix elements
		try {
			size_t i = 0;
			size_t j = 1;
			double _m = 5.0;
			cout << "trying to set migration rate from " 
					<< i << " to " << j << " to " << _m << endl;
			pop->setMigration(i,j,_m);
			double _nm = pop->getMigration(i,j);
			
			cout << "migration rate from " 
					<< i << " to " << j << " is now " << _nm << endl;
			
		}
		catch (invalid_argument& ia) {
			std::cout << "That resulted in an invalid argument exception\n" 
				+ std::string(ia.what()) << std::endl;
		}
		
		// set migration matrix elements
		try {
			size_t i = 1;
			size_t j = 0;
			double _m = 6.0;
			cout << "trying to set migration rate from " 
					<< i << " to " << j << " to " << _m << endl;
			pop->setMigration(i,j,_m);
			double _nm = pop->getMigration(i,j);
			
			cout << "migration rate from " 
					<< i << " to " << j << " is now " << _nm << endl;
			
		}
		catch (invalid_argument& ia) {
			std::cout << "That resulted in an invalid argument exception\n" 
				+ std::string(ia.what()) << std::endl;
		}
		
		// set migration matrix elements
		try {
			size_t i = 2;
			size_t j = 0;
			double _m = 6.0;
			cout << "trying to set migration rate from " 
					<< i << " to " << j << " to " << _m << endl;
			pop->setMigration(i,j,_m);
			double _nm = pop->getMigration(i,j);
			
			cout << "migration rate from " 
					<< i << " to " << j << " is now " << _nm << endl;
			
		}
		catch (invalid_argument& ia) {
			std::cout << "That resulted in an invalid argument exception\n" 
				+ std::string(ia.what()) << std::endl;
		}
		// set migration matrix elements
		try {
			size_t i = 0;
			size_t j = 2;
			double _m = 6.0;
			cout << "trying to set migration rate from " 
					<< i << " to " << j << " to " << _m << endl;
			pop->setMigration(i,j,_m);
			double _nm = pop->getMigration(i,j);
			
			cout << "migration rate from " 
					<< i << " to " << j << " is now " << _nm << endl;
			
		}
		catch (invalid_argument& ia) {
			std::cout << "That resulted in an invalid argument exception\n" 
				+ std::string(ia.what()) << std::endl;
		}
		
		// print the migration matrix
		
		cout << "Pop migration matrix is " << endl;
		cout << pop->migrationMatrixToString() << endl;
		
		return pop;	
	}		
	catch (std::exception& e) {
	
		std::cout << "std::exception in makeTestPop1:\n" + std::string(e.what()) << std::endl;
		throw;
	}
	catch (...) {
	
		std::cout << "Unknown exception type in makeTestPop1:\n" << std::endl;
		throw;
	}
	
}



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



// a population structure that will not communicate
boost::shared_ptr< hudson_ms::PopulationStructure > makePopIslandBad1(size_t Nzero)
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
		
		{
			// add a sub-pop
			size_t subNsamSize = 2;
			pop->add(subNsamSize);
		}
		
		
		{
			size_t i = 0;
			size_t j = 1;
			double _m = 0.1/(4.0*1000000); // fudge so that we get 4nZero*_m = 0.1 when Nzero = 1000000
			pop->setMigration(i,j,_m*4*Nzero);
		}
		{
			size_t i = 1;
			size_t j = 0;
			double _m = 0.1/(4.0*1000000); // fudge so that we get 4nZero*_m = 0.1 when Nzero = 1000000
			pop->setMigration(i,j,_m*4*Nzero);
		}
		
		return pop;	
	}
	catch (std::exception& e) {
	
		std::cout << "Error in makePopBad1:\n" + std::string(e.what()) << std::endl;
		throw;
	}
}

// another population structure that will not communicate
boost::shared_ptr< hudson_ms::PopulationStructure > makePopIslandBad2(size_t Nzero)
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
		
		{
			size_t i = 0;
			size_t j = 1;
			double _m = 0.1/(4.0*1000000); // fudge so that we get 4nZero*_m = 0.1 when Nzero = 1000000
			pop->setMigration(i,j,_m*4*Nzero);
		}
		{
			size_t i = 1;
			size_t j = 0;
			double _m = 0.1/(4.0*1000000); // fudge so that we get 4nZero*_m = 0.1 when Nzero = 1000000
			pop->setMigration(i,j,_m*4*Nzero);
		}
		{
			size_t i = 2;
			size_t j = 3;
			double _m = 0.1/(4.0*1000000); // fudge so that we get 4nZero*_m = 0.1 when Nzero = 1000000
			pop->setMigration(i,j,_m*4*Nzero);
		}
		{
			size_t i = 3;
			size_t j = 2;
			double _m = 0.1/(4.0*1000000); // fudge so that we get 4nZero*_m = 0.1 when Nzero = 1000000
			pop->setMigration(i,j,_m*4*Nzero);
		}
		
		return pop;	
	}
	catch (std::exception& e) {
	
		std::cout << "Error in makePopBad2:\n" + std::string(e.what()) << std::endl;
		throw;
	}
}



//example  on pg 14 ms doc
boost::shared_ptr< hudson_ms::PopulationStructure > makePopInstantSizeChange()
{
	// N0 is something or other
	// T1 = 4N0 generations back
	double t1 = 1.0; // T1/(4N0)
	// T2 = 8N0 generations back
	double t2 = 2.0; // T2/(4N0)
	boost::shared_ptr< hudson_ms::PopulationStructure > pop ( new hudson_ms::PopulationStructure );
		
	{
		// add a sub-pop
		size_t subNsamSize = 15;
		pop->add(subNsamSize);
	}
	std::vector< hudson_ms::HudsonMSDemogEventPtr > _eList;
			
	
	{
		boost::shared_ptr < const hudson_ms::HudsonMSDemogEvent > 
			ev(new hudson_ms::HudsonMSDemogEvent_N(t1, 0.1));
	
		_eList.push_back(ev);
	}
	{
		boost::shared_ptr < const hudson_ms::HudsonMSDemogEvent > 
			ev(new hudson_ms::HudsonMSDemogEvent_N(t2, 4.0));
	
		_eList.push_back(ev);
	}
	
	pop->setEventList(_eList);
			
	return pop;
}

//example  on pg 14 ms doc
boost::shared_ptr< hudson_ms::PopulationStructure > makePopOutgroupSequence()
{
	boost::shared_ptr< hudson_ms::PopulationStructure > pop ( new hudson_ms::PopulationStructure );
		
	// N0 is something or other
	// T1 = 24N0 generations back
	double t1 = 6.0; // T1/(4N0)
	
	{
		// add a sub-pop
		size_t subNsamSize = 1;
		pop->add(subNsamSize);
	}
	{
		// add a sub-pop
		size_t subNsamSize = 10;
		pop->add(subNsamSize);
	}
	
	std::vector< hudson_ms::HudsonMSDemogEventPtr > _eList;
			
	boost::shared_ptr < const hudson_ms::HudsonMSDemogEvent > 
			ev(new hudson_ms::HudsonMSDemogEvent_j(t1, 0, 1)); // pops 1,2, in ms command line terms
	
	_eList.push_back(ev);
	
	pop->setEventList(_eList);
			
	return pop;
}

//example  on pg 14 ms doc
boost::shared_ptr< hudson_ms::PopulationStructure > makePopInstantSizeChangeAndExpGrowth()
{
	int N_T2 = 10000;
	int N_T1 = 5000;
	int N_T0 = 20000;
	int N0 = N_T0;
	
	double T1 = 16000.0; //generations back
	double t1 = T1/(4*N0); // T1/(4N0)
	double T2 = 24000.0; //generations back
	double t2 = T2/(4*N0); // T1/(4N0)
	
	// increase in pop size by growth from T2 (size N_T1) to present (size N_T0)
	// N_T1 * exp(alphag * t1) = N_T0
	double alphag = (1/t1) * std::log((1.0*N_T0)/N_T1);
	
	boost::shared_ptr< hudson_ms::PopulationStructure > pop ( new hudson_ms::PopulationStructure );
		
	{
		// add a sub-pop
		size_t subNsamSize = 15;
		double growth = alphag; // should be 6.93;
		pop->add(subNsamSize, growth);
	}
		std::vector< hudson_ms::HudsonMSDemogEventPtr > _eList;
			
	{	boost::shared_ptr < const hudson_ms::HudsonMSDemogEvent > 
			ev(new hudson_ms::HudsonMSDemogEvent_G(t1, 0.0));
	
		_eList.push_back(ev);
	}
	
	{
		double ratio = (1.0*N_T2)/N0;
		boost::shared_ptr < const hudson_ms::HudsonMSDemogEvent > 
			ev(new hudson_ms::HudsonMSDemogEvent_N(t2, ratio));
	
		_eList.push_back(ev);
	}
	
	pop->setEventList(_eList);
	
	return pop;
}



// new example
boost::shared_ptr< hudson_ms::PopulationStructure > makePopBackwardsSplit()
{
	boost::shared_ptr< hudson_ms::PopulationStructure > pop ( new hudson_ms::PopulationStructure );
		
	{
		// add a sub-pop
		size_t subNsamSize = 10;
		double growth = 5.0;
		pop->add(subNsamSize, growth);
	}
	
	std::vector< hudson_ms::HudsonMSDemogEventPtr > _eList;
	
	{
		boost::shared_ptr < const hudson_ms::HudsonMSDemogEvent > 
			ev(new hudson_ms::HudsonMSDemogEvent_m(0.1, 0, 1, 0.5));  
	
	
		_eList.push_back(ev);
	}
	{
		boost::shared_ptr < const hudson_ms::HudsonMSDemogEvent > 
			ev(new hudson_ms::HudsonMSDemogEvent_m(0.1, 1, 0, 2.0));  
	
	
		_eList.push_back(ev);
	}
	{
		boost::shared_ptr < const hudson_ms::HudsonMSDemogEvent > 
			ev(new hudson_ms::HudsonMSDemogEvent_n(0.1, 1, 0.25));  // pop 2 in ms command line terms
	
	
		_eList.push_back(ev);
	}
	{
		boost::shared_ptr < const hudson_ms::HudsonMSDemogEvent > 
			ev(new hudson_ms::HudsonMSDemogEvent_s(0.1, 0, 0.5));  // pop 1 in ms command line terms
	
	
		_eList.push_back(ev);
	}
	
	
	pop->setEventList(_eList);
	
	return pop;
}

//example  on pg 15 ms doc
boost::shared_ptr< hudson_ms::PopulationStructure > makePopTwoSpeciesWithPopSizeDiff()
{
	int N_1_T1 = 10000;
	
	int N_2_T0 = 5000;
	int N_2_T2 = 2000;
	
	int N_1_T0 = 40000;
	int N0 = N_1_T0;
	
	double T1 = 5000.0; //generations back
	double t1 = T1/(4*N0); // T1/(4N0) // should be 0.03125
	double T2 = 10000.0; //generations back
	double t2 = T2/(4*N0); // T1/(4N0) // should be 0.0625
	double T3 = 15000.0; //generations back
	double t3 = T3/(4*N0); // T1/(4N0) // should be 0.09375

	// increase in first pop size by growth from T1 (size N_1_T1) to present (size N_1_T0)
	// N_1_T1 * exp(alphag1 * t1) = N_1_T0
	double alphag1 = (1/t1) * std::log((1.0*N_1_T0)/N_1_T1); // should be 44.36
	

	boost::shared_ptr< hudson_ms::PopulationStructure > pop ( new hudson_ms::PopulationStructure );
		
	{
		// add a sub-pop
		size_t subNsamSize = 3;
		pop->add(subNsamSize, alphag1);
	}
	{
		// add a sub-pop
		size_t subNsamSize = 12;
		double growth = 0.0;
		double size_ratio = (1.0*N_2_T0)/N0; // should be 0.125
		pop->add(subNsamSize, size_ratio, growth);
	}
	
	std::vector< hudson_ms::HudsonMSDemogEventPtr > _eList;
	
	{	boost::shared_ptr < const hudson_ms::HudsonMSDemogEvent > 
			ev(new hudson_ms::HudsonMSDemogEvent_g(t1, 0, 0.0)); // pop 1 in ms command line terms
		// no growth in first pop before t1
		_eList.push_back(ev);
	}
	
	{
		double size_ratio = (1.0*N_2_T2)/N0; // should be 0.05
		// second pop had instantaneous change in size at t2
		boost::shared_ptr < const hudson_ms::HudsonMSDemogEvent > 
			ev(new hudson_ms::HudsonMSDemogEvent_n(t2, 1, size_ratio));  // pop 2 in ms command line terms
		_eList.push_back(ev);
	}
	
	{
		// lineages in second pop moved to first pop at time t3
		boost::shared_ptr < const hudson_ms::HudsonMSDemogEvent > 
			ev(new hudson_ms::HudsonMSDemogEvent_j(t3, 1, 0));  // pop 2 moved to pop 1 in ms command line terms
	
		_eList.push_back(ev);
	}
	
	pop->setEventList(_eList);
	
	return pop;
}
//example  on pg 17 ms doc
boost::shared_ptr< hudson_ms::PopulationStructure > makePopSteppingStoneRecentBarrier1()
{
	// N0 is something or other
	// T1 = 8N0 generations back
	double t1 = 2.0; // T1/(4N0)
	
	boost::shared_ptr< hudson_ms::PopulationStructure > pop ( new hudson_ms::PopulationStructure );
	
	{
		// add a sub-pop
		size_t subNsamSize = 0;
		pop->add(subNsamSize);
	}
	{
		// add a sub-pop
		size_t subNsamSize = 7;
		pop->add(subNsamSize);
	}
	{
		// add a sub-pop
		size_t subNsamSize = 0;
		pop->add(subNsamSize);
	}
	{
		// add a sub-pop
		size_t subNsamSize = 0;
		pop->add(subNsamSize);
	}
	{
		// add a sub-pop
		size_t subNsamSize = 8;
		pop->add(subNsamSize);
	}
	{
		// add a sub-pop
		size_t subNsamSize = 0;
		pop->add(subNsamSize);
	}
	
	double _m = 2.5;
	
	{
		size_t i = 0;
		size_t j = 1;
		pop->setMigration(i,j,_m);
	}
	{
		size_t i = 1;
		size_t j = 0;
		pop->setMigration(i,j,_m);
	}
	{
		size_t i = 1;
		size_t j = 2;
		pop->setMigration(i,j,_m);
	}
	{
		size_t i = 2;
		size_t j = 1;
		pop->setMigration(i,j,_m);
	}
	{
		size_t i = 3;
		size_t j = 4;
		pop->setMigration(i,j,_m);
	}
	{
		size_t i = 4;
		size_t j = 3;
		pop->setMigration(i,j,_m);
	}
	{
		size_t i = 4;
		size_t j = 5;
		pop->setMigration(i,j,_m);
	}
	{
		size_t i = 5;
		size_t j = 4;
		pop->setMigration(i,j,_m);
	}
	
	std::vector< hudson_ms::HudsonMSDemogEventPtr > _eList;
	
	{	boost::shared_ptr < const hudson_ms::HudsonMSDemogEvent > 
			ev(new hudson_ms::HudsonMSDemogEvent_m(t1, 2, 3, _m));  // pops 3,4 in ms terms
	
		_eList.push_back(ev);
	}
	{	boost::shared_ptr < const hudson_ms::HudsonMSDemogEvent > 
			ev(new hudson_ms::HudsonMSDemogEvent_m(t1, 3, 2, _m)); // pops 4,3 in ms terms
	
		_eList.push_back(ev);
	}
	pop->setEventList(_eList);
	
	return pop;
}

//example  on pg 17 ms doc
boost::shared_ptr< hudson_ms::PopulationStructure > makePopSteppingStoneRecentBarrier2()
{
	// N0 is something or other
	// T1 = 8N0 generations back
	double t1 = 2.0; // T1/(4N0)

	boost::shared_ptr< hudson_ms::PopulationStructure > pop ( new hudson_ms::PopulationStructure );
		
	{
		// add a sub-pop
		size_t subNsamSize = 0;
		pop->add(subNsamSize);
	}
	{
		// add a sub-pop
		size_t subNsamSize = 7;
		pop->add(subNsamSize);
	}
	{
		// add a sub-pop
		size_t subNsamSize = 0;
		pop->add(subNsamSize);
	}
	{
		// add a sub-pop
		size_t subNsamSize = 0;
		pop->add(subNsamSize);
	}
	{
		// add a sub-pop
		size_t subNsamSize = 8;
		pop->add(subNsamSize);
	}
	{
		// add a sub-pop
		size_t subNsamSize = 0;
		pop->add(subNsamSize);
	}
	size_t n = pop->nSubPops();
	double _m = 2.5;
	std::vector < std::vector<double> > mm(n, std::vector<double>(n, 0.0));
	(mm.at(0)).at(1) = _m;
	(mm.at(1)).at(0) = _m;
	(mm.at(1)).at(2) = _m;
	(mm.at(2)).at(1) = _m;
	(mm.at(3)).at(4) = _m;
	(mm.at(4)).at(3) = _m;
	(mm.at(4)).at(5) = _m;
	(mm.at(5)).at(4) = _m;
	
	pop->setMigrationMatrix(mm);
	
	(mm.at(2)).at(3) = _m;
	(mm.at(3)).at(2) = _m;
	
	std::vector< hudson_ms::HudsonMSDemogEventPtr > _eList;
	
	{	
		boost::shared_ptr < const hudson_ms::HudsonMSDemogEvent > 
			ev(new hudson_ms::HudsonMSDemogEvent_ma(t1, mm));  
	
		_eList.push_back(ev);
	
	}
	
	pop->setEventList(_eList);
	
	return pop;
}

boost::shared_ptr< hudson_ms::PopulationStructure > makeTestPopForSS1()
{
	try {
	
		boost::shared_ptr< hudson_ms::PopulationStructure > pop ( new hudson_ms::PopulationStructure );
		
		{
			// add a sub-pop
			size_t subNsamSize = 10;
			pop->add(subNsamSize);
		}
			
		return pop;	
	}
	catch (std::exception& e) {
	
		std::cout << "Error in makeTestPopForSS1:\n" + std::string(e.what()) << std::endl;
		throw;
	}
}


