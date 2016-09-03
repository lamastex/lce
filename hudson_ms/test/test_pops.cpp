
/*
 * definition for populations testing
*/


#include "config.h"

#include "test_pops.hpp"
#include "test_pops_subs.hpp"
#include "make_test_pops.hpp"

#include <iostream>
#include <stdexcept>

using namespace std;


int main(int argc, char *argv[])
{
		//testPopStructReordering();
		//testPopStructBasic();
		//testPopStructMSstring();
		//testPopMigCheck();
		testPopStructExceptions();
		
		return 0;
	
}



// make a population and get it to make its ms command, then use ms to execute 
void testPopStructMSstring()

{ 
	try {
		
		
		cout << "\n\nTry making one-sub-pop population" << endl;
		boost::shared_ptr< hudson_ms::PopulationStructure > pop  = makeTestPop1();
		
		double g = 5.0;
		size_t seed = 1245;
		double rho = 100.0; // 4No r where r is cross proby over between sites in nsites 
		size_t nsites = 2500;
		cout << "seed = " << seed << endl;
		
		pop->setAllGrowthRates(g);
		
		std::string filename = "testPopStruct1.txt";
		
		hudson_ms::HudsonMSRecombParams recombParams(rho, nsites);
		
		testPopStructSub(recombParams, pop, filename, seed);
			
	}
	catch (exception& ee) {
		std::cout << "That resulted in exception\n" 
			+ std::string(ee.what()) << std::endl;
	}	
	
	try {
		
		cout << "\n\nTry making more complex population" << endl;
		
		boost::shared_ptr< hudson_ms::PopulationStructure > pop  = makeTestPop2();
		
		double g = 5.0;
		size_t seed = 1245;
		double rho = 100.0; // 4No r where r is cross proby over between sites in nsites 
		size_t nsites = 2500;
		cout << "seed = " << seed << endl;
		
		pop->setAllGrowthRates(g);
		
		hudson_ms::HudsonMSRecombParams recombParams(rho, nsites);
		
		std::string filename = "testPopStruct2.txt";
		
		testPopStructSub(recombParams, pop, filename, seed);
		
		cout << "\n\nEnd of testing\n\n" << endl; 	
			
	}
	catch (exception& ee) {
		std::cout << "That resulted in exception\n" 
			+ std::string(ee.what()) << std::endl;
	}
	catch (...) {
		std::cout << "That resulted in an unknown exception\n" << std::endl;
		throw;
	}	
		
}



// test the new population structures - basic testing
void testPopStructBasic()

{ 
	cout << "\n\nTest new population structure code: constructors \n" << endl;
		
	try {
		
		cout << "\n\nConstructing an empty population" << endl;
		
		hudson_ms::PopulationStructure pop;
		
		testPopConstructorSub(pop);
	}
	catch (exception& ee) {
		std::cout << "That resulted in exception\n" 
			+ std::string(ee.what()) << std::endl;
		throw;
	}	
	
	try {
		
		for (int nsam = 0; nsam < 2; ++nsam) {
			cout << "\n\nConstructing a population with one subpopulation, nsam = " << nsam << endl;
			
			hudson_ms::PopulationStructure pop(nsam);
			
			testPopConstructorSub(pop);
		}
	}
	catch (exception& ee) {
		std::cout << "That resulted in exception\n" 
			+ std::string(ee.what()) << std::endl;
		throw;
	}
	
	try {
		
		cout << "\n\nConstructing an empty population called pop" << endl;
			
		hudson_ms::PopulationStructure pop;
		
		cout << "\nAdding subpopulations:" << endl;
		
		int nsps = 4;
			
		for (int nsam = 0; nsam < nsps; ++nsam) {
			
			if (nsam < 2) {
				cout << "With pop.add(nsam) with nsam = " << nsam << endl;
				pop.add(nsam);
			}	
			else if (nsam == 2) {
				double sizeRatio = 0.5;
				double growthRate = 10.0;
				cout << "with pop.add(nsam, sizeRatio, growthRate) ";
				cout << "with nsam = " << nsam << ", sizeRatio = " << sizeRatio << ", growthRate = " << growthRate << endl;
				pop.add(nsam, sizeRatio, growthRate);
			}
			else {
				double growthRate = 3.5;
				cout << "With pop.add(nsam, growthRate) ";
				cout << "with nsam = " << nsam << ", growthRate = " << growthRate << endl;
				
				pop.add(nsam, growthRate);
			}
		}
		cout << "\nAfter adding subpopulations:" << endl;
		testPopConstructorSub(pop);
		
		
		cout << "\nSetting growth rates:" << endl;
			
		for (int i = 0; i < nsps; ++i) {
			
			double gr = i/1.0;
			std::cout << "pop.setGrowthRateAt(" << i << ", " << gr << ")" << endl; 
			pop.setGrowthRateAt(i, gr);
		}
		cout << "\nAfter setting growth rates:" << endl;
		cout << pop.toString() << endl;
				
		cout << "\n\nCopy constructor: PopulationStructure popCopy(pop) :" << endl;
		hudson_ms::PopulationStructure popCopy(pop);
		cout << popCopy.toString() << endl;
		
		double allGrowth1 = 10.4;
		cout << "\nChange all growth rates in popCopy: popCopy.setAllGrowthRates(" << allGrowth1 << ")"<< endl;
		popCopy.setAllGrowthRates(allGrowth1);
		cout << "Now we have popCopy as:" << endl;
		cout << popCopy.toString() << endl;
		
		cout << "and pop as:" << endl;
		cout << pop.toString() << endl;
		
		cout << "\n\nAssignment: PopulationStructure popAssign = pop :" << endl;
		hudson_ms::PopulationStructure popAssign = pop;
		cout << popAssign.toString() << endl;
		
		double allMigs = 5.0;
		cout << "\nSet all inter-subpop migration parameters in popAssign: popAssign.setMigrationMatrix(" << allMigs << ")"<< endl;
		popAssign.setMigrationMatrix(allMigs);
		cout << "Now we have popAssign as:" << endl;
		cout << popAssign.toString() << endl;
		
		
		cout << "\n\nstd::swap(popAssign, popCopy)" << endl;
		std::swap(popAssign, popCopy);
		cout << "Now we have popAssign as:" << endl;
		cout << popAssign.toString() << endl;
		
		cout << "and popCopy as:" << endl;
		cout << popCopy.toString() << endl;
		
		cout << "\nTest referring to subpops by order and by label" << endl;
		for (int i = 0; i < nsps; ++i) {
			
			testGetAndSetSimple(pop, i);
		}
		
		
		cout << "\nTest getting and setting migration parameters with PopCopy" << endl;
		testGetAndSetMigrationSimple(popCopy);
		cout << "\npopCopy now is:" << endl;
		cout << popCopy.toString() << endl;
		
		cout << "\nGet full migration matrix from popCopy and give it to pop:" << endl;
		cout << "vector < vector < double > > mm = popCopy.makeMigration_matrix();" << endl;
		cout << "pop.setMigrationMatrix(mm);" << endl;
		vector < vector < double > > mm = popCopy.makeMigration_matrix();
		pop.setMigrationMatrix(mm);
		cout << "pop now is:" << endl;
		cout << pop.toString() << endl;
		
		cout << "\npopAssign is currently:" << endl;
		cout << popAssign.toString() << endl;
		cout << "\nAssign PopAssign = pop" << endl;
		popAssign = pop;
		cout << "popAssign now is:" << endl;
		cout << popAssign.toString() << endl;
		
		cout << "\n\nEnd of testing\n\n" << endl; 	
		
	}
	catch (exception& ee) {
		std::cout << "That resulted in exception\n" 
			+ std::string(ee.what()) << std::endl;
		throw;
	}	
	catch (...) {
		std::cout << "That resulted in an unknown exception\n" << std::endl;
		throw;
	}	
	
	
		
}




void testPopStructExceptions()
{
	try {
		
		cout << "\n\nConstructing an empty population called pop" << endl;
			
		hudson_ms::PopulationStructure pop;
		
		cout << "pop output and mig check" << endl;
		testPopOutputWithMigCheck(pop);
		
		
		cout << "Try to access subpop characteristics that are out of range" << endl;
		{
			int index = 0;
			popOutOfRangeTestSub(pop, index);
		}
		
		
		{
			cout << "\n\nAdding subpopulations with invalid parameters:" << endl;
			{
				int nsam = -1; //x
				double sr = 0.0;
				double gr = 0.0;
				testSubPopsInvalidArgSub(nsam, sr, gr);
			}
			{
				int nsam = 0;
				double sr = 0; //x
				double gr = 0;
				testSubPopsInvalidArgSub(nsam, sr, gr);
			}
			{
				int nsam = 0;
				double sr = -0.1; //x
				double gr = 0.0;
				testSubPopsInvalidArgSub(nsam, sr, gr);
			}
			{
				int nsam = 0;
				double sr = 0.0; //x
				double gr = -0.1;
				testSubPopsInvalidArgSub(nsam, sr, gr);
			}
			{
				int nsam = -1; //x
				double sr = -0.1; //x
				double gr = -0.1;
				testSubPopsInvalidArgSub(nsam, sr, gr);
			}
			{
				int nsam = 0;
				double sr = 0.1;
				double gr = 0.1;
				testSubPopsInvalidArgSub(nsam, sr, gr);
			}
			
			
			
		}
		int nsps = 4;
		cout << "\nAdding " << nsps << " subpopulations:" << endl;
		
			
		for (int nsam = 0; nsam < nsps; ++nsam) {
			
			if (nsam < 1) {
				#if (0)
				cout << "With pop.add(nsam) with nsam = " << nsam << endl;
				#endif
				pop.add(nsam);
			}	
			else if (nsam < 3) {
				double sizeRatio = 1.0/nsam;
				double growthRate = (nsam)*1.0;
				#if (0)
				cout << "with pop.add(nsam, sizeRatio, growthRate) ";
				cout << "with nsam = " << nsam << ", sizeRatio = " << sizeRatio << ", growthRate = " << growthRate << endl;
				#endif
				pop.add(nsam, sizeRatio, growthRate);
			}
			else {
				double growthRate = (nsam)*1.0;
				#if (0)
				cout << "With pop.add(nsam, growthRate) ";
				cout << "with nsam = " << nsam << ", growthRate = " << growthRate << endl;
				#endif
				pop.add(nsam, growthRate);
			}
		}
		cout << "\nAfter adding subpopulations, pop output and mig check:" << endl;
		testPopOutputWithMigCheckFailure(pop);
		
		cout << "\n\nTry to access subpop characteristics with indices that are in range" << endl;
		{
			int index = 0;
			popOutOfRangeTestSub(pop, index);
		}
		{
			int index = nsps-1;
			popOutOfRangeTestSub(pop, index);
		}
		
		cout << "\n\nTry to access subpop characteristics with indices that are out of range" << endl;
		{
			int index = -1;
			popOutOfRangeTestSub(pop, index);
		}
		{
			int index = nsps;
			popOutOfRangeTestSub(pop, index);
		}
		
		cout << "\n\nTry to access subpop characteristics with labels that exist" << endl;
		{
			int index = 1466%nsps; //some number in the index range
			std::string label = pop.getLabelAt(index);
			popNotRecognisedTestSub(pop, label);
		}
		
		
		cout << "\n\nTry to access subpop characteristics with labels that do not exist" << endl;
		{
			int index = 1666%nsps; //some number in the index range
			std::string label = pop.getLabelAt(index);
			label = label + "a";
			popNotRecognisedTestSub(pop, label);
		}
		
		
		// setting values
		cout << "\n\nTry to set properties to invalid values" << endl;
		{
			int nsam = -1; //x
			double sr = 0.0; //x
			double gr = 0.0;
			testSubPopsSetInvalidArgSub(pop,nsam, sr, gr);
		}
		{
			int nsam = 0; 
			double sr = -0.01; //x
			double gr = -0.01;
			testSubPopsSetInvalidArgSub(pop,nsam, sr, gr);
		}
		
		// check on relabelling
		{
			int nsam = 1;
			hudson_ms::PopulationStructure popSingle(nsam);
			cout << "\n\nTry relabelling tests with a population with nSubPop = " << ( popSingle.nSubPops() ) << endl;
			testPopRelabellingInvalidArg(popSingle);
		}
		
		// check on relabelling
		{
			cout << "\n\nTry relabelling tests with a pop which has nSubPop = " << ( pop.nSubPops() ) << endl;
			testPopRelabellingInvalidArg(pop);
		}
		
		// 
		{
			cout << "\n\nTest reorder and erase with eList state" << endl; 
	
			testPopBeforeAndAfterEvents(pop);
		}
		
		cout << "\n\nEnd of testing\n\n" << endl; 
		
				
	}
	catch (exception& ee) {
		std::cout << "That resulted in a std::exception\n" 
			+ std::string(ee.what()) << std::endl;
		throw;
	}	
	catch (...) {
		std::cout << "That resulted in an unknown exception\n" << std::endl;
		throw;
	}	
}

void testPopStructReordering()
{
	try {
		
		cout << "\n\nConstructing an empty population called pop" << endl;
			
		hudson_ms::PopulationStructure pop;
		
		cout << "\nAdding subpopulations:" << endl;
		
		int nsps = 4;
			
		for (int nsam = 1; nsam <= nsps; ++nsam) {
			
			if (nsam < 2) {
				cout << "With pop.add(nsam) with nsam = " << nsam << endl;
				pop.add(nsam);
			}	
			else if (nsam < 4) {
				double sizeRatio = 1.0/nsam;
				double growthRate = (nsam-1)*1.0;
				cout << "with pop.add(nsam, sizeRatio, growthRate) ";
				cout << "with nsam = " << nsam << ", sizeRatio = " << sizeRatio << ", growthRate = " << growthRate << endl;
				pop.add(nsam, sizeRatio, growthRate);
			}
			else {
				double growthRate = (nsam-1)*1.0;
				cout << "With pop.add(nsam, growthRate) ";
				cout << "with nsam = " << nsam << ", growthRate = " << growthRate << endl;
				
				pop.add(nsam, growthRate);
			}
		}
		cout << "\nAfter adding subpopulations:" << endl;
		cout << pop.toString() << endl;
		
		vector < string > labels;
		labels.push_back("1"); labels.push_back("2"); labels.push_back("3"); labels.push_back("4");
		
		{//relabel
			int order1 = 0;
			string newlabel1("Originally1");
			testPopRelabellingSubByOrder(pop, labels, order1, newlabel1);
			
			int order3 = 2;
			string newlabel3("Originally3");
			testPopRelabellingSubByLabel(pop, labels, order3, newlabel3);
			
		}
		
		//set some migration rates
		cout << "\nSet some migration rates" << endl;
		
		{
			int order1 = 0;
			int order21 = 1;
			double m1 = 3.0;
			int order22 = 2;
			double m2 = 4.0;
			testPopMigSetSub(pop, labels, order1, order21, order22, m1, m2);
		}
		{
			int order1 = 1;
			int order21 = 2;
			double m1 = 0.1;
			int order22 = 3;
			double m2 = 1.0;
			testPopMigSetSub(pop, labels, order1, order21, order22, m1, m2);
			
		}
		{
			int order1 = 2;
			int order21 = 0;
			double m1 = 5.0;
			int order22 = 1;
			double m2 = 10.0;
			testPopMigSetSub(pop, labels, order1, order21, order22, m1, m2);
		}
		{
			int order1 = 3;
			int order21 = 0;
			double m1 = 0.6;
			int order22 = 2;
			double m2 = 0.7;
			testPopMigSetSub(pop, labels, order1, order21, order22, m1, m2);
		}
		
		testPopOutputWithMigCheck(pop);
		
		
		
		{ // reordering
			vector < string > newOrder;
			newOrder.push_back(labels[1]); newOrder.push_back(labels[0]); newOrder.push_back(labels[2]); newOrder.push_back(labels[3]);
			
			testPopOrderingSub(pop, newOrder);
			labels = pop.getLabels();
		}
		
		{	
			// resetting some migration parameters " << endl; 
			cout << endl;
			int order1 = 0;
			int order2 = 3;
			double m = 777.0;
			
			testPopMigSetSubByOrder(pop, labels, order1, order2, m);
			testPopOutputWithMigCheck(pop);
		}
		{
			// reordering
			vector < string > newOrder;
			newOrder.push_back(labels[0]); newOrder.push_back(labels[1]); newOrder.push_back(labels[3]); newOrder.push_back(labels[2]);
			
			testPopOrderingSub(pop, newOrder);
			labels = pop.getLabels();
		}
		{	// resetting some migration parameters
			cout << endl;
			int order1 = 2;
			int order2 = 1;
			double m = 999.0;
			
			testPopMigSetSubByOrder(pop, labels, order1, order2, m);
			testPopOutputWithMigCheck(pop);
		}
		
		{
			// erasing
			cout << endl;
			int orderErase = 0;
			testPopEraseByOrder(pop, labels, orderErase);
			labels = pop.getLabels();
		}
		
		{
			// relabelling in bulk
			vector < string > newlabels = labels;
			int order = 1;
			newlabels.at(order) = string("Originally4");
			
			testPopRelabellingSub(pop, newlabels);
			labels = pop.getLabels();
		}
		
		{	// resetting some migration parameters
			cout << endl;
			int order1 = 2;
			int order2 = 1;
			double m = 1999.0;
			
			testPopMigSetSubByOrder(pop, labels, order1, order2, m);
			testPopOutputWithMigCheck(pop);
		}
		
		{
			// erasing
			cout << endl;
			int orderErase = 1;
			testPopEraseByLabel(pop, labels, orderErase);
			labels = pop.getLabels();
		}
			
		cout << "\n\nEnd of testing\n\n" << endl; 
		
	}
	catch (exception& ee) {
		std::cout << "That resulted in exception\n" 
			+ std::string(ee.what()) << std::endl;
		throw;
	}
	catch (...) {
		std::cout << "That resulted in an unknown exception\n" << std::endl;
		throw;
	}		
}


void testPopMigCheck()

{ 
	size_t Nzero = 1000000;
	try {
		try {
			cout << "makePopSingle()" << endl;
			boost::shared_ptr< hudson_ms::PopulationStructure > pop  = makePopSingle();
			
			testPopMigCheckSub(pop);
			
		}
		catch (std::logic_error& le) {
			std::cout << "That resulted in a logic error exception\n" 
				+ std::string(le.what()) << std::endl;
		}	
		
		try {
			cout << "makePopIslandBasic()" << endl;
			boost::shared_ptr< hudson_ms::PopulationStructure > pop  = makePopIslandBasic(Nzero);
			
			testPopMigCheckSub(pop);
			
		}
		catch (std::logic_error& le) {
			std::cout << "That resulted in a logic error exception\n" 
				+ std::string(le.what()) << std::endl;
		}
		
		try {
			cout << "makePopIslandComplex()" << endl;
			boost::shared_ptr< hudson_ms::PopulationStructure > pop  = makePopIslandComplex(Nzero);
			
			testPopMigCheckSub(pop);
			
		}
		catch (std::logic_error& le) {
			std::cout << "That resulted in a logic error exception\n" 
				+ std::string(le.what()) << std::endl;
		}	
		
		
		try {
			cout << "Try makePopIslandBad1()" << endl;
			boost::shared_ptr< hudson_ms::PopulationStructure > pop  = makePopIslandBad1(Nzero);
			
			testPopMigCheckSub(pop);
			
		}
		catch (std::logic_error& le) {
			std::cout << "That resulted in a logic error exception\n" 
				+ std::string(le.what()) << std::endl;
		}	
		
		try {
			cout << "makePopIslandBad2()" << endl;
			boost::shared_ptr< hudson_ms::PopulationStructure > pop  = makePopIslandBad2(Nzero);
			
			testPopMigCheckSub(pop);
		}	
		catch (std::logic_error& le) {
			std::cout << "That resulted in a logic error exception\n" 
				+ std::string(le.what()) << std::endl;
		}
		
		cout << "\n\nEnd of testing\n\n" << endl; 	
	}
	catch (std::exception& e) {
	
		std::cout << "Standard exception - Error:\n" + std::string(e.what()) << std::endl;
	}
	catch (...) {
	
		std::cout << "Unknown type exception" << std::endl;
	}
}

