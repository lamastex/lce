
/*
 * definition for populations testing subroutines
*/


#include "config.h"

#include "test_pops_subs.hpp"


#include <mscplusplus/event.hpp>

#include <sstream>
#include <fstream>
#include <iostream>
#include <iterator>
#include <cassert>
#include <stdexcept>


using namespace std;





std::string getDefaultOutputDir()
{
	return std::string("../outputPopTests");
}

void makeOutputDir(const std::string& outputDir)
{
	std::ostringstream stm;
	stm << "mkdir " << outputDir;
	std::string cmd = stm.str();
	system(cmd.c_str());
	
}

/*Open a file and return the data as a string */
std::string string_from_file(const std::string name)
{
    std::ifstream is(name.c_str());
    if (!is.is_open())
        throw std::runtime_error("unable to open file '" + name + "'");
    
    return std::string((std::istreambuf_iterator<char>(is)),
                 std::istreambuf_iterator<char>());
}

void testPopStructSub(const hudson_ms::HudsonMSRecombParams& recombParams,
		const boost::shared_ptr< hudson_ms::PopulationStructure >& pop,
		const std::string& filename,
		int seed)
{
	
		
	cout << recombParams.toString() << endl;
	
	// get the ms command
	cout << "Try making the ms command using this pop";
	
	std::string thisPop = pop->toString();
		
	cout << "\nThe string description for the population is " << endl;
	cout << thisPop << "\n" << endl;
	
	std::string cmd = pop->makeMScommand(recombParams, seed);
		
	cout << "The ms command is " << endl;
	cout << cmd << endl;
	
	std::string outputDir = getDefaultOutputDir();
	
	makeOutputDir(outputDir);
	
	std::string pathAndFilename = outputDir + std::string("/") + filename;
	
	// try executing it 
	std::ostringstream stm;
	stm << cmd << "| grep ';$' > " << pathAndFilename;
	cmd = stm.str();
	cout << "***** cmd = " << cmd << endl;
	system(cmd.c_str());
	std::string tree = string_from_file(pathAndFilename);
	cout << "The tree we got back is " << endl;
	cout << tree << endl;

}


void testPopConstructorSub(const hudson_ms::PopulationStructure& pop)
{
	
	cout << "pop.empty() = ";
	bool popEmpty = pop.empty();
	cout << popEmpty << endl;
	
	cout << "pop.nSubPops() = ";
	size_t nsp = pop.nSubPops();
	cout << nsp << endl;
	
	cout << "pop.totalNsam() = ";
	size_t nsam = pop.totalNsam();
	cout << nsam << endl;
	
	cout << "pop.eListEmpty() = ";
	bool eListEmpty = pop.eListEmpty();
	cout << eListEmpty << endl;
	
	cout << "pop.toString() is" << endl;
	std::string toString = pop.toString();
	cout << toString << endl;
		
}


// assumes the subpop label will the the string version of the order
void testGetAndSetSimple(hudson_ms::PopulationStructure& pop, int order)
{
	ostringstream stm;
	stm << order+1;
	string label = stm.str();
	
	cout << "Test with index = " << order << ", label = \"" <<label << "\"" << endl;
	
	{
		size_t nsam1 = pop.getSubNsamAt(order);
		cout << "getSubNsamAt(int index) = " << nsam1 << endl;
		size_t nsam1c = nsam1;
		nsam1 += 3;
		cout << "nsam1 += 3 = " << nsam1;
		cout << " and setSubNsamAt(order, nsam1)" << endl;
		pop.setSubNsamAt(order, nsam1);
		size_t nsam2 = pop.getSubNsamAt(label);
		cout << "getSubNsamAt(label) = " << nsam2 << endl;
		assert(nsam2 == nsam1);
		nsam2 -= 3;
		cout << "nsam2 -= 0.5 = " << nsam2;
		cout << " and setSubNsamAt(label, nsam2)" << endl;
		pop.setSubNsamAt(label, nsam2);
		size_t nsam3 = pop.getSubNsamAt(label);
		cout << "getSubNsamAt(const std::string& label) = " << nsam3 << endl;
		assert(nsam3 == nsam1c);
	}
	{
		double sr1 = pop.getSizeRatioAt(order);
		cout << "getSizeRatioAt(int index) = " << sr1 << endl;
		double sr1c = sr1;
		sr1 += 0.5;
		cout << "sr1 += 0.5 = " << sr1;
		cout << " and setSizeRatioAt(order, sr1)" << endl;
		pop.setSizeRatioAt(order, sr1);
		double sr2 = pop.getSizeRatioAt(label);
		cout << "getSizeRatioAt(label) = " << sr2 << endl;
		assert(sr2 == sr1);
		sr2 -= 0.5;
		cout << "sr2 -= 0.5 = " << sr2;
		cout << " and setSizeRatioAt(label, sr2)" << endl;
		pop.setSizeRatioAt(label, sr2);
		double sr3 = pop.getSizeRatioAt(label);
		cout << "getSizeRatioAt(const std::string& label) = " << sr3 << endl;
		assert(sr3 == sr1c);
	}
	{
		double gr1 = pop.getGrowthRateAt(order);
		double gr1c = gr1;
		cout << "getGrowthRateAt(order) = " << gr1 << endl;
		gr1 += 10;
		cout << "gr1 += 10 = " << gr1;
		cout << " and setGrowthRateAt(order, gr1)" << endl;
		pop.setGrowthRateAt(order, gr1);
		double gr2 = pop.getGrowthRateAt(label);
		cout << "getGrowthRateAt(label) = " << gr2 << endl;
		gr2 -= 10;
		cout << "gr2 -= 10 = " << gr2;
		cout << " and setGrowthRateAt(label, gr2)" << endl;
		pop.setGrowthRateAt(label, gr2);
		double gr3 = pop.getGrowthRateAt(label);
		cout << "getGrowthRateAt(const std::string& label) = " << gr3 << endl;
		assert(gr3 == gr1c);
	}
		
}

// assumes the subpop label will the the string version of the order
void testGetAndSetMigrationSimple(hudson_ms::PopulationStructure& pop)
{
	
	int nsps = static_cast<int>(pop.nSubPops());
	int orderBig = 1333;
	
	int order1 = (orderBig)%nsps;
	
	ostringstream stm;
	stm << order1+1;
	string label1 = stm.str();
	
	int order2 = (order1+1)%nsps;
	ostringstream stm2;
	stm2 << order2+1;
	string label2 = stm2.str();
	
	int order3 = (order1+3)%nsps;
	ostringstream stm3;
	stm3 << order3+1;
	string label3 = stm3.str();
	
	
	cout << "Test with order1 = " << order1 << ", label1 = \"" <<label1 << "\"" << endl;
	
	{
		cout << "and order2 = " << order2 << ", label2 = \"" <<label2 << "\"" << endl;
		
		double m1 = pop.getMigration(order1, order2);
		cout << "m1 = pop.getMigration(order1, order2) = " << m1 << endl;
		
		m1 *= 2;
		cout << "m1 *= 2 = " << m1 << endl;
		cout << "pop.setMigration(order1, order2, m1)" << endl;
		pop.setMigration(order1, order2, m1);
		
		double m2 = pop.getMigration(label1, label2);
		cout << "m2 = pop.getMigration(label1, label2) = " << m2 << endl;
		
		assert(m1 == m2);
	}
	
	cout << "Test with order1 = " << order2 << ", label1 = \"" <<label2 << "\"" << endl;
	
	{
		cout << "and order2 = " << order3 << ", label2 = \"" <<label3 << "\"" << endl;
		
		double m1 = pop.getMigration(label2, label3);
		cout << "m1 = pop.getMigration(label1, label2) = " << m1 << endl;
		
		m1 /= 2;
		cout << "m1 /= 2 = " << m1 << endl;
		cout << "pop.setMigration(label1, label2, m1)" << endl;
		pop.setMigration(label2, label3, m1);
		double m2 = pop.getMigration(order2, order3);
		cout << "m2 = pop.getMigration(order1, order2) = " << m2 << endl;
		
		assert(m1 == m2);
		
	}
	cout << "Test with order1 = " << order1 << ", label1 = \"" <<label1 << "\"" << endl;
	{
		
		cout << "and order2 = " << order3 << ", label2 = \"" <<label3 << "\"" << endl;
		
		double m1 = pop.getMigration(order1, order3);
		cout << "m1 = pop.getMigration(order1, order2) = " << m1 << endl;
		
		m1 *= 4;
		cout << "m1 *= 4 = " << m1 << endl;
		cout << "pop.setMigration(label1, label2, m1)" << endl;
		pop.setMigration(label1, label3, m1);
		double m2 = pop.getMigration(order1, order3);
		cout << "m2 = pop.getMigration(order1, order2) = " << m2 << endl;
	
		assert(m1 == m2);
	}
	cout << "Test with order1 = " << order3 << ", label1 = \"" <<label3 << "\"" << endl;
	{
		
		cout << "and order2 = " << order1 << ", label2 = \"" <<label1 << "\"" << endl;
	
	
		
		double m1 = pop.getMigration(label3, label2);
		cout << "m1 = pop.getMigration(label1, label2) = " << m1 << endl;
		
		m1 /= 3;
		cout << "m1 /= 3 = " << m1 << endl;
		cout << "pop.setMigration(order1, order2, m1)" << endl;
		pop.setMigration(order3, order2, m1);
		double m2 = pop.getMigration(label3, label2);
		cout << "m2 = pop.getMigration(label1, label2) = " << m2 << endl;
		
		assert(m1 == m2);
	}
		
}




void testPopOutputWithMigCheck(const hudson_ms::PopulationStructure& pop)
{
	cout << "\npop now is " << endl;
	cout << pop.toString() << endl;
	cout << "pop.migrationMatrixToString() is " << endl;
	cout << pop.migrationMatrixToString() << endl;
	
	cout << "\npop.doMigCheck()" << endl;
	std::string migCheckResult = pop.doMigCheck();
		
	if (!migCheckResult.empty()) {
		throw std::logic_error( migCheckResult);
	}
}

// expecting a failure in migCheck
void testPopOutputWithMigCheckFailure(const hudson_ms::PopulationStructure& pop)
{
	cout << "\npop now is " << endl;
	cout << pop.toString() << endl;
	cout << "pop.migrationMatrixToString() is " << endl;
	cout << pop.migrationMatrixToString() << endl;
	
	try {
		cout << "\npop.doMigCheck()" << endl;
		std::string migCheckResult = pop.doMigCheck();
		
		if (!migCheckResult.empty()) {
			throw std::logic_error( migCheckResult);
		}
	}
	catch (std::logic_error& le) {
		cout << "std::logic_error :\n" << le.what() << endl;
	}
}

void testPopMigSetSub(hudson_ms::PopulationStructure& pop,
			const vector < string >& labels,
			int order1, int order21, int order22,
			double m1, double m2)
{
	testPopMigSetSubByLabel(pop, labels, order1, order21, m1);
	testPopMigSetSubByOrder(pop, labels, order1, order22, m2);
	
}

void testPopMigSetSubByLabel(hudson_ms::PopulationStructure& pop,
			const vector < string >& labels,
			int order1, int order21, double m1)
{
	string label1 = labels.at(order1);
	string label2 = labels.at(order21);
	
	double m_ex = pop.getMigration(order1, order21);
	if (m_ex > 0.0) {
		cout << "Change migration parameter into " << order1 << " from " << order21 << endl;
		cout << "From getMigration(" << order1 << ", " << order21 << ") = " << m_ex << " to " << m1 << endl;
	}
	else {
		cout << "Setting migration parameter into \"" << label1 << "\" from \"" << label2 << "\"" << endl;
		cout << "setMigration(\"" << label1 << "\", \"" << label2 << "\", " << m1 << ")"<< endl;
	}
	
	pop.setMigration(label1, label2, m1);
	double check = pop.getMigration(order1, order21);
	assert (m1 == check); 
	cout << "pop.getMigration(" << order1 << ", " <<  order21 << ") = " << check << endl;
}

void testPopMigSetSubByOrder(hudson_ms::PopulationStructure& pop,
			const vector < string >& labels,
			int order1, int order22, double m2)
{
	string label1 = labels.at(order1);
	string label2 = labels.at(order22);
	
	double m_ex = pop.getMigration(label1, label2);
	if (m_ex > 0.0) {
		cout << "Change migration parameter into \"" << label1 << "\" from \"" << label2 << "\"" << endl;
		cout << "From getMigration(\"" << label1 << "\", \"" << label2 << "\") = " << m_ex << " to " << m2 << endl;
	}
	else {
		cout << "Setting migration parameter into " << order1 << " from " << order22 << endl;
		cout << "setMigration(" << order1 << ", " << order22 << "\", " << m2 << ")" << endl;
	}
	
	pop.setMigration(order1, order22, m2);
	double check = pop.getMigration(label1, label2);
	assert (m2 == check); 
	cout << "pop.getMigration(\"" << label1 << "\", \"" <<  label2 << "\") = " << check << endl;
}

void testPopEraseByOrder(hudson_ms::PopulationStructure& pop,
			const vector < string >& labels,
			int orderErase)
{
	string label = labels.at(orderErase);
	
	cout << "Deleting population at index " << orderErase << " (label \"" << label << "\")" << endl;
	cout << "pop.Erase(" << orderErase << ")" << endl;
	
	pop.erase(orderErase);
	
	testPopOutputWithMigCheck(pop);
}

void testPopEraseByLabel(hudson_ms::PopulationStructure& pop,
			const vector < string >& labels,
			int orderErase)
{
	string label = labels.at(orderErase);
	
	cout << "Deleting population at index " << orderErase << " (label \"" << label << "\")" << endl;
	cout << "pop.Erase(\"" << label << "\")" << endl;
	pop.erase(label);
	
	testPopOutputWithMigCheck(pop);
}

void testPopRelabellingSubByOrder(hudson_ms::PopulationStructure& pop,
						vector < string >& labels,
						int order, const string& newlabel)
{
	//relabel
	cout << "\nRelabel(int = " << order << ", string = " << newlabel << ")" << endl;
	pop.relabel(order, newlabel);
	
	labels = pop.getLabels();
	
	cout << "After relabelling:" << endl;
	cout << pop.toString() << endl;
}

void testPopRelabellingSubByLabel(hudson_ms::PopulationStructure& pop,
						vector < string >& labels,
						int order, const string& newlabel)

{
	string oldLabel = labels[order];
	
	
	cout << "\nRelabel(string = \"" << oldLabel << "\", string = " << newlabel << ")" << endl;
	pop.relabel(oldLabel, newlabel);
	
	labels = pop.getLabels();
	
	cout << "After relabelling:" << endl;
	cout << pop.toString() << endl;
	
	
}

void testPopRelabellingSub(hudson_ms::PopulationStructure& pop,
						const vector < string >& newlabels)
{
	cout  << "\nRelabelling to newlabels =" << endl;
	ostream_iterator<string> out_it (cout,", ");
	copy ( newlabels.begin(), newlabels.end(), out_it );
	cout << "\npop.relabel(newlabels)" << endl;
	pop.relabel(newlabels);
	
	cout << "\npop now is " << endl;
	cout << pop.toString() << endl;
	cout << "pop.migrationMatrixToString() is " << endl;
	cout << pop.migrationMatrixToString() << endl;
	
}

void testPopOrderingSub(hudson_ms::PopulationStructure& pop,
						const vector < string >& newOrder)
{
	cout  << "\nReordering to newOrder =" << endl;
	ostream_iterator<string> out_it (cout,", ");
	copy ( newOrder.begin(), newOrder.end(), out_it );
	cout << "\npop.resetOrdering(newOrder)" << endl;
	pop.resetOrdering(newOrder);
	
	testPopOutputWithMigCheck(pop);
}


void testPopMigCheckSub(boost::shared_ptr< hudson_ms::PopulationStructure > pop)
{
	
	// get the ms command
	cout << "Try making the ms command using this pop";

	double g = 5.0;
	size_t seed = 1245;
	double rho = 0.0; // 4No r where r is cross proby over between sites in nsites 
	size_t nsites = 2500;
	cout << "growth " << g << ", rho " << ", nsites " << nsites << " and seed " << seed << endl;
	
	pop->setAllGrowthRates(g);
	
	hudson_ms::HudsonMSRecombParams recombParams(rho, nsites);
			
	std::string cmd = pop->makeMScommand(recombParams, seed);
			
	cout << "The ms command is " << endl;
	cout << cmd << endl;
	
	cout << "do migCheck " << endl;
	std::string migCheckResult = pop->doMigCheck();
		
	if (!migCheckResult.empty()) {
		throw std::logic_error( migCheckResult);
	}
	
	cout << "Passed migCheck\n " << endl;
	
}


void throwLogicErrorCheck(bool outOfRange)
{
	if (outOfRange) throw std::logic_error("Should not be able to do that");
	cout << "successful" << endl;
}

// test sets and gets with an index
void popOutOfRangeTestSub(const hudson_ms::PopulationStructure& pop, int index)
{
	cout << "\nTesting accessor index ranges, index = " << index << endl;
	
	int nsps = static_cast<int>(pop.nSubPops());
	
	bool outOfRange = (index < 0 || index >= nsps );
	if (outOfRange) cout << "invalid_argument exceptions are expected here"<< endl;
	else cout << "No invalid_argument exceptions are expected here"<< endl;
	
	// tests
	try {
		cout << "popCopy.getSubNsamAt(index)" << endl;
		hudson_ms::PopulationStructure popCopy(pop);
		popCopy.getSubNsamAt(index);
		throwLogicErrorCheck(outOfRange);
		
	}
	catch (std::invalid_argument& ia) {
		cout << "std::invalid_argument :\n" << ia.what() << endl;
		if (!outOfRange) throw;
	}
	try {
		cout << "popCopy.getSizeRatioAt" << endl;
		hudson_ms::PopulationStructure popCopy(pop);
		popCopy.getSizeRatioAt(index);
		
		throwLogicErrorCheck(outOfRange);
	}
	catch (std::invalid_argument& ia) {
		cout << "std::invalid_argument :\n" << ia.what() << endl;
		if (!outOfRange) throw;
	}
	try {
		cout << "popCopy.getGrowthRateAt" << endl;
		hudson_ms::PopulationStructure popCopy(pop);
		popCopy.getGrowthRateAt(index);
		
		throwLogicErrorCheck(outOfRange);
	}
	catch (std::invalid_argument& ia) {
		cout << "std::invalid_argument :\n" << ia.what() << endl;
		if (!outOfRange) throw;
	}
	try {
		cout << "popCopy.getLabelAt" << endl;
		hudson_ms::PopulationStructure popCopy(pop);
		popCopy.getLabelAt(index);
		
		throwLogicErrorCheck(outOfRange);
		
	}
	catch (std::invalid_argument& ia) {
		cout << "std::invalid_argument :\n" << ia.what() << endl;
		if (!outOfRange) throw;
	}
	try {
		cout << "popCopy.setGrowthRateAt" << endl;
		hudson_ms::PopulationStructure popCopy(pop);
		popCopy.setGrowthRateAt(index, 8.0);
		
		throwLogicErrorCheck(outOfRange);
	}
	catch (std::invalid_argument& ia) {
		cout << "std::invalid_argument :\n" << ia.what() << endl;
		if (!outOfRange) throw;
	}
	try {
		cout << "popCopy.erase" << endl;
		hudson_ms::PopulationStructure popCopy(pop);
		popCopy.erase(index);
		
		throwLogicErrorCheck(outOfRange);
	}
	catch (std::invalid_argument& ia) {
		cout << "std::invalid_argument :\n" << ia.what() << endl;
		if (!outOfRange) throw;
	}
	try {
		cout << "popCopy.relabel" << endl;
		hudson_ms::PopulationStructure popCopy(pop);
		popCopy.relabel(index, "NewLabel");
		
		throwLogicErrorCheck(outOfRange);
	}
	catch (std::invalid_argument& ia) {
		cout << "std::invalid_argument :\n" << ia.what() << endl;
		if (!outOfRange) throw;
	}
	if (nsps > 0 ) {
		int order = ( (index < 1 ? nsps-1 : index-1) )%nsps;
		double m = 5.0;
		try {
			cout << "popCopy.getMigration(" << index << ", " << order << ")" << endl;
			hudson_ms::PopulationStructure popCopy(pop);
			popCopy.getMigration(index, order);
			throwLogicErrorCheck(outOfRange);
		}
		catch (std::invalid_argument& ia) {
			cout << "std::invalid_argument :\n" << ia.what() << endl;
			if (!outOfRange) throw;
		}
		try {
			cout << "popCopy.setMigration(" << index << ", " << order << ", " << m << ")" << endl;
			hudson_ms::PopulationStructure popCopy(pop);
			popCopy.setMigration(index, order, m);
			throwLogicErrorCheck(outOfRange);
		}
		catch (std::invalid_argument& ia) {
			cout << "std::invalid_argument :\n" << ia.what() << endl;
			if (!outOfRange) throw;
		}
		try {
			cout << "popCopy.getMigration(" << order << ", " << index << ")" << endl;
			hudson_ms::PopulationStructure popCopy(pop);
			popCopy.getMigration(order, index);
			throwLogicErrorCheck(outOfRange);
		}
		catch (std::invalid_argument& ia) {
			cout << "std::invalid_argument :\n" << ia.what() << endl;
			if (!outOfRange) throw;
		}
		try {
			cout << "popCopy.setMigration(" << order << ", " << index << ", " << m << ")" << endl;
			hudson_ms::PopulationStructure popCopy(pop);
			popCopy.setMigration(order, index, m);
			throwLogicErrorCheck(outOfRange);
		}
		catch (std::invalid_argument& ia) {
			cout << "std::invalid_argument :\n" << ia.what() << endl;
			if (!outOfRange) throw;
		}
	}
	
	// all other exceptions thrown
}

// test sets and gets with a label
void popNotRecognisedTestSub(const hudson_ms::PopulationStructure& pop, 
							const std::string& label)
{
	cout << "\nTesting accessors with label, label = \"" << label << "\""<< endl;
	
	int nsps = static_cast<int>(pop.nSubPops());
	
	bool outOfRange = !(pop.hasLabel(label));
	if (outOfRange) cout << "invalid_argument exceptions are expected here"<< endl;
	else cout << "No invalid_argument exceptions are expected here"<< endl;
	
	// tests
	try {
		cout << "popCopy.getSubNsamAt(label)" << endl;
		hudson_ms::PopulationStructure popCopy(pop);
		popCopy.getSubNsamAt(label);
		throwLogicErrorCheck(outOfRange);
		
	}
	catch (std::invalid_argument& ia) {
		cout << "std::invalid_argument :\n" << ia.what() << endl;
		if (!outOfRange) throw;
	}
	try {
		cout << "popCopy.getSizeRatioAt" << endl;
		hudson_ms::PopulationStructure popCopy(pop);
		popCopy.getSizeRatioAt(label);
		
		throwLogicErrorCheck(outOfRange);
	}
	catch (std::invalid_argument& ia) {
		cout << "std::invalid_argument :\n" << ia.what() << endl;
		if (!outOfRange) throw;
	}
	try {
		cout << "popCopy.getGrowthRateAt" << endl;
		hudson_ms::PopulationStructure popCopy(pop);
		popCopy.getGrowthRateAt(label);
		
		throwLogicErrorCheck(outOfRange);
	}
	catch (std::invalid_argument& ia) {
		cout << "std::invalid_argument :\n" << ia.what() << endl;
		if (!outOfRange) throw;
	}
	
	
	try {
		cout << "popCopy.setGrowthRateAt" << endl;
		hudson_ms::PopulationStructure popCopy(pop);
		popCopy.setGrowthRateAt(label, 8.0);
		
		throwLogicErrorCheck(outOfRange);
	}
	catch (std::invalid_argument& ia) {
		cout << "std::invalid_argument :\n" << ia.what() << endl;
		if (!outOfRange) throw;
	}
	try {
		cout << "popCopy.erase" << endl;
		hudson_ms::PopulationStructure popCopy(pop);
		popCopy.erase(label);
		
		throwLogicErrorCheck(outOfRange);
	}
	catch (std::invalid_argument& ia) {
		cout << "std::invalid_argument :\n" << ia.what() << endl;
		if (!outOfRange) throw;
	}
	try {
		cout << "popCopy.relabel" << endl;
		hudson_ms::PopulationStructure popCopy(pop);
		popCopy.relabel(label, "NewLabel");
		
		throwLogicErrorCheck(outOfRange);
	}
	catch (std::invalid_argument& ia) {
		cout << "std::invalid_argument :\n" << ia.what() << endl;
		if (!outOfRange) throw;
	}
	if (nsps > 0 ) {
		std::string label2;
		for (int order = 0; order < nsps; ++order) {
			label2 = pop.getLabelAt(order);
			if (label != label2) {
				break;
			}
		}
		if (label2.empty()) throw std::logic_error("Something went wrong finding another label");
		
		double m = 5.0;
		try {
			cout << "popCopy.getMigration(\"" << label << "\", \"" << label2 << "\")" << endl;
			hudson_ms::PopulationStructure popCopy(pop);
			popCopy.getMigration(label, label2);
			throwLogicErrorCheck(outOfRange);
		}
		catch (std::invalid_argument& ia) {
			cout << "std::invalid_argument :\n" << ia.what() << endl;
			if (!outOfRange) throw;
		}
		try {
			cout << "popCopy.setMigration(\"" << label << "\", \"" << label2 << "\", " << m << ")" << endl;
			hudson_ms::PopulationStructure popCopy(pop);
			popCopy.setMigration(label, label2, m);
			throwLogicErrorCheck(outOfRange);
		}
		catch (std::invalid_argument& ia) {
			cout << "std::invalid_argument :\n" << ia.what() << endl;
			if (!outOfRange) throw;
		}
		try {
			cout << "popCopy.getMigration(\"" << label2 << "\", \"" << label << "\")" << endl;
			hudson_ms::PopulationStructure popCopy(pop);
			popCopy.getMigration(label2, label);
			throwLogicErrorCheck(outOfRange);
		}
		catch (std::invalid_argument& ia) {
			cout << "std::invalid_argument :\n" << ia.what() << endl;
			if (!outOfRange) throw;
		}
		try {
			cout << "popCopy.setMigration(\"" << label2 << "\", \"" << label << "\", " << m << ")" << endl;
			hudson_ms::PopulationStructure popCopy(pop);
			popCopy.setMigration(label2, label, m);
			throwLogicErrorCheck(outOfRange);
		}
		catch (std::invalid_argument& ia) {
			cout << "std::invalid_argument :\n" << ia.what() << endl;
			if (!outOfRange) throw;
		}
	}
	
	// all other exceptions thrown
}

void testSubPopsInvalidArgSub(int nsam, double sr, double gr)
{
	bool illegalNsam = (nsam < 0);
	bool illegalSr = (sr <= 0.0);
	
	std::cout << "\nTest PopulationStructure(...) constructors and PopulationStructure::add(...) with nsam = " << nsam << ", size ratio = " << sr << " and growth rate = " << gr << endl;
	std::cout << "nsam = " << nsam;
	if (illegalNsam) std::cout  << " is illegal ";
	else std::cout << " is legal ";
	std::cout << " and size ratio = " << sr;
	if (illegalSr) std::cout  << " is illegal" << endl;
	else std::cout << " is legal" << endl;
	
	try {	
		cout << "Try constructor pop(_sn = " << nsam << ")" << endl;
		hudson_ms::PopulationStructure pop(nsam);
		throwLogicErrorCheck(illegalNsam);
	}
	catch (std::invalid_argument& ia) {
		cout << "std::invalid_argument :\n" << ia.what() << endl;
		if (!illegalNsam) throw;
	}
	
	try {	
		cout << "Make empty pop and try add(_sn = " << nsam << ")" << endl;
		hudson_ms::PopulationStructure pop;
		pop.add(nsam);
		throwLogicErrorCheck(illegalNsam);
	}
	catch (std::invalid_argument& ia) {
		cout << "std::invalid_argument :\n" << ia.what() << endl;
		if (!illegalNsam) throw;
	}
	
	try {	
		cout << "Make empty pop and try add(_sn = " << nsam << ", _gr = " << gr << ")" << endl;
		hudson_ms::PopulationStructure pop;
		pop.add(nsam, gr);
		throwLogicErrorCheck(illegalNsam);
	}
	catch (std::invalid_argument& ia) {
		cout << "std::invalid_argument :\n" << ia.what() << endl;
		if (!illegalNsam) throw;
	}
	
	try {	
		cout << "Make empty pop and try add(_sn = " << nsam << ", _sr = " << sr << ", _gr = " << gr << ")" << endl;
		hudson_ms::PopulationStructure pop;
		pop.add(nsam, sr, gr);
		throwLogicErrorCheck(illegalNsam || illegalSr);
	}
	catch (std::invalid_argument& ia) {
		cout << "std::invalid_argument :\n" << ia.what() << endl;
		if (!illegalNsam && !illegalSr) throw;
	}
	
	if (illegalNsam && illegalSr) {
		try {
			int sr2 = 1;	
			cout << "Just illegal nsam: make empty pop and try add(_sn = " << nsam << ", _sr2 = " << sr2 << ", _gr = " << gr << ")" << endl;
			hudson_ms::PopulationStructure pop;
			pop.add(nsam, sr2, gr);
			throwLogicErrorCheck(illegalNsam);
		}
		catch (std::invalid_argument& ia) {
			cout << "std::invalid_argument :\n" << ia.what() << endl;
			if (!illegalNsam) throw;
		}
		try {
			int nsam2 = 1;	
			cout << "Just illegal size ratio: make empty pop and try add(_sn = " << nsam2 << ", _sr = " << sr << ", _gr = " << gr << ")" << endl;
			hudson_ms::PopulationStructure pop;
			pop.add(nsam2, sr, gr);
			throwLogicErrorCheck(illegalSr);
		}
		catch (std::invalid_argument& ia) {
			cout << "std::invalid_argument :\n" << ia.what() << endl;
			if (!illegalSr) throw;
		}
	}
}

void testSubPopsSetInvalidArgSub(const hudson_ms::PopulationStructure& pop,
									int nsam, double sr, double gr)
{
	int nsps = static_cast <int> (pop.nSubPops());
	
	if (nsps > 1) {
		bool illegalNsam = (nsam < 0);
		bool illegalSr = (sr <= 0.0);
		bool illegalGr = false;
		
		std::cout << "\nTest PopulationStructure::SetSubNsamAt(...) and setSizeRatioAt and setGrowthRateAt with nsam = " << nsam << ", size ratio = " << sr << " and growth rate = " << gr << endl;
		std::cout << "nsam = " << nsam;
		if (illegalNsam) std::cout  << " is illegal ";
		else std::cout << " is legal ";
		std::cout << " and size ratio = " << sr;
		if (illegalSr) std::cout  << " is illegal" << endl;
		else std::cout << " is legal" << endl;
		std::cout << "growth rate is legal" << endl;
		
		int index = 2999%nsps;
				
		std::string label = pop.getLabelAt(index);
		
		try {	
			cout << "Try setSubNsamAt(" << index << ", _sn = " << nsam << ")" << endl;
			hudson_ms::PopulationStructure popCopy(pop);
			popCopy.setSubNsamAt(index, nsam);
			throwLogicErrorCheck(illegalNsam);
		}
		catch (std::invalid_argument& ia) {
			cout << "std::invalid_argument :\n" << ia.what() << endl;
			if (!illegalNsam) throw;
		}
		try {	
			cout << "Try setSubNsamAt(\"" << index << "\", _sn = " << nsam << ")" << endl;
			hudson_ms::PopulationStructure popCopy(pop);
			popCopy.setSubNsamAt(label, nsam);
			throwLogicErrorCheck(illegalNsam);
		}
		catch (std::invalid_argument& ia) {
			cout << "std::invalid_argument :\n" << ia.what() << endl;
			if (!illegalNsam) throw;
		}
		
		try {	
			cout << "Try setSizeRatioAt(" << index << ", _sr = " << sr << ")" << endl;
			hudson_ms::PopulationStructure popCopy(pop);
			popCopy.setSizeRatioAt(index, sr);
			throwLogicErrorCheck(illegalSr);
		}
		catch (std::invalid_argument& ia) {
			cout << "std::invalid_argument :\n" << ia.what() << endl;
			if (!illegalSr) throw;
		}
		
		try {	
			cout << "Try setSizeRatioAt(\"" << label << "\", _sr = " << sr << ")" << endl;
			hudson_ms::PopulationStructure popCopy(pop);
			popCopy.setSizeRatioAt(label, sr);
			throwLogicErrorCheck(illegalSr);
		}
		catch (std::invalid_argument& ia) {
			cout << "std::invalid_argument :\n" << ia.what() << endl;
			if (!illegalSr) throw;
		}
		
		try {	
			cout << "Try setGrowthRateAt(" << index << ", _gr = " << gr << ")" << endl;
			hudson_ms::PopulationStructure popCopy(pop);
			popCopy.setGrowthRateAt(index, gr);
			throwLogicErrorCheck(illegalGr);
		}
		catch (std::invalid_argument& ia) {
			cout << "std::invalid_argument :\n" << ia.what() << endl;
			if (!illegalGr) throw;
		}
		
		try {	
			cout << "Try setGrowthRateAt(\"" << label << "\", _gr = " << gr << ")" << endl;
			hudson_ms::PopulationStructure popCopy(pop);
			popCopy.setGrowthRateAt(label, gr);
			throwLogicErrorCheck(illegalGr);
		}
		catch (std::invalid_argument& ia) {
			cout << "std::invalid_argument :\n" << ia.what() << endl;
			if (!illegalGr) throw;
		}
	}
}


void testPopRelabellingIndividualInvalidArgSub(const hudson_ms::PopulationStructure& pop,
						int index1, string oldLabel, string newLabel, bool duplicate)
{
	
	hudson_ms::PopulationStructure popCopy(pop);
	
	cout << "Try to relabel index " << index1 << " with label \"" << oldLabel << "\" to \"" << newLabel << "\"" << endl;
	
	if (duplicate) cout << "invalid_argument exceptions are expected here"<< endl;
	else cout << "No invalid_argument exceptions are expected here"<< endl;

	try {
		hudson_ms::PopulationStructure popCopy(pop);
				
		cout << "Using relabel(" << index1 << ", \"" << newLabel << "\")" << endl;
	
		popCopy.relabel(index1, newLabel);
					
		throwLogicErrorCheck(duplicate);
	}
	catch (std::invalid_argument& ia) {
		cout << "std::invalid_argument :\n" << ia.what() << endl;
		if (!duplicate) throw;
	}
	
	try {	
		hudson_ms::PopulationStructure popCopy(pop);
			
		cout << "Using relabel(\"" << index1 << "\", \"" << newLabel << "\")" << endl;
		
		popCopy.relabel(oldLabel, newLabel);
					
		throwLogicErrorCheck(duplicate);
	}
	catch (std::invalid_argument& ia) {
		cout << "std::invalid_argument :\n" << ia.what() << endl;
		if (!duplicate) throw;
	}
}
	

void testPopRelabellingIndividualInvalidArg(const hudson_ms::PopulationStructure& pop,
										int index1, int index2)
{
	
	string oldLabel = pop.getLabelAt(index1);
	string newLabel = pop.getLabelAt(index2);
	
	bool duplicate = (index1 != index2);
	
	testPopRelabellingIndividualInvalidArgSub(pop, index1, oldLabel, newLabel, duplicate);
	
}
	

void testPopRelabellingInvalidArg(const hudson_ms::PopulationStructure& pop)
{
	
	int nsps = static_cast <int> (pop.nSubPops());
	
	if (nsps > 1) {
	
		{	
			int index1 = 16767%nsps;
			int index2 = ( (index1 < 1 ? nsps-1 : index1-1) )%nsps;
			
			testPopRelabellingIndividualInvalidArg(pop, index1, index2);
		}
		
	}
	if (nsps > 0) {
		
		{	
			int index1 = 1888%nsps;
			int index2 = index1;
			
			testPopRelabellingIndividualInvalidArg(pop, index1, index2);
		}
		
		bool duplicate = false;
		try {	
			hudson_ms::PopulationStructure popCopy(pop);
			
			vector < string > labels = popCopy.getLabels();
			
			random_shuffle(labels.begin(), labels.end());
				
			cout << "Using relabel(newlabels) where newlabels is " << endl;
			ostream_iterator<string> out_it (cout,", ");
			copy ( labels.begin(), labels.end(), out_it );
			cout << endl;

			duplicate = false; // a permuation is okay
			if (duplicate) cout << "an invalid_argument exceptions is expected here"<< endl;
			else cout << "No invalid_argument exceptions are expected here"<< endl;

			
			popCopy.relabel(labels);
						
			throwLogicErrorCheck(duplicate);
		}
		
		catch (std::invalid_argument& ia) {
			cout << "std::invalid_argument :\n" << ia.what() << endl;
			if (!duplicate) throw;
		}
		try {	
			hudson_ms::PopulationStructure popCopy(pop);
			
			vector < string > labels = popCopy.getLabels();
			
			random_shuffle(labels.begin(), labels.end());
			
			int index1 = 1999%nsps;
			int index2 = (index1 + 1)%nsps;
			labels.at(index2) = labels.at(index1);
			
			cout << "Using relabel(newlabels) where newlabels is " << endl;
			ostream_iterator<string> out_it (cout,", ");
			copy ( labels.begin(), labels.end(), out_it );
			cout << endl;

			duplicate = (index1 != index2); // but not a duplication
			if (duplicate) cout << "an invalid_argument exceptions is expected here"<< endl;
			else cout << "No invalid_argument exceptions are expected here"<< endl;

			
			popCopy.relabel(labels);
						
			throwLogicErrorCheck(duplicate);
		}
		
		catch (std::invalid_argument& ia) {
			cout << "std::invalid_argument :\n" << ia.what() << endl;
			if (!duplicate) throw;
		}
		
	}
	else cout << "pop is empty" << endl;
}

void testPopBeforeAndAfterEventsSub(const hudson_ms::PopulationStructure& pop)
{
	
	bool reorderBlocked = !(pop.eListEmpty());
	
	int nsps = static_cast<int> (pop.nSubPops());
	
	if (nsps > 0) {
	
		try {	
			
			hudson_ms::PopulationStructure popCopy(pop);
					
			cout << "Using resetOrdering when eListEmpty = " << popCopy.eListEmpty() << endl;
			
			vector < string > labels = popCopy.getLabels();
			
			popCopy.resetOrdering(labels); // even using the same labels can cause it to throw
						
			throwLogicErrorCheck(reorderBlocked);
		}
		catch (std::runtime_error& re) {
			cout << "std::runtime_error :\n" << re.what() << endl;
			if (!reorderBlocked) throw;
		}
		
		int orderErase = 2222%nsps;
					
		try {	
			
			hudson_ms::PopulationStructure popCopy(pop);
			
			cout << "Using erase(int orderErase = " << orderErase << ") when eListEmpty = " << popCopy.eListEmpty() << endl;
			
			popCopy.erase(orderErase);
						
			throwLogicErrorCheck(reorderBlocked);
		}
		catch (std::runtime_error& re) {
			cout << "std::runtime_error :\n" << re.what() << endl;
			if (!reorderBlocked) throw;
		}
		
		try {	
			int index = (orderErase+1)%nsps;
			
			hudson_ms::PopulationStructure popCopy(pop);
			
			string label = popCopy.getLabelAt(index);
			
			cout << "Using erase(string label = \"" << label << "\") when eListEmpty = " << popCopy.eListEmpty() << endl;
			
			popCopy.erase(label);
						
			throwLogicErrorCheck(reorderBlocked);
		}
		catch (std::runtime_error& re) {
			cout << "std::runtime_error :\n" << re.what() << endl;
			if (!reorderBlocked) throw;
		}
	}
	else cout << "pop is empty" << endl;	
	
}

void testPopBeforeAndAfterEvents(const hudson_ms::PopulationStructure& pop)
{
	
	
	hudson_ms::PopulationStructure popCopy(pop);
	
	if (!pop.empty() && pop.eListEmpty()) {
		
		cout << "Should be able to reorder and erase when eList is empty: no runtime_error exceptions expected" << endl; 
	
		testPopBeforeAndAfterEventsSub(popCopy);
				
		cout << "Add an event to the pop's event list" << endl; 
		
		std::vector< hudson_ms::HudsonMSDemogEventPtr > _eList;
	
		
		double t = 0.0;
		double gr = 0.0;
				
		boost::shared_ptr < const hudson_ms::HudsonMSDemogEvent > 
			ev(new hudson_ms::HudsonMSDemogEvent_G(t, gr));  
	
		_eList.push_back(ev);
		
		popCopy.setEventList(_eList);
	
		cout << "Should not be able to reorder and erase when eList is not empty: runtime_errors exceptions expected" << endl; 
		testPopBeforeAndAfterEventsSub(popCopy);
	
	}
	else if (pop.empty()) cout << "pop is empty" << endl;
	else cout << "pop.eListEmpty() is not true" << endl;

	
}
