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
\brief Basic testing of Simdata sampler.

No adaptive abc, just testing the sampler,
ie sampling with one set of parameter values.

The MultiLociSimDataSampler uses libsequence
to simulate snp data for given parameters.  

*/

#include "config.h"

#include "descriptive_stats.hpp"
#include "summary_statistic_set.hpp"

#include "multi_loci_simdata_sampler.hpp"
#include "multi_loci_polytable_sampler.hpp"
#include "multi_loci_polytable_set.hpp"

#include "parameter_particle.hpp"

#include "prng.hpp"
#include "prng_gsl.hpp"

#include "utilities.hpp"

#include <mscplusplus/population_structure.hpp>

#include <boost/lambda/bind.hpp>
#include <boost/lambda/lambda.hpp>


#include <boost/shared_ptr.hpp>


#include <iostream>
#include <sstream>
#include <vector>


using namespace std;
using namespace mct;

void testSub(
	const boost::shared_ptr< const hudson_ms::PopulationStructure > pop,
	const std::string& filenamepop,
	double theta, double growth, 
	size_t nloci, size_t n_examples, size_t nsites, 
	const vector < size_t>& n_loci_nsites);
	
void doSS(boost::shared_ptr < SummaryStatisticSet > ss_set, 
			const std::string& filename);


void simdataTesting();

int main()
{
	simdataTesting();
	
	
	return 0;
	
}

	

// basic test that simdata sampling is working
void simdataTesting()
{
	try {
		int myPid = getpid();
		std::cout << "basic test for simdata " << std::endl;
		std::cout << "This process id is " << myPid << std::endl;
		
		// number to generate with 'true' parameters
		size_t n_examples = 10;
		//size_t n_examples = 10000;
		
		size_t nsam = 10; // samples
				
		// loci in each multi-loci SimData
		size_t nloci = 5; 
		
		size_t nsites = 25000; // number of sites in each locus 
		
		//vector of nloci copies of nsites
		std::vector < size_t > nloci_nsites(nloci, nsites);
		
		size_t Nzero = 1000000; // effective popn size (per sub-pop)
		double unit = 1000000000.0; // 10^9
		
		// equivalent to 50/nsites in coalescent time units of 4Nzero
		
		//mutation rates
		//mutation rates PER SITE ie 4N0(mu/site)
		double mymus[] = {5.0}; //0.001,10.0,50.0,100.0
		//double mythetas[] = {0.001,10.0,50.0,100.0}; //0.001,10.0,50.0,100.0
		std::vector <double> mus (mymus, mymus + sizeof(mymus) / sizeof(double) );
		std::vector <double> thetas(mus.size());
		
		using namespace boost::lambda;
		//mults by 4*Nzero/unit
		transform(mus.begin(), mus.end(), thetas.begin(), _1*4*Nzero/unit);
		
		//growth rates
		double mygrowths[] = {50.0};
		//double mygrowths[] = {0.0,10.0,50.0,100.0};
		std::vector <double> growths (mygrowths, mygrowths + sizeof(mygrowths) / sizeof(double) );
		
		bool doneEmptyPtr = false;
		bool doneEmptyPop = false;
		
		for (std::vector < double>::iterator theta_it = thetas.begin(); 
				theta_it < thetas.end();
				++ theta_it) {
				
				double theta = *theta_it;
					
			for (std::vector < double>::iterator growth_it = growths.begin(); 
				growth_it < growths.end();
				++ growth_it) {
				
				double growth = *growth_it;
				
				if (!doneEmptyPtr) {
					boost::shared_ptr< hudson_ms::PopulationStructure > pop;
					
					
					std::string filenamepop = "PopPtrEmpty"; 
					
					testSub(pop, filenamepop, theta, growth, nloci, n_examples, nsites, nloci_nsites);
					doneEmptyPtr = true;
				}
				if (!doneEmptyPop) {
					boost::shared_ptr< hudson_ms::PopulationStructure > 
								pop(new hudson_ms::PopulationStructure() );
					
					
					std::string filenamepop = "PopEmpty"; 
					
					testSub(pop, filenamepop, theta, growth, nloci, n_examples, nsites, nloci_nsites);
					doneEmptyPop = true;
				}
		
				//population structure
				boost::shared_ptr< hudson_ms::PopulationStructure > pop
				(new hudson_ms::PopulationStructure(nsam));
				
				std::string filenamepop = "TestSimdata";
		
				testSub(pop, filenamepop, theta, growth, 
					nloci, n_examples, nsites, nloci_nsites);
					
			}
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



void testSub(
	const boost::shared_ptr< const hudson_ms::PopulationStructure > pop,
	const std::string& filenamepop,
	double theta, double growth, 
	size_t nloci, size_t n_examples, size_t nsites, 
	const vector < size_t>& nloci_nsites)
{
	
	int seed = 2345;
	boost::shared_ptr < PRNGen > r ( new PRNGenGSL(seed) );
	boost::shared_ptr < PRNGen > r_sampler = r->spawnAnother();
		
	/*
	 * use this prng to spawn another random number generator.
	 * This ensures that if we increase the number of models the 
	 * samples will have the same seed.
	 * */
	
	cout << "\n\nnsites = " << nsites << ", nloci = " << nloci << ", n_examples = " << n_examples << endl;
	
	cout << "theta " << theta << " growth " << growth << endl;

	bool noPop = false;
	
	if (pop) {
		//population structure
		std::string popStr = pop->toString();
		std::cout << "Pop structure is " << popStr << endl;
							
		cout << "pop->totalNsam() = " << pop->totalNsam() << endl;
		if (pop->empty()) noPop = true;
	}
	else {
		cout << "pop ptr is empty " << endl;
		noPop = true;
	}
	
	cout << "Making sampler" << endl;
	
	/*make the sampler*/
	boost::shared_ptr < MultiLociPolyTableSampler > 
						sampler ( new MultiLociSimDataSampler(r_sampler) );
		
	double realvalues[] = {theta,growth}; 
	
	boost::shared_ptr < aabc::ParameterParticle > 
			p( new aabc::ParameterParticle(realvalues, 2) );

	
	cout << "Sampling" << endl;
	
	try {
		
		boost::shared_ptr < mct::MultiLociPolyTableSet > 
			refset_ptr = sampler->sample(n_examples,
			pop, nloci_nsites, p);
		
		if (noPop) throw std::logic_error("Should not be able to do that");
		
		boost::shared_ptr < SummaryStatisticSet > ss_set = 
				refset_ptr->getSummaryStatisticSet();
		
		boost::shared_ptr < SummaryStatisticSet > all_ss_set = 
				refset_ptr->getAllAvailableStatisticsSet();
	
		std::ostringstream stm;
		stm << "_ns_" << pop->totalNsam() << "_nl_" << nloci << "_n_" << n_examples;
		stm.precision(3);
		stm <<   fixed;
		stm << "_th_" << theta;
		stm << "_g_" << growth << ".txt";
		
		std::string path("../output/");

		std::string filename1 = path + filenamepop + std::string("SummaryStats") + stm.str();
		std::string filename2 = path + filenamepop + std::string("AllStats") + stm.str();
	
		cout << "\nLooking at the 'official' summary stats" << endl;
		doSS(ss_set, filename1);
		cout << "\nLooking at all the available summary stats" << endl;
		doSS(all_ss_set,filename2);
	}
	
	catch (std::invalid_argument& ia) {
		cout << "Error :\n" << ia.what() << endl;
		if(!noPop) throw;
	} 
}

void doSS(boost::shared_ptr < SummaryStatisticSet > ss_set, 
			const std::string& filename)
{
	
	try {
		
		std::cout << "Printing summary statistic set" << std::endl;
		cout << ss_set->plainString() << endl;
		
		if (ss_set->empty()) std::cout << "ss_set is empty)" << std::endl;
		
	}
	catch (std::exception& e) {
		std::cout << "Failed printing stats values, error is:\n"
		<< std::string(e.what()) << std::endl;
		throw;
	}
	
	std::vector< std::vector< double > > stats;
	
	try {	
		stats = *(ss_set->getAllValues());
		
	}
	catch (std::exception& e) {
		std::cout << "Failed to get stats values, error is:\n"
		<< std::string(e.what()) << std::endl;
		throw;
	}
		
	try {
		
		std::cout << "All stats values from getAllValues() vectors are" << std::endl;
		for (std::vector< std::vector< double > >::iterator it = stats.begin();
			it < stats.end(); ++it) {
			copy(it->begin(), it->end(), ostream_iterator<double>(std::cout, "\t"));
			std::cout << std::endl;
		}
		
		if (stats.empty()) std::cout << "(stats values is empty)" << std::endl;
		
	}
	catch (std::exception& e) {
		std::cout << "Failed printing stats values, error is:\n"
		<< std::string(e.what()) << std::endl;
		throw;
	}
	
	try {
		std::cout << "Output to file " << filename << std::endl;
		ss_set->outputToFile(filename, false); // overwrite
		std::cout << "and read in again " << std::endl;
		cout << (mct_utilities::string_from_file(filename)) << endl;
		
	}
	catch (std::exception& e) {
		std::cout << "Failed outputting or reading stats, error is:\n"
		<< std::string(e.what()) << std::endl;
		throw;
	}
	
	bool expectFailDS = (ss_set->size() < 1);
	bool expectFailStandard = ((!ss_set->empty()) && (ss_set->size() < 2));
	
	try {
		std::cout << "\nTry to get descriptive stats from summary stats" << std::endl;
		
		for (size_t i = 0; i < stats.size() ; ++i) {
				
			cout << "Descriptive stats for index " << i << endl;
			DescriptiveStats ds(ss_set, i);
			std::cout << (ds.toString()) << std::endl;
		}
	
		if (expectFailDS) throw std::logic_error("Should not be able to do that");
		
	}
	catch (std::exception& e) {
		std::cerr << "Failed on descriptive stats from summary stats, error is:\n"
		<< std::string(e.what()) << std::endl;
		if (!expectFailDS) throw;
		
	}
	
	try {	
		std::cout << "\nGet standardised stats as a summary stats set" << std::endl;
		
		boost::shared_ptr < SummaryStatisticSet > sss = 
					ss_set->makeSummaryStatisticSetStandardised();
	
		
		std::cout << "standardised stats values are:" << std::endl;
		std::cout << sss->plainString() << std::endl;

		if (sss->empty())  std::cout << "standardised stats summary stat is empty:" << std::endl;
		
		if (expectFailStandard) throw std::logic_error("Should not be able to do that");
		
		try {
			std::cout << "\nGet descriptive stats from standardised summary stats" << std::endl;
			
			for (size_t i = 0; i < stats.size() ; ++i) {
					
				cout << "Descriptive stats for index " << i << endl;
				DescriptiveStats ds(sss, i);
				std::cout << (ds.toString()) << std::endl;
			}
		
			if (expectFailDS) throw std::logic_error("Should not be able to do that");
		}
		catch (std::exception& e) {
			std::cerr << "Failed on descriptive stats from standardised summary stats test, error is:\n"
			<< std::string(e.what()) << std::endl;
			if (!expectFailDS) throw;
			
		}
	}
	catch (std::exception& e) {
		std::cout << "Failed to get all standardised stats as a summary stats set error is:\n"
		<< std::string(e.what()) << std::endl;
		if (!expectFailStandard) throw;
	}
	
	
	
}
