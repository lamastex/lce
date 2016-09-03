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
\brief Timing test with Simdata sampler.

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

#include <mscplusplus/population_structure.hpp>

#include <boost/lambda/bind.hpp>
#include <boost/lambda/lambda.hpp>

#include <boost/shared_ptr.hpp>


#include <iostream>
#include <sstream>
#include <vector>
#include <time.h>


using namespace std;
using namespace mct;

void testSub(
	const boost::shared_ptr< const hudson_ms::PopulationStructure > pop,
	double theta, double growth, 
	size_t nloci, size_t n_examples, size_t nsites, 
	const vector < size_t>& nloci_nsites);

void simdataTimingTest();

int main()
{
	simdataTimingTest();
	
	
	return 0;
	
}

	

// timing test for simdata
void simdataTimingTest()
{
	try {
		int myPid = getpid();
		std::cout << "timing test for simdata" << std::endl;
		std::cout << "This process id is " << myPid << std::endl;
		
		// number to generate with 'true' parameters
		size_t n_examples = 1000;
		//size_t n_examples = 10000;
		
		size_t nsam = 10; // samples
				
		// loci in each multi-loci SimData
		size_t nloci = 25; 
		
		
		// number of sites in each locus, we'll have the same for each locus here
		size_t nsites = 25000; 
		
		//vector of nloci copies of nsites
		std::vector < size_t > nloci_nsites(nloci, nsites);
		
		
		size_t Nzero = 1000000; // effective popn size (per sub-pop)
		double unit = 1000000000.0; // 10^9
		
		//mutation rates
		//mutation rates PER SITE ie 4N0(mu/site)
		double mymus[] = {5}; //0.001,10.0,50.0,100.0
		
		//double mythetas[] = {0.001,10.0,50.0,100.0}; //0.001,10.0,50.0,100.0
		std::vector <double> mus (mymus, mymus + sizeof(mymus) / sizeof(double) );
		std::vector <double> thetas(mus.size());
		
		using namespace boost::lambda;
		//mults by 4*Nzero/unit
		transform(mus.begin(), mus.end(), thetas.begin(), _1*4*Nzero/unit);
		
		//growth rates
		double mygrowths[] = {30.0}; // 30.0
		//double mygrowths[] = {0.0,10.0,50.0,100.0};
		std::vector <double> growths (mygrowths, mygrowths + sizeof(mygrowths) / sizeof(double) );
		
		for (std::vector < double>::iterator theta_it = thetas.begin(); 
				theta_it < thetas.end();
				++ theta_it) {
				
			double theta = *theta_it;
					
			for (std::vector < double>::iterator growth_it = growths.begin(); 
				growth_it < growths.end();
				++ growth_it) {
				
				double growth = *growth_it;
			
				//population structure
				boost::shared_ptr< hudson_ms::PopulationStructure > pop
				(new hudson_ms::PopulationStructure(nsam));
				
				testSub(pop, theta, growth, 
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
	
	boost::shared_ptr < PRNGen > r_models = r->spawnAnother();
	
	
	cout << "\n\nnsites = " << nsites << ", nloci = " << nloci << ", n_examples = " << n_examples << endl;
	
	cout << "theta " << theta << " growth " << growth << endl;

	std::string popStr = pop->toString();
	std::cout << "Pop structure is " << popStr << endl;
	cout << "pop->totalNsam() = " << pop->totalNsam() << endl;
		
	cout << "\nMaking sampler" << endl;
	
	/*make the sampler*/
	boost::shared_ptr < MultiLociPolyTableSampler > 
						sampler ( new MultiLociSimDataSampler(r_sampler) );
	
	double realvalues[] = {theta,growth}; 
	
	boost::shared_ptr < aabc::ParameterParticle > 
			p( new aabc::ParameterParticle(realvalues, 2) );
	
	// time the actual sampling
	std::cout << "Starting to sample: " << endl;
	clock_t startTime = clock();
			
	boost::shared_ptr < mct::MultiLociPolyTableSet > 
		refset_ptr = sampler->sample(n_examples,
		pop, nloci_nsites, p);
	
	clock_t endTime = clock();
	clock_t timeTaken = endTime - startTime;
	
	std::cout << "\nFinished sampling, sampling time was "
		<< ((static_cast<double>(timeTaken)) / CLOCKS_PER_SEC)
		<< "s (" << timeTaken << " clocks)" <<std::endl;
	
	// time getting samples out
	std::cout << "\nGetting sample stats into SummaryStatisticSet form: " << endl;
	startTime = clock();
	
	boost::shared_ptr < SummaryStatisticSet > all_ss_set = 
				refset_ptr->getAllAvailableStatisticsSet();
	
	endTime = clock();
	timeTaken = endTime - startTime;
	
	std::cout << "finished getting SummaryStatisticSet, time was "
		<< ((static_cast<double>(timeTaken)) / CLOCKS_PER_SEC)
		<< "s (" << timeTaken << " clocks)" <<std::endl;
	
	std::cout << "\nTry to get descriptive stats from summary stats" << std::endl;
	
	startTime = clock();
	
	if (all_ss_set->size() > 0)	{
		size_t nstats = all_ss_set->at(0)->size();
		for (size_t i = 0; i < nstats ; ++i) {
				
			cout << "Descriptive stats for index " << i << endl;
			DescriptiveStats ds(all_ss_set, i);
			std::cout << (ds.toString()) << std::endl;
		}
	}
	endTime = clock();
	timeTaken = endTime - startTime;
	
	std::cout << "finished getting DescriptiveStats, time (including printing) was "
		<< ((static_cast<double>(timeTaken)) / CLOCKS_PER_SEC)
		<< "s (" << timeTaken << " clocks)" <<std::endl;
	
				
}

