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
\brief Small test of multiloci polytable for population comparisons.


*/

#include "config.h"

#include "descriptive_stats.hpp"
#include "summary_statistic_set.hpp"

#include "multi_loci_polysites_sampler.hpp"
#include "multi_loci_simdata_sampler.hpp"
#include "multi_loci_polytable_sampler.hpp"
#include "multi_loci_polytable_set.hpp"

#include "model_set_config_builder.hpp"
#include "multi_loci_modelset_set.hpp"
#include "modelset.hpp"

#include "parameter_particle.hpp"

#include "prng.hpp"
#include "prng_gsl.hpp"

//#include "utilities.hpp"

#include "make_test_pops.hpp"


#include <mscplusplus/population_structure.hpp>


#include <boost/shared_ptr.hpp>


#include <iostream>
#include <sstream>
#include <vector>
#include <iterator>



using namespace std;
using namespace mct;


void polytableTesting();

void testSub(
	const boost::shared_ptr< const hudson_ms::PopulationStructure > pop,
	const string& filenamepop,
	double theta, double growth, 
	size_t nloci, size_t n_examples, size_t nsites, 
	const vector < size_t>& n_loci_nsites);

void testPolysitesSub(
	const boost::shared_ptr< const hudson_ms::PopulationStructure > pop,
	bool noPop,
	double theta, double growth, 
	size_t nloci, size_t n_examples, size_t nsites, 
	const vector < size_t>& n_loci_nsites);

void testSimdataSub(
	const boost::shared_ptr< const hudson_ms::PopulationStructure > pop,
	bool noPop,
	double theta, double growth, 
	size_t n_examples,  
	const vector < size_t>& nloci_nsites);
	
void doSample(boost::shared_ptr < MultiLociPolyTableSampler > sampler,
	const boost::shared_ptr< const hudson_ms::PopulationStructure > pop,
	bool noPop,
	size_t n_examples, const vector < size_t>& nloci_nsites,
	boost::shared_ptr < aabc::ParameterParticle > p);
	
void doSS(boost::shared_ptr < SummaryStatisticSet > ss_set);

int main()
{
	polytableTesting();
	
	
	return 0;
	
}

	

void polytableTesting()
{
	try {
		int myPid = getpid();
		std::cout << "test MultiLociPolyTable with subpopulation comparisons " << std::endl;
		std::cout << "This process id is " << myPid << std::endl;
		
		// number to generate with 'true' parameters
		size_t n_examples = 1;
		
		//size_t nsam = 10; // samples
		
		// loci in each multi-loci SimData
		size_t nloci = 2; 
		
		// number of sites in each locus, we'll have the same for each locus here
		size_t nsites = 5;
		
		//vector of nloci copies of nsites
		std::vector < size_t > nloci_nsites(nloci, nsites);
		
		size_t Nzero = 1000000; // effective popn size (per sub-pop)
		double unit = 1000000000.0; // 10^9
		
		
		double mu_per_site = 100.0/unit;	
		double theta = mu_per_site*4*Nzero;
		double growth = 10.0;
			
		cout << "\n\nnsites = " << nsites << ", nloci = " << nloci << ", n_examples = " << n_examples << endl;
		cout << "theta " << theta << " growth " << growth << endl;

		// do one set for each population type
				
		
		// null pop pointer
		{
			boost::shared_ptr< hudson_ms::PopulationStructure > pop;
			
			std::string filenamepop = "NullPopPtr"; // alter this one
			
			testSub(pop, filenamepop, theta, growth, nloci, n_examples, nsites, nloci_nsites);
		}
		// empty pop 
		{
			boost::shared_ptr< hudson_ms::PopulationStructure > 
								pop(new  hudson_ms::PopulationStructure);
			
			std::string filenamepop = "EmptyPop"; // alter this one
			
			testSub(pop, filenamepop, theta, growth, nloci, n_examples, nsites, nloci_nsites);
		}
		// pop with one subpop with no individuals
		{
			boost::shared_ptr< hudson_ms::PopulationStructure > pop = makePopSingleZero();
			
			std::string filenamepop = "PopIslandSingleZero"; // alter this one
			
			testSub(pop, filenamepop, theta, growth, nloci, n_examples, nsites, nloci_nsites);
		}
		// one pop 
		{
			boost::shared_ptr< hudson_ms::PopulationStructure > 
						pop = makePopSingle();
			
			
			std::string filenamepop = "PopSingle"; // alter this one
			
			testSub(pop, filenamepop, theta, growth, nloci, n_examples, nsites, nloci_nsites);
		}
		// pop with two subpops, one with no individuals
		{
			boost::shared_ptr< hudson_ms::PopulationStructure > pop = makePopIslandOneZero(Nzero);
			
			std::string filenamepop = "PopIslandOneZero"; // alter this one
			
			testSub(pop, filenamepop, theta, growth, nloci, n_examples, nsites, nloci_nsites);
		}
		// pop with 3 subpops, one with no individuals
		{
			boost::shared_ptr< hudson_ms::PopulationStructure > pop = makePopIslandThreeOneZero(Nzero);
			
			std::string filenamepop = "PopIslandThreeOneZero"; // alter this one
			
			testSub(pop, filenamepop, theta, growth, nloci, n_examples, nsites, nloci_nsites);
		}
		
		// pop with 3 subpops, one individual each
		{
			boost::shared_ptr< hudson_ms::PopulationStructure > pop = makePopIslandThreeSmallest(Nzero);
			
			std::string filenamepop = "PopIslandThreeSmallest"; // alter this one
			
			testSub(pop, filenamepop, theta, growth, nloci, n_examples, nsites, nloci_nsites);
		}
	
		// pop with 3 subpops, different numbers of individuals in each
		{
			boost::shared_ptr< hudson_ms::PopulationStructure > pop = makePopIslandPTstats(Nzero);
			
			std::string filenamepop = "PopIslandPolyTableStats"; // alter this one
			
			testSub(pop, filenamepop, theta, growth, nloci, n_examples, nsites, nloci_nsites);
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
	const string& filenamepop,
	double theta, double growth, 
	size_t nloci, size_t n_examples, size_t nsites, 
	const vector < size_t>& nloci_nsites)
{
	bool noPop = false;
	
	if (pop) {
		//population structure
		cout << "\n\n" << filenamepop << endl;
		std::string popStr = pop->toString();
		std::cout << "Pop structure is " << popStr << endl;
		if (pop->empty()) noPop = true;
		else if (!pop->totalNsam()) noPop = true;
	}
	else {
		cout << "pop ptr is empty " << endl;
		noPop = true;
	}
	
	testPolysitesSub(pop, noPop, theta, growth,
						nloci, n_examples, nsites, 
						nloci_nsites);

	testSimdataSub(pop, noPop, theta, growth, 
			n_examples, nloci_nsites);

}

void testPolysitesSub(
	const boost::shared_ptr< const hudson_ms::PopulationStructure > pop,
	bool noPop,
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
	
	cout << "\n\nTesting with Polysites sampling " << endl;

	/*use a ModelSetConfigBuilder to build the model data for one locus*/
	mct::ModelSetConfigBuilder configBuilder;
	configBuilder.addModelConfig("JC69", nsites);
	// create our models, one for each site, with default thetas
	boost::shared_ptr < mct::ModelSet > 
		mset_ptr	( 
			new mct::ModelSet(configBuilder.getModelSetConfigs()
									, r_models)
					);
	
	/* Make a multilocus model set with nloci copies of the same modelset*/
	boost::shared_ptr < mct::MultiLociModelSetSet > 
		models_ptr	( 
			new mct::MultiLociModelSetSet(mset_ptr, nloci) );
	
	
	// configure sampler with multilocus model set
	boost::shared_ptr < MultiLociPolyTableSampler > 
				sampler ( new MultiLociPolySitesSampler(models_ptr, r_sampler) );
	
	double realvalues[] = {theta,growth}; 
	
	boost::shared_ptr < aabc::ParameterParticle > 
			p( new aabc::ParameterParticle(realvalues, 2) );

	doSample(sampler, pop, noPop, n_examples, nloci_nsites, p);
}

void testSimdataSub(
	const boost::shared_ptr< const hudson_ms::PopulationStructure > pop,
	bool noPop,
	double theta, double growth, 
	size_t n_examples,  
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
	
	cout << "\n\nTesting with Simdata sampling " << endl;

	/*make the sampler*/
	boost::shared_ptr < MultiLociPolyTableSampler > 
						sampler ( new MultiLociSimDataSampler(r_sampler) );
	double realvalues[] = {theta,growth}; 
	
	boost::shared_ptr < aabc::ParameterParticle > 
			p( new aabc::ParameterParticle(realvalues, 2) );

	doSample(sampler, pop, noPop, n_examples, nloci_nsites, p);
}

void doSample(boost::shared_ptr < MultiLociPolyTableSampler > sampler,
	const boost::shared_ptr< const hudson_ms::PopulationStructure > pop,
	bool noPop,
	size_t n_examples, const vector < size_t>& nloci_nsites,
	boost::shared_ptr < aabc::ParameterParticle > p)
{	
	cout << "Sampling\n" << endl;
	
	try {
		
		boost::shared_ptr < mct::MultiLociPolyTableSet > 
			refset_ptr = sampler->sample(n_examples,
			pop, nloci_nsites, p);
		
		if (noPop) throw std::logic_error("Should not be able to do that");
		
		cout << "Have sample, size "; 
		cout << refset_ptr->size() << endl;
		
		//for each element in the sample
		for (size_t i = 0; i < refset_ptr->size(); ++i) {
			
			boost::shared_ptr < const MultiLociPolyTable > mlpt = refset_ptr->at(i);
			
			// print out the results
			mlpt->printSubpopData(cout);
			
			boost::shared_ptr < const SummaryStatistic  > withinPopAvs 
								= mlpt->withinPopAverageDifferences();
			
			std::cout << "\nwithin pop average differences are " <<std::endl;
			cout << withinPopAvs->plainString() << endl;
			
			boost::shared_ptr < const SummaryStatistic  > betweenPopAvs 
								= mlpt->betweenPopAverageDifferences();
			std::cout << "\nbetween pop average differences are " <<std::endl;
			cout << betweenPopAvs->plainString() << endl;
			
			boost::shared_ptr < const SummaryStatistic  > correctedAvs 
								= mlpt->correctedBetweenPopAverageDifferences();
			std::cout << "\ncorrected between pop average differences are " <<std::endl;
			cout << correctedAvs->plainString() << endl;
			
			boost::shared_ptr < const SummaryStatistic  > pairwiseAvs 
								= mlpt->pairwiseAverageDifferences();
			std::cout << "\npairwise pop average differences are " <<std::endl;
			cout << pairwiseAvs->plainString() << endl;
			
			boost::shared_ptr < const SummaryStatistic  > Fst 
													= mlpt->fst();
			cout << "\nFst is " << Fst->plainString() << endl;
			
			boost::shared_ptr < const SummaryStatistic  > pairwiseFsts 
										= mlpt->pairwiseFst();
			std::cout << "\npairwiseFsts are " <<std::endl;
			cout << pairwiseFsts->plainString() << endl;

		}
		
		#if(0)
			boost::shared_ptr < SummaryStatisticSet > all_ss_set = 
					refset_ptr->getAllAvailableStatisticsSet();
		
			cout << "\nLooking at all the available summary stats" << endl;
			doSS(all_ss_set);
		#endif
	}
	
	catch (std::invalid_argument& ia) {
		cout << "Error :\n" << ia.what() << endl;
		if(!noPop) throw;
	} 
}

void doSS(boost::shared_ptr < SummaryStatisticSet > ss_set)
{
	
	if (!ss_set->empty()) {
	
		size_t nstats = ss_set->at(0)->size();
		for (size_t i = 0; i < nstats ; ++i) {
				
			cout << "Descriptive stats for index " << i << endl;
			DescriptiveStats ds(ss_set, i);
			std::cout << (ds.toString()) << std::endl;
		}
	
		
	}
	
	
}
