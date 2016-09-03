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
\brief Small test of polysites sampler.

No adaptive abc, just testing the sampler,
ie sampling with one set of parameter values.

The MultiLociPolySitesSampler uses our MCT models to
generate sequences (which it can then convert to
snp representation).  
* 
Run this test with debugging on in the sampler to see the 
samples and converstion to the libsequence SNP structures.

*/

#include "config.h"

#include "descriptive_stats.hpp"
#include "summary_statistic_set.hpp"

#include "multi_loci_polysites_sampler.hpp"
#include "multi_loci_polytable_sampler.hpp"
#include "multi_loci_polytable_set.hpp"

#include "model_set_config_builder.hpp"
#include "multi_loci_modelset_set.hpp"
#include "modelset.hpp"

#include "parameter_particle.hpp"

#include "prng.hpp"
#include "prng_gsl.hpp"

#include "utilities.hpp"

#include <mscplusplus/population_structure.hpp>

#include <boost/shared_ptr.hpp>


#include <iostream>
#include <sstream>
#include <vector>




using namespace std;
using namespace mct;


void polysitesTesting();

void testSub(
	const boost::shared_ptr< const hudson_ms::PopulationStructure > pop,
	double theta, double growth, 
	size_t nloci, size_t n_examples, size_t nsites, 
	const vector < size_t>& n_loci_nsites);
	
void doSS(boost::shared_ptr < SummaryStatisticSet > ss_set);

int main()
{
	polysitesTesting();
	
	
	return 0;
	
}

	

// small test that polysites sampling is working
// run with debugging on in the polysites sampler
void polysitesTesting()
{
	try {
		int myPid = getpid();
		std::cout << "basic testing  for polysites " << std::endl;
		std::cout << "This process id is " << myPid << std::endl;
		
		// number to generate with 'true' parameters
		size_t n_examples = 1;
		
		size_t nsam = 10; // samples
		
		// loci in each multi-loci SimData
		size_t nloci = 1; 
		
		// number of sites in each locus, we'll have the same for each locus here
		size_t nsites = 10;
		
		//vector of nloci copies of nsites
		std::vector < size_t > nloci_nsites(nloci, nsites);
		
		size_t Nzero = 1000000; // effective popn size (per sub-pop)
		double unit = 1000000000.0; // 10^9
		
		
		double mu_per_site = 200.0/unit;	
		double theta = mu_per_site*4*Nzero;
		double growth = 50.0;
			
		//population structure
		boost::shared_ptr< hudson_ms::PopulationStructure > pop
		(new hudson_ms::PopulationStructure(nsam));
		
		testSub(pop, theta, growth, 
			nloci, n_examples, nsites, nloci_nsites);
					
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

	
	cout << "Sampling" << endl;
	
	try {
		
		boost::shared_ptr < mct::MultiLociPolyTableSet > 
			refset_ptr = sampler->sample(n_examples,
			pop, nloci_nsites, p);
		
		if (noPop) throw std::logic_error("Should not be able to do that");
		
		boost::shared_ptr < SummaryStatisticSet > all_ss_set = 
				refset_ptr->getAllAvailableStatisticsSet();
	
		cout << "\nLooking at all the available summary stats" << endl;
		doSS(all_ss_set);
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
