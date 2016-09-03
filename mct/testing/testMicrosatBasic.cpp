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
\brief Basic testing of microsat sampler.

No adaptive abc, just testing the sampler,
ie sampling with one set of parameter values.

The MultiLociMicrosatSampler uses our MCT models to
generate microsate repeat data.  

*/

#include "config.h"

#include "descriptive_stats.hpp"
#include "summary_statistic_set.hpp"

#include "multi_loci_microsat_sampler.hpp"
#include "multi_loci_microsat_set.hpp"

#include "model_set_config_builder.hpp"
#include "modelset.hpp"

#include "parameter_particle.hpp"

#include "prng.hpp"
#include "prng_gsl.hpp"

#include "utilities.hpp"

#include "make_test_pops.hpp"

#include <mscplusplus/population_structure.hpp>

#include <boost/lambda/bind.hpp>
#include <boost/lambda/lambda.hpp>


#include <boost/shared_ptr.hpp>


#include <iostream>
#include <sstream>
#include <vector>


using namespace std;
using namespace mct;


void microsatTesting();

void testSub(
	const boost::shared_ptr< const hudson_ms::PopulationStructure > pop,
	const std::string& filenamepop,
	double theta, double growth, 
	size_t nloci, size_t n_examples, size_t _nsites, 
	const vector < size_t>& nsites);


int main()
{
	microsatTesting();
	
	
	return 0;
	
}

	
// basic test that microsats sampling is working
void microsatTesting()
{
	try {
		int myPid = getpid();
		std::cout << "basic testing  for microsats " << std::endl;
		std::cout << "This process id is " << myPid << std::endl;
		
		// loci in each multi-loci SimData
		size_t myloci[] = {1, 5}; 
		std::vector <size_t> nls (myloci, myloci + sizeof(myloci) / sizeof(size_t) );
		
		size_t mysites[] = {1, 3}; // number of sites in each locus 
		std::vector <size_t> nss (mysites, mysites + sizeof(mysites) / sizeof(size_t) );
		
		size_t Nzero = 1000000; // effective popn size (per sub-pop)
		double unit = 100000.0; // 10^5
		
		//mutation rates
		//mutation rates PER SITE ie 4N0(mu/site)
		double mymus[] = {1.0}; //0.001,10.0,50.0,100.0
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
		
		// number to generate with 'true' parameters
		//size_t n_examples = 1000;
		size_t n_examples = 1;
		
		bool doneEmptyPtr = false;
		bool doneEmptyPop = false;
				
		for (std::vector < size_t>::iterator ns_it = nss.begin(); 
				ns_it < nss.end();
				++ ns_it) {
				
				size_t nsites = *ns_it; 
					
			for (std::vector < size_t>::iterator nl_it = nls.begin(); 
				nl_it < nls.end();
				++ nl_it) {
			
				size_t nloci = *nl_it; 
				
				//vector of nloci copies of nsites
				std::vector < size_t > nloci_nsites(nloci, nsites);
			
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
							
							
							std::string filenamepop = "PopPtrEmpty"; // alter this one
							
							testSub(pop, filenamepop, theta, growth, nloci, n_examples, nsites, nloci_nsites);
							doneEmptyPtr = true;
						}
						if (!doneEmptyPop) {
							boost::shared_ptr< hudson_ms::PopulationStructure > 
										pop(new hudson_ms::PopulationStructure() );
							
							
							std::string filenamepop = "PopEmpty"; // alter this one
							
							testSub(pop, filenamepop, theta, growth, nloci, n_examples, nsites, nloci_nsites);
							doneEmptyPop = true;
						}
						// do one set for each population type
						{
							boost::shared_ptr< hudson_ms::PopulationStructure > 
										pop = makePopSingle();
							
							
							std::string filenamepop = "PopSingle"; // alter this one
							
							testSub(pop, filenamepop, theta, growth, nloci, n_examples, nsites, nloci_nsites);
						}
						// more complex pop
						{
							boost::shared_ptr< hudson_ms::PopulationStructure > pop = makePopIslandComplex(Nzero);
							
							std::string filenamepop = "PopIslandComplex"; // alter this one
							
							testSub(pop, filenamepop, theta, growth, nloci, n_examples, nsites, nloci_nsites);
						}
					}
				}
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
	
	cout << "\n\n*** " << filenamepop << " ***" << endl;
		
	int seed = 2345;
	boost::shared_ptr < PRNGen > r ( new PRNGenGSL(seed) );
	
	/*
	 * use this prng to spawn others for sampler and models.
	 * This ensures that results are independent of how
	 * internals of either work, and we can increase samples
	 * and models and first results will be the same as before
	 * the increase
	 * */
	
	boost::shared_ptr < PRNGen > r_sampler = r->spawnAnother();

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
	
	/*use a ModelSetConfigBuilder to build the model data */
	mct::ModelSetConfigBuilder configBuilder;
	configBuilder.addModelConfig("SMM", nsites); // SMM model
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
	
	cout << "about to make sampler" << endl;
				
	// configure a sampler with modelset
	boost::shared_ptr < MultiLociMicrosatSampler > 
				sampler ( new MultiLociMicrosatSampler(models_ptr, r_sampler) );
	
	double realvalues[] = {theta,growth}; 
	
	boost::shared_ptr < aabc::ParameterParticle > p( new aabc::ParameterParticle(realvalues, 2) );

	cout << "about to use sampler" << endl;
	
	try {
		boost::shared_ptr < mct::MultiLociMicrosatSet > 
			refset_ptr = sampler->sample(n_examples,
			pop, nloci_nsites, p);
		
		if (noPop) throw std::logic_error("Should not be able to do that");
		
		cout << "The multi-locus microsat set is\n" << *refset_ptr << endl;
		
		cout << "\nLooking at the summary stats which should be the means for each individual over the reps" << endl;
		
		boost::shared_ptr < SummaryStatisticSet > ss_set = refset_ptr->getSummaryStatisticSet();
		
		cout << ss_set->toString() << endl;
	}
	catch (std::invalid_argument& ia) {
		cout << "Error :\n" << ia.what() << endl;
		if(!noPop) throw;
	} 
	
}

