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
\brief Testing of microsat sampler with a population of subpopulations,
looking at the inter-population (population comparison) statistics.

Sampling with one set of parameter values.

The MultiLociMicrosatSampler uses our MCT models to
generate microsate repeat data.  

*/

#include "config.h"

#include "descriptive_stats.hpp"
#include "summary_statistic_set.hpp"

#include "multi_loci_microsat_sampler.hpp"
#include "multi_loci_microsat_set.hpp"
#include "multi_loci_microsat.hpp"

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


void microsatPopStatsTesting();

void testSub(
	const boost::shared_ptr< const hudson_ms::PopulationStructure > pop,
	const std::string& filenamepop,
	double theta, double growth, 
	size_t nloci, size_t n_examples, size_t nsites, 
	const vector < size_t>& nloci_nsites);


int main()
{
	microsatPopStatsTesting();
	
	
	return 0;
	
}

	
/* test that microsats population stats are working
 * NOTE this does not say the stats are right, just that 
 * they are happening when they are supposed to! */
void microsatPopStatsTesting()
{
	try {
		int myPid = getpid();
		std::cout << "basic testing  for microsats " << std::endl;
		std::cout << "This process id is " << myPid << std::endl;
		
		// loci in each multi-loci SimData
		size_t nloci = 5; 
		
		size_t nsites = 1; // number of sites in each locus 
		vector < size_t> nloci_nsites(nloci, nsites);
		
		size_t Nzero = 1000000; // effective popn size (per sub-pop)
		double unit = 100000.0; // 10^5
		
		//mutation rates
		//mutation rates PER SITE ie 4N0(mu/site)
		double mu = 1.0; 
		
		
		//mults by 4*Nzero/unit
		double theta = mu*4*Nzero/unit;
		
		//growth rates
		double growth = 50.0;
		
		// number to generate with 'true' parameters
		size_t n_examples = 1;
		
				
		// do one set for each population type
		{
			boost::shared_ptr< hudson_ms::PopulationStructure > 
						pop = makePopSingle();
			
			
			std::string filenamepop = "PopSingle"; // alter this one
			
			testSub(pop, filenamepop, theta, growth, nloci, n_examples, nsites, nloci_nsites);
		}
		// pop with 2 sub populations, one subpop with no individuals
		{
			boost::shared_ptr< hudson_ms::PopulationStructure > pop = makePopIslandOneZero(Nzero);
			
			std::string filenamepop = "PopIslandOneZero"; // alter this one
			
			testSub(pop, filenamepop, theta, growth, nloci, n_examples, nsites, nloci_nsites);
		}
		// pop with 3 subpops, one individual each
		{
			boost::shared_ptr< hudson_ms::PopulationStructure > pop = makePopIslandThreeSmallest(Nzero);
			
			std::string filenamepop = "PopIslandThreeSmallest"; // alter this one
			
			testSub(pop, filenamepop, theta, growth, nloci, n_examples, nsites, nloci_nsites);
		}
		// pop with 3 subpops, one with no samples
		{
			boost::shared_ptr< hudson_ms::PopulationStructure > pop = makePopIslandThreeOneZero(Nzero);
			
			std::string filenamepop = "PopIslandThreeOneZero"; // alter this one
			
			testSub(pop, filenamepop, theta, growth, nloci, n_examples, nsites, nloci_nsites);
		}
		// pop with 4 subpops, two with no samples
		{
			boost::shared_ptr< hudson_ms::PopulationStructure > pop = makePopIslandFourTwoZero(Nzero);
			
			std::string filenamepop = "PopIslandFourTwoZero"; // alter this one
			
			testSub(pop, filenamepop, theta, growth, nloci, n_examples, nsites, nloci_nsites);
		}
		
		// pop with 3 subpops, two individuals each
		{
			boost::shared_ptr< hudson_ms::PopulationStructure > pop = makePopIslandThreeSmall(Nzero);
			
			std::string filenamepop = "PopIslandThreeSmall"; // alter this one
			
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
	
	cout << "theta " << theta << " growth " << growth << endl;
	
	//population structure
	std::string popStr = pop->toString();
	std::cout << "Pop structure is " << popStr << endl;
	
	
	
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
	
	cout << "about to use sampler" << endl;
				
	// configure a sampler with modelset
	boost::shared_ptr < MultiLociMicrosatSampler > 
				sampler ( new MultiLociMicrosatSampler(models_ptr, r_sampler) );
	
	double realvalues[] = {theta,growth}; 
	
	boost::shared_ptr < aabc::ParameterParticle > p( new aabc::ParameterParticle(realvalues, 2) );

	boost::shared_ptr < mct::MultiLociMicrosatSet > 
		refset_ptr = sampler->sample(n_examples,
		pop, nloci_nsites, p);
	
	//for each element in the sample
	for (size_t i = 0; i < refset_ptr->size(); ++i) {
		
		boost::shared_ptr < const MultiLociMicrosat > mlms = refset_ptr->at(i);
	
		// print out the results
		mlms->printSubpopData(cout);
		
		
		int prec = 5;
		
		std::cout << "\nindividual average repeat sizes are" << std::endl;
		boost::shared_ptr < const SummaryStatistic  > means
			= mlms->getAverageRepeats();
		cout << means->plainString(prec) << endl;

		std::cout <<  "\nindividual repeat size covariances are" << std::endl;
		boost::shared_ptr < const SummaryStatistic  > covars
			= mlms->getCovariances();
		cout << covars->plainString(prec) << endl;


		boost::shared_ptr < const SummaryStatistic  > withinPopAvs 
							= mlms->withinPopAverageDifferences();
			
		std::cout << "\nwithin pop average differences are " <<std::endl;
		cout << withinPopAvs->plainString(prec) << endl;
		
		boost::shared_ptr < const SummaryStatistic  > betweenPopAvs 
							= mlms->betweenPopAverageDifferences();
		std::cout << "\nbetween pop average differences are " <<std::endl;
		cout << betweenPopAvs->plainString(prec) << endl;
		
		boost::shared_ptr < const SummaryStatistic  > correctedAvs 
							= mlms->correctedBetweenPopAverageDifferences();
		std::cout << "\ncorrected between pop average differences are " <<std::endl;
		cout << correctedAvs->plainString(prec) << endl;
		
		boost::shared_ptr < const SummaryStatistic  > pairwiseAvs 
							= mlms->pairwiseAverageDifferences();
		std::cout << "\npairwise pop average differences are " <<std::endl;
		cout << pairwiseAvs->plainString(prec) << endl;
		
		boost::shared_ptr < const SummaryStatistic  > Fst 
												= mlms->fst();
		cout << "\nFst is " << Fst->plainString(prec) << endl;
		
		boost::shared_ptr < const SummaryStatistic  > pairwiseFsts 
									= mlms->pairwiseFst();
		std::cout << "\npairwiseFsts are " <<std::endl;
		cout << pairwiseFsts->plainString(prec) << endl;

		boost::shared_ptr < const SummaryStatistic  > rst 
									= mlms->rst();
		std::cout << "\nRst is " <<std::endl;
		cout << rst->plainString(prec) << endl;
		
		boost::shared_ptr < const SummaryStatistic  > betweenRsts 
									= mlms->distanceSlatkinRstBetweenPop();
		std::cout << "\nbetweenPopRsts are " <<std::endl;
		cout << betweenRsts->plainString(prec) << endl;

		std::cout << "\nbetween populations DeltaMu" << std::endl;
		boost::shared_ptr < const SummaryStatistic  > ss
			= mlms->distanceGoldsteinDeltaMuSqBetweenPop();
		cout << ss->plainString(prec) << endl;
		
		boost::shared_ptr < const SummaryStatistic  > pairwiseAverageSumSqDiff 
									= mlms->pairwiseAverageSumSqDiffs();
		std::cout << "\npairwiseAverageSumSqDiff are " <<std::endl;
		cout << pairwiseAverageSumSqDiff->plainString(prec) << endl;

		
	}
}

