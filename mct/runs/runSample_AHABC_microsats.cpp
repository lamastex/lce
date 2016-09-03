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
\brief Sample mct::MultiLociMicrosat objects using a 
mct::MultiLociMicrosatSampler.
 
mct::MultiLociMicrosat objects represent multiple loci of microsat data
summarised as microsat repeats data.  
  
The aim here is to use some particular parameter values and to get a large 
number of mct::MultiLociMicrosat reps for any particular
set or combination of them (ie consider these a sample from the population
of possible mct::MultiLociPolyTable objects given those values), and get
the summary statistics for the reps.  

The summary statistics can be
further analysed using AHABC techniques.

*/

#include "config.h"

#include "multi_loci_microsat_sampler.hpp"
#include "multi_loci_microsat_set.hpp"
#include "multi_loci_microsat.hpp"

#include "model_set_config_builder.hpp"
#include "modelset.hpp"

#include "descriptive_stats.hpp"
#include "summary_statistic_set.hpp"

#include "parameter_particle.hpp"

#include "prng.hpp"
#include "prng_gsl.hpp"

#include "utilities.hpp"

#include "make_run_pops.hpp"

#include <mscplusplus/population_structure.hpp>

#include <boost/lambda/bind.hpp>
#include <boost/lambda/lambda.hpp>

#include <boost/shared_ptr.hpp>


#include <iostream>
#include <sstream>
#include <vector>
#include <stack>

using namespace std;
using namespace mct;
using namespace aabc;

void microsatsRaaz();

void microsatsRaazInnerLoop(double theta, double growth, double alpha, double beta,
					size_t subNsamSizeEach, size_t nloci, size_t _nsites,
					vector < size_t>& nsites, size_t n_examples );

int main()
{
	microsatsRaaz();
	
	
	return 0;
	
}

// microsat results for various theta and growth and migration rates
// for Raaz, October 2011
void microsatsRaaz()
{
	try {
		int myPid = getpid();
		std::cout << "Microsat runs for Raaz " << std::endl;
		std::cout << "This process id is " << myPid << std::endl;
		
		// size of sample from each of the (two) subpops, ie total sample is twice this
		size_t subNsamSizeEach = 15;
		
		// loci in each multi-loci microsat
		size_t nloci = 10; 
		
		size_t nsites = 1; // number of sites in each locus 
		vector < size_t> nloci_nsites(nloci,nsites);
		
		size_t Nzero = 10000; // effective popn size (per sub-pop) 10^4
		double unit = 100000.0; // 10^5
		// unit just helps with specifying mu_per_site
		// what is given here is really mu/site * unit 
		
		//mutation rates PER SITE
		double mymus[] = {1.0}; 
		// 'real' mu/site is mu_per_site divided by unit
		
		std::vector <double> mus (mymus, mymus + sizeof(mymus) / sizeof(double) );
		std::vector <double> thetas(mus.size());
		
		using namespace boost::lambda;
		
		transform(mus.begin(), mus.end(), thetas.begin(), _1*4*Nzero/unit);
		
		//growth rates
		double mygrowths[] = {0.0};
		std::vector <double> growths (mygrowths, mygrowths + sizeof(mygrowths) / sizeof(double) );
		
		//migration rates into pop0 from pop1
		// want to end up with 0.01, 0.1, 1.0, 10
		// need to put in _/(4*Nzero)
		double myalphas[] = {2.5e-7, 2.5e-6, 2.5e-5, 2.5e-4}; 
									// done 0.1 for beta 0.01, 0.1
		std::vector <double> alphas (myalphas, myalphas + sizeof(myalphas) / sizeof(double) );
		// what goes to MS is in terms of 4*N_0*m where m is fraction each subpop made up of new migrants each generation
		transform(alphas.begin(), alphas.end(), alphas.begin(), _1*4*Nzero); 
		
		
		//migration rates into pop1 from pop0
		double mybetas[] = {2.5e-7, 2.5e-6, 2.5e-5, 2.5e-4}; 
		std::vector <double> betas (mybetas, mybetas + sizeof(mybetas) / sizeof(double) );
		// what goes to MS is in terms of 4*N_0*m where m is fraction each subpop made up of new migrants each generation
		transform(betas.begin(), betas.end(), betas.begin(), _1*4*Nzero); 
		
		// number to generate with 'true' parameters
		//size_t n_examples = 1000000;
		size_t n_examples = 10000;
		
		for (std::vector < double>::iterator theta_it = thetas.begin(); 
				theta_it < thetas.end();
				++ theta_it) {
					
			for (std::vector < double>::iterator growth_it = growths.begin(); 
				growth_it < growths.end();
				++ growth_it) {
					
				for (std::vector < double>::iterator alpha_it = alphas.begin(); 
					alpha_it < alphas.end();
					++ alpha_it) {
						
					for (std::vector < double>::iterator beta_it = betas.begin(); 
						beta_it < betas.end();
						++ beta_it) {	
			
						double theta = *theta_it;
						double growth = *growth_it;
						double alpha = *alpha_it;
						double beta = *beta_it;
						
						microsatsRaazInnerLoop(theta, growth, alpha, beta,
								subNsamSizeEach, nloci, nsites, 
								nloci_nsites, n_examples );
										
					}
				}
										
			}
		}
		
	}
	catch (std::exception& e) {
	
		std::cout << "Error:\n" + std::string(e.what()) << std::endl;
	}
}


void microsatsRaazInnerLoop(double theta, double growth, double alpha, double beta,
					size_t subNsamSizeEach, size_t nloci, size_t nsites,
					vector < size_t>& nloci_nsites, size_t n_examples )
{
	
	int seed = 2345;
	boost::shared_ptr < PRNGen > r ( new PRNGenGSL(seed) );
	
	boost::shared_ptr < PRNGen > r_sampler = r->spawnAnother();

	boost::shared_ptr < PRNGen > r_models = r->spawnAnother();
	
	
	cout << "theta " << theta << " growth " << growth
		<< " alpha " << alpha << " beta " << beta << endl;
	
	//population structure
	boost::shared_ptr< hudson_ms::PopulationStructure > pop 
			= makePopRaazRuns(subNsamSizeEach,
								alpha, beta);

	std::string popStr = pop->toString();
	std::cout << "Pop structure is " << popStr << endl;
	
	// fit the precision to the magnitudes of the data
	int precTheta = mct_utilities::fitPrecision(theta);
	int precGrowth = mct_utilities::fitPrecision(growth);
	int precAlpha = mct_utilities::fitPrecision(alpha);
	int precBeta = mct_utilities::fitPrecision(beta);

	std::ostringstream stm;
	stm << "_ns_" << pop->totalNsam() << "_nl_" << nloci << "_n_" << n_examples;
	stm.precision(precTheta);
	stm <<   fixed;
	stm << "_th_" << theta;
	stm.precision(precGrowth);
	stm << "_g_" << growth;
	stm.precision(precAlpha);
	stm << "_a_" << alpha;
	stm.precision(precBeta);
	stm << "_b_" << beta << ".txt";
	
	std::string filenamestart = "RaazMicrosatStats";
	std::string filenametype1 = "DeltaMu";
	std::string filenametype2 = "Fst";
	std::string filenametype3 = "Rst";
	
	
	std::string filename1 = filenamestart + filenametype1 + stm.str();
	std::string filename2 = filenamestart + filenametype2 + stm.str();
	std::string filename3 = filenamestart + filenametype3 + stm.str();
	
	std::string path("../RaazOutput/");
			
	filename1 = path + filename1;
	filename2 = path + filename2;
	filename2 = path + filename3;
					
	//use a ModelSetConfigBuilder to build the model data 
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

	{
		// configure a sampler with modelset
		boost::shared_ptr < MultiLociMicrosatSampler > 
					sampler ( new MultiLociMicrosatSampler(models_ptr, r_sampler) );
		
		double realvalues[] = {theta,growth}; 
		
		boost::shared_ptr < ParameterParticle > 
			p( new ParameterParticle(realvalues, 2) );

		boost::shared_ptr < mct::MultiLociMicrosatSet > 
			refset = sampler->sample(n_examples,
			pop, nloci_nsites, p);
			
		cout << "finished sampling, now getting stats ..." << endl;
		
		
		//deltamus
		{
			boost::shared_ptr < SummaryStatisticSet > 
					stats = refset->getDistanceDeltaMuSqBetweenPop();
			
			stats->outputToFile(filename1, true); // append
			
			if (!stats->empty() ) {
				for ( size_t index = 0; index < stats->at(0)->size(); ++index ) {
					cout << "\tAt index " << index << "\t";
					DescriptiveStats ds(stats, index);
					std::cout << (ds.toString()) << std::endl;
				}
				std::cout << std::endl;
			}
			
		}
		
		// Fst
		{
			boost::shared_ptr < SummaryStatisticSet > 
					stats = refset->getBetweenPopFst();
			stats->outputToFile(filename2, true); // append
			
			if (!stats->empty() ) {
				for ( size_t index = 0; index < stats->at(0)->size(); ++index ) {
					cout << "\tAt index " << index << "\t";
					DescriptiveStats ds(stats, index);
					std::cout << (ds.toString()) << std::endl;
				}
				std::cout << std::endl;
			}
		}
		
		// Rst
		{
			boost::shared_ptr < SummaryStatisticSet > 
					stats = refset->getDistanceSlatkinRstBetweenPop();
			stats->outputToFile(filename3, true); // append
			
			if (!stats->empty() ) {
				for ( size_t index = 0; index < stats->at(0)->size(); ++index ) {
					cout << "\tAt index " << index << "\t";
					DescriptiveStats ds(stats, index);
					std::cout << (ds.toString()) << std::endl;
				}
				std::cout << std::endl;
			}
		}
		cout << "\n\n" << endl;
	}			
}

