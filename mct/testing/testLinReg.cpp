/* MCT - Markov Chains on Trees.


   Copyright (C) 2011, 2012 Jennifer Harlow

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
\brief Basic testing of SimdData sampler.

No adaptive abc, just testing the sampler,
ie sampling with one set of parameter values.

The MultiLociSimDataSampler uses libsequence
to simulate snp data for given parameters.  

*/

#include "config.h"

#include "pop_mc_abc_simulator_adaptive.hpp"
#include "descriptive_stats.hpp"
#include "summary_statistic_set.hpp"

#include "parameter_particle_validator_multi_loci_polytable.hpp"
#include "multi_loci_simdata_sampler.hpp"
#include "multi_loci_polytable_sampler.hpp"
#include "multi_loci_polytable_set.hpp"

#include "parameter_particle.hpp"

#include "parameter_prior_uniform.hpp"
#include "parameter_prior.hpp"
#include "parameter_prior_set.hpp"

#include "prng.hpp"
#include "prng_gsl.hpp"

#include "reg_linreg_step.hpp"
#include "reg_transformer_null.hpp"
#include "reg_transformer_log.hpp"
#include "reg_transformer_tan.hpp"

#include <mscplusplus/population_structure.hpp>

#include <boost/lambda/bind.hpp>
#include <boost/lambda/lambda.hpp>

#include <boost/shared_ptr.hpp>

#include <iostream>
#include <sstream>
#include <vector>
#include <stack>


using namespace std;
using namespace aabc;
using namespace reg;

/* test lin reg with SimData objects
 * use lin reg to do all distance tolerance checks*/
void linregTesting1();

/* test lin reg with SimData objects
 * use simulator to do all distance tolerance checks
 * does not produce output for Kevin's system */
void linregTesting2();


int main()
{
	linregTesting1();

	#if(0)
	linregTesting2();
	#endif
	
	return 0;
	
}

	

// test lin reg with SimData objects
// use lin reg to do all distance tolerance checks
void linregTesting1()
{
	try {
		
		int myPid = getpid();
		std::cout << "Testing linear regression " << myPid << std::endl;
		std::cout << "This process id is " << myPid << std::endl;
		
		int seed = 2345;
		boost::shared_ptr < mct::PRNGen > r ( new mct::PRNGenGSL(seed) );

		size_t nsam = 10; // samples
		// done (100, 10), (100,0), (100,100), (10,0) (10,10), (10,100) (50,0), (50,10), (50,100)
		
		// loci in each multi-loci SimData
		size_t nloci = 5; 
		
		size_t nsites = 25000; // number of sites in each locus 
		
		//vector of nloci copies of nsites
		std::vector < size_t > nloci_nsites(nloci, nsites);
		
		// theta is scaled mutation rate per locus
		// ie 4Ne*mu, mu = mut. rate per site * nsites// 
		double theta = 50.0; 
		
		double growth = 100.0; // growth rate, scaled per locus
		
		// number of repetitions in each comparison sample generated for each parameter particle
		size_t n_reps = 1; 
		
		//size_t iterations = 10;
		size_t iterations = 100000;
		//size_t iterations = 1000;
		
		// the percentile distance on which to base calculations of epsilon
		double baseEpsilonPercentile = 1.0;
		
		/*
		 * use this prng to spawn others for priors and simulator
		 * and sampler and models.
		 * This ensures that results are independent of how
		 * internals of any of them work, and we can increase samples
		 * and models tc and first results will be the same as before
		 * the increase
		 * */
		
		boost::shared_ptr < mct::PRNGen > r_priors = r->spawnAnother();

		boost::shared_ptr < mct::PRNGen > r_simulator = r->spawnAnother();
		
		boost::shared_ptr < mct::PRNGen > r_sampler = r->spawnAnother();
	
		boost::shared_ptr < mct::PRNGen > r_true_sampler = r->spawnAnother();
	
		
		boost::shared_ptr < mct::MultiLociPolyTableSampler > sampler ( new mct::MultiLociSimDataSampler(r_sampler) );
		
		boost::shared_ptr < mct::MultiLociPolyTableSampler > true_sampler ( new mct::MultiLociSimDataSampler(r_true_sampler) );
		
		double realvalues[] = {theta,growth}; 
		
		boost::shared_ptr < ParameterParticle > true_p( new ParameterParticle(realvalues, 2) );
		
		// number to generate with 'true' parameters
		size_t n_examples = 5;
		//size_t n_examples = 1;
		
		//population structure
		boost::shared_ptr< hudson_ms::PopulationStructure > pop
		(new hudson_ms::PopulationStructure(nsam));
		
		//get a reference set using the known parameters
		boost::shared_ptr < mct::MultiLociPolyTableSet > 
				refset = true_sampler->sample(n_examples, pop, nloci_nsites, true_p);
		
		const boost::shared_ptr < const mct::MultiLociPolyTable > 
						ref = refset->at(0);
		
		boost::shared_ptr < mct::SummaryStatisticSet > 
						ref_stats = refset->getSummaryStatisticSet();
		
		std::string ref_filename = "linregTestRefStats.txt";
		
		std::string path("../output/");
			
		ref_filename = path + ref_filename;
		
		//output the reference stats files so that we can use the original ABCreg on them
		ref_stats->outputToFile(ref_filename, false); // overwrite the file
		
		ParameterPriorSet pset(2);
		
		// need mins and maxes for priors
		double mymins[] = {0.0001,0.0}; 
		double mymaxs[] = {100.0,200.0};
		std::vector < double > mins (mymins, mymins + sizeof(mymins) / sizeof(double) );
		std::vector < double > maxs (mymaxs, mymaxs + sizeof(mymaxs) / sizeof(double) );
		
		// theta priors
		pset.at(0) = boost::shared_ptr < ParameterPrior > 
					( new ParameterPriorUniform(mins.at(0),maxs.at(0),r_priors)); // theta
		
		//growth priors
		pset.at(1) = boost::shared_ptr < ParameterPrior > 
						( new ParameterPriorUniform(mins.at(1),maxs.at(1),r_priors)); //growth
		
		boost::shared_ptr < ParameterParticleValidatorMultiLociPolyTable > 
					validator( new ParameterParticleValidatorMultiLociPolyTable(
									ref, sampler, n_reps) );
		
		PopMCabcSimulatorAdaptive simulator(pset,
								validator,
								baseEpsilonPercentile,
								r_simulator);
								
		double mydoubles[] = {1.0}; //0.5, 0.75, 0.75, 1.0
		std::stack < double > epsilonMults (std::deque <double>(mydoubles, mydoubles + sizeof(mydoubles) / sizeof(double) ));
		
		
		std::string filenameStart = "linregTest";
		
		//iterate
		while (!epsilonMults.empty()) {
			
			simulator.iterate(epsilonMults.top(), iterations);
			
			double epsilon = simulator.getCurrentEpsilon();
			std::cout << "epsilon = " << epsilon << std::endl;
			std::cout << std::endl;
			
			boost::shared_ptr < ParameterParticleSet >
					particles = simulator.getCurrentParticles();
			
			boost::shared_ptr < mct::SummaryStatisticSet >
					stats = simulator.getCurrentSummaryStatistics();
			/*
			std::ostringstream stm1;
			stm1.precision(3);
			stm1 <<   fixed;
			stm1 << filenameStart << "_eM_" << epsilonMults.top();
			std::string filename = stm1.str();
			*/
			
			
			std::string filenameP = filenameStart + "_P.txt";
			std::string filenameS = filenameStart + "_S.txt";
			
			std::string path("../output/");
			
			filenameP = path + filenameP;
			filenameS = path + filenameS;
			
			//output the files so that we can use the original ABCreg on them
			particles->outputToFile(filenameP, false); // overwrite file
			stats->outputToFile(filenameS, false); // overwrite file
			
			
			double tol = 0.01;
			//double tol = 0.8;
			
			boost::shared_ptr < reg::Transformer > trans(new reg::TanTransformer(mins,maxs) );
			reg::LinearRegressionStep lr(trans);
			std::string basefile = "../output/output";
			std::ostringstream stm2;
			stm2.precision(3);
			stm2 <<   fixed;
			stm2 << basefile << "_th_" << theta;
			stm2 << "_g_" << growth;
			stm2 << "_bEP_" << baseEpsilonPercentile;
			stm2 << "_eM_" << epsilonMults.top();
			stm2 << "_tol_" << tol;	
			basefile = stm2.str();
						
			lr.linreg(particles, stats, ref_stats, tol, basefile);

			epsilonMults.pop();
			
		}
		std::cout << "\n\nEnd of testing 1\n" << std::endl;
	
	}
	catch (std::exception& e) {
	
		std::cout << "std::exception::\n" + std::string(e.what()) << std::endl;
	}
	catch (...) {
	
		std::cout << "Unknown error\n" << std::endl;
	}
		
}

// test lin reg with SimData objects
// use simulator to do all distance tolerance checks
// does not produce output for Kevin's system

/*
 * Note that for tan transformations, Thornton uses mins and maxs 
 * from data referred to as 'prior', ie
 * the parameter values to be transformed before being regressed.
 * This is not what Hamilton et al use - they use (or say they do) the 
 * actual min and max bounds on the prior as the distribution
 * from which possible parameter values are generated. 
 * See Hamilton et al. 2005 PNAS 7476 .
 * 
 * I use the same method as Hamilton.  See the docs for the 
 * reg::TransformerTan class.  I think that Kevin is doing
 * something wrong here.  My output will differ to his if
 * the tan transformation is used
 */
	
void linregTesting2()
{
	try {
		
		int myPid = getpid();
		std::cout << "Testing linear regression " << myPid << std::endl;
		std::cout << "This process id is " << myPid << std::endl;
		
		int seed = 2345;
		boost::shared_ptr < mct::PRNGen > r ( new mct::PRNGenGSL(seed) );
		
		size_t nsam = 10; // samples
		
		// loci in each multi-loci SimData
		size_t nloci = 5; 
		
		size_t nsites = 25000; // number of sites in each locus 
		
		//vector of nloci copies of nsites
		std::vector < size_t > nloci_nsites(nloci, nsites);
		
		// theta is scaled mutation rate per locus
		// ie 4Ne*mu, mu = mut. rate per site * nsites// 
		double theta = 50.0; 
		
		double growth = 100.0; // growth rate, scaled per locus
		
		// number of repetitions in each comparison sample generated for each parameter particle
		size_t n_reps = 1; 
		
		int priors_seed = 1 + r->uniform_int( r->getMax() );
		
		int simulator_seed = 1 + r->uniform_int( r->getMax() );
		
		int sampler_seed = 1 + r->uniform_int( r-> getMax() );
				
		boost::shared_ptr < mct::PRNGen > r_true_sampler = r->spawnAnother();
		
		boost::shared_ptr < mct::MultiLociPolyTableSampler > 
				true_sampler ( new mct::MultiLociSimDataSampler(r_true_sampler) );
		
		double realvalues[] = {theta,growth}; 
		
		boost::shared_ptr < ParameterParticle > true_p( new ParameterParticle(realvalues, 2) );
		
		// number to generate with 'true' parameters
		size_t n_examples = 5;
		
		//population structure
		boost::shared_ptr< hudson_ms::PopulationStructure > pop
		(new hudson_ms::PopulationStructure(nsam));
		
		//get a reference set using the known parameters
		boost::shared_ptr < mct::MultiLociPolyTableSet > 
				refset = true_sampler->sample(n_examples, pop, nloci_nsites, true_p);
		
		std::string ref_filename = "linregTestNewRefStats.txt";
		
		std::string path("../output/");
			
		ref_filename = path + ref_filename;
		
		//output the reference stats files so that we can use the original ABCreg on them
		{
			boost::shared_ptr < mct::SummaryStatisticSet > 
						ref_stats = refset->getSummaryStatisticSet();
			ref_stats->outputToFile(ref_filename, false); // overwrite the file
		}
		
		// do for each set of ref stats
		for (size_t refindex = 0; refindex < refset->size(); ++refindex) {
			std::cout << "refindex = " << refindex << std::endl;
			
			boost::shared_ptr < mct::PRNGen > r_simulator( new mct::PRNGenGSL (simulator_seed) );
		
			boost::shared_ptr < mct::PRNGen > r_sampler( new mct::PRNGenGSL (sampler_seed) );
			
			boost::shared_ptr < mct::PRNGen > r_priors( new mct::PRNGenGSL (priors_seed) );
			
			ParameterPriorSet pset(2);
			
			// theta priors
			pset.at(0) = boost::shared_ptr < ParameterPrior > 
						( new ParameterPriorUniform(0.0001,200.0,r_priors)); // theta
			
			//growth priors
			pset.at(1) = boost::shared_ptr < ParameterPrior > 
							( new ParameterPriorUniform(0.0,100.0,r_priors)); //growth
							
			// sampler
			boost::shared_ptr < mct::MultiLociPolyTableSampler > 
					sampler ( new mct::MultiLociSimDataSampler(r_sampler) );
					
			// ref stat and a set just of this ref stat for this loop
			const boost::shared_ptr < const mct::MultiLociPolyTable > 
						ref = refset->at(refindex);
			
			boost::shared_ptr < mct::SummaryStatisticSet > 
						this_ref_stats(new mct::SummaryStatisticSet());
			this_ref_stats->add(ref->getSummaryStatistic());

			// validator
			boost::shared_ptr < ParameterParticleValidatorMultiLociPolyTable > 
					validator( new ParameterParticleValidatorMultiLociPolyTable(
									ref, sampler, n_reps) );
			
			// the percentile distance on which to base calculations of epsilon
			double baseEpsilonPercentile = 0.01;
		
			// simulator
			PopMCabcSimulatorAdaptive simulator(pset,
									validator,
									baseEpsilonPercentile,
									r_simulator);
			
			double mydoubles[] = {1.0}; //0.5, 0.75, 0.75, 1.0
			std::stack < double > epsilonMults (std::deque <double>(mydoubles, mydoubles + sizeof(mydoubles) / sizeof(double) ));
			
			size_t iterations = 1000;
			//size_t iterations = 1000;
		
			//iterate
			while (!epsilonMults.empty()) {
				
				std::cout << "starting iterations " << std::endl;
				
				simulator.iterate(epsilonMults.top(), iterations);
				
				double epsilon = simulator.getCurrentEpsilon();
				std::cout << "epsilon = " << epsilon << std::endl;
				std::cout << std::endl;
							
				boost::shared_ptr < ParameterParticleSet >
					particles = simulator.getCurrentParticles();
			
				boost::shared_ptr < mct::SummaryStatisticSet >
					stats = simulator.getCurrentSummaryStatistics();
				
				//double tol = 0.01;
				double tol = 1.0;
				
				boost::shared_ptr < reg::Transformer > trans(new reg::LogTransformer() );
				reg::LinearRegressionStep lr(trans);
				std::string basefile = "../output/output_new";
				std::ostringstream stm2;
				stm2.precision(3);
				stm2 <<   fixed;
				stm2 << basefile << "_th_" << theta;
				stm2 << "_g_" << growth;
				stm2 << "_bEP_" << baseEpsilonPercentile;
				stm2 << "_eM_" << epsilonMults.top();
				stm2 << "_tol_" << tol << "." << refindex;;	
			
				basefile = stm2.str();
							
				lr.linreg(particles, stats, this_ref_stats, tol, basefile);
				
				epsilonMults.pop();
				
			}// end while through epsilon mults
			
		} // end loop through ref stats
		
		// now do it again but to produce Kevin's numbers
		{
			std::cout << "doing Kevin's files " << std::endl;
			
			boost::shared_ptr < mct::PRNGen > r_simulator( new mct::PRNGenGSL (simulator_seed) );
			
			boost::shared_ptr < mct::PRNGen > r_sampler( new mct::PRNGenGSL (sampler_seed) );
			
			boost::shared_ptr < mct::PRNGen > r_priors( new mct::PRNGenGSL (priors_seed) );
			
			ParameterPriorSet pset(2);
			
			// theta priors
			pset.at(0) = boost::shared_ptr < ParameterPrior > 
						( new ParameterPriorUniform(0.0001,100.0,r_priors)); // theta
			
			//growth priors
			pset.at(1) = boost::shared_ptr < ParameterPrior > 
							( new ParameterPriorUniform(0.0,200.0,r_priors)); //growth
							
			// sampler
			boost::shared_ptr < mct::MultiLociPolyTableSampler > 
					sampler ( new mct::MultiLociSimDataSampler(r_sampler) );
					
			// ref stat and a set just of this ref stat for this loop
			// ref stat should not matter - we are going to take all of them
			const boost::shared_ptr < const mct::MultiLociPolyTable > 
						ref = refset->at(0);
			
			boost::shared_ptr < mct::SummaryStatisticSet > 
						this_ref_stats(new mct::SummaryStatisticSet());
			this_ref_stats->add(ref->getSummaryStatistic());

			// validator
			boost::shared_ptr < ParameterParticleValidatorMultiLociPolyTable > 
					validator( new ParameterParticleValidatorMultiLociPolyTable(
									ref, sampler, n_reps) );
			
			// the percentile distance on which to base calculations of epsilon
			double baseEpsilonPercentile = 1.0;
		
			// simulator
			PopMCabcSimulatorAdaptive simulator(pset,
									validator,
									baseEpsilonPercentile,
									r_simulator);
			
			double mydoubles[] = {1.0}; //0.5, 0.75, 0.75, 1.0
			std::stack < double > epsilonMults (std::deque <double>(mydoubles, mydoubles + sizeof(mydoubles) / sizeof(double) ));
			
			
			std::string filenameStart = "linregTest";
			
			size_t iterations = 100000;
			//size_t iterations = 1000;
		
			//iterate
			while (!epsilonMults.empty()) {
				
				std::cout << "starting iterations " << std::endl;
				
				simulator.iterate(epsilonMults.top(), iterations);
				
				double epsilon = simulator.getCurrentEpsilon();
				std::cout << "epsilon = " << epsilon << std::endl;
				std::cout << std::endl;
							
				boost::shared_ptr < ParameterParticleSet >
					particles = simulator.getCurrentParticles();
			
				boost::shared_ptr < mct::SummaryStatisticSet >
					stats = simulator.getCurrentSummaryStatistics();
				/*	
				std::ostringstream stm1;
				stm1.precision(3);
				stm1 <<   fixed;
				stm1 << filenameStart << "_eM_" << epsilonMults.top();
				std::string filename = stm1.str();
				* */
				std::string filenameP = filenameStart + "_P.txt";
				std::string filenameS = filenameStart + "_S.txt";
				
				std::string path("../output/");
			
				filenameP = path + filenameP;
				filenameS = path + filenameS;
					
				//output the files so that we can use the original ABCreg on them
				particles->outputToFile(filenameP, false); // overwrite file
				stats->outputToFile(filenameS, false); // overwrite file
							
				epsilonMults.pop();
				
				
			}// end while through epsilon mults
			
		} // end producing Kevin's stats
		
		std::cout << "\n\nEnd of testing 2\n" << std::endl;
	
	}
	catch (std::exception& e) {
	
		std::cout << "std::exception::\n" + std::string(e.what()) << std::endl;
	}
	catch (...) {
	
		std::cout << "Unknown error\n" << std::endl;
	}
		
}

