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
\brief Run AABC using a mct::MultiLociSimDataSampler but fix
the process so that I use a 'good' piece of observed data, ie
one in a high density area of the distribution of data given
the true parameter values. In this run I use a 
mct::MultiLociSimDataSampler to simulate possible observed data
objects and set bounds for what summary statistics it can have
to be acceptable (and use the first one that meets all criteria).
* 
The aim of this is to contrast what happens with 'good' data
compared to what can happen with 'bad' data (from a low density 
part of the distribution).

AABC is aiming to create a population of parameter particles which seem 
'close to' the parameter values that might have created some observed
data.  The type of the observed data here is a mct::MultiLociPolyTable object
(multi-locus nucleotide data summarised as snp information).

For any parameter particle, 'close to' is here assessed by simulating other
mct::MultiLociPolyTable objects using that particle and calculating
a distance between the simulations and the observed data (using
summary statistics etc etc). 

mct::MultiLociPolyTable objects represent multiple loci of nucleotide data
summarised as snp data.  
 
mct::MultiLociPolyTable objects can be simulated in various ways. The
'simdata' sampler uses libsequence directly to do the simulation
of the snp data. We could also use our own mct models (by using
a mct::MultiLociPolySitesSampler). For a simple
one-population population structure, the main difference is
that simdata sampling is much quicker.  If you wanted to do this for
a population of subpopulations, you would have to adapt this to 
use a mct::MultiLociPolySitesSampler instead

*/

#include "config.h"

#include "pop_mc_abc_simulator_adaptive.hpp"
#include "parameter_particle_validator_multi_loci_polytable.hpp"
#include "multi_loci_simdata_sampler.hpp"
#include "multi_loci_polytable_sampler.hpp"
#include "multi_loci_polytable_set.hpp"
#include "multi_loci_polytable.hpp"

#include "parameter_prior_set.hpp"
#include "parameter_prior_uniform.hpp"
#include "parameter_prior.hpp"

#include "parameter_particle_set_analysis.hpp"
#include "parameter_particle_set.hpp"
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
#include <stack>

using namespace std;
using namespace mct;
using namespace aabc;

void runAABC_SimdataFixObserved();

void outputStandardisedParticles(
	boost::shared_ptr < const ParameterParticleSet > particles,
	const std::string& filename,
	const std::string& filenameInfo,
	boost::shared_ptr < const ParameterParticle > true_p);
	
int main()
{
	runAABC_SimdataFixObserved();
	
	
	return 0;
	
}

void runAABC_SimdataFixObserved()
{
	try {
		
		int myPid = getpid();
		std::cout << "Sim data aabc using a 'good' summary stat" << std::endl;
		std::cout << "This process id is " << myPid << std::endl;
		
		int seed = 2345;
		boost::shared_ptr < PRNGen > r ( new PRNGenGSL(seed) );

		/* How many iterations we want in each generation - ie how many members of the 
		 * population we are aiming for in each generation */
		size_t iterations = 10000;
				
		size_t nsam = 10; // samples
		size_t nsites = 25000; // number of sites in each locus 
		
		// loci in each multi-loci SimData
		size_t nloci = 5; 
		
		//vector of nloci copies of nsites
		std::vector < size_t > nloci_nsites(nloci, nsites);
		
		size_t Nzero = 1000000; // effective popn size (per sub-pop)
		
		double unit = 1000000000.0; // 10^9
		// unit just helps with specifying mu_per_site
		// what is given here is mu/site * unit
		double mu_per_site = 0.5; 
		// 'real' mu/site is mu_per_site divided by unit
		
		/* theta
		 * theta PER SITE ie 4N0(mu/site) */
		double theta = mu_per_site*4*Nzero/unit;  
		
		#if(0)
			double mythetas[] = {1.0, 2.5, 5.0, 7.5, 10.0}; //0.001,10.0,50.0,100.0
			std::vector <double> mus (mymus, mymus + sizeof(mymus) / sizeof(double) );
			std::vector <double> thetas(mus.size());
		
			using namespace boost::lambda;
			mults by 4*Nzero/lambda
			transform(mus.begin(), mus.end(), thetas.begin(), _1*4*Nzero/unit);
		#endif
		
		double growth = 10.0; // growth rate, scaled per locus
		
		/*
		 * use this prng to spawn others for priors and simulator
		 * and sampler and models.
		 * This ensures that results are independent of how
		 * internals of any of them work, and we can increase samples
		 * and models etc and first results will be the same as before
		 * the increase (I hope anyway ...)
		 * */
		
		boost::shared_ptr < PRNGen > r_priors = r->spawnAnother();

		boost::shared_ptr < PRNGen > r_simulator = r->spawnAnother();
		
		boost::shared_ptr < PRNGen > r_sampler = r->spawnAnother();
	
		boost::shared_ptr < PRNGen > r_true_sampler = r->spawnAnother();

		// make the sampler we will use with the aabc stuff
		boost::shared_ptr < MultiLociPolyTableSampler > 
				sampler ( new MultiLociSimDataSampler(r_sampler) );
		
		/* also want one that works on the true parameter values, so
		 * that we can use this to generate a possible observed value*/
		boost::shared_ptr < MultiLociPolyTableSampler > 
				true_sampler ( new MultiLociSimDataSampler(r_true_sampler) );
		
		double realvalues[] = {theta,growth}; 
		
		// a paramter particle of these true values = note the order = theta and then growth
		boost::shared_ptr < ParameterParticle > true_p( new ParameterParticle(realvalues, 2) );
		
		// number to generate with 'true' parameters
		size_t n_examples = 1;
		
		/* Create a simple one-subpop population structure for nsam samples */
		boost::shared_ptr< hudson_ms::PopulationStructure > pop
		(new hudson_ms::PopulationStructure(nsam));
		
		/* Use the true sampler with various parameters including the 
		 * number of examples to get, the population structure, loci, sites
		 * and the parameter particle with the true values in, to
		 * get a reference object  
		 * using the true parameters (this would be from 'observed' data)*/
		boost::shared_ptr < mct::MultiLociPolyTableSet > 
				refset = true_sampler->sample(n_examples, pop, nloci_nsites, true_p);
		
		/* The reference set should contain one object (assuming n_example = 1)
		 * and we need to get that individual object itself */
		boost::shared_ptr < const mct::MultiLociPolyTable > ref = refset->at(0);
		
		/* BUT I want to fix it so that this is a 'good' bit of observed data, ie
		 * in a high density part of the distribution of thingummies given
		 * these true parameters. I must have looked at the distribution of summary
		 * stats for these parameters from another run and decided to hunt 
		 * for something that is in the bands specified below.  I get the
		 * summmary stats for the observed data I am looking at now ...*/
		boost::shared_ptr < const SummaryStatistic > refstats = ref->getSummaryStatistic();
		
		// and compare to the min and max bounds (seg sites and heterozygosity)
		
		/* Need to take care setting these - make sure they are appropriate for 
		 * the 'true' parameters - this means having done some previous run
		 * really to look at the shape of the distribution */
		double ss_min = 44.5;
		double ss_max = 45.5;
		double h_min = 14.0;
		double h_max = 14.5;
		// and keep sampling while not in the bounds		
		while ( refstats->at(0) < ss_min || refstats->at(0) > ss_max
				|| refstats->at(1) < h_min || refstats->at(1) > h_max) { // this is not what we want
				
			refset = true_sampler->sample(n_examples, pop, nloci_nsites, true_p);
			ref = refset->at(0);
			refstats = ref->getSummaryStatistic();
		
		}
		// I'll drop out of the loop when I have a 'good' piece of observed data
		
		/* Note - would be better to put some limit on how many tries here - if the
		 * bounds are inappropriate it could just keep going for ever ... */
		cout << "Have reference object, SummaryStats = " 
			<< *(ref->getSummaryStatistic()) << std::endl;	
		
		// Now make an object to hold the priors for the parameter values
		// It will hold two priors (one for theta, one for growth, in that order)
		ParameterPriorSet pset(2);
		
		// make and add theta priors (uniform) - remember theta is per site
		double theta_min = 0.00000001; // dont use 0.0!
		double theta_max = 0.005;
		// the theta parameter prior also gets a prng, r_priors
		pset.at(0) = boost::shared_ptr < ParameterPrior > 
				( new ParameterPriorUniform(theta_min, theta_max, r_priors)); // theta
		
		//growth priors
		double g_min = 0.0;
		double g_max = 100.0;
		// the growth parameter prior also gets the same prng, r_priors
		pset.at(1) = boost::shared_ptr < ParameterPrior > 
				( new ParameterPriorUniform(g_min, g_max, r_priors)); //growth
		
		/* Setting up the validator.  This is the thing that says how
		 * close a given particle is to what wanted.  It can do this 
		 * lots of ways but the one we use for the multi-loci stuff
		 * uses the given particle to generate a sample or some samples and then
		 * calculates a distance using summary statistics from the samples
		 * compared to summary statistics from the reference object*/
		
		/* n_reps is the number of repetitions in each comparison sample 
		 * generated for each parameter particle
		 * This would probably always be 1.  
		 * I set it up so that the validator might be told to simulate multiple objects with
		 * the same particle and measure 'close enough' using the means of the statistics from
		 *  all of them, but I am not sure whether this might violate some of the 
		 * theory behind the sampling etc etc. Stick to 1 I  think! */
		size_t n_reps = 1; 
		
		/* Now create the validator that the aabc simulator will use to check how close
		 * some simulated particle is to what it wants to home in on.
		 * The validator we use for the multi-loci stuff needs
		 * to know about the reference object, how many simulations it uses
		 * to assess distance to the reference object 
		 * (n_rep, usually 1, as discussed above), and
		 * it needs the sampler so that it can actually generate samples
		 * from using each parameter particle it is testing. */
		boost::shared_ptr < ParameterParticleValidatorMultiLociPolyTable > 
					validator( new ParameterParticleValidatorMultiLociPolyTable(
									ref, sampler, n_reps) );
		
		
		/* Now the actual aabc parameter particle simulator. This needs
		 * the validator and a prng, but it also needs to some idea of 
		 * 'close enough'.  This is based on a percentile value (which
		 * takes away the problem of judging an appropriate epsilon
		 * before having any idea of the range of distances that
		 * we might get between
		 * simulated multi-loci objects and the reference object,
		 * but gets a little complicated in the actual simulations - see
		 * docs for the simulator and hope that I've written something
		 * intelligible there.  Meanwhile, here I am specifying
		 * an initial percentile of 10%. */
		
		// the percentile distance on which to base calculations of epsilon
		double baseEpsilonPercentile = 0.1;
		
		// Finally create the simulator itself
		PopMCabcSimulatorAdaptive simulator(pset,
								validator,
								baseEpsilonPercentile,
								r_simulator);

		/* In each generation, we get tighter and tighter in our 
		 * idea of 'close enough'.  These multipliers are applied 
		 * progressively, ie 0.5 * the value we get from the base 
		 * epsilon criteria for the first loop, then 0.75 * that, then
		 * 0.75 * that, then 1.0 * that */
		double mydoubles[] = {0.5, 0.75, 0.75, 1.0}; //2.0, 3.0, 4.0, 5.0, 10.0
		std::stack < double > epsilonMults (std::deque <double>(
				mydoubles, mydoubles + sizeof(mydoubles) / sizeof(double) ));
				
		/* The number of epsilon mults we give controls how many loops
		 * (generations) we do - each loop pops a value from the stack */
		
		// create a string that summarises all our different values, to label output
		std::ostringstream stm;
		stm << "_s_" << seed << "_sam_" << nsam << "_ns_" << nsites << "_nl_" << nloci 
				<< "_nr_" << n_reps << "_ni_" << iterations;
		stm.precision(2);
		stm <<   fixed;
		stm <<  "_bE_" << baseEpsilonPercentile;
		stm.precision(1);
		stm <<   fixed;
		stm.precision(5);
		stm << "_th_" << theta;
		stm.precision(3);
		stm << "_g_" << growth;
		std::string filenameStart = "SimdataFix" + stm.str();
		
		stm << "\n" << "Priors: {" << pset.toString() << "}\n";	
		stm << "Reference object SummaryStats = " 
			<< *(ref->getSummaryStatistic()) << "\n"	
			<< "Seed = " << seed << "\n";	
		std::string fileline = "SimdataFix" + stm.str();
		
		/* I am going to use mrs histograms on the outputs of
		 * parameter populations, but because of the
		 * difference in scale between our theta-per-site
		 * and our growth, if we use the histograms on a
		 * simulated population of parameter values directly 
		 * we will find that the growth dimension is chopped and the 
		 * theta-per-site dimension is not.  I implemented a
		 * way of standardising (centering, rescaling by standard
		 * deviation) but I am also doing a quick and dirty 
		 * rescaling to get both parameters into the same rough
		 * order of magnitude. Note the the mults are supplied in
		 * the same order as the parameter particles hold their
		 * values, ie theta first then growth.*/
		
		std::vector < double > mults;
		mults.push_back(nsites); // multiplier for theta output
		mults.push_back(1.0); // multiplier for growth output
		
		// specifiy a path - that directory had better exist or
		// there will be an error outputting the file. 
		std::string path("../output/");
		
		/* send the whole of the set of summary stats simulated with the
		 * 'true' parameter values to file for reference later */
		std::string filename_true = path + filenameStart + "RefStatSample.txt";
		boost::shared_ptr < SummaryStatisticSet > ref_ss_set = 
				refset->getSummaryStatisticSet();
		ref_ss_set->outputToFile(filename_true, false); // overwrite
		
		std::cout << "\nstarting AABC iterations, true parameters are " 
					<< true_p->toString() << endl;
	
		
		//iterate
		while (!epsilonMults.empty()) {
			
			std::ostringstream stm1;
			stm1.precision(3);
			stm1 <<   fixed;
			stm1 << fileline << "epsilonMult= " << epsilonMults.top();
			std::string thisfileline = stm1.str();
			
			std::cout << "\n" << thisfileline << std::endl;
			
			//This is the line that actually does the simulation
			simulator.iterate(epsilonMults.top(), iterations);
	
			// Not essential - just output to help us to know what's happening
			double epsilon = simulator.getCurrentEpsilon();
			std::cout << "epsilon = " << epsilon << std::endl;
			std::cout << std::endl;
			
			
			/* Get the current particles out of the simulator
			 * - I do this so that I can histogram each
			 * generation and show how the simulated population
			 * changes with each loop */
			boost::shared_ptr < ParameterParticleSet >
					particles = simulator.getCurrentParticles();
			
			/* I'm not actually interested in the summary statistics
			themselves - they are just a means to an end here */
			#if(0)
				boost::shared_ptr < SummaryStatisticSet >
					stats = simulator.getCurrentSummaryStatistics();
			#endif
			
			/* I want to send the particles to a file so that we can 
			 * input them later to an mrs histogram (or do something else
			 * with them ...) */
			
			std::string filename1 = path;
			std::string filename2 = path;
			std::string filename3 = path;
			std::string filename4 = path;
			
			// Construct filenames that tells me what this epsilon actually is
			std::ostringstream stm2;
			stm2.precision(3);
			stm2 <<   fixed;
			stm2 <<  "_eM_" << epsilonMults.top() 
				<< "_e_" << epsilon;
			filename1 = filename1 + filenameStart + stm2.str() +  ".txt";
			filename2 = filename2 + filenameStart 
					+ "Rescaled" + stm2.str()  +  ".txt";
			filename3 = filename3 + filenameStart 
					+ "Stndrd" + stm2.str() +  ".txt";
			filename4 = filename4 + filenameStart 
					+ "Stndrd" + stm2.str() +  "_Info.txt";
					
			// put a first line in the files
			mct_utilities::outputToFile(thisfileline, filename1, false); // overwrite file
			mct_utilities::outputToFile(thisfileline, filename2, false); // overwrite file
			mct_utilities::outputToFile(thisfileline, filename3, false); // overwrite file
			mct_utilities::outputToFile(thisfileline, filename4, false); // overwrite file
			
			
			/* and ask the current particles to output itself */
			particles->outputToFile(filename1, true); // append to the file
			/* and scaled - supply the mults here to rescale the output */
			particles->outputToFile(mults, filename2, true); // append to the file
			
			// or we can get the current particles standardised
			outputStandardisedParticles(particles, filename3, filename4, true_p);
			
			// pop the stack			
			epsilonMults.pop();
			
		}
		
	}
	catch (std::exception& e) {
		std::cout << "Exception " << e.what() << std::endl; 
	}
		
}


void outputStandardisedParticles(
	boost::shared_ptr < const ParameterParticleSet > particles,
	const std::string& filename,
	const std::string& filenameInfo,
	boost::shared_ptr < const ParameterParticle > true_p)
{
	/* make an analysis of the particles
	 * this just enables us to get the standardised values etc easily*/
	ParameterParticleSetAnalysis panalysis(particles);
	
	/* have a look at the means and sample sds while we are at it */
	boost::shared_ptr < const mct::SummaryStatistic > 
							means = panalysis.getMeans();
	std::vector<double> meansVec = means->getValues();
	
	/* have a look at the means and sample sds while we are at it */
	boost::shared_ptr < const mct::SummaryStatistic > 
							sds = panalysis.getSampleSDs();
	std::vector<double>sdsVec = sds->getValues();
	
	std::cout << "means are:" << std::endl;
			std::cout << means->plainString() << endl;
	std::cout << "sds are:" << std::endl;
			std::cout << sds->plainString() << endl;
	
	/* and the rmses against the true parameter value  */
	boost::shared_ptr < const mct::SummaryStatistic > 
							rmses = panalysis.getRMSEs(true_p);
	std::cout << "RMSEs are:" << std::endl;
			std::cout <<rmses->plainString() << endl;
	
	/* and get the standardised values */
	boost::shared_ptr < mct::SummaryStatisticSet > sss =
		panalysis.makeSummaryStatisticSetStandardisedParameterValues();
	
	/* and ask standardised values to output themselves */
	sss->outputToFile(filename, true); // append to the file
	
	// add to file filenameInfo for some info for the standardised values
	std::ostringstream stm1;
	std::ostringstream stm2;
	std::ostringstream stm3;
	for (size_t i = 0; i < meansVec.size(); ++i) {
			
		stm1 << "\t" << i;
		stm2 << "\t" << meansVec[i];
		stm3 << "\t" << sdsVec[i];
		
	}
	std::string headersLine("index");
	headersLine += stm1.str();
	std::string meansLine("means");
	meansLine += stm2.str();
	std::string sdsLine("sds");
	sdsLine += stm3.str();
	
	mct_utilities::outputToFile(headersLine, filenameInfo, true); // append to file
	mct_utilities::outputToFile(meansLine, filenameInfo, true); // append to file
	mct_utilities::outputToFile(sdsLine, filenameInfo, true); // append to file
					
} 

