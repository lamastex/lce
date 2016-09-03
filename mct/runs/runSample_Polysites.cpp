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
\brief Sample mct::MultiLociPolyTable objects using a 
mct::MultiLociPolySitesSampler.
 
mct::MultiLociPolyTable objects represent multiple loci of nucleotide data
summarised as snp data.  
 
The aim here is to use some particular parameter values and to get a large 
number of mct::MultiLociPolyTable reps for any particular
set or combination of them (ie consider these a sample from the population
of possible mct::MultiLociPolyTable objects given those values), and get
the summary statistics for the reps.  The summary statistics can be
further analysed to give an idea of the shape of their distribution.

mct::MultiLociPolyTable objects can be simulated in various ways. The
'polysites' sampler uses our mct models to simulate nucleotide 
sequences which are then turned into the snp data.  We could also use
libsequence simulate snp data directly.  For a simple
one-population population structure, the main difference is
that simdata sampling is much quicker (using mct
to get to the snp data is very time-consuming compared to libsequence 
if we do not need to take advantage of the flexibility that the mct 
method gives for using different models for different sites), 
but if you wanted to have
a population of subpopulations, you would have to use the polysites
sampler.  
* 
I have done an equivalent
of this program for the mct::MultiLociSimDataSampler in a similar
simple one-population case.

*/

#include "config.h"

#include "multi_loci_polysites_sampler.hpp"
#include "multi_loci_polytable_sampler.hpp"
#include "multi_loci_polytable_set.hpp"
#include "multi_loci_polytable.hpp"

#include "model_set_config_builder.hpp"
#include "modelset.hpp"

#include "descriptive_stats.hpp"
#include "summary_statistic_set.hpp"

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

void runSample_Polysites();

void doOutput(const std::string& filename1,
			const std::string& filename2,
			const std::string& filename3,
			boost::shared_ptr < const SummaryStatisticSet > ss_set);

int main()
{
	runSample_Polysites();
	
	
	return 0;
	
}

	

void runSample_Polysites()
{
	try {
		
		int myPid = getpid();
		std::cout << "Range of summary stats for MultiLociPolySitesSampler sampler " << std::endl;
		std::cout << "This process id is " << myPid << std::endl;
		
		/* number of reps to generate for each combination of possible parameter value of interest
		 * be careful - the polysites sampler takes a lot longer than the simdata sampler
		 * I have used a smaller number of reps here because of that */
		size_t n_examples = 10000;
		
		size_t nsam = 10; // samples
		
		size_t nsites = 25000; // number of sites in each locus 
		
		size_t Nzero = 1000000; // effective popn size (per sub-pop) 10^6
		
		/* Setting up to do the simulations for a number of different 
		 * number of loci, as well as different theta and growth 
		 */
		 
		size_t myloci[] = {5, 10, 20, 25}; 
		std::vector <size_t> nls (myloci, myloci + sizeof(myloci) / sizeof(size_t) );
		
		double unit = 1000000000.0; // 10^9
		// unit just helps with specifying mu_per_site
		// what is given here is mu/site * unit
		//mutation rates per site
		double mymus[] = {1.0, 2.5, 5.0, 7.5, 10.0}; 
		// 'real' mu/site is mu_per_site divided by unit
		
		std::vector <double> mus (mymus, mymus + sizeof(mymus) / sizeof(double) );
		std::vector <double> thetas(mus.size()); // going to put theta's in next...
		
		//mutation parameter theta PER SITE ie 4N0(mu/site)
		using namespace boost::lambda;
		//mults by 4*Nzero/unit
		transform(mus.begin(), mus.end(), thetas.begin(), _1*4*(Nzero/unit));
		
		//growth rates
		double mygrowths[] = {0.0,10.0,50.0,100.0};
		
		std::vector <double> growths (mygrowths, mygrowths + sizeof(mygrowths) / sizeof(double) );
		
		// loop through combinations of nloci, theta and growth ...
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
					
					// a new prng for each inner loop, each with the same seed	
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
					
					// extract the parameter value we are using for this inner loop
					double growth = *growth_it;
					
					cout << "\n\nnloci " << nloci << " theta " << theta << " growth " << growth << endl;
					
					/* Create a simple one-subpop population structure for nsam samples */
					boost::shared_ptr< hudson_ms::PopulationStructure > pop
					(new hudson_ms::PopulationStructure(nsam));
					
					/* if you had a more complex population you might want to check the migration
					 * matrix was okay - see the hudson_ms::PopulationStructure documentation */
					#if(0)
						// check the migration matrix
						std::string migCheckResult = pop.doMigCheck();
		
						if (!migCheckResult.empty()) {
							throw std::logic_error( migCheckResult);
						}
					#endif
					
					// create a string that summarises all our different values, to label output
					std::ostringstream stm;
					stm << "_ns_" << nsam << "_nl_" << nloci << "_n_" << n_examples;
					stm.precision(3);
					stm <<   fixed;
					stm << "_th_" << theta;
					stm << "_g_" << growth;
					
					std::string filenameStart("PolySitesStats");
					
					// specifiy a path - that directory had better exist or
					// there will be an error outputting the file. 
					std::string path("../output/");
					
					// filenames for output
					std::string filename1 = path + filenameStart + stm.str() +  ".txt";
					std::string filename2 = path + filenameStart + "Stndrd" + stm.str() +  ".txt";
					std::string filename3 = path + filenameStart + "Stndrd" + stm.str() +  "_Info.txt";
					
					/*use a ModelSetConfigBuilder to build the model data */
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
					
					/* Make the sampler for this inner loop: the sampler
					 * will create MultiLociPolytable objects from which
					 * we get the summary statistics.  The sampler
					 * is configured with the models and a prng  */
					boost::shared_ptr < MultiLociPolyTableSampler > 
								sampler ( new MultiLociPolySitesSampler(models_ptr, r_sampler) );
					
					/* make a parameter particle with the theta and growth values
					 * being used for this inner loop.  */
					double realvalues[] = {theta,growth}; 
					
					boost::shared_ptr < ParameterParticle > 
								p( new ParameterParticle(realvalues, 2) );
			
					/* Make a set of reps as a sample of mct::MultiLociPolyTable
					 * that can be generated using these parameter values.*/
					boost::shared_ptr < mct::MultiLociPolyTableSet > 
						refset = sampler->sample(n_examples, pop, nloci_nsites, p);
					
					/* and get the summary statistics for this set */
					boost::shared_ptr < SummaryStatisticSet > ss_set = 
							refset->getSummaryStatisticSet();
					
					//output
					doOutput(filename1, filename2, filename3, ss_set);
				
				} // end growth loop
			} // end theta loop
		} // end nloci loop
			
	} // end try (main routine)
	catch (std::exception& e) {
	
		std::cout << "Error:\n" + std::string(e.what()) << std::endl;
	}
}

void doOutput(const std::string& filename1,
				const std::string& filename2,
				const std::string& filename3,
				boost::shared_ptr < const SummaryStatisticSet > ss_set)
{
	std::vector < double > means;
	std::vector < double > sds;
	
	// output to a file (remember that filename includes path)
	ss_set->outputToFile(filename1, false); // overwrite
	
	/* The summary statistics set contains the particular
	 * statistics I was originally using for the aabc
	 * stuff - segsites and heterozygosity.  We can get
	 * them out directly but it's a weak point in using 
	 * the samplers for sampling the actual data and 
	 * summary stats directly that one has to just
	 * know what the order is (segsites first then heterozygosity)
	 * - this should be fixed up if we do more of this
	 * but there are also methods to get the stats
	 * by name as well so it's not that important */
					
	
	/* Use the DescriptiveStats class to analyse the summary stats 
	 * - see documentation for that class - it can 
	 * produce a quick summary with min,max, mean, median etc etc
	 * which we send directly to standard output here 
	 * - this can be very useful to get an idea of what
	 * is happening without having to use the histograms
	 * as well.  */				
	if (ss_set->size() > 0)	{
		std::cout << "Descriptive statistics for all summary stats" << endl;
		size_t nstats = ss_set->at(0)->size();
		for (size_t i = 0; i < nstats ; ++i) {
				
			cout << "For index " << i << endl;
			DescriptiveStats ds(ss_set, i);
			std::cout << (ds.toString()) << std::endl;
			// keep the means and sds in the vectors
			means.push_back(ds.mean());
			sds.push_back(ds.sampleSD());
		}
	}

	/* make the standardised values:
	 * standardising makes it easier to do the histograms
	 * because we don't get problems with scale 
	 * of seg sites being different to heterozygosity */
	boost::shared_ptr < SummaryStatisticSet > ss_set_stndrd = 
		ss_set->makeSummaryStatisticSetStandardised();
	
	// output standardised stats to file
	ss_set_stndrd->outputToFile(filename2, false); // overwrite
	
	// make a file that contains some info for the standardised values
	mct_utilities::outputToFile(filename3, filename3, false); // overwrite file
	std::ostringstream stm1;
	std::ostringstream stm2;
	std::ostringstream stm3;
	for (size_t i = 0; i < means.size(); ++i) {
			
		stm1 << "\t" << i;
		stm2 << "\t" << means[i];
		stm3 << "\t" << sds[i];
		
	}
	std::string headersLine("index");
	headersLine += stm1.str();
	std::string meansLine("means");
	meansLine += stm2.str();
	std::string sdsLine("sds");
	sdsLine += stm3.str();
	
	mct_utilities::outputToFile(headersLine, filename3, true); // append to file
	mct_utilities::outputToFile(meansLine, filename3, true); // append to file
	mct_utilities::outputToFile(sdsLine, filename3, true); // append to file
	
}
