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
mct::MultiLociSimDataSampler.

mct::MultiLociPolyTable objects represent multiple loci of nucleotide data
summarised as snp data.  
 
The aim here is to use some particular parameter values and to get a large 
number of mct::MultiLociPolyTable reps for any particular
set or combination of them (ie consider these a sample from the population
of possible mct::MultiLociPolyTable objects given those values), and get
the summary statistics for the reps.  The summary statistics can be
further analysed to give an idea of the shape of their distribution.

mct::MultiLociPolyTable objects can be simulated in various ways. The
'simdata' sampler uses libsequence directly to do the simulation
of the snp data.  We could also use our own mct models (by using
a mct::MultiLociPolySitesSampler). For a simple
one-population population structure, the main difference is
that simdata sampling is much quicker.  If you wanted to do this for
a population of subpopulations, you would have to adapt this to 
use a mct::MultiLociPolySitesSampler instead.  
* 
I have done an equivalent
of this program for the mct::MultiLociPolySitesSampler in a similar
simple one-population case.

*/

#include "config.h"

#include "mct/multi_loci_simdata_sampler.hpp"
#include "mct/multi_loci_polytable_sampler.hpp"
#include "mct/multi_loci_polytable_set.hpp"
#include "mct/multi_loci_polytable.hpp"

#include "mct/descriptive_stats.hpp"
#include "mct/summary_statistic_set.hpp"

#include "mct/prng.hpp"
#include "mct/prng_gsl.hpp"

#include "mct/utilities.hpp"

#include "aabc/parameter_particle.hpp"


#include "make_run_pops.hpp"

#include <mscplusplus/population_structure.hpp>
#include <mscplusplus/recomb_params.hpp>
#include <mscplusplus/segsite_params.hpp>

#include <boost/lambda/bind.hpp>
#include <boost/lambda/lambda.hpp>

#include <boost/shared_ptr.hpp>


#include <iostream>
#include <sstream>
#include <vector>
#include <stack>
#include <iomanip>

using namespace std;
using namespace mct;
using namespace aabc;

void runNucScenarios();

void doOutput(const std::string& filename1,
			const std::string& filename2,
			const std::string& filename3,
			boost::shared_ptr < const SummaryStatisticSet > ss_set);

int main()
{
	runNucScenarios();
	
	
	return 0;
	
}

	

void runNucScenarios()
{
	try {
		
		int myPid = getpid();
		std::cout << "Sharyn runs: nucleotide simulations: population scenarios" << std::endl;
		std::cout << "This process id is " << myPid << std::endl;

		// number of reps to generate for each combination of possible parameter value of interest
		size_t n_examples = 1000000;
		//size_t n_examples = 10;

		size_t nsites = 600; // number of sites in each locus 
		
		/* Population size N0 for each sub-population is 1000 haploids, 
		 * ie 500 diploids */ 
		size_t nZero = 500; // effective diploid popn size
		
		/* Setting up to do the simulations for a number of different 
		 * number of loci, as well as different theta and growth 
		 */
		 
		size_t myloci[] = {1}; 
		std::vector <size_t> nls (myloci, myloci + sizeof(myloci) / sizeof(size_t) );
		
		double unit = 1000000.0; // 10^6
		// unit just helps with specifying mu_per_site
		// what is given below is mu/site * unit
		//mutation rates per site
		double mymus[] = {2.0}; 
		// 'real' mu/site is mu_per_site divided by unit
		
		std::vector <double> mus (mymus, mymus + sizeof(mymus) / sizeof(double) );
		std::vector <double> thetas(mus.size()); // going to put theta's in next...
		
		//mutation parameter theta PER SITE ie 4N0(mu/site)
		using namespace boost::lambda;
		//mults by 4*nZero/unit
		transform(mus.begin(), mus.end(), thetas.begin(), _1*4*(nZero/unit));
		
		/* set up populations.
		 * will need a pointer to a population and a name for each one
		 */
		
		std::vector < boost::shared_ptr< hudson_ms::PopulationStructure > > pops;
		std::vector < std::string > popnames;
		
		pops.push_back ( makePopScenario1(nZero) );
		popnames.push_back("Scenario1");
		pops.push_back ( makePopScenario2(nZero) );
		popnames.push_back("Scenario2");
		pops.push_back ( makePopScenario3(nZero) );
		popnames.push_back("Scenario3");
		pops.push_back ( makePopScenario4(nZero) );
		popnames.push_back("Scenario4");
		pops.push_back ( makePopScenario4contrast(nZero) );
		popnames.push_back("Scenario4contrast");
		pops.push_back ( makePopScenario5(nZero) );
		popnames.push_back("Scenario5");

		
		// specifiy a path - that directory had better exist or
		// there will be an error outputting the file. 
		std::string path("../sharyn_output/");
					
		// loop through pop and pop names together
		std::vector < std::string >::iterator name_it = popnames.begin();
		for (std::vector < boost::shared_ptr< hudson_ms::PopulationStructure > >::iterator
					pop_it = pops.begin(); 
					pop_it < pops.end();
					++pop_it, ++name_it) {
			
			/* population for this loop */
			boost::shared_ptr< hudson_ms::PopulationStructure > pop_ptr = *pop_it;
					
			std::string filenameStart = *name_it;
			
			cout << "\n\npopulation " << filenameStart << endl;
			// Get the pop to give description of itself */
			std::cout << "Pop structure is " << ( pop_ptr->toString() ) << endl;
						
			// loop through combinations of nloci, theta ...
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
					
					cout << "\nnloci " << nloci << " theta per site " << theta << endl;
					
					// a new prng for each inner loop, each with the same seed	
					int seed = 2345;
					
					hudson_ms::HudsonMSRecombParams recombParams;
					hudson_ms::HudsonMSSegSiteParams ssParams(theta*nsites);
					std::cout << "Equivalent ms command (with some particular seed) is " << endl;
					std::cout << pop_ptr->makeMScommand(recombParams, ssParams, seed) << std::endl;
			
					boost::shared_ptr < PRNGen > r ( new PRNGenGSL(seed) );
					/*
					 * use this prng to spawn another for sampler .
					 * This ensures that results are independent of how
					 * internals work, and we can increase samples
					 * and first results will be the same as before
					 * the increase
					 * */
					boost::shared_ptr < PRNGen > r_sampler = r->spawnAnother();

					
					// create a string that summarises all our different values, to label output
					std::ostringstream stm;
					stm << "_N0_" << nZero << "_si_" << nsites << "_nl_" << nloci;
					stm.precision(4);
					stm <<   fixed;
					stm << "_th_" << theta  << "_n_" << n_examples;
					
					
					// filenames for output
					std::string filename1 = path + filenameStart + stm.str() +  ".txt";
					std::string filename2 = path + filenameStart + "Stndrd" + stm.str() +  ".txt";
					std::string filename3 = path + filenameStart + "Stndrd" + stm.str() +  "_Info.txt";
					
					/* Make the sampler for this inner loop: the sampler
					 * will create MultiLociPolytable objects from which
					 * we get the summary statistics.  The sampler
					 * is configured with its prng */
					boost::shared_ptr < MultiLociPolyTableSampler > 
							sampler ( new MultiLociSimDataSampler(r_sampler) );
			
					/* make a parameter particle with the theta and growth values
					 * being used for this inner loop.  */
					double realvalues[] = {theta}; 
					boost::shared_ptr < ParameterParticle > p( new ParameterParticle(realvalues, 1) );
			
					/* Make a set of reps as a sample of mct::MultiLociPolyTable
					 * that can be generated using these parameter values.*/
					boost::shared_ptr < mct::MultiLociPolyTableSet > 
						refset = sampler->sample(n_examples, pop_ptr, nloci_nsites, p);
					
					/* and get the summary statistics for this set */
					boost::shared_ptr < SummaryStatisticSet > ss_set = 
							refset->getSummaryStatisticSet();
					
					//output
					doOutput(filename1, filename2, filename3, ss_set);
									
				} // end theta loop
			} // end nloci loop
		} // end pops loop
			
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
		std::cout << "\nDescriptive statistics for all summary stats" << endl;
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
	
	//use more precision for outputting the means and sds info
	int prec = 6;
	stm2 << fixed << setprecision(prec);
	stm3 << fixed << setprecision(prec);
	
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
