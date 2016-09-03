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
mct::MultiLociSimDataSampler and get out summary stats to use in AHABC.

mct::MultiLociPolyTable objects represent multiple loci of nucleotide data
summarised as snp data.  
 
The aim here is to use some particular parameter values and to get a large 
number of mct::MultiLociPolyTable reps for any particular
set or combination of them (ie consider these a sample from the population
of possible mct::MultiLociPolyTable objects given those values), and get
the summary statistics for the reps.  

The summary statistics can be
further analysed using AHABC techniques.

mct::MultiLociPolyTable objects can be simulated in various ways. The
'simdata' sampler uses libsequence directly to do the simulation
of the snp data.  We could also use our own mct models (by using
a mct::MultiLociPolySitesSampler). For a simple
one-population population structure, the main difference is
that simdata sampling is much quicker.  If you wanted to do this for
a population of subpopulations, you would have to adapt this to 
use a mct::MultiLociPolySitesSampler instead.  

*/

#include "config.h"

#include "multi_loci_simdata_sampler.hpp"
#include "multi_loci_polytable_sampler.hpp"
#include "multi_loci_polytable_set.hpp"
#include "multi_loci_polytable.hpp"

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

void nucsRaaz();

void doMoreOutput(const std::string& filename2,
				const std::string& filename3,
				boost::shared_ptr < const SummaryStatisticSet > ss_set);

int main()
{
	nucsRaaz();
	
	
	return 0;
	
}

void nucsRaaz()
{
	try {
		int myPid = getpid();
		std::cout << "Nucleotide runs for Raaz, using simdata " << std::endl;
		std::cout << "This process id is " << myPid << std::endl;
		
		size_t nsam = 30; // samples
		
		size_t nsites = 25000; // number of sites in each locus 
		
		// number to generate with 'true' parameters
		size_t n_examples = 1000000;
		
		//number of loci in a multilocus
		size_t myloci[] = {10}; 
		
		std::vector <size_t> nls (myloci, myloci + sizeof(myloci) / sizeof(size_t) );
		
		double Nzero = 10000.0; // effective popn size (per sub-pop) 10^4
		
		double unit = 1000000000.0; // 10^9
		
		//mutation rates PER SITE before scaling by unit
		double mymus[] = {0.01, 0.1, 1.0, 10.0}; // 
		
		std::vector <double> mus (mymus, mymus + sizeof(mymus) / sizeof(double) );
		// convert the mus into thetas
		std::vector <double> thetas(mus.size());
		
		using namespace boost::lambda;
		//mults by 4*Nzero/unit
		transform(mus.begin(), mus.end(), thetas.begin(), _1*4*(Nzero/unit));
		
		//growth rates
		double mygrowths[] = {0.0,0.1,1.0,10.0,100.0};
		std::vector <double> growths (mygrowths, mygrowths + sizeof(mygrowths) / sizeof(double) );
		
		// overall locus sizes
		for (std::vector < size_t>::iterator nl_it = nls.begin(); 
				nl_it < nls.end();
				++ nl_it) {
					
			size_t nloci = *nl_it; 
			//vector of nloci copies of nsites
			std::vector < size_t > nloci_nsites(nloci, nsites);
		
			// over all thetas	
			for (std::vector < double>::iterator theta_it = thetas.begin(); 
					theta_it < thetas.end();
					++ theta_it) {
				
				double theta = *theta_it;
					
				//over all growths		
				for (std::vector < double>::iterator growth_it = growths.begin(); 
					growth_it < growths.end();
					++ growth_it) {
						
					int seed = 2345;
					boost::shared_ptr < PRNGen > r ( new PRNGenGSL(seed) );
					
					boost::shared_ptr < PRNGen > r_sampler = r->spawnAnother();

					double growth = *growth_it;
					
					// fit the precision to the magnitudes of the data
					int precTheta = mct_utilities::fitPrecision(theta);
					int precGrowth = mct_utilities::fitPrecision(growth);
					
					cout << "\n\nnloci " << nloci << " theta " << theta << " growth " << growth << endl;
					
					std::ostringstream stm;
					stm << "_si_" << nsites << "_ns_" << nsam << "_nl_" << nloci << "_n_" << n_examples;
					stm.precision(precTheta);
					stm <<   fixed;
					stm << "_th_" << theta;
					stm.precision(precGrowth);
					stm << "_g_" << growth;
					
					std::string path("../RaazOutput/");
					std::string filenameStart("RaazSimDataStats");
			
					std::string filename1 = path + filenameStart + stm.str() +  ".txt";
					std::string filename2 = path + filenameStart + "Stndrd" + stm.str() +  ".txt";
					std::string filename3 = path + filenameStart + "Stndrd" + stm.str() +  "_Info.txt";
					
					// sampler
					boost::shared_ptr < MultiLociPolyTableSampler > 
							sampler ( new MultiLociSimDataSampler(r_sampler) );
			
					double realvalues[] = {theta,growth}; 
					
					boost::shared_ptr < ParameterParticle > p( new ParameterParticle(realvalues, 2) );
			
					//population structure - just a single vanilla population
					boost::shared_ptr< hudson_ms::PopulationStructure > pop
					(new hudson_ms::PopulationStructure(nsam));
					
					// check the migration matrix (a bit unnecessary here!)
					std::string migCheckResult = pop->doMigCheck();
		
					if (!migCheckResult.empty()) {
						throw std::logic_error( migCheckResult);
					}
					
					/* can have memory problems - do in sections so that
					 * each section is written to file as we go */
					
					// do runs in sections of sectionSize at most
					size_t sectionSize = 500;
					size_t stillToDo = n_examples;
					bool firstRun = true;
					
					// overall container for all the available statistics for the dataset
					boost::shared_ptr < SummaryStatisticSet > all_ss_set (new SummaryStatisticSet());
					
					while (stillToDo > 0) {
						
						size_t thisRun = std::min(sectionSize, stillToDo);
						
						// sample
						boost::shared_ptr < mct::MultiLociPolyTableSet > 
							dataset = sampler->sample(thisRun, pop, nloci_nsites, p);
						
						// get all the available statistics for the dataset
						boost::shared_ptr < SummaryStatisticSet > ss_set = 
								dataset->getAllAvailableStatisticsSet();
						
						if (firstRun) {
							ss_set->outputToFile(filename1, false); // overwrite
						}
						else {
							ss_set->outputToFile(filename1, true); // append
						}
						
						// add to the overall container
						all_ss_set->add(ss_set);
						
						firstRun = false;
						stillToDo -= thisRun;
					} // end loop through sections
					
					/* the overall container should have collection of 
					 * all the summary statistics
					 * use a subroutine to do some more output
					 * for descriptive stats, standardised output etc */
					doMoreOutput(filename2, filename3, all_ss_set);
						
				}
			}
		}
			
	}
	catch (std::exception& e) {
	
		std::cout << "Error:\n" + std::string(e.what()) << std::endl;
	}
}

void doMoreOutput(const std::string& filename2,
				const std::string& filename3,
				boost::shared_ptr < const SummaryStatisticSet > ss_set)
{
	std::vector < double > means;
	std::vector < double > sds;
	
					
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

	// make the standardised values
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
