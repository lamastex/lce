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

std::string makeFilename(	const std::string& path,
							const std::string& filenamestart,
							const std::string& filenamepop,
							const std::string& filenametype,
							const std::string& filenamedetails);

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
		size_t n_examples = 3;
		
		for (std::vector < double>::iterator theta_it = thetas.begin(); 
				theta_it < thetas.end();
				++ theta_it) {
			
			double theta = *theta_it;
					
			for (std::vector < double>::iterator growth_it = growths.begin(); 
				growth_it < growths.end();
				++ growth_it) {
					
				
				double growth = *growth_it;
				
				// do one set for each population type
				{
					boost::shared_ptr< hudson_ms::PopulationStructure > 
								pop = makePopSingle();
					
					
					std::string filenamepop = "PopSingle"; // alter this one
					
					testSub(pop, filenamepop, theta, growth, nloci, n_examples, nsites, nloci_nsites);
				}
				// pop with one subpop with no individuals
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
				// pop with 3 subpops, two individuals each
				{
					boost::shared_ptr< hudson_ms::PopulationStructure > pop = makePopIslandThreeSmall(Nzero);
					
					std::string filenamepop = "PopIslandThreeSmall"; // alter this one
					
					testSub(pop, filenamepop, theta, growth, nloci, n_examples, nsites, nloci_nsites);
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
	
	cout << "theta " << theta << " growth " << growth << endl;
	
	//population structure
	std::string popStr = pop->toString();
	std::cout << "Pop structure is " << popStr << endl;
	
	std::ostringstream stm;
	stm << "_ns_" << pop->totalNsam() << "_nl_" << nloci << "_n_" << n_examples;
	stm.precision(3);
	stm <<   fixed;
	stm << "_th_" << theta;
	stm << "_g_" << growth;
	
	std::string filenamedetails = stm.str();
	
	std::string filenamestart = "MicrosatTest";
	
	std::string filenametype1 = "Stats";
	std::string filenametype2 = "Means";
	std::string filenametype3 = "Covar";
	std::string filenametype4 = "DeltaMu";
	std::string filenametype5 = "RstsBetween";
	std::string filenametype6 = "PairwiseAvSumSqDiff";
	std::string filenametype7 = "FstsBetween";
	std::string filenametype8 = "Rsts";
	std::string filenametype9 = "Fsts";
	
	std::string path("../output/");


	std::string filename1 = makeFilename(path, filenamestart, filenamepop,
							filenametype1, filenamedetails);
	std::string filename2 = makeFilename(path, filenamestart, filenamepop,
							filenametype2, filenamedetails);
	std::string filename3 = makeFilename(path, filenamestart, filenamepop,
							filenametype3, filenamedetails);
	std::string filename4 = makeFilename(path, filenamestart, filenamepop,
							filenametype4, filenamedetails);
	std::string filename5 = makeFilename(path, filenamestart, filenamepop,
							filenametype5, filenamedetails);
	std::string filename6 = makeFilename(path, filenamestart, filenamepop,
							filenametype6, filenamedetails);
	std::string filename7 = makeFilename(path, filenamestart, filenamepop,
							filenametype7, filenamedetails);
	std::string filename8 = makeFilename(path, filenamestart, filenamepop,
							filenametype8, filenamedetails);
	std::string filename9 = makeFilename(path, filenamestart, filenamepop,
							filenametype9, filenamedetails);
	
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
	
	cout << "The multi-locus microsat set is\n" << *refset_ptr << endl;
	
	cout << "\nGetting population stats" << endl;
	
	
	std::vector < boost::shared_ptr < SummaryStatisticSet > > ss_sets;
	std::vector < std::string > filenames;
	
	try {
		ss_sets.push_back(	refset_ptr->getSummaryStatisticSet() );
		filenames.push_back(filename1);
	}
	catch (std::runtime_error& re) {
		cout << "Runtime error:\n" << re.what() << endl;
	}
	try {
		ss_sets.push_back(	refset_ptr->getAverageRepeats() );
		filenames.push_back(filename2);
	}
	catch (std::runtime_error& re) {
		cout << "Runtime error:\n" << re.what() << endl;
	}
	try {
		ss_sets.push_back(	refset_ptr->getCovariances() );
		filenames.push_back(filename3);
	}
	catch (std::runtime_error& re) {
		cout << "Runtime error:\n" << re.what() << endl;
	}
	
	try {
		ss_sets.push_back(	refset_ptr->getDistanceDeltaMuSqBetweenPop() );
		filenames.push_back(filename4);
	}
	catch (std::runtime_error& re) {
		cout << "Runtime error:\n" << re.what() << endl;
		throw;
	}
	try {
		ss_sets.push_back(	refset_ptr->getDistanceSlatkinRstBetweenPop() );
		filenames.push_back(filename5);
	}
	catch (std::runtime_error& re) {
		cout << "Runtime error:\n" << re.what() << endl;
		throw;
	}
	
	try {
		ss_sets.push_back(	refset_ptr->getPairwiseAverageSumSqDifferences() );
		filenames.push_back(filename6);
	}
	catch (std::runtime_error& re) {
		cout << "Runtime error:\n" << re.what() << endl;
		throw;
	}
	
	try {
		ss_sets.push_back(	refset_ptr->getBetweenPopFst() );
		filenames.push_back(filename7);
	}
	catch (std::runtime_error& re) {
		cout << "Runtime error:\n" << re.what() << endl;
		throw;
	}
	
	try {
		ss_sets.push_back(	refset_ptr->getRst() );
		filenames.push_back(filename8);
	}
	catch (std::runtime_error& re) {
		cout << "Runtime error:\n" << re.what() << endl;
		throw;
	}
	
	try {
		ss_sets.push_back(	refset_ptr->getFst() );
		filenames.push_back(filename9);
	}
	catch (std::runtime_error& re) {
		cout << "Runtime error:\n" << re.what() << endl;
		throw;
	}
	
	
	assert( ss_sets.size() == filenames.size() );
	
	cout << "Sending output to files" << endl;
	
	for (size_t i = 0; i < filenames.size(); ++i) {
		ss_sets[i]->outputToFile(filenames[i], false); // overwrite
	}
	
	cout << "\nRead back in from the files and send to standard output: " << endl;
	
	std::vector < boost::shared_ptr < SummaryStatisticSet > >::iterator sit = ss_sets.begin();
	for (std::vector < std::string >::iterator it = filenames.begin();
		it < filenames.end();
		++it, ++sit) {
		cout << "\nfile " << (*it) << endl;
		cout << (mct_utilities::string_from_file(*it)) << endl;
		cout << "\nand get descriptive statistics from population stats: " << endl;
		if (!(*sit)->empty() ) {
			for ( size_t index = 0; index < (*sit)->at(0)->size(); ++index ) {
				cout << "\tAt index " << index << "\t";
				DescriptiveStats ds((*sit), index);
				std::cout << (ds.toString()) << std::endl;
			}
		}
	}

}

std::string makeFilename(	const std::string& path,
							const std::string& filenamestart,
							const std::string& filenamepop,
							const std::string& filenametype,
							const std::string& filenamedetails)
{
	return path + filenamestart + filenamepop + filenametype + filenamedetails + std::string(".txt");
}
