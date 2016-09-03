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
\brief Test our partitioned arg sequence generator with various models.


*/

#include "config.h"

#include "partitioned_arg_sequence_generator.hpp"
#include "mct_arg_factory.hpp"
#include "population_labeled_arg.hpp"

#include "model_set_config_builder.hpp"
#include "modelset.hpp"

#include "prng.hpp"
#include "prng_gsl.hpp"

#include "make_test_pops.hpp"

#include <mscplusplus/population_structure.hpp>
#include <Sequence/PolySites.hpp>

#include <boost/shared_ptr.hpp>

#include <iostream>
#include <iterator>

using namespace std;
using namespace mct;

void subGenSeqs(PartitionedARGSequenceGenerator& partitioned_arg,
		boost::shared_ptr < mct::ModelSet > models_ptr,
		double theta, size_t nsites,
		bool msatOnly, bool nucOnly);

void subPartitionedArgSequenceGenerator(
	const boost::shared_ptr< const hudson_ms::PopulationStructure > pop,
	size_t nsites, size_t nsam, double theta, double growth, double rho);
		
void testPartitionedArgSequenceGenerator();

int main()
{
	testPartitionedArgSequenceGenerator();
	
	
	return 0;
	
}

	
void testPartitionedArgSequenceGenerator()

{ 
	try {
		
		size_t nsites = 10;
		size_t nsam = 5; // will split this between 2 populations
		
		double theta = 10.0;
		double growth = 5.0;
		
		size_t Nzero = 1000000;
		
		// make a population structure 
		
		cout << "\n\nSetting up population structure" << endl;
		boost::shared_ptr< hudson_ms::PopulationStructure > 
						pop = makePopIslandBasic(Nzero);
		
		double rho = 0.0; // 4No r where r is cross proby over between sites in nsites 
		
		pop->setAllGrowthRates(growth);

		subPartitionedArgSequenceGenerator( pop, nsites, nsam, theta, growth, rho);

		std::cout << "\nEnd of testing\n" << endl;
	}
	catch (std::exception& e) {
	
		std::cout << "std::exception:\n" + std::string(e.what()) << std::endl;
	}
	catch (...) {
	
		std::cout << "Unknown error" << std::endl;
	}
	
	
}


void subPartitionedArgSequenceGenerator(
	const boost::shared_ptr< const hudson_ms::PopulationStructure > pop,
	size_t nsites, size_t nsam, double theta, double growth, double rho)

{ 
	cout << "nsites = " << nsites << " nsam = " << nsam << endl;
	
	int seed = 2345;
	boost::shared_ptr < PRNGen > r ( new PRNGenGSL(seed) );

	boost::shared_ptr < PRNGen > r_models = r->spawnAnother();
	
	// make one for the args
	boost::shared_ptr < PRNGen > r_args = r->spawnAnother();
	
	cout << "theta " << theta << " growth " << growth << ", rho " << ", nsites " << nsites << endl;
	
	// make a population structure 
	
	assert(pop->totalNsam() == nsam);
		
	// print migration matrix
	cout << "Pop migration matrix is " << endl;
	cout << pop->migrationMatrixToString() << endl;
	
	
	try {
		cout << "\n\nTry using empty population labeled arg pointer to make partitioned arg sequence generator " << endl;
		
		boost::shared_ptr < PopulationLabeledARG > empty_ptr;
	
		PartitionedARGSequenceGenerator partitioned_arg(nsites, empty_ptr);
		
	}
	catch (std::invalid_argument& ee) {
		std::cout << "That resulted in an invalid argument exception\n" 
			+ std::string(ee.what()) << std::endl;
	}
	
	cout << "\nTry using the population to make a population labeled arg pointer using mct::ARGFactory" << endl;
	
	boost::shared_ptr < PRNGen > rr = r_args->spawnAnother();
	
	boost::shared_ptr < PopulationLabeledARG > labeled_hist = 
		mct::ARGFactory::makeHistFromMS(pop, nsites, rho, rr);
	
	std::string labeled_arg = labeled_hist->toString();
	cout << "The labeled_arg is\n" << labeled_arg << endl; 
	
	cout << "\nTry making a PartitionedARGSequenceGenerator" << endl;
		
	PartitionedARGSequenceGenerator partitioned_arg(nsites, labeled_hist);
	
	// we can use the same PartitionedARGSequenceGenerator with lots of different models
	
	try {
		size_t wrong_nsites = nsites - 1;
		cout << "\n\nTry using this with models with wrong number of sites " << wrong_nsites << endl;
		
		// make some models
		/*use a ModelSetConfigBuilder to build the model data */
		mct::ModelSetConfigBuilder configBuilder;
		
		configBuilder.addModelConfig("SMM", wrong_nsites); // model_factory.cpp for options
		// create our models, one for each site, with default thetas
		boost::shared_ptr < mct::ModelSet > 
			models_ptr	( 
				new mct::ModelSet(configBuilder.getModelSetConfigs()
										, r_models)
						);
		
		bool msatOnly = true;
		bool nucOnly = false;
		
		subGenSeqs(partitioned_arg, models_ptr, theta, wrong_nsites,
				msatOnly, nucOnly);
		
	}
	catch (std::invalid_argument& ee) {
		std::cout << "That resulted in an invalid argument exception\n" 
			+ std::string(ee.what()) << std::endl;
	}
	
	try {
		cout << "\n\nTry using this with a microsat model with correct number of sites " << endl;
		
		// make some models
		/*use a ModelSetConfigBuilder to build the model data */
		mct::ModelSetConfigBuilder configBuilder;
		
		size_t half = nsites/2;
		if (half > 0) {
			
			configBuilder.addModelConfig("SMM", half); // model_factory.cpp for options
			configBuilder.addModelConfig("REV", nsites - half);
		}
		else {
			configBuilder.addModelConfig("SMM", nsites); // model_factory.cpp for options
			// create our models, one for each site, with default thetas
		}
		
		boost::shared_ptr < mct::ModelSet > 
			models_ptr	( 
				new mct::ModelSet(configBuilder.getModelSetConfigs()
										, r_models)
						);
		
		bool msatOnly = true;
		bool nucOnly = false;
		
		subGenSeqs(partitioned_arg, models_ptr, theta, nsites,
				msatOnly, nucOnly);
		
	}
	catch (std::invalid_argument& ee) {
		std::cout << "That resulted in an invalid argument exception\n" 
			+ std::string(ee.what()) << std::endl;
	}
	
	try {
		cout << "\n\nTry using this with a nucleotide model with correct number of sites " <<  endl;
		
		// make some models
		/*use a ModelSetConfigBuilder to build the model data */
		mct::ModelSetConfigBuilder configBuilder;
		
		
		configBuilder.addModelConfig("JC69", nsites); // model_factory.cpp for options
			// create our models, one for each site, with default thetas
				
		boost::shared_ptr < mct::ModelSet > 
			models_ptr	( 
				new mct::ModelSet(configBuilder.getModelSetConfigs()
										, r_models)
						);
		
		bool msatOnly = false;
		bool nucOnly = true;
		
		subGenSeqs(partitioned_arg, models_ptr, theta, nsites,
				msatOnly, nucOnly);
		
	}
	catch (std::invalid_argument& ee) {
		std::cout << "That resulted in an invalid argument exception\n" 
			+ std::string(ee.what()) << std::endl;
	}
	
	try {
		cout << "\n\nTry using this with a mixed nucleotide and microsat model with correct number of sites " <<  endl;
		
		// make some models
		/*use a ModelSetConfigBuilder to build the model data */
		mct::ModelSetConfigBuilder configBuilder;
		
		size_t half = nsites/2;
		if (half > 0) {
			
			size_t quarter = half/2;
			if (quarter > 0) {
				
				configBuilder.addModelConfig("SMM", quarter); // model_factory.cpp for options
				configBuilder.addModelConfig("JC69", half - quarter);
				configBuilder.addModelConfig("REV", quarter); // model_factory.cpp for options
				configBuilder.addModelConfig("JC69", nsites - half- quarter);
			}
			
			else {
				configBuilder.addModelConfig("SMM", half); // model_factory.cpp for options
				configBuilder.addModelConfig("JC69", half); // model_factory.cpp for options
			}
		
		
		
			boost::shared_ptr < mct::ModelSet > 
				models_ptr	( 
					new mct::ModelSet(configBuilder.getModelSetConfigs()
											, r_models)
							);
			
			bool msatOnly = false;
			bool nucOnly = false;
			
			subGenSeqs(partitioned_arg, models_ptr, theta, nsites,
					msatOnly, nucOnly);
		}
		// do nothing if we cant do multiple models	
	}
	catch (std::invalid_argument& ee) {
		std::cout << "That resulted in an invalid argument exception\n" 
			+ std::string(ee.what()) << std::endl;
	}
	
}


void subGenSeqs(PartitionedARGSequenceGenerator& partitioned_arg,
		boost::shared_ptr < mct::ModelSet > models_ptr,
		double theta, size_t nsites,
		bool msatOnly, bool nucOnly)
{

	//reset thetas in the models to be per site theta
	std::vector < double > new_thetas(models_ptr->size(),theta/nsites);
	models_ptr->resetThetas(new_thetas);

	boost::shared_ptr < mct::MCTSampleSequence > sample_seqs_ptr
					= partitioned_arg.createSampleSeqs(models_ptr);

	/* print the sequences of the sample in full form */
	std::cout << "\nsequence in raw format (with labels) is" << std::endl;
	std::cout << *sample_seqs_ptr << std::endl;


	/* print the sequences of the sample population in the phylip format */
	std::cout << "\nsequence in stringPhylipFormat() is" << std::endl;
	std::cout << sample_seqs_ptr->stringPhylipFormat();
	
	try {
		std::cout << "\ngetting sequence in Sequence::PolySites polyNucTableFormat() format" << std::endl;
		Sequence::PolySites polysites = 
								sample_seqs_ptr->polyNucTableFormat();
								
		if (!nucOnly) {
			throw std::logic_error("Should not be able to do that");
		}
		
		std::cout << "polysites format polysites.print(...) gives" << std::endl;
		polysites.print(cout);

	}
	catch (std::domain_error& de) {
		std::cout << "That resulted in a std::domain_error\n" 
			+ std::string(de.what()) << std::endl;
		if (nucOnly) throw;
	}
	
	try {
		std::cout << "\ngetting sequence in microsat repeats format" << std::endl;
		std::vector < site_seq > reps = 
								sample_seqs_ptr->microsatRepeatsFormat();
								
		if (!msatOnly) {
			throw std::logic_error("Should not be able to do that");
		}
		
		std::cout << "repeats format is" << std::endl;
		for (std::vector < site_seq >::iterator it
			= reps.begin();
			it < reps.end();
			++it) {
	
			cout << "( ";
			
			std::ostream_iterator<int> out_it (cout, "\t");
			copy ( it->begin(), it->end(), out_it );
		  
			cout << ")\n";
			cout << endl;
		}

	}
	catch (std::domain_error& de) {
		std::cout << "That resulted in a std::domain_error\n" 
			+ std::string(de.what()) << std::endl;
		if (msatOnly) throw;
	}
	
}
