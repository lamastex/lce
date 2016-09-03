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
\brief Testing ParameterParticleAnalysis.

*  
*/

#include "config.h"

#include "parameter_particle_set_analysis.hpp"
#include "parameter_particle_set.hpp"
#include "parameter_particle.hpp"

#include <iostream>
#include <vector>
#include <iterator>



using namespace std;
using namespace aabc;


void subTest(boost::shared_ptr < ParameterParticleSet > pset, 
				boost::shared_ptr < ParameterParticle > tru,
				bool consistent);

void paramAnalysisTest();

int main()
{
	paramAnalysisTest();
	
	
	return 0;
	
}

	
void paramAnalysisTest()
{
	try {
		int myPid = getpid();
		std::cout << "ParameterParticleAnalysis test " << std::endl;
		std::cout << "This process id is " << myPid << std::endl;
		
		double true_values[] = {2.5, 12.5}; 
		
		std::vector <double> values (true_values, true_values + sizeof(true_values) / sizeof(double) );
		
		boost::shared_ptr < ParameterParticle > tru
			( new ParameterParticle(values) );
			
			
		
		boost::shared_ptr < ParameterParticleSet > pset( new ParameterParticleSet() ); 
		
		
		{
			std::cout << "\n\n********** Testing with parameter particles analysis with empty parameter particles set *********:" << std::endl;
		
			bool consistent = true;
			
			subTest(pset, tru, consistent);
		}
			
		std::cout << "\n\n************ Testing with parameter particles analysis with consistent size ************:" << std::endl;
				
		//mutation rates
		double mythetas[] = {1.0, 2.0, 3.0, 4.0}; //0.001,10.0,50.0,100.0
		
		
		//growth rates
		double mygrowths[] = {11.0, 12.0, 13.0, 14.0};
		
		assert( sizeof(mythetas) / sizeof(double) == sizeof(mygrowths) / sizeof(double) );
		assert( sizeof(mythetas) / sizeof(double) > 2 );
		
		
		for (size_t i = 0; i < 1 ; ++i ) {
		
	
				double theta = mythetas[i];
				double growth = mygrowths[i];
				
				double values[] = {theta,growth}; 
				
				boost::shared_ptr < ParameterParticle > new_p
					( new ParameterParticle(values, sizeof(values) / sizeof(double)) );
				pset-> add(new_p);
				
		}
		
		
		{
			std::cout << "\nOne particle to analyse" << std::endl;
		
			bool consistent = true;
			
			subTest(pset, tru, consistent);
		}
		
		for (size_t i = 1; i < 2 ; ++i ) {
		
	
				double theta = mythetas[i];
				double growth = mygrowths[i];
				
				double values[] = {theta,growth}; 
				
				boost::shared_ptr < ParameterParticle > new_p
					( new ParameterParticle(values, sizeof(values) / sizeof(double)) );
				pset-> add(new_p);
				
		}
		
		
		{
			std::cout << "\nTwo particles to analyse" << std::endl;
		
			bool consistent = true;
			
			subTest(pset, tru, consistent);
		}
		
		for (size_t i = 2; i < (sizeof(mythetas) / sizeof(double)) ; ++i ) {
		
	
				double theta = mythetas[i];
				double growth = mygrowths[i];
				
				double values[] = {theta,growth}; 
				
				boost::shared_ptr < ParameterParticle > new_p
					( new ParameterParticle(values, sizeof(values) / sizeof(double)) );
				pset-> add(new_p);
				
		}
		
		
		{
			std::cout << "\nAll particles to analyse" << std::endl;
		
			bool consistent = true;
			
			subTest(pset, tru, consistent);
		}
		
		
		std::cout << "\n\n************ Testing parameter particle set with parameter particles of inconsistent size ************:" << std::endl;
		
		{
			double true_values[] = {5.0, 21.0, 100};
			boost::shared_ptr < ParameterParticle > new_p
					( new ParameterParticle(true_values, sizeof(true_values) / sizeof(double)) );
			pset-> add(new_p); 
		}
		
		{
			double true_values[] = {6.0}; 
			boost::shared_ptr < ParameterParticle > new_p
					( new ParameterParticle(true_values, sizeof(true_values) / sizeof(double)) );
			pset-> add(new_p);
		}
		
		{
			double true_values[] = {7.0, 31.0, 101, 200.0}; 
			boost::shared_ptr < ParameterParticle > new_p
					( new ParameterParticle(true_values, sizeof(true_values) / sizeof(double)) );
			pset-> add(new_p); 
		}
		
		{
			double true_values[] = {8.0, 41.0}; 
			boost::shared_ptr < ParameterParticle > new_p
					( new ParameterParticle(true_values, sizeof(true_values) / sizeof(double)) );
			pset-> add(new_p); 
		}
				
		
		std::vector< std::vector< double > > parameters;
		try {
			parameters = *(pset->getAllValues());
		}
		catch (std::exception& e) {
			std::cerr << "Failed in pset->getAllValues(), error is:\n"
			<< std::string(e.what()) << std::endl;
		}
		
		try {	
			std::cout << "parameter values are:" << std::endl;
			for (std::vector< std::vector< double > >::iterator it = parameters.begin();
			it < parameters.end(); ++it) {
				copy(it->begin(), it->end(), ostream_iterator<double>(std::cout, "\t"));
				std::cout << std::endl;
			}
		}
		catch (std::exception& e) {
			std::cerr << "Failed printing parameter values, error is:\n"
			<< std::string(e.what()) << std::endl;
		}
		
		
		std::cout << "\n\n************ Testing parameter particle set analysis with parameter particles of inconsistent size ************:" << std::endl;
		
		
		
		{
			cout << "The constructor call should fail" << endl;
		
			bool consistent = false;
			
			subTest(pset, tru, consistent);
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
	

void subTest(boost::shared_ptr < ParameterParticleSet > pset, 
				boost::shared_ptr < ParameterParticle > tru,
				bool consistent)
{
	std::cout << "about to make analysis" << std::endl;
	
	bool canMakePPA = true;
	
	try {		
		ParameterParticleSetAnalysis tmp(pset);
		
		if (!consistent) throw std::logic_error("Should not be able to do that");
	
	}
	catch (std::exception& e) {
		std::cout << "Failed to make analysis, error is:\n"
		<< std::string(e.what()) << std::endl;
		if (consistent && (!pset->empty())) throw;
		else canMakePPA = false;
	}
	
	if (canMakePPA) {
	
		ParameterParticleSetAnalysis panalysis(pset);
		
		std::vector< std::vector< double > > parameters;
		
		try {
			parameters = *(panalysis.getAllParameters());
		
		}
		catch (std::exception& e) {
			std::cout << "Failed to get all values, error is:\n"
			<< std::string(e.what()) << std::endl;
			throw;
		}
		
		try {	
			if (!parameters.empty()) {	
				std::cout << "parameter values are:" << std::endl;
				for (std::vector< std::vector< double > >::iterator it = parameters.begin();
				it < parameters.end(); ++it) {
					copy(it->begin(), it->end(), ostream_iterator<double>(std::cout, "\t"));
					std::cout << std::endl;
				}
			}
			else std::cout << "all values is empty:" << std::endl;
		}
		catch (std::exception& e) {
			std::cerr << "Failed printing parameter values, error is:\n"
			<< e.what() << std::endl;
			throw;
		}
		
		
		bool expectFail = (panalysis.nParticles() < 1);
		
		try {	
			boost::shared_ptr < const mct::SummaryStatistic > means = panalysis.getMeans();
			
			if (expectFail) throw std::logic_error("Should not be able to do that");
			
			std::cout << "means are:" << std::endl;
			std::cout << means->plainString() << endl;
		}
		catch (std::exception& e) {
			std::cerr << "Failed on means test, error is:\n" << 
			e.what() << std::endl;
			if (!expectFail) throw;
			
		}
		
		try {	
			boost::shared_ptr < const mct::SummaryStatistic > vars = panalysis.getSampleVars();
			
			if (expectFail) throw std::logic_error("Should not be able to do that");
			
			std::cout << "vars are:" << std::endl;
			std::cout << vars->plainString() << endl;
		}
		catch (std::exception& e) {
			std::cerr << "Failed on vars test, error is:\n" << 
			e.what() << std::endl;
			if (!expectFail) throw;
			
		}		
		
		try {
			boost::shared_ptr < const mct::SummaryStatistic > sds = 
							panalysis.getSampleSDs();
			
			if (expectFail) throw std::logic_error("Should not be able to do that");
			
			std::cout << "sds are:" << std::endl;
			std::cout << sds->plainString() << endl;
		}
		catch (std::exception& e) {
		std::cerr << "Failed on sds test, error is:\n" << 
			e.what() << std::endl;
			if (!expectFail) throw;
		}
		
		try {	
			boost::shared_ptr < const mct::SummaryStatistic > RMSEs = 
						panalysis.getRMSEs(tru); 
			
			if (expectFail) throw std::logic_error("Should not be able to do that");
			
			std::cout << "true values is " << (*tru) << endl;
			std::cout << "RMSEs against true values are:" << std::endl;
			std::cout << RMSEs->plainString() << endl;
			
		}
		catch (std::exception& e) {
			std::cerr << "Failed on RMSEs test, error is:\n" << 
			e.what() << std::endl;
			if (!expectFail) throw;
		}
		
		bool expectFailStan = (panalysis.nParticles() < 2);
		try {	
			boost::shared_ptr < std::vector < std::vector < double > > >
					s_vals =	panalysis.getAllParametersStandardised();
			
			if (expectFailStan) throw std::logic_error("Should not be able to do that");
			
			std::cout << "All standardised values as vectors are:" << std::endl;
			for (std::vector< std::vector< double > >::iterator it = s_vals->begin();
			it < s_vals->end(); ++it) {
				copy(it->begin(), it->end(), ostream_iterator<double>(std::cout, "\t"));
				std::cout << std::endl;
			}
			
		}
		catch (std::exception& e) {
			std::cerr << "Failed on standardised values as vector test, error is:\n" << 
			e.what() << std::endl;
			if (!expectFailStan) throw;
		}
		try {	
			boost::shared_ptr < mct::SummaryStatisticSet > sss =
				panalysis.makeSummaryStatisticSetStandardisedParameterValues();
			
			if (expectFailStan) throw std::logic_error("Should not be able to do that");
			
			std::cout << "All standardised values as summary statistic set are:" << std::endl;
			std::cout << (*sss) << std::endl;
			
		}
		catch (std::exception& e) {
			std::cerr << "Failed on standardised values as summary statistic set test, error is:\n" << 
			e.what() << std::endl;
			if (!expectFailStan) throw;
		}
	}
}
