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
\brief Testing SummaryStatistic.

*  
*/

#include "config.h"

#include "summary_statistic_set.hpp"
#include "summary_statistic.hpp"
#include "descriptive_stats.hpp"

#include <iostream>
#include <vector>
#include <iterator>
#include <stdexcept>



using namespace std;
using namespace mct;


void summaryStatsTest();

int main()
{
	summaryStatsTest();
	
	
	return 0;
	
}

void testSub(boost::shared_ptr < SummaryStatisticSet > ss_set, bool consistent);	

void summaryStatsTest()
{
	try {
		int myPid = getpid();
		std::cout << "summary stats test " << std::endl;
		std::cout << "This process id is " << myPid << std::endl;
		
		boost::shared_ptr < SummaryStatisticSet > ss_set( new SummaryStatisticSet() ); 
			
		try {
			
			std::cout << "\n\n************ Testing with summary statistics of consistent size ************:" << std::endl;
			
			
			bool consistent = true;
			
			std::cout << "\nEmpty summary statistics set" << std::endl;
			
			testSub(ss_set, consistent);
			
					
			//mutation rates
			double mythetas[] = {1.0, 2.0, 3.0, 4.0}; //0.001,10.0,50.0,100.0
			
			
			//growth rates
			double mygrowths[] = {5.0, 10.0, 15.0, 20.0};
			
			assert( sizeof(mythetas) / sizeof(double) == (sizeof(mygrowths) / sizeof(double)) );
			assert( sizeof(mythetas) / sizeof(double) > 2 );
			
			std::cout << "\nOne summary statisic in set" << std::endl;
			
			for (size_t i = 0; i < 1 ; ++i ) {
			
		
					double theta = mythetas[i];
					double growth = mygrowths[i];
					
					double values[] = {theta,growth}; 
					
					boost::shared_ptr < SummaryStatistic > new_s
						( new SummaryStatistic(values, sizeof(values) / sizeof(double)) );
					ss_set-> add(new_s);
					
			}
			
			testSub(ss_set, consistent);
			
			
			std::cout << "\n Two summary statisics in set" << std::endl;
			
			for (size_t i = 1; i < 2; ++i ) {
			
		
					double theta = mythetas[i];
					double growth = mygrowths[i];
					
					double values[] = {theta,growth}; 
					
					boost::shared_ptr < SummaryStatistic > new_s
						( new SummaryStatistic(values, sizeof(values) / sizeof(double)) );
					ss_set-> add(new_s);
					
			}
			
			testSub(ss_set, consistent);
			
			std::cout << "\n" << (sizeof(mythetas) / sizeof(double)) << " summary statisics in set" << std::endl;
			
			
			for (size_t i = 2; i < sizeof(mythetas) / sizeof(double) ; ++i ) {
			
		
					double theta = mythetas[i];
					double growth = mygrowths[i];
					
					double values[] = {theta,growth}; 
					
					boost::shared_ptr < SummaryStatistic > new_s
						( new SummaryStatistic(values, sizeof(values) / sizeof(double)) );
					ss_set-> add(new_s);
					
			}
			
			testSub(ss_set, consistent);
			
			
		}
		catch (std::exception& e) {
				std::cerr << "Failed on consistent size tests, error is:\n"
				<< std::string(e.what()) << std::endl;
				throw;
		}
		
		try {
			std::cout << "\n\n************ Testing summary statistic set with summary statistics of inconsistent size ************:" << std::endl;
			
			boost::shared_ptr < SummaryStatisticSet > ss_set_new( new SummaryStatisticSet() ); 
		
			{
				double true_values[] = {5.0, 21.0, 100};
				boost::shared_ptr < SummaryStatistic > new_s
					( new SummaryStatistic(true_values, sizeof(true_values) / sizeof(double)) );
				ss_set_new->add(new_s);
			}
			
			{
				double true_values[] = {6.0}; 
				boost::shared_ptr < SummaryStatistic > new_s
					( new SummaryStatistic(true_values, sizeof(true_values) / sizeof(double)) );
				ss_set_new->add(new_s);
			}
			
			{
				double true_values[] = {7.0, 31.0, 101, 200.0}; 
				boost::shared_ptr < SummaryStatistic > new_s
					( new SummaryStatistic(true_values, sizeof(true_values) / sizeof(double)) ); 
				ss_set_new->add(new_s);
			}
			
			{
				double true_values[] = {8.0, 41.0}; 
				boost::shared_ptr < SummaryStatistic > new_s
					( new SummaryStatistic(true_values, sizeof(true_values) / sizeof(double)) ); 
				ss_set_new->add(new_s);
			}
					
			ss_set->add(ss_set_new);
			
			bool consistent = false;
			
			testSub(ss_set, consistent);
			
		}
		catch (std::exception& e) {
			std::cerr << "Failed on inconsistent size tests, error is:\n"
			<< std::string(e.what()) << std::endl;
			throw;
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


void testSub(boost::shared_ptr < SummaryStatisticSet > ss_set, bool consistent)
{
	
	{
		std::vector< std::vector< double > > stats;
	
		try {	
			std::cout << "Try to get all stats values" << std::endl;
			stats = *(ss_set->getAllValues());
			
		}
		catch (std::exception& e) {
			std::cout << "Failed to get stats values, error is:\n"
			<< std::string(e.what()) << std::endl;
			throw;
		}
			
		try {
			if (!stats.empty()) {	
				std::cout << "Try to print all stats values" << std::endl;
				for (std::vector< std::vector< double > >::iterator it = stats.begin();
				it < stats.end(); ++it) {
					copy(it->begin(), it->end(), ostream_iterator<double>(std::cout, "\t"));
					std::cout << std::endl;
				}
			}
			else std::cout << "stats values is empty:" << std::endl;
			
			
		}
		catch (std::exception& e) {
			std::cout << "Failed printing stats values, error is:\n"
			<< std::string(e.what()) << std::endl;
			throw;
		}
	}
	
	{
		std::vector< std::vector< double > > stats;
	
		try {	
			std::cout << "Try to get all stats values in summary stats layout" 
					<< std::endl;
			stats = *(ss_set->getAllValuesSumStatLayout());
			
		}
		catch (std::exception& e) {
			std::cout << "Failed to get stats values in summary stats layout, error is:\n"
			<< std::string(e.what()) << std::endl;
			throw;
		}
			
		try {
			if (!stats.empty()) {	
				std::cout << "Try to print all stats values in summary stats layout" << std::endl;
				for (std::vector< std::vector< double > >::iterator it = stats.begin();
				it < stats.end(); ++it) {
					copy(it->begin(), it->end(), ostream_iterator<double>(std::cout, "\t"));
					std::cout << std::endl;
				}
			}
			else std::cout << "stats values is empty:" << std::endl;
			
			
		}
		catch (std::exception& e) {
			std::cout << "Failed printing stats values in summary stats layout, error is:\n"
			<< std::string(e.what()) << std::endl;
			throw;
		}
	}
	
	std::vector< std::vector< double > > stats = *(ss_set->getAllValues());
	
	bool expectFailDS = (ss_set->size() < 1);
	bool expectFailStats = (ss_set->size() < 1 || !consistent);
	bool expectFailStandard = (!ss_set->empty()) && (ss_set->size() < 2 || !consistent);
	
	try {
		std::cout << "\nTry to get descriptive stats from summary stats" << std::endl;
		
		for (size_t i = 0; i < stats.size() ; ++i) {
				
			cout << "Descriptive stats for index " << i << endl;
			DescriptiveStats ds(ss_set, i);
			std::cout << (ds.toString()) << std::endl;
		}
	
		if (expectFailDS) throw std::logic_error("Should not be able to do that");
		
		try {
			int index = -1;
			cout << "Descriptive stats for index " << index << endl;
			DescriptiveStats ds(ss_set, index);
			
			throw std::logic_error("Should not be able to do that");
		}
		catch (std::exception& e) {
			std::cerr << "Failed on descriptive stats from summary stats, error is:\n"
			<< std::string(e.what()) << std::endl;
			
		}
		try {
			int index = static_cast<int> (stats.size());
			cout << "Descriptive stats for too-big index " << index << endl;
			DescriptiveStats ds(ss_set, index);
			
			throw std::logic_error("Should not be able to do that");
		}
		catch (std::exception& e) {
			std::cerr << "Failed on descriptive stats from summary stats, error is:\n"
			<< std::string(e.what()) << std::endl;
			
		}
	}
	catch (std::exception& e) {
		std::cerr << "Failed on descriptive stats from summary stats, error is:\n"
		<< std::string(e.what()) << std::endl;
		if (!expectFailDS) throw;
		
	}
	
	try {
		std::cout << "\nTry to get descriptive stats from values vector" << std::endl;
		
		for (size_t i = 0; i < stats.size() ; ++i) {
				
			cout << "Descriptive stats for index " << i << endl;
			DescriptiveStats ds(stats, i);
			std::cout << (ds.toString()) << std::endl;
		}
	
		if (expectFailDS) throw std::logic_error("Should not be able to do that");
		
		try {
			int index = -1;
			cout << "Descriptive stats for index " << index << endl;
			DescriptiveStats ds(stats, index);
			
			throw std::logic_error("Should not be able to do that");
		}
		catch (std::exception& e) {
			std::cerr << "Failed on descriptive stats from all values, error is:\n"
			<< std::string(e.what()) << std::endl;
			
		}
		try {
			int index = static_cast<int> (stats.size());
			cout << "Descriptive stats for too-big index " << index << endl;
			DescriptiveStats ds(stats, index);
			
			throw std::logic_error("Should not be able to do that");
		}
		catch (std::exception& e) {
			std::cerr << "Failed on descriptive stats from all values, error is:\n"
			<< std::string(e.what()) << std::endl;
			
		}
	}
	catch (std::exception& e) {
		std::cerr << "Failed on descriptive stats from values vector test, error is:\n"
		<< std::string(e.what()) << std::endl;
		if (!expectFailDS) throw;
		
	}
	
	try {
		std::cout << "\nTry to get means" << std::endl;
			
		boost::shared_ptr < const SummaryStatistic > means = ss_set->getSummaryStatisticMeans();
		std::cout << "means are:" << std::endl;
		std::cout << means->plainString() << endl;
		
		if (expectFailStats) throw std::logic_error("Should not be able to do that");
	}
	catch (std::exception& e) {
		std::cerr << "Failed on means test, error is:\n"
		<< std::string(e.what()) << std::endl;
		if (!expectFailStats) throw;
	}
	
	try {	
		std::cout << "Try to get sample variances" << std::endl;
		
		boost::shared_ptr < const SummaryStatistic > vars = ss_set->getSummaryStatisticVars();
		std::cout << "vars are:" << std::endl;
		std::cout << vars->plainString() << endl;
		
		if (expectFailStats) throw std::logic_error("Should not be able to do that");
		
		
	}
	catch (std::exception& e) {
		std::cerr << "Failed on vars test, error is:\n"
		<< std::string(e.what()) << std::endl;
		if (!expectFailStats) throw;
	}		
	
	try {
		std::cout << "Try to get all sample standard deviations" << std::endl;
		
		boost::shared_ptr < const SummaryStatistic > sds = 
						ss_set->getSummaryStatisticSDs();
		std::cout << "sds are:" << std::endl;
		std::cout << sds->plainString() << endl;
		
		if (expectFailStats) throw std::logic_error("Should not be able to do that");
	}
	catch (std::exception& e) {
		std::cerr << "Failed on sds test, error is:\n"
		<< std::string(e.what()) << std::endl;
		if (!expectFailStats) throw;
		
	}
	
	{
		std::vector< std::vector< double > > s_stats;
		
		try {	
			std::cout << "\nTry to get all standardised stats values in vectors" << std::endl;
			
			s_stats = *(ss_set->getAllValuesStandardised());
		
			if (!s_stats.empty()) {
				std::cout << "standardised stats values are:" << std::endl;
				for (std::vector< std::vector< double > >::iterator it = s_stats.begin();
				it < s_stats.end(); ++it) {
					copy(it->begin(), it->end(), ostream_iterator<double>(std::cout, "\t"));
					std::cout << std::endl;
				}
			}
			else std::cout << "standardised stats values is empty:" << std::endl;
			
			if (expectFailStandard) throw std::logic_error("Should not be able to do that");
			
			try {
				std::cout << "\nTry to get descriptive stats from standardised values vector" << std::endl;
				
				for (size_t i = 0; i < s_stats.size() ; ++i) {
						
					cout << "Descriptive stats for index " << i << endl;
					DescriptiveStats ds(s_stats, i);
					std::cout << (ds.toString()) << std::endl;
				}
			
				if (expectFailDS) throw std::logic_error("Should not be able to do that");
			}
			catch (std::exception& e) {
				std::cerr << "Failed on descriptive stats from standardised values vector test, error is:\n"
				<< std::string(e.what()) << std::endl;
				if (!expectFailDS) throw;
				
			}
			
		}
		catch (std::exception& e) {
			std::cout << "Failed to get standardised values, error is:\n"
			<< std::string(e.what()) << std::endl;
			if (!expectFailStandard) throw;
		}
		
		try {	
			std::cout << "Try to get all standardised stats as a summary stats set" << std::endl;
			
			boost::shared_ptr < SummaryStatisticSet > sss = 
						ss_set->makeSummaryStatisticSetStandardised();
		
			if (!sss->empty()) {
				std::cout << "standardised stats values from their summary stats are:" << std::endl;
				std::cout << (*sss) << std::endl;

			}
			else std::cout << "standardised stats summary stat is empty:" << std::endl;
			
			if (expectFailStandard) throw std::logic_error("Should not be able to do that");
			
			try {
				std::cout << "\nTry to get descriptive stats from standardised summary stats" << std::endl;
				
				for (size_t i = 0; i < s_stats.size() ; ++i) {
						
					cout << "Descriptive stats for index " << i << endl;
					DescriptiveStats ds(sss, i);
					std::cout << (ds.toString()) << std::endl;
				}
			
				if (expectFailDS) throw std::logic_error("Should not be able to do that");
			}
			catch (std::exception& e) {
				std::cerr << "Failed on descriptive stats from standardised summary stats test, error is:\n"
				<< std::string(e.what()) << std::endl;
				if (!expectFailDS) throw;
				
			}
		}
		
		catch (std::exception& e) {
			std::cout << "Failed to get standardised values, error is:\n"
			<< std::string(e.what()) << std::endl;
			if (!expectFailStandard) throw;
		}
	}
	
	{
		std::vector< std::vector< double > > s_stats;
		
		try {	
			std::cout << "\nTry to get all standardised stats values in vectors" 
				<< " in summary stats layout"  << std::endl;
			
			s_stats = *(ss_set->getAllValuesStandardisedSumStatLayout());
		
			if (!s_stats.empty()) {
				std::cout << "standardised stats values are:" << std::endl;
				for (std::vector< std::vector< double > >::iterator it = s_stats.begin();
				it < s_stats.end(); ++it) {
					copy(it->begin(), it->end(), ostream_iterator<double>(std::cout, "\t"));
					std::cout << std::endl;
				}
			}
			else std::cout << "standardised stats values is empty:" << std::endl;
			
			if (expectFailStandard) throw std::logic_error("Should not be able to do that");
			
		}
		catch (std::exception& e) {
			std::cout << "Failed to get standardised values in summary stats layout, error is:\n"
			<< std::string(e.what()) << std::endl;
			if (!expectFailStandard) throw;
		}
		
		
	}
	
	
	
}
