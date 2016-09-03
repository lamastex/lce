/* MCT - Markov Chains on Trees.


   Copyright (C) 2012 Jennifer Harlow

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
\brief Testing of prior lims routine.


*/

#include "testing_toolsAHABC.hpp"

#include "config.h" // mct

#include "utilities.hpp"

#include <mscplusplus/population_structure.hpp>

#include <boost/shared_ptr.hpp>

#include <iostream>
#include <vector>
#include <utility>
#include <cassert>
#include <ctime>

using namespace std;
using namespace mct;
using namespace aabc;
	

int main()
{
	int myPid = getpid();
	std::cout << "test prior lims" << std::endl;
	std::cout << "This process id is " << myPid << std::endl;
	
	size_t n_obs = 50; // number of 'observed values' to have
	
	size_t nsam = 10; // number of individuals in a sample
			
	// loci in each multi-loci SimData
	size_t nloci = 1; 
	
	size_t nsites = 25000; // number of sites in each locus 
	
	size_t Nzero = 1000000; // effective popn size (per sub-pop)
	double unit = 1000000000.0; // 10^9
	
	int seed = 2345;
	
	//true mutation rate
	//mutation rates PER SITE ie 4N0(mu/site)
	double mu = 5.0; //before scaling by dividing by unit
	double theta_per_site = mu*4*Nzero/unit;
	
	std::vector < double > trueParams (1, theta_per_site);
	
	//true growth rates
	double growth = 50.0;
	
	/* Create an empty population structure */
	boost::shared_ptr< hudson_ms::PopulationStructure > pop
	(new hudson_ms::PopulationStructure());
	/* Add a  simple one-subpop population structure for nsam samples, growth */
	pop->add(nsam, growth);
	
	
	std::vector < size_t > n_reps_vec; 
	{
		//size_t tmp[] = {10};
		size_t tmp[] = {10, 50, 100};
		n_reps_vec.insert (n_reps_vec.begin(), tmp, tmp+3);
		//n_reps_vec.insert (n_reps_vec.begin(), tmp, tmp+1);
	}
	std::vector < size_t > gridSegments_vec; 
	{
		size_t tmp[] = {10};
		//size_t tmp[] = {10, 20, 50};
		//gridSegments_vec.insert (gridSegments_vec.begin(), tmp, tmp+3);
		gridSegments_vec.insert (gridSegments_vec.begin(), tmp, tmp+1);
	}
	std::vector < double > percentClosest_vec; 
	{
		double tmp[] = {0.2};
		//double tmp[] = {0.1, 0.2, 0.3};
		//percentClosest_vec.insert (percentClosest_vec.begin(), tmp, tmp+3);
		percentClosest_vec.insert (percentClosest_vec.begin(), tmp, tmp+1);
	}
	
	for (size_t i = 0; i < percentClosest_vec.size(); ++i) {
		
		double percentClosest = percentClosest_vec[i];
		
		for (size_t j = 0; j < n_reps_vec.size(); ++j) {
		
			size_t n_reps_test = n_reps_vec[j];
			
			for (size_t k = 0; k < gridSegments_vec.size(); ++k) {
		
				size_t gridSegments = gridSegments_vec[k];
	
				{
					cout << "\n\nSmall heterozygosity only example, n_reps_test = " << n_reps_test
							<< ", gridSegments = " << gridSegments 
							<< ", percentClosest = " << percentClosest << endl;
							
					std::vector< std::pair<double, double> > priorLims;
					// prior limits for theta_per_site
					priorLims.push_back( pair<double, double>
						//I think I should use 0.50*4 etc as min here
						(0.50*4*Nzero/unit, 50.0*4*Nzero/unit));
				
				
					// only want to use heterozygosity stat
					std::vector < size_t > reqStats(1, 2);
					
					clock_t startTime = clock();
					
					getPriorLims(
						n_obs,
						nloci,
						nsites,
						pop,
						n_reps_test,	
						trueParams,
						priorLims,
						reqStats,
						seed,
						gridSegments,
						percentClosest);
					
					clock_t endTime = clock();
	
					double timing = ((static_cast<double>(endTime - startTime)) / CLOCKS_PER_SEC);
					
					cout << "\nTiming\t" << timing;
					
				}
				
				{
				
					cout << "\n\nSmall heterozygosity and segsites example, n_reps_test = " << n_reps_test
							<< ", gridSegments = " << gridSegments 
							<< ", percentClosest = " << percentClosest << endl;
						
					std::vector< std::pair<double, double> > priorLims;
					// prior limits for theta_per_site
					priorLims.push_back( pair<double, double>
							//I think I should use 0.50*4 etc as min here
							(0.50*4*Nzero/unit, 50.0*4*Nzero/unit));
					// prior limits for growth
					priorLims.push_back( pair<double, double>
						(30.00, 70.0));
					// all stats
					std::vector < size_t > reqStats;
					
					clock_t startTime = clock();
					
					getPriorLims(
						n_obs,
						nloci,
						nsites,
						pop,
						n_reps_test,	
						trueParams,
						priorLims,
						reqStats,
						seed,
						gridSegments,
						percentClosest);
					
					clock_t endTime = clock();
	
					double timing = ((static_cast<double>(endTime - startTime)) / CLOCKS_PER_SEC);
					
					cout << "\nTiming\t" << timing;
				}
				
				#if(1)
				{
					cout << "\n\nLarger heterozygosity and segsites example, n_reps_test = " << n_reps_test
							<< ", gridSegments = " << gridSegments 
							<< ", percentClosest = " << percentClosest << endl;
					
					std::vector< std::pair<double, double> > priorLims;
					priorLims.push_back( pair<double, double>
						//I think I should use 0.50*4 etc as min here
						(0.50*4*Nzero/unit, 50.0*4*Nzero/unit));
					// prior limits for growth
					priorLims.push_back( pair<double, double>
						(0.00, 100.0));
						
					// all stats
					std::vector < size_t > reqStats;
					
					clock_t startTime = clock();
					
					getPriorLims(
						n_obs,
						nloci,
						nsites,
						pop,
						n_reps_test,	
						trueParams,
						priorLims,
						reqStats,
						seed,
						gridSegments,
						percentClosest);
					
					clock_t endTime = clock();
	
					double timing = ((static_cast<double>(endTime - startTime)) / CLOCKS_PER_SEC);
					
					cout << "\nTiming\t" << timing;
				}
				#endif
			}
		}
	}
				
	
	
	return 0;
	
}


