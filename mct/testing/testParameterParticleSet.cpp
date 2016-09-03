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
\brief Testing ParameterParticleSet.

*  
*/

#include "config.h"

#include "parameter_particle_set.hpp"

#include <iostream>



using namespace std;
using namespace aabc;


void testParameterParticleSet();

int main()
{
	testParameterParticleSet();
	
	
	return 0;
	
}

	

void testParameterParticleSet() 
{
	gsl_matrix* cov = NULL;
	
	try {
		ParameterParticleSet pset;
		
		std::cout << "Try getting covariance (should fail)" << endl;
		try {
				cov = pset.covariance();
				
			}
		catch (exception& ee) {
			try {
				if (cov != NULL) gsl_matrix_free(cov);
				cov = NULL;
			}
			catch (std::exception& ee) {} // catch and swallow
			std::cout << "That resulted in exception\n" 
				+ std::string(ee.what()) << std::endl;
		}
			
			
	
		vector < boost::shared_ptr < const ParameterParticle > > toAdd;
		
		{
			double realvalues[] = {30.0, 60., 90.}; 
			boost::shared_ptr < ParameterParticle > new_p( new ParameterParticle(realvalues, 3));
			toAdd.push_back( new_p );
		}
		
		{
			double realvalues[] = {20.,30.,25.}; 
			boost::shared_ptr < ParameterParticle > new_p( new ParameterParticle(realvalues, 3));
			toAdd.push_back( new_p );
		}
		
		{
			double realvalues[] = {35.,30.,40.}; 
			boost::shared_ptr < ParameterParticle > new_p( new ParameterParticle(realvalues, 3));
			toAdd.push_back( new_p );
		}
		
		{
			double realvalues[] = {10.,15.,10.}; 
			boost::shared_ptr < ParameterParticle > new_p( new ParameterParticle(realvalues, 3));
			toAdd.push_back( new_p );
		}
		
		std::cout << "Try adding parameter particles" << endl;
		try {
				pset = ParameterParticleSet(toAdd);
		
				
		}
		catch (exception& ee) {
			std::cout << "That resulted in exception\n" 
				+ std::string(ee.what()) << std::endl;
				throw;
		}
		
		std::cout << "Try getting covariance (should work)" << endl;
		try {
				cov = pset.covariance();
				std::cout << "have cov, size1 = " << cov->size1 << endl;
		
		
				for (size_t i = 0; i < cov->size1; ++i) {
					for (size_t j = 0; j < cov->size2 ; ++j) {
						cout << gsl_matrix_get(cov, i, j) << "\t";
					}
					cout <<  endl;
				}
				cout <<  endl;
				try {
					if (cov != NULL) gsl_matrix_free(cov);
					cov = NULL;
				}
				catch (std::exception& ee) {} // catch and swallow
	
				
		}
		catch (exception& ee) {
			try {
				if (cov != NULL) gsl_matrix_free(cov);
				cov = NULL;
			}
			catch (std::exception& ee) {} // catch and swallow
			std::cout << "That resulted in exception\n" 
				+ std::string(ee.what()) << std::endl;
		}
	
	
		std::cout << "\n\nEnd of testing\n" << std::endl;
	
	}
	catch (std::exception& e) {
		
		try {
			if (cov != NULL) gsl_matrix_free(cov);
			cov = NULL;
		}
		catch (std::exception& ee) {} // catch and swallow
	
		std::cout << "std::exception::\n" + std::string(e.what()) << std::endl;
	}
	catch (...) {
		
		try {
			if (cov != NULL) gsl_matrix_free(cov);
			cov = NULL;
		}
		catch (std::exception& ee) {} // catch and swallow
	
		std::cout << "Unknown error\n" << std::endl;
	}
	
}		
