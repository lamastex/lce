/* MCT - Markov Chains on Trees.

   Copyright (C) 2009 Brendan Bycroft <brb44@student.canterbury.ac.nz>
   Copyright (C) 2011 Jenny Harlow

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
\brief ARGFactory definitions.
*/

#include "mct_arg_factory.hpp"
#include "arg_reader.hpp"
#include "arg_ms_converter.hpp"
#include "libseq.hpp"

#include <stdexcept>
#include <sstream>
#include <limits> //for max double

#define MAXSEED 999999998 // max seed ms can be given
// anything larger and 10th pos is ignored, at least when args are read from cmd line

//#define MYDEBUG
#if defined (MYDEBUG)
	#include <iostream>
	
#endif

using namespace mct;

   
boost::shared_ptr < PopulationLabeledARG > 
				ARGFactory::makeHistFromLibseq(
				const boost::shared_ptr< hudson_ms::PopulationStructure > pop,
				const size_t nsites,
				const double growth,
				const double rho,
				const boost::shared_ptr < const PRNGen >& r_ptr)
{
	if (!(pop)) {
		throw std::invalid_argument("ARGFactory::makeHistFromLibseq(...) : pop pointer null");
	}
	
	if (pop->nSubPops() > 1) {
		throw std::invalid_argument("ARGFactory::makeHistFromLibseq(...) : > 1 sub-population");
	}

	boost::shared_ptr < Sequence::arg >
		hist(new Sequence::arg( LibSeq::initialise_hist(
						pop->totalNsam(),
						static_cast<int> (nsites),
						r_ptr,
						growth, // change in pop^n size
						0.0, // t_begin
						std::numeric_limits<double>::max(), // t_end
						rho, //rho
						-1)  //size_at_end
						) );
	
	boost::shared_ptr < PopulationLabeledARG > result 
		(new PopulationLabeledARG( hist, pop->makePopulationLabels()) );
	
	return result;
}

boost::shared_ptr < PopulationLabeledARG > 
				ARGFactory::makeHistFromMS(
				const boost::shared_ptr< const hudson_ms::PopulationStructure > pop,
				const size_t nsites,
				//const double growth,
				const double rho,
				const boost::shared_ptr < const PRNGen >& r_ptr)
{
	
	if (!(pop)) {
		throw std::invalid_argument("ARGFactory::makeHistFromMS(...) : pop pointer empty");
	}
	
	size_t max_seed = MAXSEED;
	
	//int seed = 1 + gsl_rng_uniform_int(r, max_seed);
	int seed = 1 + r_ptr->uniform_int( max_seed );
	
	boost::shared_ptr < Sequence::arg >  hist = ARGmsConverter::makeArg(
							nsites,
							rho,
							pop,
							seed);
	
	#ifdef MYDEBUG
		std::string cmd = pop->makeMScommand(rho, nsites, seed, growth);
		std::cout << "The ms command line would be" << std::endl;
		std::cout << cmd << std::endl;
	#endif
								
	boost::shared_ptr < PopulationLabeledARG > result
		(new PopulationLabeledARG( hist, pop->makePopulationLabels()) );
		
	return result;	

}


boost::shared_ptr < PopulationLabeledARG > 
				ARGFactory::makeHistFromNewickString(
				const std::string& newicktreestring,
				const boost::shared_ptr< hudson_ms::PopulationStructure > pop,
				const size_t nsites)
{
	if (!(pop)) {
		throw std::invalid_argument(
			"ARGFactory::makeHistFromNewickString(...) : pop pointer empty");
	}
	// make an arg reader
	ARGReader argReader;
	
	boost::shared_ptr < Sequence::arg > hist 
					= argReader.makeTreefileIntoArg(
										newicktreestring, nsites); 
	
	return boost::shared_ptr < PopulationLabeledARG >
		(new PopulationLabeledARG( 
							hist, pop->makePopulationLabels()) );

}

