/* MCT - Markov Chains on Trees.

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
\brief PopMCabcSimulatorAdaptive definitions.
*/


#include "pop_mc_abc_simulator_adaptive.hpp"
#include "parameter_prior.hpp"
#include "parameter_particle_validation_result.hpp"
#include "summary_statistic.hpp"

#include "rmv.h" // multivariate normal density and random values

#include <boost/function.hpp>
#include <boost/bind.hpp>

#include <algorithm> // for partition and for_each
#include <cmath>
#include <ctime> // for timestamp
//#include <iostream>
#include <fstream> //file stream
#include <iomanip> //change io precision etc
#include <stdexcept>
#include <cassert>


#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>


#define DEBUG	
//#define MYDEBUG
//#define MYDEBUG_EXTRA
//#define MYDEBUG_COMPARE

#if defined (DEBUG) || defined (MYDEBUG) || defined (MYDEBUG_EXTRA) || defined (MYDEBUG_COMPARE)
	#include <iostream>
	
#endif


using namespace aabc;

PopMCabcSimulatorAdaptive::PopMCabcSimulatorAdaptive(
			const ParameterPriorSet& ps,
			const boost::shared_ptr < ParameterParticleValidator >& v,
			const double bEP,
			const boost::shared_ptr < mct::PRNGen >& _r_ptr) :
				priors(ps),
				n_params(ps.size()),
				validator(v),
				baseEpsilonPercentile(bEP),
				epsilon(0.0), // initialise to 0.0
				lastAcceptanceRatio(bEP) // initialise to base epsilon percentile
{
	particles = boost::shared_ptr < ParameterParticleSet >(
							new ParameterParticleSet() );
	summaryStatistics = boost::shared_ptr < mct::SummaryStatisticSet >(
							new mct::SummaryStatisticSet() );

	if (0 == n_params) {
		throw std::invalid_argument(
			"PopMCabcSimulatorAdaptive::PopMCabcSimulatorAdaptive(...) : priors empty");
	}
	if (!v) {
		throw std::invalid_argument(
			"PopMCabcSimulatorAdaptive::PopMCabcSimulatorAdaptive(...) : validator null");
	}
	if (bEP <= 0) {
		throw std::invalid_argument(
			"PopMCabcSimulatorAdaptive::PopMCabcSimulatorAdaptive(...) : base epsilon percentile <= 0" );
	}
	if (bEP > 1) {
		throw std::invalid_argument(
			"PopMCabcSimulatorAdaptive::PopMCabcSimulatorAdaptive(...) : base epsilon percentile > 1" );
	}

	initialisePRNGs(_r_ptr);
}

		
PopMCabcSimulatorAdaptive::~PopMCabcSimulatorAdaptive() 
{
		//std::cout << "PopMCabcSimulatorAdaptive destructor called" << std::endl;
}


		
void PopMCabcSimulatorAdaptive::iterate(const double epsilonMult,
								const size_t n_reps,
								const std::string& logfilename)
{
	gsl_matrix *tau = NULL;
	
	try {
		
		if (epsilonMult <= 0.0) {
			throw std::invalid_argument(
				"PopMCabcSimulatorAdaptive::PopMCabcSimulatorAdaptive(...) : epsilonMult is <= 0.0");
		}
		
		if (n_reps < 3) {
			throw std::invalid_argument(
				"PopMCabcSimulatorAdaptive::PopMCabcSimulatorAdaptive(...) : n_reps < 3");
		}
		
		clock_t startIt = clock();
		
		boost::shared_ptr< ParameterParticleValidationResultSet > vrs;
		
		size_t attempts = 0; // track how many attempts we make
		
		// are we initialising or iterating from existing particles?
		bool initialising = particles->empty();	
		
		#ifdef DEBUG
			if (initialising) 
				std::cout << "Running with DEBUG defined\n" << std::endl;
		#endif

		/* If we are not initialising, need to calculate tau*/
		 if (!initialising) {
			tau = calcTau();
			#ifdef DEBUG
				std::cout << "tau is: " << std::endl;
				for (size_t j = 0; j < tau->size1; ++j)  {
					for (size_t i = 0; i < tau->size2; ++i) {
						std::cout << gsl_matrix_get (tau, i, j) << "\t";
					}
					std::cout << std::endl;
				}
				
			#endif
			
			/* and also to calculate the epsilon for this iteration*/
			updateEpsilon(epsilonMult);
			
			#ifdef DEBUG
				std::cout << "Using epsilon = " << epsilon << std::endl;
			#endif
		}
		
		size_t current_size = 0;
		
		/* constructor sets lastAcceptanceRatio
		to the baseEpsilonPercentage so if we are initialising, 
		and epsilonMult = 1, the acceptance ratio will be the 
		percentage that the baseEpsilonPercentage says we should
		pass first time around.  If we are not initialsing we 
		make an initial guess of the acceptance ratio based on the last
		one and the epsilonMult.*/
		double acceptanceRatio = lastAcceptanceRatio * epsilonMult; 
		
		while (current_size < n_reps) {
			
			/* we keep track of the acceptance ratio so that we
			 can estimate how many particles we need to generate
			 to get the required number of acceptances, but we keep 
			 looping and updating the acceptance ratio until we
			 have the required number of acceptances.
			 In the initialisation step, /*/
			size_t to_add = static_cast<size_t> 
					(std::ceil((n_reps - current_size)
						/ acceptanceRatio) );
			#ifdef MYDEBUG
				std::cout << "current size = " << current_size 
				<< " and acc_ratio = " << acceptanceRatio 
				<< " and to_add = " << to_add 
				<< std::endl;
			#endif
			
			
			boost::shared_ptr < ParameterParticleSet > 
						newParticles;
			
			/* If this is the first iteration, we need to do the
			 initial particle generation, otherwise we get
			 new proposals from the existing particles.*/			
			if (initialising) {
				#ifdef MYDEBUG
					std::cout << "initialising " << to_add << std::endl;
				#endif
						
				newParticles  = initialParticlePtrs(to_add);
				#ifdef MYDEBUG_COMPARE
					{
						std::cout << "First five parameter particles:" << std::endl;
						for (size_t i = 0; i < 5; ++i) {
							std:: cout << *(newParticles->at(i)) << std::endl;
						}
						std::cout << std::endl;
		
					}
					#endif
			}
			else {
				#ifdef MYDEBUG
					std::cout << "new particles " << to_add << std::endl;
				#endif
				newParticles = newParticlePtrs(to_add, tau);
			}
			
			
			/* if this is the first time through for this iteration
			 then set the normaliser in the validator (normalises
			 distance calculations), otherwise just get results
			 with distances normalised with existing normaliser.*/
			if (vrs) { // just adding to results for this iteration
				vrs = validator->updateValidationResultSet(
													vrs, newParticles);
			}
			else { // first run through loop for this iteration
			
				#ifdef DDDEBUG
					std::cout << "setting up the validator ..." << std::endl;
				#endif
			
				vrs = validator->setNormaliserAndGetValidationResultSet(
														newParticles);
														
				#ifdef DDDEBUG
					std::cout << "                    ... done" << std::endl;
				#endif
														
				/* note that if initialising and the epsilonMult for 
				 * the intial run is 1, then the initial epsilon 
				 * will allow exactly baseEpsilonPercentile of the 
				 * generated particles to pass, and since we aim to 
				 * generate n_reps/baseEpsilonPercentile of them,
				 * we will get exactly the right number on the initial run.*/
				if (initialising) {
					setInitialEpsilon(epsilonMult, vrs);
					
				}
			}
			
			/* keep track of how many attempts (proposals) - used 
			 for the acceptance ratio and just for interest.*/
			attempts += to_add;
			
			
			
			#ifdef MYDEBUG_EXTRA
			{
				std::cout << "attempts " << attempts << std::endl;
				std::cout << "results " << std::endl;
				ParameterParticleValidationResultSet::iterator it;
				for (it = vrs->begin();
					it < vrs->end();
					++it) {
						std::cout << (*it)->toString() << std::endl;
					}
			}				
			#endif
			
			
			/* discard from the results any with distance > epsilon. 
			* note not >= epsilon
			* otherwise we will discard too many if epsilon falls in the middle of a
			* a stretch of consecutive distances all == e
			* use discardByN later to take care of the rest if necessary */

			vrs->discardByDistance(epsilon);
			// vrs now only has compliant results
			
			#ifdef MYDEBUG
			{
				std::cout << "after discardByDistance(" << epsilon << "), results size = " << vrs->size() << std::endl;
			}			
			#endif
			
			// update current size 
			current_size = vrs->size();
			//update acceptance ratio
			if ( !current_size )
				throw std::runtime_error(
					"PopMCabcSimulatorAdaptive::PopMCabcSimulatorAdaptive(...) : no particles accepted");

			#if(0)
				acceptanceRatio = static_cast<double>(ceil(current_size)/attempts); // old version
			#endif
			acceptanceRatio = static_cast<double>(current_size)/attempts; // double division
			#ifdef MYDEBUG
			{
				std::cout << "acceptance ratio is now " << acceptanceRatio << std::endl;
			}			
			#endif
			
		} // end while loop
		
		lastAcceptanceRatio = acceptanceRatio; // save it for another iteration
		
		/* vrs should now contain at least n_reps results,
		 * but it might contain more - need to take out the extras
		 */
		if ( current_size > n_reps) {
			vrs->discardN(current_size - n_reps);
			#ifdef MYDEBUG
			{
				std::cout << "after discardN, results size = " << vrs->size() << std::endl;
			}			
			#endif
		}
		
		#ifdef DEBUG
			if (initialising) {
				std::cout << "final initial epsilon = " << epsilon << std::endl;
			}
		#endif
		
		assert(vrs->size() == n_reps);
		/*
		#ifdef MYDEBUG
		{
			std::cout << "after discardN(" << (current_size - n_reps) << "), results " << std::endl;
			ParameterParticleValidationResultSet::iterator it;
			for (it = vrs->begin();
				it < vrs->end();
				++it) {
					std::cout << (*it)->toString() << std::endl;
				}
		}				
		#endif
		*/
		
		// should now have n_reps results
		boost::shared_ptr < ParameterParticleSet > 
						acceptedNewParticles = vrs->getParameterParticleSet();
		
		assert(acceptedNewParticles->size() == n_reps);
		
		#ifdef DEBUG
			std::cout << "\nTotal attempts (including successes) = " 
					<< attempts << std::endl;
		#endif
		
		//calculate the new weights
		if( initialising ) {
			initialiseWeights(acceptedNewParticles);
		}
		else {
			setNewWeights(acceptedNewParticles, tau);
		}
		if (tau != NULL) gsl_matrix_free (tau);	
		
		
		// copy in the new particles and get the new summary statistics
		particles = acceptedNewParticles;
		assert(particles->size() == n_reps);
		
		summaryStatistics = vrs->getSummaryStatisticSet();
		assert(summaryStatistics->size() == n_reps);
		
		clock_t endIt = clock();
		/* if initialising, start log file, and in all cases, log summary */
		if (initialising) startLogFile(logfilename);
		logIteration(logfilename, 
			n_reps, epsilonMult, attempts, 
			static_cast<double>(n_reps)/attempts,
			((static_cast<double>(endIt - startIt)) / CLOCKS_PER_SEC));
		
	}
	catch (...) {
		try {
			if (tau != NULL) gsl_matrix_free (tau);	
		}
		catch (...) {}
		throw;	
	}
}
			
const boost::shared_ptr < ParameterParticleSet >& 
				PopMCabcSimulatorAdaptive::getCurrentParticles() const
{
	return particles;
}

const boost::shared_ptr < mct::SummaryStatisticSet >& 
				PopMCabcSimulatorAdaptive::getCurrentSummaryStatistics() const
{
	return summaryStatistics;
}						
			
double PopMCabcSimulatorAdaptive::getCurrentEpsilon() const
{
	return epsilon;
}

//protected



boost::shared_ptr< ParameterParticle > 
			PopMCabcSimulatorAdaptive::newParticlePtr(
			const gsl_matrix * const tau) const
{
	gsl_vector* mean = NULL;
			
	gsl_vector* result = NULL;
	
	int maxtries = MAXTRIES;
	
	#ifdef MYDEBUG
		// keep the particles we try with
		ParameterParticleSet pset(maxtries);
	#endif
	
	try {
		bool success = false;
		
		boost::shared_ptr < ParameterParticle > newPtr;
		
		while ((!success) & (maxtries > 0)) {
			
			size_t picked = particles->pickOne(unif_func);
			
			mean = particles->at(picked)->get_gsl_vector();
			
			result = gsl_vector_calloc(n_params);
			
			// puts a randomly generated values into result
			simulateVector(mean, tau, result);
			
			gsl_vector_free(mean);
			mean = NULL;
			
			newPtr = boost::shared_ptr < ParameterParticle > 
								 ( new ParameterParticle(result) ); 
			
			gsl_vector_free(result);
			result = NULL;
			
			
			double wt = priors.pdf(newPtr);
						
			if ( wt > 0.0 ) {
				success = true;
			}
			else {
				#ifdef MYDEBUG
					// keep the particles we try with using pdf as weight
					newPtr->setWeight(wt);
					pset.add(newPtr);
					
					
				#endif
			}
			
			maxtries--;
		}
		if (!success) {
			
			#ifdef MYDEBUG
				// output the particles we tried
				pset.outputToFile(std::string("RejectedParticles"), false);
			#endif
			
			throw std::runtime_error(
				std::string("PopMCabcSimulatorAdaptive::newParticlePtr(const gsl_matrix * const)")
				+ std::string(":\nCannot find new proposal particle within prior support") );
		}
		
		return newPtr;
	}
	catch (...) {
		try {
			if (mean != NULL) gsl_vector_free(mean); 
		}
		catch (...) {}
		try {
			if (result != NULL) gsl_vector_free(result); 
		}
		catch (...) {}
		throw;
	}
}


boost::shared_ptr< ParameterParticleSet > 
			PopMCabcSimulatorAdaptive::newParticlePtrs(
			const size_t n_reps,
			const gsl_matrix * const tau) const
{
	boost::shared_ptr< ParameterParticleSet > 
						pset( new ParameterParticleSet(n_reps) );
						
	for (size_t i = 0; i < n_reps; ++i) {
		
		#ifdef DEBUG
			if ( (n_reps <= 10) || ( i % (n_reps/10) == 0) ) {
				std::cout << ".";
				std::cout.flush();
			}
		#endif
		
		
		
		pset->add( newParticlePtr(tau) );
	
	}
	
	return pset;
}

#if(0) // old method before parameter prior set generated its own parameter particles
boost::shared_ptr< ParameterParticle > 
			PopMCabcSimulatorAdaptive::initialParticlePtr() const
{
	// make a parameter particle using the priors
	boost::shared_ptr < ParameterParticle > 
								p ( new ParameterParticle(n_params) );
	for (size_t j = 0; j < n_params; ++j)
	{
		p->at(j) = priors.randFrom(j);
	}
	return p;
}
#endif

boost::shared_ptr < ParameterParticleSet >
			PopMCabcSimulatorAdaptive::initialParticlePtrs(
					const size_t n_reps) const
{
	return priors.generateParameterParticleSet(n_reps); 
}


void PopMCabcSimulatorAdaptive::simulateVector(const gsl_vector * const mean,
						const gsl_matrix * const tau,
						gsl_vector * const result) const
{
	// puts a randomly generated multivariate normal vector into result
	rmvnorm(ugaussian_func , mean, tau, result);
}
		


void PopMCabcSimulatorAdaptive::initialisePRNGs(const boost::shared_ptr < mct::PRNGen >& _r_ptr)
{
	r_rand_ptr = _r_ptr->spawnAnother();

	//use the current random number generator to get a seed
	// seed will be between 1 and gsl_rng_max(r) inclusive
	/*size_t seed = 1 + 
				gsl_rng_uniform_int(r, gsl_rng_max(rand));		
	gsl_rng_set(rand, seed);
	*/
	unif_func =
		boost::bind<double>(&mct::PRNGen::uniform01_double, r_rand_ptr);
		//boost::bind<double>(gsl_ran_flat, rand, 0.0, 1.0);
	
	ugaussian_func =
		boost::bind<double>(&mct::PRNGen::ugaussian_double, r_rand_ptr);
		//boost::bind<double>(gsl_ran_ugaussian, rand);
}


gsl_matrix * PopMCabcSimulatorAdaptive::calcTau() const
{
	gsl_matrix * cov = NULL;
	gsl_matrix * tau = NULL;
	try {
		cov = particles->weightedCovariance();
		tau = gsl_matrix_calloc(cov->size1, cov->size2);
		
		gsl_matrix_memcpy (tau, cov);
		
		gsl_matrix_free(cov);
		cov = NULL;
		
		double scalar = 2.0; // scalar multiplier tau = scalar * cov
		gsl_matrix_scale (tau, scalar);
		
		return tau;
	}
	catch (...) {
		try {
			if (cov != NULL) gsl_matrix_free(cov); 
		}
		catch (...) {}
		try {
			if (tau != NULL) gsl_matrix_free(tau); 
		}
		catch (...) {}
		throw;
	}
}		
		
// give the particle their default weights 
void PopMCabcSimulatorAdaptive::initialiseWeights(
			const boost::shared_ptr < ParameterParticleSet >& 
			newParticles) const
{
	size_t N = newParticles->size();
	double defaultWeight = 1.0/N; // note 1.0 to ensure double
	for (size_t i = 0; i < N; ++i) {
		newParticles->at(i)->setWeight( defaultWeight );
	}
}

void PopMCabcSimulatorAdaptive::setNewWeights(
			const boost::shared_ptr < ParameterParticleSet >& newParticles,
			const gsl_matrix * const tau) const
{
	std:: size_t N = newParticles->size();
	std:: size_t presentN = particles->size();
		
	std::vector < gsl_vector * > newVectors( N, NULL );
	std::vector < gsl_vector * > means( presentN, NULL );
	
	try {
		
		/* set up containers of the gsl_vectors corresponding to new 
		 * and old particles once, to save redoing this (for the old particles
		 * at least) multiple times - not great but probably worth it
		 * */
		for (size_t i = 0; i < N; ++i) {
			newVectors.at(i) = newParticles->at(i)->get_gsl_vector();
		}
		
		for (size_t j = 0; j < presentN; ++j) {
			means.at(j) = particles->at(j)->get_gsl_vector();
		}
			
		std::vector < double > newWeights( N );
		double accNewWeights = 0.0;
		
		// for each new particle
		for (size_t i = 0; i < N; ++i) {
			
			//get the prior on this new particle
			double priorDensity = priors.pdf(newParticles->at(i));
			
			// sum over all old particles
			double divisor = 0.0;
			for (size_t j = 0; j < presentN; ++j) {
			
				divisor += ( (particles->at(j))->getWeight() 
						* calcWeightAdjustment(i,
										j,
										newVectors, means,
										tau)
						 );
			}
			double newWeight = priorDensity/divisor;
			newWeights.at(i) = newWeight;
			accNewWeights += newWeight;
			
		}
		//get rid of the gsl_vectors
		for (std::vector< gsl_vector * >::iterator it = newVectors.begin();
				it < newVectors.end();
				++it) {
			gsl_vector_free(*it);
			(*it) = NULL;
		}
		for (std::vector< gsl_vector * >::iterator it = means.begin();
				it < means.end();
				++it) {
			gsl_vector_free(*it);
			(*it) = NULL;
		}
		
		//set the new weights, normalising by the total new weight
		for (size_t i = 0; i < N; ++i) {
			newParticles->at(i)->setWeight( newWeights.at(i)/accNewWeights );
		}
	}
	catch (...) {
		
		//get rid of the gsl_vectors
		for (std::vector< gsl_vector * >::iterator it = newVectors.begin();
				it < newVectors.end();
				++it) {
			try {
				if ((*it) != NULL) gsl_vector_free(*it);
				(*it) = NULL;
			}
			catch (...) {}
		}
		for (std::vector< gsl_vector * >::iterator it = means.begin();
				it < means.end();
				++it) {
			try {
				if ((*it) != NULL) gsl_vector_free(*it);
				(*it) = NULL;
			}
			catch (...) {}
		}
		throw;
	}
}

double PopMCabcSimulatorAdaptive::calcWeightAdjustment(
			size_t i, size_t j,
			const std::vector < gsl_vector * >& newVectors,
			const std::vector < gsl_vector * >& means,
			const gsl_matrix * const tau) const 
{
	return dmvnorm(newVectors[i], means[j], tau); // from rmv.cpp
}

void PopMCabcSimulatorAdaptive::setInitialEpsilon(
	const double epsilonMult,
	const boost::shared_ptr < ParameterParticleValidationResultSet >& vrs)
{
	
	// calculate the epsilon
	epsilon = epsilonMult * 
			vrs->getDistancePercentile(baseEpsilonPercentile);
	// debugx
	#ifdef MYDEBUG
		std::cout << "setting initial epsilon: epsilonMult = " << epsilonMult << "and bEP = "<< baseEpsilonPercentile << std::endl;
		std::cout << "epsilon = " << epsilon << std::endl;
	#endif
	
}

void PopMCabcSimulatorAdaptive::updateEpsilon(const double epsilonMult)
{
	epsilon *= epsilonMult;
}

// start a log file
// s is the file name
void PopMCabcSimulatorAdaptive::startLogFile(const std::string& s) const
{
    // Make a string with filename and timestamp to start log file
    time_t rawtime;
    struct tm * timeinfo;
    time ( &rawtime );
    timeinfo = localtime ( &rawtime );
    std::ofstream os(s.c_str());         // replace data
    if (os.is_open()) {
        os << "File " << s << " created " <<  asctime (timeinfo) << std::endl;
		
		// add the description of the priors
		os << "\n" << priors.toString() << std::endl;
		
		// add a new line and the base epsilon
		os << "baseEpsilonPercentile\t" << baseEpsilonPercentile << std::endl;
		
		// add headers
		os << "\nreps\temult\teps\tatmpt\tratio\tsecs";
		
        os.close();
    }
    else {
        std::cerr << "Error: could not open file named "
            << s << std::endl << std::endl;
    }

}

void PopMCabcSimulatorAdaptive::logIteration(const std::string& s,
										size_t n_reps,
										double epsilonMult,
										size_t attempts,
										double acceptanceRatio,
										double secs,
										int prec) const
{
	std::ofstream os(s.c_str(), std::ios::app);         // append
    if (os.is_open()) {
		//set output precision for doubles
		os << std::setprecision(prec);
		
		// add line
		os << "\n" << n_reps
			<< "\t" << epsilonMult 
			<< "\t" << epsilon 
			<< "\t" << attempts 
			<< "\t" << acceptanceRatio
			<< "\t" << secs;
		
		os.close();
	}
    else {
        std::cerr << "Error: could not open file named "
            << s << std::endl << std::endl;
    }
}

