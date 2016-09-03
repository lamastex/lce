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
\brief PopMCabcSimulatorAdaptive declarations.
*/

#ifndef _INC_POPMCABC_SIMULATOR_ADAPTIVE_H
#define _INC_POPMCABC_SIMULATOR_ADAPTIVE_H


#include "parameter_prior_set.hpp"

#include "parameter_particle_validation_result_set.hpp"
#include "parameter_particle_validator.hpp"
#include "parameter_particle_set.hpp"
#include "parameter_particle.hpp"
#include "summary_statistic_set.hpp"

#include "prng.hpp"

#include <gsl/gsl_matrix.h>

#include <boost/smart_ptr.hpp>
#include <boost/function.hpp>

#include <cstdlib>
#include <vector>
#include <string>

#define MAXTRIES 1000000 // max tries for generating one acceptable parameter particle


namespace aabc {
	
	/** @brief A class for an adaptive (sequential)
	 population Monte Carlo approximate
	 Bayesian computation simulator.
	 
	 The simulator uses a sequential sampling approach as described in
	 Beaumont, M.A., Cornuet, J., Marin, J., and Robert, C.P (2009), 
	 "Adaptive approximate Bayesian computation", <em>Biometrika</em>,
	 96(4), pp. 983-990.  
	 
	 A PopMCabcSimulatorAdaptive uses an adaptive population Monte Carlo 
	 approximate Bayesian computation algorithm to be able to:
	 <ul>
	 <li>initialise a ParameterParticleSet from a 
	 given ParameterPriorSet,</li>
	 <li>iterate from one ParameterParticleSet to another,</li>
	 <li>return the current ParameterParticleSet when asked.</li>
	 <li>return the mct::SummaryStatisticSet for the current 
	 ParameterParticleSet when asked.</li>
	 </ul>
	 
	 The PopMCabcSimulatorAdaptive is configured with a 
	 ParameterParticleValidator which it uses to validate proposed 
	 new ParameterParticles.  The %ParameterParticleValidator is
	 responsible for providing a mct::SummaryStatistic for
	 each propoposed ParameterParticle and for calculating a
	 distance from that SummaryStatistic to some target
	 SummaryStatistic which describe the suitability of a
	 proposed %ParameterParticle.
	 
	 The PopMCabcSimulatorAdaptive uses a <b>rejection sampling approach</b>, 
	 comparing the distances calculated
	 for each proposed %ParameterParicle to some tolerance or 
	 \f$ \varepsilon \f$ to decide whether to accept or reject each
	 proposal.  The PopMCabcSimulatorAdaptive approach is 
	 <b>sequential</b>, using a number of iterations of the rejection 
	 sampling process to
	 produce successively better proposals for a given target 
	 distribution.
	 
	 The tolerance or \f$ \varepsilon \f$ used in each iteration is
	 some fraction of that used in the previous iteration, with the
	 \f$ \varepsilon \f$ used in the first iteration being calculated
	 as some percentile of the empirical distribution of
	 distances calculated for ParameterParticles generated in that
	 first iteration. 
	 */
	 

	class PopMCabcSimulatorAdaptive {
		public:
		
			/*! \brief Constructor. 
			 
			 The simulator is configured with a ParameterPriorSet, a
			 ParameterParticleValidator, and a 
			 \link aabc::PopMCabcSimulatorAdaptive::baseEpsilonPercentile baseEpsilonPercentile\endlink
			 value.  The baseEpsilonPercentile is used to calculate the 
			 \f$ \varepsilon_1 \f$ used in the first 
			 \link aabc::PopMCabcSimulatorAdaptive::iterate iteration\endlink.
			 
			 For example, if \a bEP = 0.1, then the \f$ \varepsilon \f$
			 will be based on the 0.1 percentile distance over 
			 a number of replications in the first iteration.
			 
			 \param ps The ParameterPriorSet to use.
			 \param v The ParameterParticleValidator to use.
			 \param bEP The percentile of all sample distances to use as the base
			 for calculating the 'close enough' distance \f$ \varepsilon_1 \f$
			 (0 < bEP <= 1).
			 \param _r_ptr A pointer to a random number generator.
			\pre v is not a null pointer, and 0 < bEP <= 1.
			 */
			PopMCabcSimulatorAdaptive(
				const ParameterPriorSet& ps,
				const boost::shared_ptr < ParameterParticleValidator >& v,
				const double bEP,
				const boost::shared_ptr < mct::PRNGen >& _r_ptr);
		
			virtual ~PopMCabcSimulatorAdaptive();
			
			/*! \brief Iterate by generating a new ParameterParticleSet.
			
			If the current ParameterParticleSet is empty,
			then initialise, otherwise generate from the current set.
			
			The simulation method used is as described in 
			Beaumont, M.A., Cornuet, J., Marin, J., and Robert, C.P.,
			(2009), "Adaptive approximate Bayesian computation", 
			<em>Biometrika</em>, 96(4), pp. 983-990. 
			
			Some adaptations and interpretations of that method are
			included here.
			
			For a parameter of interest \f$ \theta \f$ and 
			\f$ N \f$ as the number of parameter values to generate in
			each iteration, and iterations \f$ t=1 \ldots T \f$, with
			a corresponding sequence of tolerance thresholds 
			\f$ \varepsilon_1 \geq \ldots \geq \varepsilon_T \f$, and a 
			prior \f$ \pi(\theta) \f$, and a 
			reference sample \f$ y \f$, and a measure of the distance
			between \f$ y \f$ and some other sample \f$ x \f$ 
			\f$ \varrho(x, y) \f$:
			
			<ol>
			<li><b>Initialisation step: \f$ t=1 \f$</b></li>
			<ul>
			<li>let \f$ i = 1 \f$</li>
			<li>while \f$ i \leq N \f$</li>
			<ul>
			<li>Simulate \f$ \theta_{\star}^{(1)} \sim \pi(\theta) \f$ and \f$ x \sim f(x | \theta_{\star}^{(1)}) \f$</li>
			<li>If \f$ \varrho(x, y) < \varepsilon_1 \f$ set \f$ \theta_i^{(1)} = \theta_{\star}^{(1)} \f$, \f$i=i+1 \f$</li>
			</ul>
			<li>For each \f$i = 1\ldots N \f$:</li>
			<ul>
			<li>Set the weight associated with \f$ \theta_i^{(1)}, \,\,\omega_i^{(1)} = \frac{1}{N} \f$</li>
			</ul>
			</ul>
			<li><b>Iteration step: for \f$ t=2, \ldots , T \f$</b></li>
			<ul>
			<li>Take \f$ \tau^2_{t} \f$ as twice the weighted empirical variance of the \f$ \theta_j^{(t-1)} \f$</li>
			<li>let \f$i = 1 \f$</li>
			<li>while \f$i \leq N \f$</li>
			<ul>
			<li>Pick \f$ \theta_{\star}^{(t-1)} \f$ from the \f$ \theta_j^{(t-1)} \f$ with probabilities \f$ \omega_j^{(t-1)} \f$</li>
			<li>Repeat</li>
			<ul> 
			<li>Generate \f$ \theta_{\star}^{(t)} | \theta_{\star}^{(t-1)} \sim \mathcal{N}(\theta_{\star}^{(t-1)},\, \tau^2_t) \f$</li>
			</ul>
			<li>Until \f$ \pi(\theta_i^{(t)}) > 0 \f$</li>
			<li>\f$ x \sim f(x | \theta_{\star}^{(t)}) \f$</li>
			<li>If \f$ \varrho(x, y) < \varepsilon_t \f$ set \f$ \theta_i^{(t)} = \theta_{\star}^{(t)} \f$, \f$ i=i+1 \f$</li>
			</ul>
			<li>For each \f$ \theta_i^{(t)} \f$</li>
			<ul>
			<li>Set \f$ \omega_i^{(t)} \propto \frac{\pi(\theta_i^{(t)})}{ \sum_{j=1}^N \phi \{ \tau_t^{-1}( \theta_i^{(t)}- \theta_j^{(t-1)}) \} } \f$</li>
			</ul>
			</ul>
			</ol>
			
			The generation of test data 
			\f$ x \sim f(x | \theta_{\star}^{(t)}) \f$
			and the measure of the distance \f$ \varrho(x, y) \f$ are 
			the responsibility of the simulator's ParameterParticleValidator. 
			
			The tolerance \f$ varepsilon_1 \f$ used in 
			the first iteration \f$ t=1 \f$
			is calculated using a specified percentile value from the
			empirical distribution of distances \f$ \varrho(x, y) \f$ from
			that first iteration
			(see 
			\link aabc::PopMCabcSimulatorAdaptive::baseEpsilonPercentile baseEpsilonPercentile\endlink,
			a configuration parameter for the simulator), 
			and a specified multiple of that percentile, \a epsilonMult 
			for the first iteration.
			
			For example, if baseEpsilonPercentile = 0.1 then the 
			tolerance for the first iteration 
			will be based on the \f$ d_{[0.1 percentile]} \f$,
			the 0.1 percentile value of the distances 
			\f$ \varrho(x, y) \f$, and if \a epsilonMult for
			that first iteration is  = 0.75, 
			then for this iteration 
			\f$ \varepsilon = 0.75d_{[0.1 percentile]} \f$/
			
			\note The number of tries to get one acceptable particle should 
			be limited (to prevent the program 'hanging')
			and if this limit is exceeded, the method
			should abort leaving the object as it was at the start of
			the iteration.
						
			\param epsilonMult The multiple of the current epsilon
			\f$ \varepsilon_{t-1} \f$ to use for this iteration.
			\note \f$ \varepsilon_0 \f$ or the 'base' epsilon
			is set using the 
			\link aabc::PopMCabcSimulatorAdaptive::baseEpsilonPercentile baseEpsilonPercentile\endlink
			and replications from the first iteration.  In the first
			iteration, \f$ t=1 \f$, \f$ \varepsilon_1 = \f$ \a epsilonMult \f$ \varepsilon_0 \f$
			using the value for \a epsilonMult supplied for that iteration.
			\param n_reps The number of repetitions to generate.
			\param logfilename The filename to use for logging.  If
			no filename is supplied the name defaults to the empty string "".
			If \a logfilename = "", no logging will take place. Otherwise
			a summary of the results will be put into the file
			identified by a logfilename.  On the first iteration the
			contents of any file of the same name will be overwritten;
			on subsequent iterations the log will be appended to the 
			existing content if a file named \a logfilename exists.
			\pre \a epsilonMult > 0.0.
			\pre \a n_reps >= 3. */
			virtual void iterate(const double epsilonMult,
									const size_t n_reps,
									const std::string& logfilename = "");
			
			/*! \brief Return a copy of the current
			ParameterParticleSet. */
			virtual const boost::shared_ptr < ParameterParticleSet >&  
						getCurrentParticles() const;
			
			/*! \brief Return a copy of the current
			mct::SummaryStatisticSet. */
			virtual const boost::shared_ptr < mct::SummaryStatisticSet >& 
				getCurrentSummaryStatistics() const;
			
			/*! \brief Return the value of the current epsilon. */
			virtual double getCurrentEpsilon() const;
								
		protected:
		
			PopMCabcSimulatorAdaptive();
									
			PopMCabcSimulatorAdaptive(const PopMCabcSimulatorAdaptive& other);
			
			PopMCabcSimulatorAdaptive& operator=(PopMCabcSimulatorAdaptive tmp);
			
			#if(0)
				virtual boost::shared_ptr< ParameterParticle > 
								initialParticlePtr() const;
			#endif
								
			virtual boost::shared_ptr < ParameterParticleSet >
				initialParticlePtrs(
					const size_t n_reps) const;
								
			virtual boost::shared_ptr< ParameterParticle > newParticlePtr(
								const gsl_matrix * const tau) const;	
			
			virtual boost::shared_ptr< ParameterParticleSet > 
				newParticlePtrs(
					const size_t n_reps,
					const gsl_matrix * const tau) const;
			
			virtual void simulateVector(const gsl_vector * const mean,
						const gsl_matrix * const tau,
						gsl_vector * const result) const;
			
			virtual void initialisePRNGs(const boost::shared_ptr < mct::PRNGen >& _r_ptr);
			
			// get the weighted covariance matrix of the current parameter particles and scale it
			virtual gsl_matrix * calcTau() const;
			
			// give the particle their default weights 
			virtual void initialiseWeights(
					const boost::shared_ptr < ParameterParticleSet >&
											newParticles) const;
						
			// set weights for the particles in newParticles
			virtual void setNewWeights(
					const boost::shared_ptr < ParameterParticleSet >&
											newParticles,
					const gsl_matrix * const tau) const;
			
			
			virtual double calcWeightAdjustment(
						size_t i, size_t j,
						const std::vector < gsl_vector * >& newVectors,
						const std::vector < gsl_vector * >& means,
						const gsl_matrix * const tau) const;
			
			/* set the initial epsilon based replicates with distance
			 * in baseEpsilonPercentile*/
			virtual void setInitialEpsilon(
				const double epsilonMult,
				const boost::shared_ptr < 
						ParameterParticleValidationResultSet >& vrs);
			
			/* update the epsilon for a new iteration */
			virtual void updateEpsilon( const double epsilonMult );
		
			/* start a log file */
			virtual void startLogFile(const std::string& s) const;
			
			/* log the summary from one iteration */
			virtual void logIteration(const std::string& s,
										size_t n_reps,
										double epsilonMult,
										size_t attempts,
										double acceptanceRatio,
										double secs,
										int prec = 3) const;
										//int prec = 3) const;
										
			//------------- data members -----------------
																	
			/*! \brief The priors. */
			const ParameterPriorSet priors;
			
			/* \brief The size of the parameter vectors we want. */
			const size_t n_params;
			
			/*! \brief The strategy for validating a 
			proposed ParameterParticle. */
			const boost::shared_ptr < ParameterParticleValidator > 
															validator;
			
			/* \brief The current set of 
			\link aabc::ParameterParticle ParameterParticles\endlink. */
			boost::shared_ptr < ParameterParticleSet > particles;
			
			/* \brief The current set of 
			\link SummaryStatistic SummaryStatistics\endlink. */
			boost::shared_ptr < mct::SummaryStatisticSet > summaryStatistics;
			
			/*! \brief The baseEpsilonPercentile is used to calculate the 
			 base \f$ \varepsilon \f$ (see 
			 \link aabc::PopMCabcSimulatorAdaptive::epsilon epsilon\endlink).*/
			 const double baseEpsilonPercentile;
			
			/*! \brief The current
			 \f$ \varepsilon \f$ used in the current
			 \link aabc::PopMCabcSimulatorAdaptive::iterate iteration\endlink.
			 
			 After the first iteration has been 
			 completed the 'base epsilon'
			 \link aabc::PopMCabcSimulatorAdaptive::epsilon epsilon\endlink
			 = \f$ \varepsilon_0 \f$ can be calculated from the \a baseEpsilonPercentile of the 
			 distances for the replicates generated in that first iteration.  
			
			 \f$ \varepsilon_1 \f$, the \f$ \varepsilon \f$ used for
			 rejection sampling in the first iteration, is calculated 
			 as some fraction (\a epsilonMult for that iteration) of
			 \f$ \varepsilon_0 \f$, 
			 \f$ \varepsilon_2 \f$ as some fraction of \f$ \varepsilon_1 \f$, etc
			 (see the \link aabc::PopMCabcSimulatorAdaptive::iterate iteration\endlink 
			 process).*/
			double epsilon;
			
			// the acceptance ratio for the last completed iteration
			double lastAcceptanceRatio;
			
			// prng
			boost::shared_ptr < mct::PRNGen > r_rand_ptr;
			
			// uniform 0,1 prng
			boost::function<double ()> unif_func;
			
			// standard normal prng
			boost::function<double ()> ugaussian_func;
			
		private:											
			
	};

} // end namespace aabc

#endif /* _INC_POPMCABC_SIMULATOR_H */

