/*
 * Copyright (C) 2005--2009 Raazesh Sainudiin and Kevin Thornton
 *
 * This file is part of lce, a C++ class library for lumped coalescent experiments.
 *
 * lce is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */
/*! \file ebc_params.hpp
        \brief prototype for command line input arguments management
*/
#ifndef __PARAMS_HPP__
#define __PARAMS_HPP__
#include <string>

/*! a struct for all parameters of the program

    The argc and argv to the main are processed by process_options to return an object of params
*/
struct params
{
   std::string MovesFile; /*!< The file with the set of moves from Markov basis*/
   std::string DataFile; /*!< The file with the sfs data*/
   std::string OutFile; /*!< The file with output*/
   std::string LklOutFile; /*!< The file with likelihood output*/
   std::string PrmOutFile; /*!< The file with parameter (grid) output*/
   std::string TsPsFile; /*!< The file with Ts:= Total tree size and Ps:=relative lengths of lineages leading to singleton,...,(n-1)-ton mutations*/

   int      MoveType; /*!< Type of move MoveType=2 ==> Markov basic of (S,PI) and MoveType=3 ==> Markov basis of (S,PI,Folded Singletons)*/
   int      Num_Multinomial_Starts;
   int      sieve_Construct_Type; /*!< Type_Of_Sieve_Construction (0=Exhaustive, 1=Uniform Sampler, 2=Dirichlet Sampler)*/
   int      Num_Sfs_Integrations; /*!< number of summations over set of sfs or sfs sieve*/
   int      MAX_RS_NUM; /*!< Maximum number of rejection samples*/
   int      MAX_ARGTRYNUM; /*!< Maximum number of ARGs to try or propose*/
   int      BURNIN_ARG; /*!< burin-in period during MCMC-based dependent proposals from ARG sapce*/
   int      BURNIN_SFS; /*!< burin-in period during MCMC-based dependent proposals from sfs sieve*/
   int      THIN_OUT_SFS; /*!< thinning rate to reduce dependency between samples in MCMC over sfs sieve*/
   int      length; /*!< length of locus in bp*/

   bool     boolImportanceSampling; /*!< bool for whether we are importance sampling*/

   unsigned     NumOfSfsData; /*!< number of sfs data to analyse*/
   unsigned    numImportanceSamples; /*! number of importance samples to generate if importance sampling*/
   unsigned NOISY; /*!< type of rubbish on screen*/
   unsigned sieve_NeededSamples; /*!< number of sfs samples needed in the sfs sieve*/
   unsigned sieve_MaxTrials;  /*!< maximum number of trials during integration over sfs sieve*/
   unsigned seed; /*!< seed of the random number generator*/
   unsigned SfsSz; /*!< size of the sfs (should equal number of samples - 1)*/
   unsigned sieve_temp_NeededSamples; /*!< sfs samples needed in the temporary/priming sfs sieve*/
   unsigned sieve_temp_MaxTrials; /*!< maximum number of trials in the temporary/priming sfs sieve*/
   unsigned Num_Of_Replicates; /*!< number of replications of the experiment*/
   unsigned BUFFER_SIZE; /*!< to print out contents of output in bursts of BUFFER_SIZE*/
   unsigned DATA_id; /*!< id number of the data being processed*/

   bool     USER_SEED; /*!< whether user-specified seed for RNG is used*/
   bool     USE_AVG; /*!< whether we want to use the average tree or ARG as well*/
   bool     Arg_CondlOnSfs; /*! whether we want to only use the coalescent trees (ARGS w/ rho=0) that are conditional on the observed sfs*/

   double   prior_alpha, theta_per_locus_min, theta_per_locus_max, theta_per_locus_incr,
        growth_rate_min, growth_rate_max, growth_rate_incr, rho_min, rho_max, rho_incr, heat,
        theta_number_of_points, growth_number_of_points;


  params();
  //bool valid() const;
};

#endif
