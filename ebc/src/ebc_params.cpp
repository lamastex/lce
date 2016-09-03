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
/*! \file ebc_params.cpp
        \brief command line input arguments management
*/
#include <ebc_params.hpp>

/*! constructor with default initial parameter values for params*/
params::params() :

   MovesFile( std::string() ), DataFile( std::string() ), OutFile( std::string() ),
   LklOutFile( std::string() ), PrmOutFile( std::string() ), TsPsFile( std::string() ),

   MoveType(2), Num_Multinomial_Starts(1), sieve_Construct_Type(2),
   Num_Sfs_Integrations(100), MAX_RS_NUM(1), MAX_ARGTRYNUM(1),
   BURNIN_ARG(0), BURNIN_SFS(0), THIN_OUT_SFS(1), length(1000),

   boolImportanceSampling(false),

   NumOfSfsData(1), numImportanceSamples(0),
   NOISY(0),
   sieve_NeededSamples(1000), sieve_MaxTrials(100000),
   seed(1010101), SfsSz(0), sieve_temp_NeededSamples(1000), sieve_temp_MaxTrials(10000),
   Num_Of_Replicates(1), BUFFER_SIZE(10000),  DATA_id(0),

   USER_SEED(false), USE_AVG(false), Arg_CondlOnSfs(false),

  //rho=4Nr is recombination rate across region with l base pairs
   prior_alpha(1.0), theta_per_locus_min(0.001), theta_per_locus_max(10.0), theta_per_locus_incr(0),
   growth_rate_min(0.0), growth_rate_max(10.0), growth_rate_incr(0), rho_min(0.0), rho_max(0.0), rho_incr(1.0),
   heat(10.0), theta_number_of_points(10.), growth_number_of_points(10.)

{
}

/*
bool params::valid() const
{
  return ( (!priorfile.empty()) &&
       (!datafile.empty()) &&
       (!basename.empty()) &&
       nsumm > 0 &&
       nparams > 0 &&
       tolerance > 0. );
}
*/
