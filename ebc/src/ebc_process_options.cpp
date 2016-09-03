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
/*! \file ebc_process_options.cpp
    \brief methods to process options from command line input arguments
*/
#include <ebc_process_options.hpp>
#include <cstdlib>
#include <getoptFix.h>

/*! help on usage of the program*/
void usage (std::string& prog)
{
  std::cerr << "usage: " << prog << " -m Movesfile -M MoveType (2 or 3) -d SfsDataFile -D NumOfSfsData -Y Data_identitY -q VerbosityLevel" << std::endl;
  std::cerr << "        -n Sieve_Size -N MaxTrials -C Type_Of_Sieve_Construction (0=Exhaustive, 1=Uniform Sampler, 2=Dirichlet Sampler)" << std::endl;
  std::cerr << "        -t theta_per_locus_min -T theta_per_locus_max -l length of locus in bp" << std::endl;
  std::cerr << "        -z random_seed -b No_Of_Sfs_bins(sample size - 1) -g Growth_Rate_Min -G Growth_Rate_Max" << std::endl;
  std::cerr << "        -x Recombinatio_Rate_Across_Locus_Min -X Recombinatio_Rate_Across_Locus_Max " << std::endl;
  std::cerr << "        -J Multinomial_Jumbles_EBC[0=approximate MLE used] -h sieve_temp_NeededSamples -i sieve_temp_MaxTrials" << std::endl;
  std::cerr << "        -H Heat_for_Proposed_Multinomial_in_EBC -P theta_number_of_points -Q growth_number_of_points" << std::endl;
  std::cerr << "        -s MAX_RS_NUM = number of rejection samples based on observed Ss as a Poisson of the total tree length * theta" << std::endl;
  std::cerr << "        -S MAX_ARGTRYNUM_for_Indep_M-H_chain_On_Simplex_with_PosteriorDirichletWeights -p prior_alpha_for_Dirichlet"<<std::endl;
  std::cerr << "        -R No_Replicates_over_args(i.i.d. OR MCMC over Dirich Posterior on simplex | growth, rho) -I No_Sfs_To_Integrate_EBC " << std::endl;
  std::cerr << "        -B burn in time on arg space - W wasted (burn in) time on conditional SFS space -K keep every Kth sample on SFS space" << std::endl;
  std::cerr << "        -A use average tree also -c simulate trees conditional on SFS -E TsPsFile (stores E(T), E(P_i's) for efficient simulations)" << std::endl;
}

/*! process options from command line input arguments and return params with all parameetr settings*/
params process_options(int argc, char **argv)
{
  params p;
  //extern int optind;
  int c;
  std::string prog = std::string (argv[0]);

  while ((c = getopt (argc, argv, "?:m:M:d:D:Y:q:C:n:N:z:b:p:t:T:g:G:x:X:J:h:i:H:R:c:s:S:B:I:W:K:P:Q:A:o:E:l:F:f:")) != -1)
    {
      switch (c)
        {
        case 'm':
          p.MovesFile = std::string (optarg);
          break;
        case 'M':
          p.MoveType = atoi (optarg);
          break;
        case 'd':
          p.DataFile = std::string (optarg);
          break;
        case 'D':
          p.NumOfSfsData = unsigned(atoi (optarg));
          break;
        case 'Y':
          p.DATA_id = unsigned(atoi (optarg));//data identitY
          break;
        case 'q':
          p.NOISY = unsigned(atoi(optarg)); // for quiet output -- level of rubbish
          break;
        case 'n':
          p.sieve_NeededSamples = unsigned(atoi (optarg));
          break;
        case 'N':
          p.sieve_MaxTrials = unsigned(atoi (optarg));
          break;
        case 'C':
          p.sieve_Construct_Type = (atoi (optarg));
      if(p.sieve_Construct_Type < -1 || p.sieve_Construct_Type > 2)
        {std::cerr << "ERROR : sieve_Construct_Type = 0,1, 2 or -1  BUT it is " << p.sieve_Construct_Type << std::endl; exit(0);}
          break;
    case 'z':
      p.seed = unsigned(atoi(optarg));
      p.USER_SEED = true;
      break;
    case 'b':// pass -b nbins
      p.SfsSz = unsigned(atoi(optarg));
      break;
        case 'p': //prior for Dirichlet
          p.prior_alpha = double(atof (optarg));
          break;
        case 't':
          p.theta_per_locus_min = double(atof (optarg));
          break;
        case 'T':
          p.theta_per_locus_max = double(atof (optarg));
          break;
        case 'g':
          p.growth_rate_min = double(atof (optarg));
          break;
        case 'G':
          p.growth_rate_max = double(atof (optarg));
          break;
        case 'x':
          p.rho_min = double(atof (optarg));
          break;
        case 'X':
          p.rho_max = double(atof (optarg));
          break;
        case 'J':
          p.Num_Multinomial_Starts = atoi (optarg);
          break;
        case 'h':
          p.sieve_temp_NeededSamples = unsigned(atoi (optarg));
          break;
        case 'i':
          p.sieve_temp_MaxTrials = unsigned(atoi (optarg));
          break;
        case 'H'://Local_Proposal heat for proposed multinomial during integration over sfs -- x/100
        // same value is used to get the alpha_i's for the Dirichlet prior during sieve building
          p.heat = double(atof (optarg));
          break;
        case 'R':
          p.Num_Of_Replicates = unsigned(atoi (optarg));
          break;
        case 'A':
          p.USE_AVG = true;
          break;
        case 'c':
      if(unsigned (atoi (optarg)) == 0)
        p.Arg_CondlOnSfs = false;
      else p.Arg_CondlOnSfs = true; // when rho=0 we simulate coalescent genealogies conditional on observed SFS
          break;
        case 'P':
          p.theta_number_of_points = double(atoi (optarg));
      if(p.theta_number_of_points < 2) {std::cerr << "-P # should be an integer > 1" << std::endl; exit(101); }
          break;
        case 'Q':
          p.growth_number_of_points = double(atoi (optarg));
      if(p.growth_number_of_points < 1) {std::cerr << "-P # should be an integer > 0" << std::endl; exit(101); }
          break;
        case 'I':
          p.Num_Sfs_Integrations = atoi (optarg);
          break;
        case 'l'://length of locus in bp
          p.length = atoi (optarg);
          if(p.length < 0) {std::cerr << "-s # should be an integer > 0 " << std::endl; exit(101);}
          break;
        case 's'://MAX_RS_NUM
          p.MAX_RS_NUM = atoi (optarg);
          if(p.MAX_RS_NUM < 0) {std::cerr << "-s # should be an integer >= 0 " << std::endl; exit(101);}
          break;
        case 'S'://MAX_ARGTRYNUM = maximum number of samplesfrom an independent M-H chain sampling args in propr. (pi's | posterior Dirichlet)*Poisson(T | S)
          p.MAX_ARGTRYNUM = atoi (optarg);
          if(p.MAX_ARGTRYNUM <= 0) {std::cerr << "-S # should be an integer > 0 " << std::endl; exit(101);}
          break;
        case 'B'://BURN IN NUMBER for MCMC in arg-compressed-into-simplex space
          p.BURNIN_ARG = atoi (optarg);
          if(!(p.BURNIN_ARG >= 0)) {std::cerr << "-B # should be an integer > -1 " << std::endl; exit(101);}
          break;
        case 'W'://BURN IN NUMBER for MCMC on conditional sfs space
          p.BURNIN_SFS = atoi (optarg);
          if(p.BURNIN_SFS < 0) {std::cerr << "-W # should be an integer > -1 " << std::endl; exit(101);}
          break;
        case 'K'://THIN_OUT rate for MCMC on conditional sfs space
          p.THIN_OUT_SFS = atoi (optarg);
          if(p.THIN_OUT_SFS < 0) {std::cerr << "-K # should be an integer > 0 " << std::endl; exit(101);}
          break;
    case 'o':
      p.OutFile = std::string(optarg);
      break;
        case 'E':
          p.TsPsFile = std::string (optarg);//Name of the File that will store the precomputed Expectations of Ts and Ps
          break;
          case 'F':
          if (atoi (optarg) == 1)
            p.boolImportanceSampling = true; //boolean for whether we are importance sampling
                                              // defaults to false
          break;
          case 'f':
          p.numImportanceSamples = atoi (optarg); //how many importance samples to do
          break;

        case '?':
        default:
          usage (prog);
      exit(1001);
          break;
        }
    }

  //HERE WE MAKE OTHER PARAMETER PROCESSING
  if(p.OutFile.empty())
  {
    std::cerr << "No output file specified SO output in file " << p.DataFile << "_out" << std::endl;
    p.OutFile = p.DataFile;
    p.OutFile += "_out";
  }
  p.LklOutFile = p.OutFile ; //p.LklOutFile += "_lkl";
  p.PrmOutFile = p.TsPsFile; p.PrmOutFile += "_prms";

  // setting the increments in grid-based computations
  if(p.theta_per_locus_max != p.theta_per_locus_min )
    p.theta_per_locus_incr = (p.theta_per_locus_max - p.theta_per_locus_min) / (p.theta_number_of_points-1.0);
  else p.theta_per_locus_incr = p.theta_per_locus_max;

  if( (p.growth_rate_max != p.growth_rate_min) && (p.growth_number_of_points > 1) )
    p.growth_rate_incr = (p.growth_rate_max - p.growth_rate_min) / (p.growth_number_of_points - 1.0);
  else p.growth_rate_incr = p.growth_rate_max;

  // just making the increment large for the case when only sampling from growth_rate_min once
  if(p.growth_number_of_points == 1) p.growth_rate_incr = p.growth_rate_max + 1.0;
  if(p.theta_number_of_points == 1) p.theta_per_locus_incr = p.theta_per_locus_max + 1.0;

  //setting rho increment as:
  p.rho_incr = 0.1;

  return p;
}
