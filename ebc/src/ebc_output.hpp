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
/*! \file ebc_output.hpp
        \brief prototypes of methods and templatized classes for output management
*/

#ifndef __OUTPUT_HPP__
#define __OUTPUT_HPP__
#include <ebc_sfstypes.hpp>
#include <ebc_params.hpp>
#include <iomanip>
#include <fstream>

#include <ebc_graphtypes.hpp>


using namespace std;


/*! a vanilla class for calls to output parameters in bursts*/
struct output_params
{
  double rho_per_locus;/*!< recombination rate per locus*/
  double growth_rate;/*!< exponential growth rate*/
  double theta_per_locus;/*!< scaled mutation rate theta per locus*/

  output_params( double r, double gr, double t);
  std::ostream & output(std::ostream& o) const;
};

/*! overloaded operator for outputting output_params*/
inline std::ostream & operator<< (std::ostream & o, const output_params& oparams)
{
  return oparams.output(o);
}

/*! a vanilla class for calls to output likelihoods in bursts*/
struct output_Lkl
{
  double lnl;/*!< log likelihood*/

  output_Lkl( double LNL );
  std::ostream & output(std::ostream& o);
};

/*! overloaded operator for outputting output_Lkl*/
inline std::ostream & operator<< (std::ostream & o, output_Lkl & ores)
{
  return ores.output(o);
}

/*! a vanilla class for calls to output likelihoods and Shannon's information in bursts*/
struct output_LklShi
{
  double lnl;/*!< log likelihood*/
  double shi;/*!< Shannon's information*/

  output_LklShi( double LNL, double SHI );
  std::ostream & output(std::ostream& o);

};

/*! overloaded operator for outputting output_LklShi*/
inline std::ostream & operator<< (std::ostream & o, output_LklShi & ores)
{
  return ores.output(o);
}

/*! a vanilla class for calls to output parameters and likelihoods in bursts*/
struct output_results
{
  double rho_per_locus;/*!< recombination rate per locus*/
  double growth_rate;/*!< exponential growth rate*/
  double theta_per_locus;/*!< scaled mutation rate theta per locus*/
  double lnl;/*!< log likelihood*/

  output_results( double r, double gr, double t, double LNL );
  std::ostream & output(std::ostream& o);
};

/*! overloaded operator for outputting output_results*/
inline std::ostream & operator<< (std::ostream & o, output_results & ores)
{
  return ores.output(o);
}

/*! a vanilla class for calls to output parameters, likelihoods and Shannon's information in bursts*/
struct output_results2
{
  double rho_per_locus;/*!< recombination rate per locus*/
  double growth_rate;/*!< exponential growth rate*/
  double theta_per_locus;/*!< scaled mutation rate theta per locus*/
  double shi;/*!< Shannon's information*/
  double lkl;/*!< likelihood*/

  output_results2( double r, double gr, double t, double SHI, double LKL );
  std::ostream & output(std::ostream& o);
};

/*! overloaded operator for outputting output_results2*/
inline std::ostream & operator<< (std::ostream & o, output_results2 & ores)
{
  return ores.output(o);
}

/*! print each element of type T in std::vector X*/
template<class T>
void printvector(const std::vector<T>& X)
{
  for(unsigned i = 0; i < X.size()-1; i++) {
    std::cerr << ( (i>0) ? std::setw(3) : std::setw(1) ) << X[i] << " ";
  }
  std::cerr << std::setw(3) << X[X.size()-1] << std::endl;
}

/*! print each element of type T in std::vector X after multiplying by S of type T*/
template<class T>
void printMultiProb(const std::vector<T>& X, T S)
{
  for(unsigned i = 0; i < X.size()-1; i++) {
    std::cerr << ( (i>0) ? std::setw(3) : std::setw(1) ) << X[i]*S << " ";
  }
  std::cerr << std::setw(3) << X[X.size()-1]*S << std::endl;
}

/*! print one dimensional valarray of type T to std::out*/
template<class T>
void printVA1D (const std::valarray<T> & va)
{
  for(size_t i=0; i<va.size()-1; i++) {
    std::cout << ( (i>0) ? std::setw(3) : std::setw(1) ) << va[i] << " ";
  }
  std::cout << std::setw(3) << va[va.size()-1] << std::endl;
}

/*! print one dimensional valarray of type T to std::out as a matrix with rows many Rows*/
template<class T>
void printVA2D (const std::valarray<T> & va, size_t rows)
{
  size_t k =0;
  size_t R = va.size()/rows;
  for(size_t i=0; i<R; i++) {
    for(size_t j=0; j<rows-1; j++) {
      k = i*rows+j;
      std::cout << ( (j>0) ? std::setw(3) : std::setw(2) ) << std::right << va[k] << " ";
    }
    std::cout << std::setw(3) << va[k+1] << std::right << std::endl;
  }
}

void printValarray1Dcout (const sfs_array_type& va);
std::ostream & printValarray1Dostream (std::ostream & o, const sfs_array_type& va);
void printValarray1D (const sfs_array_type& va);
void PrintLoadGridParams(params& p, std::vector<output_params> & output_GridParamsPTP);

/*! start a file for the likelihood output - makes header lines. */
std::string startLikelihoodOutputFile(ofstream& outfile,
                                const int myPid,
                                const bool importanceSampling,
                                const size_t importanceSampleSize,
                                const params& p,
                                const double coverageCriteria,
                                const double& true_rho_per_locus,
                                const double& true_growth_rate,
                                const double& true_theta_per_locus);


/*! Send a map from xstars to map from ns to fsequence counts to an outfile.
*/
void sendToOutfile(const PaddedXStarToN_FSequenceCountMap& paddedXStarMap,
                            const size_t min_n,
                            const size_t max_n,
                            const std::string& filename);

/*! Output results for one SFS to likelihood output file.*/
void outputResultsForOneSFS(const string& filename,
                            const unsigned indexSFS,
                            const unsigned actual_Num_Of_Replicates_per_theta,
                            const vector<unsigned>& maxIndexX,
                            const vector<output_params>& output_GridParams,
                            const vector<double>& maxX,
                            const vector<double>& shanInfoX,
                            const vector<double>& coverageProb,
                            const double coverageCriteria,
                            const unsigned sfsSieveSize,
                            const double true_rho_per_locus,
                            const double true_growth_rate,
                            const double true_theta_per_locus,
                            const double divergenceKL);

/*! Output a summary of all results to likelihood output file.*/
void outputSummaryResults(const string& filename,
                        const vector<output_params>& resultsMaxLikelihoods,
                        vector<double>& resultsKLDivergence,
                        const size_t resultsGrowthCoverageAccumulation,
                        const size_t resultsThetaCoverageAccumulation,
                        const size_t resultsJointGrowthThetaCoverageAccumulation,
                        const double true_rho_per_locus,
                        const double true_growth_rate,
                        const double true_theta_per_locus,
                        const double coverageCriteria);

#endif
