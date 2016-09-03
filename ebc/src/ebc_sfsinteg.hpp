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
/*! \file ebc_sfsinteg.hpp
    \brief prototypes of methods for constructing an sfs sieve and integrating over it
*/
#ifndef __SFSINTEG_HPP__ 
#define __SFSINTEG_HPP__ 

#include<ebc_moves.hpp>
#include<ebc_summaries.hpp>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_statistics.h>
#include<functional>
#include<algorithm>
#include<numeric>
#include<iterator>
#include<limits>

#define MOST_NEGATIVE_DOUBLE (- ( std::numeric_limits<double>::max() ) )
#define MOST_POSITIVE_DOUBLE (  ( std::numeric_limits<double>::max() ) )
#define MOST_NEGATIVE_DOUBLE_ALLOWED -1.0e+308
#define MOST_POSITIVE_DOUBLE_ALLOWED +1.0e+308

void CheckUserDefinedDoubleTruncation();

void PisToMultinomialPs(std::vector<double>& Multinom_Ps, 
		const std::valarray<double>& Pis, double heat, int PRINT );

double Multinomial_Prob_Ln(const size_t K, const double p[], const unsigned int n[]);
double Multinomial_Unordering_Factor_Ln(const size_t K, const unsigned int n[]);
double Multinomial_Ordered_Lnpdf(const size_t K, const double p[], const unsigned int n[]);

void MultinomialPsToSfs(sfs_array_type& va, int S, const std::valarray<double>& Pis, int PRINT );

void SfsToMultinomialPs(const sfs_array_type& va, std::vector<double>& Multinom_Ps, int PRINT );

void MeanPost_UnifDirichletPrior(const sfs_array_type& va, std::vector<double>& Multinom_Ps, double alpha_i, int PRINT );

double build_sieve(gsl_rng* rgsl, 
		   sfs_sieve & SfsSieve,
		   sfs_array_type& vaInit, 
		   std::vector<OneMove> & MovesVector, 
		   unsigned NumberOfMoves,
		   const std::vector<double> & Ps, 
		   double * prob_min_sieve, 
		   const params& p, unsigned NeededSamples, unsigned MaxTrials);

double expand_uniformly_sieve(gsl_rng* rgsl, 
		   sfs_sieve & SfsSieve,
		   sfs_array_type& vaInit, 
		   std::vector<OneMove> & MovesVector, 
		   unsigned NumberOfMoves,
		   const std::vector<double> & Ps, 
		   double * prob_min_sieve, 
		   const params& p, unsigned NeededSamples, unsigned MaxTrials);

double graph_sieve(sfs_sieve & SfsSieve,
		   sfs_array_type& vaInit, 
		   std::vector<OneMove> & MovesVector, 
		   unsigned NumberOfMoves,
		   const std::vector<double> & Ps, 
		   double * prob_min_sieve, 
		   const params& p, unsigned NeededSamples, unsigned MaxTrials);

bool integrate_over_sieve(sfs_sieve & SfsSieve, const std::valarray<double>& Pis, 
		const params& p, double & sfs_sieve_integral, double & sfs_sieve_ShanInfo);

double integrate_over_conditional_sfs(gsl_rng* rgsl, const sfs_array_type& vaInit,
		std::vector<OneMove> & MovesVector, unsigned NumberOfMoves, 
		const std::valarray<double>& Pis, const params& p);

template<typename T>
struct kahan_sum
{
  T s,c,y,t;
  kahan_sum() : s(0.),c(0.),y(0.),t(0.){}
  T & operator()(  T & v,const T & i )
  {
    y=i-c;	// c is zero or close to it
    t=s+y;	//If s is big and y small, then low-order digits of y are lost during the summation into t
    c=(t-s)-y;	//(t - s) recovers the high-order part of y; subtracting y recovers -(low part of y)
    s=t;	//Algebriacally, c should always be zero. Beware eagerly optimising compilers!
    return s;	//Next time around, the lost low part will be added to y in a fresh attempt.
  }
};

void SummarizeDensityGrid(std::vector<double> & X, const double volume, std::vector<unsigned> & MaxIndexX, std::vector<double> & MaxX, 
		std::vector<double> & CoverageProb, const unsigned TrueIndexX, std::vector<double> & ShanInfoX, params& p);

void CoverageProbShannonInfo(const std::vector<double> & X, std::vector<unsigned> & MaxIndexX, std::vector<double> & MaxX, 
		std::vector<double> & CoverageProb, const unsigned TrueIndexX, std::vector<double> & ShanInfoX, params& p);

#endif
