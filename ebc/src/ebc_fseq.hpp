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
/*! \file ebc_fseq.hpp
    \brief prototype of methods for unvintaged (sized) n-coalescent and associated shape statistics
*/
#ifndef __ARGSMAKE_HPP__
#define __ARGSMAKE_HPP__
#include<ebc_output.hpp>
#include<ebc_params.hpp>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
//libsequence headers



std::vector<unsigned> SFStoDefiniteItons(const sfs_array_type & X);

void SFStoDefiniteItonsDbl(std::valarray<double> & DefiniteItons, const sfs_array_type & X);

unsigned SampleBins(std::valarray<double> & p, double r);

void GenerateFSequence(gsl_rng* rgsl, fseq_type& fs, double& Probfs, double& PropProbfs, size_t nsam);

std::valarray<double>& make_StdNeutral_EpochTimes( gsl_rng* rgsl, params& p, std::valarray<double>& EpochTimes);

void make_set_of_pees_tees(gsl_rng* rgsl, Pees & PEES, Tees & TEES, const double & G, params& p);

std::valarray<double> EpochTimesProdFseq(std::valarray<double> & EpochTimes, fseq_type & fs, size_t nsam);


#endif
