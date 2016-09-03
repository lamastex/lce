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
/*! \file ebc_precomp_TsPs.hpp
	\brief prototypes for precomputing parameter-specific Tee and Pee in ParamsTeePee
*/
//Copyright(c) 2005, Raazesh Sainudiin and Kevin Thornton --  Wed Aug 17 23:44:58 EDT 2005
#include<ebc_sfstypes.hpp>
#include<ebc_moves.hpp>
#include<ebc_summaries.hpp>
#include<ebc_output.hpp>
#include<ebc_sfsinteg.hpp>
#include<ebc_sfsparticles.hpp>
#include<ebc_sfsdata.hpp>
#include<ebc_params.hpp>
#include<ebc_process_options.hpp>
#include<ebc_summaries.hpp> 

/*! struct for calls to output Expected Ts and Ps in bursts*/
struct ParamsTeePee
{
  double rho_per_locus;/*!< per-locus recombination rate*/
  double growth_rate;/*!< exponential growth rate*/
  double Tee;
  std::valarray<double> Pee;
  unsigned size;

  ParamsTeePee( const size_t & Psz );
  ParamsTeePee( double r, double gr, double T, const std::valarray<double>& P, unsigned Psz );
  std::ostream & output(std::ostream& o);
  std::istream & input(std::istream& i);
};

/*! overloaded operator for outputting ParamsTeePee opt to std::ostream o*/
inline std::ostream & operator<< (std::ostream & o, ParamsTeePee & otp)
{
  return otp.output(o);
}  

/*! overloaded operator for inputting ParamsTeePee opt to std::istream i*/
inline std::istream & operator>> (std::istream & i, ParamsTeePee & itp)
{
  return itp.input(i);
}  

bool ReadPTPs(std::vector<ParamsTeePee>& PTPVector, params& p);

void PrintAverageSFS(sfs_particlesV& SFS_Particles);
	
void PrintGridParamsPTP(std::vector< ParamsTeePee >&  PTPVector, params& p, std::vector<output_params> & output_GridParamsPTP);
