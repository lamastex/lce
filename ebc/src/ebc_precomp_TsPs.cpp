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
/*! \file ebc_precomp_TsPs.cpp
    \brief methods for precomputing parameter-specific Tee and Pee in ParamsTeePee
*/
//Copyright(c) 2005, Raazesh Sainudiin and Kevin Thornton --  Wed Aug 17 23:44:58 EDT 2005
#include<ebc_precomp_TsPs.hpp>

/*! constructor for ParamsTeePee*/
ParamsTeePee::ParamsTeePee( const size_t & Psz = 0 ) : rho_per_locus(0.), growth_rate(0.), Tee(0.), Pee( std::valarray<double>(Psz) ), size(Psz)
{
}

/*! constructor for ParamsTeePee*/
ParamsTeePee::ParamsTeePee( double r, double gr, double T, const std::valarray<double>& P, unsigned Psz ) : rho_per_locus(r), growth_rate(gr), Tee(T), Pee(P), size(Psz)
{
}

/*! member function for outputting ParamsTeePee to ostream o*/
std::ostream & ParamsTeePee::output(std::ostream& o)
{
  o << rho_per_locus << "\t" << growth_rate << "\t" << Tee << "\t";
  for(unsigned i=0; i< (size-1); i++)
  {
    o << ( (i>0) ? std::setw(3) : std::setw(1) ) << Pee[i] << " ";
  }
  o << std::setw(3) << Pee[size-1];
  return o;
}

/*! member function for inputting ParamsTeePee from istream I*/
std::istream & ParamsTeePee::input(std::istream& I)
{
  I >> rho_per_locus;
  I >> growth_rate;
  I >> Tee;
  for(unsigned i=0; i< size; i++)
  {
    I >> Pee[i];
  }
  return I;
}

/*! read paramter values and Tee and Pee from the file TsPsFile and store in a std::vector<ParamsTeePee>*/
bool ReadPTPs(std::vector<ParamsTeePee>& PTPVector, params& p)
{
  std::ifstream  TsPsFileStream;
  if(p.TsPsFile.empty())
  {
    std::cerr << "No file specified for Having Stored Ts and Ps... assuming file is MeanTsPs" << std::endl;
    p.TsPsFile = "MeanTsPs";
  }
  TsPsFileStream.open(p.TsPsFile.c_str(), std::ios::in); // opening TsPsFile in read mode -- file must exist
  if(!TsPsFileStream)
  {
    std::cerr << "Unable to open TsPsFile " << p.TsPsFile << " for reading :(" << std::endl;
    return false;
    //exit(1001);
  }
  ParamsTeePee PTP_Temp(p.SfsSz);       // temp variable of the struct ParamsTeePee
  do
  {
    TsPsFileStream >> PTP_Temp;
    PTPVector.push_back(PTP_Temp);
    TsPsFileStream >> std::ws;
  }
  while(!TsPsFileStream.eof());
  TsPsFileStream.close();
  return true;
}

/*! print the average site frequency spectrum in the container SFS_Particles*/
void PrintAverageSFS(sfs_particlesV& SFS_Particles)
{

  if (SFS_Particles.empty()) {std::cerr << "SFS_Particles is empty !" << std::endl; exit (11);}
  std::valarray<int> temp_sfs_array((SFS_Particles[0]).size());
  temp_sfs_array=0;
  for(unsigned particle_index=0; particle_index<SFS_Particles.size(); ++particle_index)
  {
    //std::cout << SFS_Particles[particle_index] << "\t";
    //printValarray1D(SFS_Particles[particle_index]);
    temp_sfs_array += SFS_Particles[particle_index];
  }
  std::cout << std::endl << "Average = " << std::endl;
  printValarray1D((temp_sfs_array/int(SFS_Particles.size())));
}

/*! Printing Grid of model-parameters from std::vector< ParamsTeePee > PTPVector to model-paramter outfile PrmOutFile*/
void PrintGridParamsPTP(std::vector< ParamsTeePee >&  PTPVector, params& p, std::vector<output_params> & output_GridParamsPTP)
{
  std::ofstream PrmOutFileStream;
  PrmOutFileStream.open(p.PrmOutFile.c_str(), std::ios::out); // opening OutFile in append mode -- clean start!!!
  if(!PrmOutFileStream)
  {
    std::cerr << "Unable to open OutFile " << p.PrmOutFile << " for writing :(" << std::endl;
    exit(1001);
  }
  PrmOutFileStream.close(); // closing for now
  std::cout << std::endl << "Printing Grid of Parameters from PTP Vector" << std::endl;
  for(unsigned PTP_index = 0; PTP_index < PTPVector.size(); PTP_index++ )
  {
    //std::cout << std::endl << PTPVector[PTP_index] << std::endl;

    //Now we loop over theta
    for( int theta_i = 0 ; theta_i < p.theta_number_of_points ; ++theta_i)
    {
      double theta_per_locus = ( double(theta_i) * p.theta_per_locus_incr) + p.theta_per_locus_min;
      //double theta_per_site = theta_per_locus/(double(p.length));
      // now push the output_params into std::vector<output_params> output_GridParamsPTP
      output_GridParamsPTP.push_back
      (
        output_params( PTPVector[PTP_index].rho_per_locus, PTPVector[PTP_index].growth_rate, theta_per_locus)
      );
    // output parameter values
      if ( output_GridParamsPTP.size() == p.BUFFER_SIZE )
      {
        //print out contents of output in bursts of BUFFER_SIZE
        PrmOutFileStream.open(p.PrmOutFile.c_str(),std::ios::out|std::ios::app);
        PrmOutFileStream << std::scientific << std::showpoint << std::setprecision(6);
        for(unsigned output_index=0; output_index < output_GridParamsPTP.size(); output_index++)
        {
          PrmOutFileStream << output_GridParamsPTP[output_index] << std::endl;
        }
        PrmOutFileStream.close();
        output_GridParamsPTP.clear(); //erase the buffer.
      }
    }// theta loop
  }
  //print out the remaining crap
  if (!output_GridParamsPTP.empty())
  {
    PrmOutFileStream.open(p.PrmOutFile.c_str(),std::ios::out|std::ios::app);
    PrmOutFileStream << std::scientific << std::showpoint << std::setprecision(6);
    for(unsigned output_index=0; output_index < output_GridParamsPTP.size(); output_index++)
    {
      PrmOutFileStream << output_GridParamsPTP[output_index] << std::endl;
    }
    PrmOutFileStream.close();
    output_GridParamsPTP.clear(); //erase the buffer.
  }
}

