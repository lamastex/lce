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
/*! \file ebc_summaries.cpp
	\brief methods for some population genetic summary statistics
*/
#include<ebc_summaries.hpp>
#include<iomanip>

/*! constructor for SfsInfo*/
SfsInfo::SfsInfo( const size_t & sfsSz = 0 ) : S(0), uS(0), Pi(0), D(0.), SFS( sfs_array_type(sfsSz) ), uSFS( std::vector<unsigned>(sfsSz) ) 
{
}	

/*! get the various summaries from the sfs*/
void SfsInfo::Inform( const sfs_array_type & sfs ) 
{
  SFS=sfs;
  for(unsigned i=0; i < sfs.size(); i++) uSFS[i]= unsigned(sfs[i]);
  S= Ss(sfs);
  uS = unsigned(S);
  Pi= PI(sfs);
  D= tajd_SFS(uS, sfs);
}	

/*! output the summaries in SfsInfo*/
std::ostream & SfsInfo::output(std::ostream& o)
{
  o << "# S =  "<< S << std::endl;
  o << "# PI = "<< Pi << std::endl;
  o << "# D  = "<< D << std::endl;
  o << "# The Observed Site Freq Spec. is :" << std::endl;
  o << "# "; 
  for(unsigned i=0; i< SFS.size()-1; i++) 
  {
    o << ( (i>0) ? std::setw(3) : std::setw(1) ) << SFS[i] << " ";
  }
  o << std::setw(3) << SFS[SFS.size()-1] << std::endl;    
  return o;
}

/*! an intermediate quantity needed by tajd_SFS*/ 
double a_sub_n(const unsigned & nsam)
{
   double a = 0.;
   for(unsigned i=1;i<nsam;++i)
     a += 1./double(i);
   return a;
}

/*! an intermediate quantity needed by tajd_SFS*/ 
double b_sub_n(const unsigned & nsam)
{
   double b = 0.;
   for(unsigned i=1;i<nsam;++i)
     b += 1./pow(double(i),2.);
   return b;
}

//double tajd(const unsigned & S, const vector<unsigned> & sfs)

/*! return Tajima's D statistic for sfs of type sfs_array_type with S segregating sites*/
double tajd_SFS(const unsigned S, const sfs_array_type& sfs)
{
   const unsigned nsam = sfs.size()+1;
   double denom = nsam*(nsam-1);
   double pi = 0.;
   for(unsigned i=0;i<sfs.size();++i)
     {
       if(sfs[i]>0)
     {
       unsigned bin = i+1;
       //      double ssh = 1.-double(bin)*double(bin-1)/denom;
       //cerr << i << ' ' << sfs[i] << ' ' << bin << ' ' << ssh << ' ';
       pi += 2.*double(sfs[i])*double(bin)*double(nsam-bin)/denom;
       //cerr << pi << '\n';
     }
     }
   //  cerr << "pi = " << pi << '\n';
   double a1, a2, b1, b2, c1, c2, e1, e2;

   a1 = a_sub_n (nsam);
   a2 = b_sub_n (nsam);
   b1 = (double(nsam) + 1.0) / (3.0 * (double(nsam) - 1.0));
   b2 = (2.0 * (pow (double(nsam), 2.0)
            + double(nsam) + 3.0)) / (9.0 * double(nsam) *
                      (double(nsam) - 1.0));
   c1 = b1 - 1.0 / a1;
   c2 = b2 - (double(nsam) + 2.0) / (a1 * double(nsam)) + a2 / pow (a1, 2.0);
   e1 = c1 / a1;
   e2 = c2 / (pow (a1, 2.0) + a2);
   //  cerr << pi << ' ' << double(S)/a1 << ' ' << ( pi - double(S)/ a1 ) / pow ((e1 * S + e2 * S * (S - 1.0)), 0.5) <<'\n';
   return ( pi - double(S)/a1 ) / pow ((e1 * S + e2 * S * (S - 1.0)),  0.5);
}

/*! return Tajima's D statistic for sfs of type std::vector<unsigned> with S segregating sites*/
double tajd(const unsigned S, const std::vector<unsigned> & sfs)
{
   const unsigned nsam = sfs.size()+1;
   double denom = nsam*(nsam-1);
   double pi = 0.;
   for(unsigned i=0;i<sfs.size();++i)
     {
       if(sfs[i]>0)
     {
       unsigned bin = i+1;
       //      double ssh = 1.-double(bin)*double(bin-1)/denom;
       //cerr << i << ' ' << sfs[i] << ' ' << bin << ' ' << ssh << ' ';
       pi += 2.*double(sfs[i])*double(bin)*double(nsam-bin)/denom;
       //cerr << pi << '\n';
     }
     }
   //  cerr << "pi = " << pi << '\n';
   double a1, a2, b1, b2, c1, c2, e1, e2;

   a1 = a_sub_n (nsam);
   a2 = b_sub_n (nsam);
   b1 = (double(nsam) + 1.0) / (3.0 * (double(nsam) - 1.0));
   b2 = (2.0 * (pow (double(nsam), 2.0)
            + double(nsam) + 3.0)) / (9.0 * double(nsam) *
                      (double(nsam) - 1.0));
   c1 = b1 - 1.0 / a1;
   c2 = b2 - (double(nsam) + 2.0) / (a1 * double(nsam)) + a2 / pow (a1, 2.0);
   e1 = c1 / a1;
   e2 = c2 / (pow (a1, 2.0) + a2);
   //  cerr << pi << ' ' << double(S)/a1 << ' ' << ( pi - double(S)/ a1 ) / pow ((e1 * S + e2 * S * (S - 1.0)), 0.5) <<'\n';
   return ( pi - double(S)/a1 ) / pow ((e1 * S + e2 * S * (S - 1.0)),  0.5);
}

/*! return S = the number of Segregating sites*/
int Ss (const sfs_array_type& va)
{
  return va.sum();
  //return accumulate(va.begin(),va.end(),0);
}

/*! return PI -- the unaveraged heterozygosity*/
int PI (const sfs_array_type& va)
{
  int z = va.size(), pi=0;
  for(int i=0; i < z; i++) 
    {
      assert( va[i]>=0 );
      pi += ( (i+1) * (z-i) * va[i]);
    }
  return pi;
}

/*! return first+last or folded singletons*/
int First_Plus_Last(const sfs_array_type& va)
{
  int z = va.size()-1;
  assert( va[0]>=0 );
  assert( va[z]>=0 );
  return (va[0]+va[z]);
}
