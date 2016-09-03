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
/*! \file ebc_summaries.hpp
	\brief prototypes and structs for some population genetic summary statistics
*/
#ifndef __SUMMARIES_HPP__
#define __SUMMARIES_HPP__

#include<ebc_sfstypes.hpp>
#include <cmath>
#include <cassert>

/*! struct for a site frequency spectrum and its various statistics*/
struct SfsInfo{
    int S; /*!< number of segregating sites as an int*/
    unsigned uS;/*!< number of segregating sites as an unsigned int*/
    int Pi;/*!< unaveraged heterozygosity as an int */
    double D;/*!< Tajima's D statistic -- function of S, Pi and n = sample size*/ 
    sfs_array_type SFS;/*!< site frequency spectrum or sfs of type sfs_array_type*/
    std::vector<unsigned> uSFS;/*!< site frequency spectrum or sfs as vector of unsigned int`s*/

    SfsInfo( const size_t & sfsSz);
    void Inform( const sfs_array_type & sfs );
    std::ostream & output(std::ostream& o);
};

/*! overloaded operator to output SfsInfo*/
inline std::ostream & operator<< (std::ostream & o, SfsInfo & sfsi)
{
	  return sfsi.output(o);
}
//prototypes of various functions associated with computing S, PI and folded singletons
double a_sub_n(const unsigned & nsam);
double b_sub_n(const unsigned & nsam);
//double tajd(const unsigned & S, const vector<unsigned> & sfs)
double tajd_SFS(const unsigned S, const sfs_array_type& sfs);
double tajd(const unsigned S, const std::vector<unsigned> & sfs);
int Ss (const sfs_array_type& va);
int PI (const sfs_array_type& va);
int First_Plus_Last(const sfs_array_type& va);

#endif
