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
#ifndef __SFSTYPES_HPP__
#define __SFSTYPES_HPP__
/*! \file ebc_sfstypes.hpp
	\brief structs and typedefs of basic entities
*/
#include<valarray>
#include<vector>
#include<map>
#include<list>
#include <iostream>

/*! templatized function object for lexicographical sorting of vectors whose elements have total ordering
*/
template <class T> 
class LexicoSorting 
{
  public:
    bool operator() (const T& t1, const T& t2) const {
      return std::lexicographical_compare(&t1[0], &t1[t1.size()-1], &t2[0], &t2[t2.size()-1]);
      //return lexicographical_compare(t1.begin(), t1.end(), t2.begin(), t2.end());
    }
};

/*! structure for tracking counts and weights during importance sampling
*/			
struct Count_PRw_ISw
{
	unsigned Count; /*!< counts = typically number of visits*/
	double PRw; /*! probability/weight */ 
	double ISw; /*!< importance sample weight*/

	Count_PRw_ISw(unsigned cnt, double prw, double isw); /*! constructor*/
};

/*! sfs_array_type is site frequency spectrum or sfs as an STL valarray of int`s*/
typedef std::valarray<int> sfs_array_type;
/*! sfs_sieve is a sorted associative container (STL map) of site frequency spectra and the corresponding probabilities*/
typedef std::map<sfs_array_type, double, LexicoSorting<sfs_array_type> > sfs_sieve;
/*! sfs_filter is a sorted associative container (STL map) of site frequency spectra and the corresponding Count_PRw_ISw*/
typedef std::map<sfs_array_type, Count_PRw_ISw, LexicoSorting<sfs_array_type> > sfs_filter;
/*! sfs_particlesV is an STL vector of sfs_array_type`s*/ 
typedef std::vector< sfs_array_type > sfs_particlesV;
/*! Tees is an STL vector of double`s for the total size (sum of all branch-lengths) of each n-coalescent tree in a collection*/
typedef std::vector<double> Tees;
/*! Pees is an STL vector of the EpochTimes of each n-coalescent tree in a collection and EpochTimes is a valarray of n-1 double`s 

	Typically the EpochTimes valarray has the indexing convention: 
	EpochTimes[nsam-2] = t_n = Time to 1st Coal event, ..., EpochTimes[0] = t_2 = Time to last Coal event
*/ 
typedef std::vector< std::valarray<double> > Pees;
/*! fseq_type is a valarray of double`s for an f-sequence matrix*/
typedef std::valarray<double> fseq_type;
#endif
