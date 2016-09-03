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
/*! \file ebc_moves.hpp 
	\brief prototypes of methods for moves in a Markov basis
*/
#ifndef __MOVES_HPP__
#define __MOVES_HPP__

#include<ebc_output.hpp> // contains -- #include<ebc_sfstypes.hpp>
#include <gsl/gsl_rng.h>
#include<ebc_params.hpp>

/*! struct representing one move or element of a Markov basis

	mostly used for reading in the Markov Basis from a text file
*/
struct OneMove{
	
  int MoveId;/*!< Id for the move*/
  sfs_array_type move;/*!< the move as sfs_array_type but both positive and negative entries are allowed*/
  unsigned size;/*size of move*/

  OneMove(const size_t & SpcSz);/*!< constructor with a default sample size of 0*/

  std::istream & read(std::istream& I);/*!< to read in a move*/
};

/*! inut operator for reading in a move*/
inline std::istream& operator>> (std::istream& I, OneMove& M)
{
  return M.read(I);
}  

void Make_A_Valid_Move(gsl_rng* rgsl, sfs_array_type& va, sfs_array_type& vaProp,
		std::vector<OneMove> & MovesVector, unsigned NumberOfMoves);

void Make_Nbhd_Valid_Moves(sfs_array_type& va, std::vector<sfs_array_type>& vaNbhd,
		std::vector<OneMove> & MovesVector, unsigned NumberOfMoves);

/*! read in all the moves from the MovesFile in p and store them in MovesVector*/
void ReadMoves(std::vector<OneMove>& MovesVector, params& p);
#endif
