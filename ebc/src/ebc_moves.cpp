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
/*! \file ebc_moves.cpp 
	\brief methods for moves in a Markov basis
*/
#include<ebc_moves.hpp>


OneMove::OneMove(const size_t & SpcSz = 0) : MoveId(0), move( sfs_array_type(SpcSz) ),size(SpcSz)
{
}

std::istream & OneMove::read(std::istream& I)
{
  for (unsigned i=0; i< size; i++) 
    {
      I >> move[i];
    }
  return I;
}

/*! proose an sfs vaProp that remains in the non-negative integer lattice about the 1-Markov-Basis-Neighborhood of current sfs va 
  
	the new sfs proposed is vaProp and we get this by adding or subtracting a move that is chosen uniformly 
	at random chosen from the MovesVector to the current sfs va
*/
void Make_A_Valid_Move(gsl_rng* rgsl, sfs_array_type& va, sfs_array_type& vaProp,
		std::vector<OneMove> & MovesVector, unsigned NumberOfMoves)
{
//  unsigned NumberOfMoves = MovesVector.size();
  do
  { // propose an sfs
    // this unbiased mirroring of +/- is added so that moves file can be smaller with the -1 * moves dropped
    if( gsl_rng_uniform_pos (rgsl) > 0.5 ) 
         vaProp = (va + MovesVector[gsl_rng_uniform_int (rgsl, NumberOfMoves)].move);
    else vaProp = (va - MovesVector[gsl_rng_uniform_int (rgsl, NumberOfMoves)].move);
//vaProp = ( ( gsl_rng_uniform_pos (rgsl) > 0.5 ) ? 
//(va + MovesVector[gsl_rng_uniform_int (rgsl, NumberOfMoves)].move) : 
//(va - MovesVector[gsl_rng_uniform_int (rgsl, NumberOfMoves)].move) );
  }
  //confirming that vaProp is in the reference set, in terms of, being non-negative
  while( vaProp.min() < 0 );
}

/*! produce the valid 1-Markov-Basis-Neighborhood of site frequency spectra as vaNbhd about the current sfs va 

	we obtain the 1-Markov-Basis-Neighborhood by simply adding and subtracting each move in the MovesVector to 
	the current sfs va and further ensuring that they remain in the non-negative integer lattice
*/
void Make_Nbhd_Valid_Moves(sfs_array_type& va, std::vector<sfs_array_type>& vaNbhd,
		std::vector<OneMove> & MovesVector, unsigned NumberOfMoves)
{
//  unsigned NumberOfMoves = MovesVector.size();
  sfs_array_type vaTmp(va);
  for(unsigned i=1; i<NumberOfMoves; i++) // i > 0 as we want to avoid self-loops
  { // propose an sfs via + move i
    vaTmp = (va + MovesVector[i].move);
    if(vaTmp.min() >= 0) vaNbhd.push_back(vaTmp);
    // propose an sfs via - move i
    vaTmp = (va - MovesVector[i].move);
    if(vaTmp.min() >= 0) vaNbhd.push_back(vaTmp);
  }
}


void ReadMoves(std::vector<OneMove>& MovesVector, params& p)
{
  std::ifstream  MovesFileStrm;
  OneMove OneMoveTemp(p.SfsSz);		// temp variable of the struct OneMove
  //printValarray1D(OneMoveTemp.move); // default constructor makes OneMoveTemp = (0,0,...,0) of length SfsSz
  MovesVector.push_back(OneMoveTemp); // inserting self-loop "staying put" move of (0,0,...,0)
  if(p.MovesFile.empty() ) {
    std::cerr << "Gotto give the Moves through a file < -m MovesFile> when < -t 2 > is called" << std::endl;
    exit(3);
  }
  else {
      MovesFileStrm.open(p.MovesFile.c_str()); //----------------------opening Moves file
      if (!MovesFileStrm) {
	      std::cerr << "Unable to open input file of moves :" << std::endl <<  p.MovesFile << std::endl;
	  exit(1);   // call system to stop
      } //--------------------------------
    unsigned counter = 1;
    do  
      {
        MovesFileStrm >> OneMoveTemp;//std::cerr << counter++ << std::endl;
        OneMoveTemp.MoveId = counter;
        counter++;
        MovesVector.push_back(OneMoveTemp);
        MovesFileStrm >> std::ws;
      }
    while(!MovesFileStrm.eof());
    //std::cerr << counter << '\n';
    MovesFileStrm.close();
  }//----------------end of reading in the moves and storing them in MovesVector

  //sfs_array_type tmpmove(p.SfsSz);
  std::vector<OneMove>::iterator OneMove_Iter;
  if(p.NOISY > 89) {
    for ( OneMove_Iter = MovesVector.begin( ) ; OneMove_Iter != MovesVector.end( ) ; OneMove_Iter++ ) {
      OneMoveTemp = *OneMove_Iter;
      printValarray1D(OneMoveTemp.move);
    }
  }
  //std::cerr << "size  = " << MovesVector.size() << std::endl;
  //-----------------------------end of checking that the vector MovesVector does contain the right stuff
  MovesFileStrm.close();
}
