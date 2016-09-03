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
/*! \file ebc_sfsdata.cpp
	\brief method for reading sfs data
*/

#include<ebc_moves.hpp>

/*! read in the sfs data from the file DataFile into container DataSfspectra of type std::vector<sfs_array_type> */
void ReadSfsData(std::vector<sfs_array_type>& DataSfspectra, params& p)
{
  std::vector<sfs_array_type>::iterator Sfspectra_Iter;
  std::ifstream DataFileStrm;
  OneMove OneMoveTemp(p.SfsSz);
  sfs_array_type va(p.SfsSz);

  if(p.DataFile.empty()) {
    std::cerr << "Gotto give the Data through a file < -d DataFile> " << std::endl;
    exit(3);
  }
  else {
    DataFileStrm.open(p.DataFile.c_str()); //----------------------opening Data file
    if (!DataFileStrm) 
    {
      std::cerr << "Unable to open input file of data :" << std::endl <<  p.DataFile << std::endl;
      exit(1);   // call system to stop
    } //--------------------------------
    //unsigned counter=1;
    do 
    {
      DataFileStrm >> OneMoveTemp;
      //std::cerr << counter++ << std::endl;
      DataSfspectra.push_back(OneMoveTemp.move); // remember a move is basically a valarray<int>
      DataFileStrm >> std::ws;
    }
    while(!DataFileStrm.eof() );

    if(p.NOISY > 99) {
      for ( Sfspectra_Iter = DataSfspectra.begin( ) ; Sfspectra_Iter != DataSfspectra.end( ) ; Sfspectra_Iter++ ) {
        va = *Sfspectra_Iter;
        printValarray1D(va);
      }
    }
  DataFileStrm.close();
  }
}
