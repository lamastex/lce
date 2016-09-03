/*! \file multiarray.cpp
	\brief ignore [experimental] -- to possibly replace valarray with boost::multiarray for arithmetics with f-sequences

	\todo Consider replacing std::valarray with boost::multiarray for the representation of f-sequences and subsequent
	arithmetics with them.  This should be efficient when sample size gets large and the f-sequences become 
	inefficiently represented by std::valarray`s due to their sparseness.
*/
#include "boost/multi_array.hpp"
#include <cassert>
#include <iostream>

int main () {

	typedef boost::multi_array<double, 3> array_type;
	  typedef array_type::index index;
	    array_type A(boost::extents[3][4][2]);

	      // Assign values to the elements

int values = 0;
  for(index i = 0; i != 3; ++i) 
	      for(index j = 0; j != 4; ++j)
		            for(index k = 0; k != 2; ++k)
				            A[i][j][k] = values++;

    // Verify values
int verify = 0;
  for(index i = 0; i != 3; ++i) 
	      for(index j = 0; j != 4; ++j)
		            for(index k = 0; k != 2; ++k)
				    std::cout << A[i][j][k] << ' ' ;
//				            assert(A[i][j][k] == verify++);

    return 0;
}

