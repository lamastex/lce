#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>

int main () {
    using namespace boost::numeric::ublas;
    compressed_matrix<double> m (5, 5, 5 * 5);
    compressed_matrix<double> n (5, 5, 5 * 5);
    //if(m==n) std:cout << "m==n" << std::endl;
    //if(m!=n) std:cout << "m!=n" << std::endl;
    //if(m<n) std:cout << "m<n" << std::endl;
    for (unsigned i = 0; i < m.size1 (); ++ i)
      for (unsigned j = 0; j < m.size2 (); ++ j)
        if(i==j) m (i, j) = 3 * i + j;
    std::cout << m << std::endl;
    for (unsigned i = 0; i < m.size1 (); ++ i)
      for (unsigned j = 0; j < m.size2 (); ++ j)
        if(m(i,j) != 0) std::cout << i <<  ' ' << j << ' ' << m(i,j) << std::endl;
}
