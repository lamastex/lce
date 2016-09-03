//c++ -o msfreq msfreq.cc -lsequence -O2
//c++ -o msfreq msfreq.cc -lsequence -O2 -I/Users/raazesh/include/ -L/Users/raazesh/lib/ -L/usr/local/lib -I/usr/local/include/
#include <Sequence/SimData.hpp>
#include <Sequence/SimParams.hpp>
#include <iostream>
#include <algorithm>
#include <functional>
#include <iterator>
int main(int argc, char **argv)
{
  Sequence::SimParams p;
  Sequence::SimData d;
  int rv;  
  rv = p.fromfile(stdin);
  if(rv == EOF) exit(0);
  while ( (rv=d.fromfile(stdin)) != EOF )
    {
      std::vector<double> freq(p.totsam()-1,0);
      for(Sequence::SimData::const_site_iterator i = d.sbegin();
	  i < d.send() ; ++i)
	{
	  unsigned nder = std::count(i->second.begin(),i->second.end(),'1');
	  freq[nder-1]++;
	}
      std::copy(freq.begin(),freq.end()-1,std::ostream_iterator<double>(std::cout," "));
      std::cout << *(freq.end()-1) << '\n';
    }
}
