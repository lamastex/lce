
/*
 * reading nucleotide data from file
*/


#include "config.h"
#include "poly_data_sequence_reader.hpp"
#include "multi_loci_polytable.hpp"


#include <boost/smart_ptr.hpp>
#include "cxsc.hpp"
#include <iomanip>

#include <iostream>
#include <vector>
#include <stdexcept>

using namespace std;
using namespace mct;


void doStats(boost::shared_ptr < const MultiLociPolyTable > mlpt);


void test1();
void test2();
			
int main()
{
	//test1();
	//test2();
	
	double d1 = 0.3;
	double d2 = 2.0;
	double d3 = d1/d2;
	
	cout << d3 << endl;
	cout << fixed << setprecision(1) << d3 << endl;
	
	cxsc::real xLow;
	cxsc::real xUp;
	"0.3" >> RndDown >> xLow;
	"0.3" >> RndUp >> xUp;
	
	cout << XLow << endl;
	cout << XUp << endl;
	
	return 0;
	
}

void test1()
{
	std::vector < std::string> filenames;
	
	filenames.push_back("testInputNucs.txt");
	filenames.push_back("testInputNucs1.txt");
	filenames.push_back("testInputNucs2.txt");
	filenames.push_back("testInputNucs3.txt");
	filenames.push_back("testInputNucs4.txt");
	filenames.push_back("testInputNucs5.txt");
	filenames.push_back("testInputNucs6.txt");
	filenames.push_back("testInputNucs7.txt");
	
	for (size_t i = 0; i < filenames.size(); ++i) {
		std::string s = filenames[i];
		
		cout << "\nReading " << s << endl;
		
		try {
			
			boost::shared_ptr < mct::MultiLociPolyTable > mlpt = 
			PolyDataSequenceReader::readMultiLociPolytableFromArlequinDNAfileSingleLocus
					(s);
	
			cout << "\nResult successful " << s << endl;
			
			cout << "The final polytable is " << endl;
			mlpt->printSubpopData(cout);

			doStats(mlpt);
	
		}
		catch (std::runtime_error& re) {
			cout << "\nResult not successful " << s << endl;
		}	
		
	}	
		
	
}

void test2()
{
	std::vector < std::string> filenames;
	
	filenames.push_back("testInputNucs8.txt"); // all upper case
	filenames.push_back("testInputNucs9.txt"); // some lower case
	
	for (size_t i = 0; i < filenames.size(); ++i) {
		std::string s = filenames[i];
		
		cout << "\nReading " << s << endl;
		
		try {
			
			boost::shared_ptr < mct::MultiLociPolyTable > mlpt = 
			PolyDataSequenceReader::readMultiLociPolytableFromArlequinDNAfileSingleLocus
					(s);
	
			cout << "\nResult successful " << s << endl;
			
			cout << "The final polytable is " << endl;
			mlpt->printSubpopData(cout);

			doStats(mlpt);
	
		}
		catch (std::runtime_error& re) {
			cout << "\nResult not successful " << s << endl;
		}	
		
	}	
		
	
}

void doStats(boost::shared_ptr < const MultiLociPolyTable > mlpt)
{
			
	boost::shared_ptr < const SummaryStatistic  > withinPopAvs 
						= mlpt->withinPopAverageDifferences();
	
	std::cout << "\nwithin pop average differences are " <<std::endl;
	cout << withinPopAvs->plainString() << endl;
	
	boost::shared_ptr < const SummaryStatistic  > betweenPopAvs 
						= mlpt->betweenPopAverageDifferences();
	std::cout << "\nbetween pop average differences are " <<std::endl;
	cout << betweenPopAvs->plainString() << endl;
	
	boost::shared_ptr < const SummaryStatistic  > correctedAvs 
						= mlpt->correctedBetweenPopAverageDifferences();
	std::cout << "\ncorrected between pop average differences are " <<std::endl;
	cout << correctedAvs->plainString() << endl;
	
	boost::shared_ptr < const SummaryStatistic  > pairwiseAvs 
						= mlpt->pairwiseAverageDifferences();
	std::cout << "\npairwise pop average differences are " <<std::endl;
	cout << pairwiseAvs->plainString() << endl;
	
	boost::shared_ptr < const SummaryStatistic  > Fst 
											= mlpt->fst();
	cout << "\nFst is " << Fst->plainString() << endl;
	
	boost::shared_ptr < const SummaryStatistic  > pairwiseFsts 
								= mlpt->pairwiseFst();
	std::cout << "\npairwiseFsts are " <<std::endl;
	cout << pairwiseFsts->plainString() << endl;

}

