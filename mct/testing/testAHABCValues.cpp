/* MCT - Markov Chains on Trees.


   Copyright (C) 2011, 2012 Jennifer Harow

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor,
   Boston, MA 02110-1301, USA.  */

/*! \file
\brief Testing ParameterParticleSet.

*  
*/

#include "config.h"

#include "ahabc_values.hpp"
#include "parameter_particle_set.hpp"
#include "summary_statistic_set.hpp"

#include <iostream>
#include <iterator>



using namespace std;
using namespace aabc;
using namespace mct;

void testAHABCValues();
void testAHABCValuesSelected() ;

boost::shared_ptr < AHABCValues > getAHABCValuesLarger();
boost::shared_ptr < AHABCValues > getAHABCValuesLargerSVD();
boost::shared_ptr < AHABCValues > getAHABCValuesSmall();
boost::shared_ptr < AHABCValues > getAHABCValuesSmallScaled() ;
boost::shared_ptr < AHABCValues > getAHABCValuesAligned() ;
boost::shared_ptr < AHABCValues > getAHABCValuesSelected(
		const std::vector < size_t >& reqParams,
		const std::vector < size_t >& reqStats);


int main()
{

	//testAHABCValues();
	testAHABCValuesSelected() ;
	
	return 0;
	
}



void testAHABCValues() 
{
		
	std::vector < boost::shared_ptr < AHABCValues > > vals;
	
	vals.push_back(getAHABCValuesLargerSVD());
	//vals.push_back(getAHABCValuesSmall());
	vals.push_back(getAHABCValuesSmallScaled());
	vals.push_back(getAHABCValuesAligned()); 	
	
	for (size_t i = 0; i < vals.size(); ++i) {
		
		
		cout << "\n\nExamples\n" << endl;
		
		boost::shared_ptr < AHABCValues > AHABCValuesPtr = vals[i];
		
		size_t np = AHABCValuesPtr->getNumberParameters();
		size_t ns = AHABCValuesPtr->getNumberSummaryStatistics();

		cout << "AHABCValuesPtr->getNumberParameters() = " 
							<< np << endl;
		cout << "AHABCValuesPtr->getNumberSummaryStatistics() = " 
							<< ns << endl;
		
		cout << "AHABCValuesPtr->size() = " 
							<< AHABCValuesPtr->size() << endl;
		
		//some pretend reference values - a vector of 1's
		std::vector < double > refvec(np+ns, 1.0);
		
		std::vector<std::vector<double> > data;
		AHABCValuesPtr->fillDataContainer(data);
		
		cout << "\nData is = " << endl;
		for (std::vector< std::vector<double> >::iterator it = data.begin();
						it < data.end();
						++it) {
			ostream_iterator<double> out_it (cout,"\t");
			copy ( it->begin(), it->end(), out_it );
			cout << endl;
		}
		cout << endl;
		
		AHABCValuesPtr->dataSVD();
		
	
		{
			std::vector<std::vector<double> > transformedData;
			AHABCValuesPtr->fillTransformedDataContainer(transformedData);
			
			cout << "\nTransformed data is = " << endl;
			for (std::vector< std::vector<double> >::iterator it = transformedData.begin();
							it < transformedData.end();
							++it) {
				ostream_iterator<double> out_it (cout,"\t");
				copy ( it->begin(), it->end(), out_it );
				cout << endl;
			}
			cout << endl;
				
					
			std::vector < std::vector < double > > hm
								= AHABCValuesPtr->getHouseholderMatrix();
			
			cout << "\nHouseholder matrix is = " << endl;
			for (std::vector< std::vector<double> >::iterator it = hm.begin();
							it < hm.end();
							++it) {
				ostream_iterator<double> out_it (cout,"\t");
				copy ( it->begin(), it->end(), out_it );
				cout << endl;
			}
			cout << endl;
			
			cout << "\nReference vector is " << endl;
			{
				ostream_iterator<double> out_it (cout,"\t");
				copy ( refvec.begin(), refvec.end(), out_it );
				cout << endl;
			}
						
			std::vector < double > htRefvec 
				= AHABCValuesPtr->householderTransform(refvec);
				
			cout << "\nTransformed reference vector is " << endl;
			{
				ostream_iterator<double> out_it (cout,"\t");
				copy ( htRefvec.begin(), htRefvec.end(), out_it );
				cout << endl;
			}
			cout << endl;
		
		}
		{
			std::vector < double > means
								= AHABCValuesPtr->getMeans();
			
			cout << "\nMeans are " << endl;
			{
				ostream_iterator<double> out_it (cout,"\t");
				copy ( means.begin(), means.end(), out_it );
				cout << endl;
			}
			cout << endl;
			
			std::vector < double > sds
								= AHABCValuesPtr->getSDs();
			
			cout << "\nSDs are " << endl;
			{
				ostream_iterator<double> out_it (cout,"\t");
				copy ( sds.begin(), sds.end(), out_it );
				cout << endl;
			}
			cout << endl;
								
			std::vector<std::vector<double> > standardisedData;
			AHABCValuesPtr->fillStandardisedDataContainer(standardisedData);
			
			cout << "\nStandardised data is = " << endl;
			for (std::vector< std::vector<double> >::iterator it = standardisedData.begin();
							it < standardisedData.end();
							++it) {
				ostream_iterator<double> out_it (cout,"\t");
				copy ( it->begin(), it->end(), out_it );
				cout << endl;
			}
			cout << endl;
			
			std::vector<std::vector<double> > transformedStandardisedData;
			AHABCValuesPtr->fillTransformedStandardisedDataContainer(
									transformedStandardisedData);
			
			cout << "\nTransformed standardised data is = " << endl;
			for (std::vector< std::vector<double> >::iterator it 
								= transformedStandardisedData.begin();
							it < transformedStandardisedData.end();
							++it) {
				ostream_iterator<double> out_it (cout,"\t");
				copy ( it->begin(), it->end(), out_it );
				cout << endl;
			}
			cout << endl;
				
					
			std::vector < std::vector < double > > hm
					= AHABCValuesPtr->getHouseholderMatrixStandardised();
			
			cout << "\nHouseholder matrix for standardised data is = " << endl;
			for (std::vector< std::vector<double> >::iterator it = hm.begin();
							it < hm.end();
							++it) {
				ostream_iterator<double> out_it (cout,"\t");
				copy ( it->begin(), it->end(), out_it );
				cout << endl;
			}
			cout << endl;
			
			cout << "\nReference vector is " << endl;
			{
				ostream_iterator<double> out_it (cout,"\t");
				copy ( refvec.begin(), refvec.end(), out_it );
				cout << endl;
			}
						
			std::vector < double > htRefvec 
				= AHABCValuesPtr->standardiseAndHouseholderTransform(refvec);
				
			cout << "\nStandardised and transformed reference vector is " << endl;
			{
				ostream_iterator<double> out_it (cout,"\t");
				copy ( htRefvec.begin(), htRefvec.end(), out_it );
				cout << endl;
			}
			cout << endl;
		}

	}

	/* Output for largerSVD test should be 
	AHABCValuesPtr->getNumberParameters() = 2
	AHABCValuesPtr->getNumberSummaryStatistics() = 2
	AHABCValuesPtr->size() = 5                      

	Data is = 
	30      60      10      15
	20      30      12      11
	35      30      13      14.5
	25      30      12      12  
	35      20      13      15  

	sigma is : 
	30.4539 13.4817 0.84809 0.28898

	v is = 
	-0.0386156      0.9631  0.18081 0.195587
	0.995875        0.0323135       0.0775752       -0.0342107
	-0.0760904      0.0411766       0.613506        -0.784936 
	0.0308467       0.264004        -0.764788       -0.586899 


	Leading eigenvector (as a row vector) is = 
	0.0386156       -0.995875       0.0760904       -0.0308467

	Transformed data is = 
	-58.2958        -31.4636        16.9883 12.167
	-28.5302        -20.2713        15.841  9.44288
	-27.9828        -35.2424        17.9849 12.4792
	-28.3679        -25.2826        16.2239 10.2877
	-18.0395        -34.9423        17.1979 13.2982


	Householder matrix is = 
	0.0386156       -0.995875       0.0760904       -0.0308467
	-0.995875       -0.0316036      0.0788203       -0.0319533
	0.0760904       0.0788203       0.993978        0.00244141
	-0.0308467      -0.0319533      0.00244141      0.99901   


	Reference vector is 
	1       1       1       1

	Transformed reference vector is 
	-0.912016       -0.980612       1.15133 0.938652


	Means are 
	29      34      12      13.5


	SDs are 
	6.5192  15.1658 1.22474 1.87083


	Standardised data is = 
	0.153393        1.71439 -1.63299        0.801784
	-1.38054        -0.263752       0       -1.33631
	0.920358        -0.263752       0.816497        0.534522
	-0.613572       -0.263752       0       -0.801784       
	0.920358        -0.923133       0.816497        0.801784


	Leading eigenvector (as a row vector) is = 
	0.41096 -0.596409       0.665335        0.180937

	Transformed standardised data is = 
	-1.90085        -0.365556       0.687327        1.43279
	-0.651827       0.474073        -0.823094       -1.56015
	1.17549 -0.00542612     0.528316        0.456152        
	-0.239921       0.114573        -0.422047       -0.916559
	1.61711 -0.217664       0.0294987       0.587761         


	Householder matrix for standardised data is = 
	0.41096 -0.596409       0.665335        0.180937
	-0.596409       0.39613 0.673658        0.1832  
	0.665335        0.673658        0.248489        -0.204372
	0.180937        0.1832  -0.204372       0.944421         


	Reference vector is 
	1       1       1       1

	Standardised and transformed reference vector is 
	-7.65193        -5.57487        -5.18974        -5.65038



	*/
		
	
}		


void testAHABCValuesSelected() 
{
		
	cout << "\n\nExamples\n" << endl;
		
	std::vector < size_t > reqParams;
	std::vector < size_t > reqStats;
		
	cout << "Try with empty required positions" << endl;	
	try {	
		boost::shared_ptr < AHABCValues > 
			AHABCValuesPtr = getAHABCValuesSelected(reqParams, reqStats);
	}
	catch (std::invalid_argument& ia) {
		cout << "Failed:\n" << ia.what() << endl;
	}
	reqParams.push_back(0);
	try {	
		boost::shared_ptr < AHABCValues > 
			AHABCValuesPtr = getAHABCValuesSelected(reqParams, reqStats);
	}
	catch (std::invalid_argument& ia) {
		cout << "Failed:\n" << ia.what() << endl;
	}
	cout << "Try with invalid positions" << endl;	
	reqStats.push_back(0);
	try {	
		boost::shared_ptr < AHABCValues > 
			AHABCValuesPtr = getAHABCValuesSelected(reqParams, reqStats);
	}
	catch (std::invalid_argument& ia) {
		cout << "Failed:\n" << ia.what() << endl;
	}
	reqStats[0] = 1;
	try {	
		boost::shared_ptr < AHABCValues > 
			AHABCValuesPtr = getAHABCValuesSelected(reqParams, reqStats);
	}
	catch (std::invalid_argument& ia) {
		cout << "Failed:\n" << ia.what() << endl;
	}
	reqParams[0] = 1;
	reqStats.push_back(4);
	try {	
		boost::shared_ptr < AHABCValues > 
			AHABCValuesPtr = getAHABCValuesSelected(reqParams, reqStats);
	}
	catch (std::invalid_argument& ia) {
		cout << "Failed:\n" << ia.what() << endl;
	}
	reqParams.push_back(3);
	reqStats[1] = 3;
	try {	
		boost::shared_ptr < AHABCValues > 
			AHABCValuesPtr = getAHABCValuesSelected(reqParams, reqStats);
	}
	catch (std::invalid_argument& ia) {
		cout << "Failed:\n" << ia.what() << endl;
	}
	
	reqParams.resize(1);
	reqParams[0] = 2;
	cout << "Valid positions:" << endl;
	{	
		cout << "reqParams:" << endl;
		{
			ostream_iterator<size_t> out_it (cout,"\t");
			copy ( reqParams.begin(), reqParams.end(), out_it );
			cout << endl;
		}
		cout << "reqStats:" << endl;
		{
			ostream_iterator<size_t> out_it (cout,"\t");
			copy ( reqStats.begin(), reqStats.end(), out_it );
			cout << endl;
		}	
			
		boost::shared_ptr < AHABCValues > 
			AHABCValuesPtr = getAHABCValuesSelected(reqParams, reqStats);
		
		size_t np = AHABCValuesPtr->getNumberParameters();
		size_t ns = AHABCValuesPtr->getNumberSummaryStatistics();

		cout << "AHABCValuesPtr->getNumberParameters() = " 
						<< np << endl;
		cout << "AHABCValuesPtr->getNumberSummaryStatistics() = " 
							<< ns << endl;
		
		cout << "AHABCValuesPtr->size() = " 
							<< AHABCValuesPtr->size() << endl;
		
		assert(np = reqParams.size());
		assert(ns = reqStats.size());
		
		std::vector<std::vector<double> > data;
		AHABCValuesPtr->fillDataContainer(data);
		
		cout << "\nSelected Data is = " << endl;
		for (std::vector< std::vector<double> >::iterator it = data.begin();
						it < data.end();
						++it) {
			ostream_iterator<double> out_it (cout,"\t");
			copy ( it->begin(), it->end(), out_it );
			cout << endl;
		}
		cout << endl;
		
		std::vector < double > means
							= AHABCValuesPtr->getMeans();
		
		cout << "\nMeans are " << endl;
		{
			ostream_iterator<double> out_it (cout,"\t");
			copy ( means.begin(), means.end(), out_it );
			cout << endl;
		}
		cout << endl;
		
		std::vector < double > sds
							= AHABCValuesPtr->getSDs();
		
		cout << "\nSDs are " << endl;
		{
			ostream_iterator<double> out_it (cout,"\t");
			copy ( sds.begin(), sds.end(), out_it );
			cout << endl;
		}
		cout << endl;
	
	}
	
	cout << "compare to full data:" << endl;
	{	
		boost::shared_ptr < AHABCValues > 
			AHABCValuesPtr = getAHABCValuesLarger();
		
		size_t np = AHABCValuesPtr->getNumberParameters();
		size_t ns = AHABCValuesPtr->getNumberSummaryStatistics();

		cout << "AHABCValuesPtr->getNumberParameters() = " 
						<< np << endl;
		cout << "AHABCValuesPtr->getNumberSummaryStatistics() = " 
							<< ns << endl;
		
		cout << "AHABCValuesPtr->size() = " 
							<< AHABCValuesPtr->size() << endl;
		
		std::vector<std::vector<double> > data;
		AHABCValuesPtr->fillDataContainer(data);
		
		cout << "\nFull Data is = " << endl;
		for (std::vector< std::vector<double> >::iterator it = data.begin();
						it < data.end();
						++it) {
			ostream_iterator<double> out_it (cout,"\t");
			copy ( it->begin(), it->end(), out_it );
			cout << endl;
		}
		cout << endl;
		
		std::vector < double > means
							= AHABCValuesPtr->getMeans();
		
		cout << "\nMeans are " << endl;
		{
			ostream_iterator<double> out_it (cout,"\t");
			copy ( means.begin(), means.end(), out_it );
			cout << endl;
		}
		cout << endl;
		
		std::vector < double > sds
							= AHABCValuesPtr->getSDs();
		
		cout << "\nSDs are " << endl;
		{
			ostream_iterator<double> out_it (cout,"\t");
			copy ( sds.begin(), sds.end(), out_it );
			cout << endl;
		}
		cout << endl;
	
	}
	
	
	
}		


boost::shared_ptr < AHABCValues > getAHABCValuesLarger() 
{
		
	ParameterParticleSet pTmp;
	
	{
		double realvalues[] = {30.0, 60.}; 
		boost::shared_ptr < ParameterParticle > new_p( new ParameterParticle(realvalues, 2));
		pTmp.add( new_p );
	}
	
	{
		double realvalues[] = {20.,30.}; 
		boost::shared_ptr < ParameterParticle > new_p( new ParameterParticle(realvalues, 2));
		pTmp.add( new_p );
	}
	
	{
		double realvalues[] = {35.,30.}; 
		boost::shared_ptr < ParameterParticle > new_p( new ParameterParticle(realvalues, 2));
		pTmp.add( new_p );
	}
	
	SummaryStatisticSet sTmp;
	{
		double realvalues[] = {10.,15.,10.}; 
		boost::shared_ptr < SummaryStatistic > new_s( new SummaryStatistic(realvalues, 3));
		sTmp.add( new_s );
	}
	{
		double realvalues[] = {12.,11.,13.5}; 
		boost::shared_ptr < SummaryStatistic > new_s( new SummaryStatistic(realvalues, 3));
		sTmp.add( new_s );
	}
	{
		double realvalues[] = {13.,14.5,15.}; 
		boost::shared_ptr < SummaryStatistic > new_s( new SummaryStatistic(realvalues, 3));
		sTmp.add( new_s );
	}
		
	boost::shared_ptr< const aabc::ParameterParticleSet > p_ptr
		= boost::shared_ptr< const ParameterParticleSet >
			( new ParameterParticleSet(pTmp) );
		
	
	boost::shared_ptr < const SummaryStatisticSet > s_ptr
		= boost::shared_ptr< const SummaryStatisticSet >
			( new SummaryStatisticSet(sTmp) );
	
	boost::shared_ptr < AHABCValues > abc(new AHABCValues(p_ptr, s_ptr));
	
	return abc;
	
}		

// more values to ensure that n >= d
boost::shared_ptr < AHABCValues > getAHABCValuesLargerSVD() 
{
		
	ParameterParticleSet pTmp;
	
	{
		double realvalues[] = {30.0, 60.}; 
		boost::shared_ptr < ParameterParticle > new_p( new ParameterParticle(realvalues, 2));
		pTmp.add( new_p );
	}
	
	{
		double realvalues[] = {20.,30.}; 
		boost::shared_ptr < ParameterParticle > new_p( new ParameterParticle(realvalues, 2));
		pTmp.add( new_p );
	}
	
	{
		double realvalues[] = {35.,30.}; 
		boost::shared_ptr < ParameterParticle > new_p( new ParameterParticle(realvalues, 2));
		pTmp.add( new_p );
	}
	{
		double realvalues[] = {25.,30.}; 
		boost::shared_ptr < ParameterParticle > new_p( new ParameterParticle(realvalues, 2));
		pTmp.add( new_p );
	}
	{
		double realvalues[] = {35.,20.}; 
		boost::shared_ptr < ParameterParticle > new_p( new ParameterParticle(realvalues, 2));
		pTmp.add( new_p );
	}
	
	SummaryStatisticSet sTmp;
	{
		double realvalues[] = {10.,15.}; 
		boost::shared_ptr < SummaryStatistic > new_s( new SummaryStatistic(realvalues, 2));
		sTmp.add( new_s );
	}
	{
		double realvalues[] = {12.,11.}; 
		boost::shared_ptr < SummaryStatistic > new_s( new SummaryStatistic(realvalues, 2));
		sTmp.add( new_s );
	}
	{
		double realvalues[] = {13.,14.5}; 
		boost::shared_ptr < SummaryStatistic > new_s( new SummaryStatistic(realvalues, 2));
		sTmp.add( new_s );
	}
	{
		double realvalues[] = {12.,12.}; 
		boost::shared_ptr < SummaryStatistic > new_s( new SummaryStatistic(realvalues, 2));
		sTmp.add( new_s );
	}
	{
		double realvalues[] = {13.,15.}; 
		boost::shared_ptr < SummaryStatistic > new_s( new SummaryStatistic(realvalues, 2));
		sTmp.add( new_s );
	}	
	boost::shared_ptr< const aabc::ParameterParticleSet > p_ptr
		= boost::shared_ptr< const ParameterParticleSet >
			( new ParameterParticleSet(pTmp) );
		
	
	boost::shared_ptr < const SummaryStatisticSet > s_ptr
		= boost::shared_ptr< const SummaryStatisticSet >
			( new SummaryStatisticSet(sTmp) );
	
	boost::shared_ptr < AHABCValues > abc(new AHABCValues(p_ptr, s_ptr));
	
	return abc;
	
}		


boost::shared_ptr < AHABCValues > getAHABCValuesSmallScaled() 
{
		
		

	ParameterParticleSet pTmp;
	
	{
		double realvalues[] = {1000.0}; 
		boost::shared_ptr < ParameterParticle > new_p( new ParameterParticle(realvalues, 1));
		pTmp.add( new_p );
	}
	
	{
		double realvalues[] = {1500.}; 
		boost::shared_ptr < ParameterParticle > new_p( new ParameterParticle(realvalues, 1));
		pTmp.add( new_p );
	}
	
	{
		double realvalues[] = {2000.}; 
		boost::shared_ptr < ParameterParticle > new_p( new ParameterParticle(realvalues, 1));
		pTmp.add( new_p );
	}
	
	SummaryStatisticSet sTmp;
	{
		double realvalues[] = {11.}; 
		boost::shared_ptr < SummaryStatistic > new_s( new SummaryStatistic(realvalues, 1));
		sTmp.add( new_s );
	}
	{
		double realvalues[] = {16.}; 
		boost::shared_ptr < SummaryStatistic > new_s( new SummaryStatistic(realvalues, 1));
		sTmp.add( new_s );
	}
	{
		double realvalues[] = {19.}; 
		boost::shared_ptr < SummaryStatistic > new_s( new SummaryStatistic(realvalues, 1));
		sTmp.add( new_s );
	}
		
	boost::shared_ptr< const aabc::ParameterParticleSet > p_ptr
		= boost::shared_ptr< const ParameterParticleSet >
			( new ParameterParticleSet(pTmp) );
		
	
	boost::shared_ptr < const SummaryStatisticSet > s_ptr
		= boost::shared_ptr< const SummaryStatisticSet >
			( new SummaryStatisticSet(sTmp) );
	
	boost::shared_ptr < AHABCValues > abc(new AHABCValues(p_ptr, s_ptr));
	
	return abc;
}		
	
	boost::shared_ptr < AHABCValues > getAHABCValuesSmall() 
{
		
		

	ParameterParticleSet pTmp;
	
	{
		double realvalues[] = {10.0}; 
		boost::shared_ptr < ParameterParticle > new_p( new ParameterParticle(realvalues, 1));
		pTmp.add( new_p );
	}
	
	{
		double realvalues[] = {15.}; 
		boost::shared_ptr < ParameterParticle > new_p( new ParameterParticle(realvalues, 1));
		pTmp.add( new_p );
	}
	
	{
		double realvalues[] = {20.}; 
		boost::shared_ptr < ParameterParticle > new_p( new ParameterParticle(realvalues, 1));
		pTmp.add( new_p );
	}
	
	SummaryStatisticSet sTmp;
	{
		double realvalues[] = {11.}; 
		boost::shared_ptr < SummaryStatistic > new_s( new SummaryStatistic(realvalues, 1));
		sTmp.add( new_s );
	}
	{
		double realvalues[] = {16.}; 
		boost::shared_ptr < SummaryStatistic > new_s( new SummaryStatistic(realvalues, 1));
		sTmp.add( new_s );
	}
	{
		double realvalues[] = {19.}; 
		boost::shared_ptr < SummaryStatistic > new_s( new SummaryStatistic(realvalues, 1));
		sTmp.add( new_s );
	}
		
	boost::shared_ptr< const aabc::ParameterParticleSet > p_ptr
		= boost::shared_ptr< const ParameterParticleSet >
			( new ParameterParticleSet(pTmp) );
		
	
	boost::shared_ptr < const SummaryStatisticSet > s_ptr
		= boost::shared_ptr< const SummaryStatisticSet >
			( new SummaryStatisticSet(sTmp) );
	
	boost::shared_ptr < AHABCValues > abc(new AHABCValues(p_ptr, s_ptr));
	
	return abc;
}		
boost::shared_ptr < AHABCValues > getAHABCValuesAligned() 
{
		
	ParameterParticleSet pTmp;
	
	{
		double realvalues[] = {-5.}; 
		boost::shared_ptr < ParameterParticle > new_p( new ParameterParticle(realvalues, 1));
		pTmp.add( new_p );
	}
	
	{
		double realvalues[] = {-2.0}; 
		boost::shared_ptr < ParameterParticle > new_p( new ParameterParticle(realvalues, 1));
		pTmp.add( new_p );
	}
	
	{
		double realvalues[] = {1.0}; 
		boost::shared_ptr < ParameterParticle > new_p( new ParameterParticle(realvalues, 1));
		pTmp.add( new_p );
	}
	
	{
		double realvalues[] = {3.0}; 
		boost::shared_ptr < ParameterParticle > new_p( new ParameterParticle(realvalues, 1));
		pTmp.add( new_p );
	}
	
	{
		double realvalues[] = {4.5}; 
		boost::shared_ptr < ParameterParticle > new_p( new ParameterParticle(realvalues, 1));
		pTmp.add( new_p );
	}
	
	SummaryStatisticSet sTmp;
	{
		double realvalues[] = {-6.0, 5.5}; 
		boost::shared_ptr < SummaryStatistic > new_s( new SummaryStatistic(realvalues, 2));
		sTmp.add( new_s );
	}
	{
		double realvalues[] = {-2.5, 3.0}; 
		boost::shared_ptr < SummaryStatistic > new_s( new SummaryStatistic(realvalues, 2));
		sTmp.add( new_s );
	}
	{
		double realvalues[] = {0.5, -1.5}; 
		boost::shared_ptr < SummaryStatistic > new_s( new SummaryStatistic(realvalues, 2));
		sTmp.add( new_s );
	}
	
	{
		double realvalues[] = {2.5, -2.5}; 
		boost::shared_ptr < SummaryStatistic > new_s( new SummaryStatistic(realvalues, 2));
		sTmp.add( new_s );
	}
	
	{
		double realvalues[] = {5.0, -4.0}; 
		boost::shared_ptr < SummaryStatistic > new_s( new SummaryStatistic(realvalues, 2));
		sTmp.add( new_s );
	}
		
	boost::shared_ptr< const aabc::ParameterParticleSet > p_ptr
		= boost::shared_ptr< const ParameterParticleSet >
			( new ParameterParticleSet(pTmp) );
		
	
	boost::shared_ptr < const SummaryStatisticSet > s_ptr
		= boost::shared_ptr< const SummaryStatisticSet >
			( new SummaryStatisticSet(sTmp) );
	
	boost::shared_ptr < AHABCValues > abc(new AHABCValues(p_ptr, s_ptr));
	
	return abc;
	
}		

boost::shared_ptr < AHABCValues > getAHABCValuesSelected(
		const std::vector < size_t >& reqParams,
		const std::vector < size_t >& reqStats) 
{
		
	ParameterParticleSet pTmp;
	
	{
		double realvalues[] = {30.0, 60.}; 
		boost::shared_ptr < ParameterParticle > new_p( new ParameterParticle(realvalues, 2));
		pTmp.add( new_p );
	}
	
	{
		double realvalues[] = {20.,30.}; 
		boost::shared_ptr < ParameterParticle > new_p( new ParameterParticle(realvalues, 2));
		pTmp.add( new_p );
	}
	
	{
		double realvalues[] = {35.,30.}; 
		boost::shared_ptr < ParameterParticle > new_p( new ParameterParticle(realvalues, 2));
		pTmp.add( new_p );
	}
	
	SummaryStatisticSet sTmp;
	{
		double realvalues[] = {10.,15.,10.}; 
		boost::shared_ptr < SummaryStatistic > new_s( new SummaryStatistic(realvalues, 3));
		sTmp.add( new_s );
	}
	{
		double realvalues[] = {12.,11.,13.5}; 
		boost::shared_ptr < SummaryStatistic > new_s( new SummaryStatistic(realvalues, 3));
		sTmp.add( new_s );
	}
	{
		double realvalues[] = {13.,14.5,15.}; 
		boost::shared_ptr < SummaryStatistic > new_s( new SummaryStatistic(realvalues, 3));
		sTmp.add( new_s );
	}
		
	boost::shared_ptr< const aabc::ParameterParticleSet > p_ptr
		= boost::shared_ptr< const ParameterParticleSet >
			( new ParameterParticleSet(pTmp) );
		
	
	boost::shared_ptr < const SummaryStatisticSet > s_ptr
		= boost::shared_ptr< const SummaryStatisticSet >
			( new SummaryStatisticSet(sTmp) );
	
	boost::shared_ptr < AHABCValues > abc(
			new AHABCValues(p_ptr, reqParams, s_ptr, reqStats));
	
	return abc;
	
}		
