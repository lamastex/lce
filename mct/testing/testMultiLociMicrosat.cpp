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
\brief Basic testing of MultiLociMicrosat class.

*/

#include "config.h"

#include "multi_loci_microsat.hpp"
#include "mct_microsat_locus.hpp"

#include "site_seq.hpp"

#include <mscplusplus/population_structure.hpp>

#include <boost/shared_ptr.hpp>


#include <iostream>
#include <iterator>
#include <stdexcept>



using namespace std;
using namespace mct;


void testMultiLociMicrosat();



int main()
{
	testMultiLociMicrosat();
	
	
	return 0;
	
}

	
void testMultiLociMicrosat()
{
	try {
		{
			size_t nsam = 2;
			// make an empty multiloci microsat table for nsam = 2
			boost::shared_ptr< hudson_ms::PopulationStructure > pop
				(new hudson_ms::PopulationStructure(nsam));
			MultiLociMicrosat mlms(pop);
			
			cout << "mlms has nsam " << mlms.nsam() << " and totalNsites " << mlms.totalNsites() << std::endl;
			cout << "mlms empty =  " << mlms.empty() << std::endl;
					
			try {
				cout << "toString is" << endl;
				cout << mlms.toString() << endl;
			}
			catch (std::invalid_argument& ia) {
				std::cout << "That resulted in exception\n" 
					+ std::string(ia.what()) << std::endl;
				throw;
			}
			std::cout << "Try getting averages (should give nans)" << endl;
			try {
				boost::shared_ptr<const mct::SummaryStatistic> avs = mlms.getAverageRepeats();
				cout << "Averages are" << endl;
				std::cout << (*avs) << std::endl;
					
			}
			catch (std::invalid_argument& ia) {
				std::cout << "That resulted in exception\n" 
					+ std::string(ia.what()) << std::endl;
				throw;
			}
			std::cout << "Try getting covariance " << endl;
			try {
				boost::shared_ptr<const mct::SummaryStatistic> cov = mlms.getCovariances();
				cout << "Covariances are" << endl;
				std::cout << (*cov) << std::endl;
					
			}
			catch (std::invalid_argument& ia) {
				std::cout << "That resulted in exception\n" 
					+ std::string(ia.what()) << std::endl;
				throw;
			}
			
			std::cout << "Try getting delta mus (this should be okay)" << endl;
			try {
				boost::shared_ptr<const mct::SummaryStatistic> 
					deltaMus = mlms.distanceGoldsteinDeltaMuSqBetweenPop();
				cout << "Delta mus are" << endl;
				std::cout << (*deltaMus) << std::endl;
			}
			catch (std::runtime_error& re) {
				std::cout << "That resulted in exception\n" 
					+ std::string(re.what()) << std::endl;
				throw;
			}
			
			// this should fail because nsam is wrong
			std::cout << "Try adding a site_seq (this should fail): adding" << endl; 
			try {
				vector < site_seq > seqs;
				{
					int myints[] = {30}; 
					mct::site_seq seq (myints, myints + sizeof(myints) / sizeof(int) );
					copy(seq.begin(), seq.end(), ostream_iterator<int>(std::cout, "\t"));
					std::cout << std::endl;
					
					seqs.push_back(seq);
				}
				mlms.add(seqs);
				throw std::logic_error("Should not be able to do that");
				
			
			}
			catch (std::invalid_argument& ia) {
				std::cout << "That resulted in exception\n" 
					+ std::string(ia.what()) << std::endl;
			}		
			
			std::cout << "Try adding a site_seq (this should be okay): adding" << endl; 
			try {
				vector < site_seq > seqs;
				{
					int myints[] = {30}; 
					mct::site_seq seq (myints, myints + sizeof(myints) / sizeof(int) );
					copy(seq.begin(), seq.end(), ostream_iterator<int>(std::cout, "\t"));
					std::cout << std::endl;
					
					seqs.push_back(seq);
				}
				{
					int myints[] = {60}; 
					mct::site_seq seq (myints, myints + sizeof(myints) / sizeof(int) );
					copy(seq.begin(), seq.end(), ostream_iterator<int>(std::cout, "\t"));
					std::cout << std::endl;
					
					seqs.push_back(seq);
				}
				mlms.add(seqs);
				cout << "mlms has nsam " << mlms.nsam() << " and totalNsites " << mlms.totalNsites() << std::endl;
				cout << "mlms empty =  " << mlms.empty() << std::endl;
				
				cout << "mlms.nsites() is " << std::endl;
				vector < size_t > nsites = mlms.nsites();
				copy(nsites.begin(), nsites.end(), ostream_iterator<size_t>(std::cout, "\t"));
					std::cout << std::endl;
				
				cout << "toString is" << endl;
				cout << mlms.toString() << endl;
			
			}
			catch (std::invalid_argument& ia) {
				std::cout << "That resulted in exception\n" 
					+ std::string(ia.what()) << std::endl;
				throw;
			}		
			
			std::cout << "Try adding a site_seq (this should fail): adding" << endl; 
			try {
				vector < site_seq > seqs;
				{
					int myints[] = {30}; 
					mct::site_seq seq (myints, myints + sizeof(myints) / sizeof(int) );
					copy(seq.begin(), seq.end(), ostream_iterator<int>(std::cout, "\t"));
					std::cout << std::endl;
					
					seqs.push_back(seq);
				}
				{
					int myints[] = {60}; 
					mct::site_seq seq (myints, myints + sizeof(myints) / sizeof(int) );
					copy(seq.begin(), seq.end(), ostream_iterator<int>(std::cout, "\t"));
					std::cout << std::endl;
					
					seqs.push_back(seq);
				}
				{
					int myints[] = {90}; 
					mct::site_seq seq (myints, myints + sizeof(myints) / sizeof(int) );
					copy(seq.begin(), seq.end(), ostream_iterator<int>(std::cout, "\t"));
					std::cout << std::endl;
					
					seqs.push_back(seq);
				}
				mlms.add(seqs);
				throw std::logic_error("Should not be able to do that");
			
			}
			catch (std::invalid_argument& ia) {
				std::cout << "That resulted in exception\n" 
					+ std::string(ia.what()) << std::endl;
			}		
			
			std::cout << "Try adding a site_seq (this should be okay): adding" << endl; 
			try {
				vector < site_seq > seqs;
				{
					int myints[] = {50, 45}; 
					mct::site_seq seq (myints, myints + sizeof(myints) / sizeof(int) );
					copy(seq.begin(), seq.end(), ostream_iterator<int>(std::cout, "\t"));
					std::cout << std::endl;
					
					seqs.push_back(seq);
				}
				{
					int myints[] = {40, 35}; 
					mct::site_seq seq (myints, myints + sizeof(myints) / sizeof(int) );
					copy(seq.begin(), seq.end(), ostream_iterator<int>(std::cout, "\t"));
					std::cout << std::endl;
					
					seqs.push_back(seq);
				}
				
				mlms.add(seqs);
				cout << "mlms has nsam " << mlms.nsam() << " and totalNsites " << mlms.totalNsites() << std::endl;
				cout << "mlms empty =  " << mlms.empty() << std::endl;
				
				cout << "mlms.nsites() is " << std::endl;
				vector < size_t > nsites = mlms.nsites();
				copy(nsites.begin(), nsites.end(), ostream_iterator<size_t>(std::cout, "\t"));
					std::cout << std::endl;
				
				cout << "toString is" << endl;
				cout << mlms.toString() << endl;
			
			}
			catch (std::invalid_argument& ia) {
				std::cout << "That resulted in exception\n" 
					+ std::string(ia.what()) << std::endl;
				throw;
			}
			
			std::cout << "Try adding a site_seq (this should be okay): adding" << endl; 
			try {
				vector < site_seq > seqs;
				{
					int myints[] = {10}; 
					mct::site_seq seq (myints, myints + sizeof(myints) / sizeof(int) );
					copy(seq.begin(), seq.end(), ostream_iterator<int>(std::cout, "\t"));
					std::cout << std::endl;
					
					seqs.push_back(seq);
				}
				{
					int myints[] = {20}; 
					mct::site_seq seq (myints, myints + sizeof(myints) / sizeof(int) );
					copy(seq.begin(), seq.end(), ostream_iterator<int>(std::cout, "\t"));
					std::cout << std::endl;
					
					seqs.push_back(seq);
				}
				
				mlms.add(seqs);
				cout << "mlms has nsam " << mlms.nsam() << " and totalNsites " << mlms.totalNsites() << std::endl;
				cout << "mlms.nloci() =  " << mlms.nloci() << std::endl;
				
				cout << "mlms.nsites() is " << std::endl;
				vector < size_t > nsites = mlms.nsites();
				copy(nsites.begin(), nsites.end(), ostream_iterator<size_t>(std::cout, "\t"));
					std::cout << std::endl;
				
				cout << "toString is" << endl;
				cout << mlms.toString() << endl;
			
			}
			catch (std::invalid_argument& ia) {
				std::cout << "That resulted in exception\n" 
					+ std::string(ia.what()) << std::endl;
				throw;
			}
			
			std::cout << "and another the same: adding" << endl; 
			try {
				vector < site_seq > seqs;
				{
					int myints[] = {10}; 
					mct::site_seq seq (myints, myints + sizeof(myints) / sizeof(int) );
					copy(seq.begin(), seq.end(), ostream_iterator<int>(std::cout, "\t"));
					std::cout << std::endl;
					
					seqs.push_back(seq);
				}
				{
					int myints[] = {20}; 
					mct::site_seq seq (myints, myints + sizeof(myints) / sizeof(int) );
					copy(seq.begin(), seq.end(), ostream_iterator<int>(std::cout, "\t"));
					std::cout << std::endl;
					
					seqs.push_back(seq);
				}
				
				mlms.add(seqs);
				cout << "mlms has nsam " << mlms.nsam() << " and totalNsites " << mlms.totalNsites() << std::endl;
				cout << "mlms.nloci() =  " << mlms.nloci() << std::endl;
				
				cout << "mlms.nsites() is " << std::endl;
				vector < size_t > nsites = mlms.nsites();
				copy(nsites.begin(), nsites.end(), ostream_iterator<size_t>(std::cout, "\t"));
					std::cout << std::endl;
				
				cout << "toString is" << endl;
				cout << mlms.toString() << endl;
			
			}
			catch (std::invalid_argument& ia) {
				std::cout << "That resulted in exception\n" 
					+ std::string(ia.what()) << std::endl;
				throw;
			}
			std::cout << "Try getting averages" << endl;
			try {
				boost::shared_ptr<const mct::SummaryStatistic> avs = mlms.getAverageRepeats();
				cout << "averages are" << endl;
				std::cout << (*avs) << std::endl;
					
			}
			catch (std::invalid_argument& ia) {
				std::cout << "That resulted in exception\n" 
					+ std::string(ia.what()) << std::endl;
				throw;
			}
			std::cout << "Try getting covariance " << endl;
			try {
				boost::shared_ptr<const mct::SummaryStatistic> cov = mlms.getCovariances();
				cout << "Covariances are" << endl;
				std::cout << (*cov) << std::endl;
					
			}
			catch (std::invalid_argument& ia) {
				std::cout << "That resulted in exception\n" 
					+ std::string(ia.what()) << std::endl;
				throw;
			}
			
			std::cout << "Try getting delta mus (this should be okay)" << endl;
			try {
				boost::shared_ptr<const mct::SummaryStatistic> 
					deltaMus = mlms.distanceGoldsteinDeltaMuSqBetweenPop();
				cout << "Delta mus are" << endl;
				std::cout << (*deltaMus) << std::endl;
			}
			catch (std::runtime_error& re) {
				std::cout << "That resulted in exception\n" 
					+ std::string(re.what()) << std::endl;
				throw;
			}
		}
		cout << "Now try a larger one" << endl;
		
		{
			size_t nsam = 3;
			
			boost::shared_ptr< hudson_ms::PopulationStructure > pop
				(new hudson_ms::PopulationStructure(nsam));
			
			MultiLociMicrosat mlms( pop );
			
			cout << "mlms has nsam " << mlms.nsam() << " and totalNsites " << mlms.totalNsites() << std::endl;
			cout << "mlms empty =  " << mlms.empty() << std::endl;
					
			std::cout << "Try getting averages (should give nans)" << endl;
			try {
				boost::shared_ptr<const mct::SummaryStatistic> avs = mlms.getAverageRepeats();
				cout << "Averages are" << endl;
				std::cout << (*avs) << std::endl;
					
			}
			catch (std::invalid_argument& ia) {
				std::cout << "That resulted in exception\n" 
					+ std::string(ia.what()) << std::endl;
				throw;
			}
			std::cout << "Try getting covariance " << endl;
			try {
				boost::shared_ptr<const mct::SummaryStatistic> cov = mlms.getCovariances();
				cout << "Covariances are" << endl;
				std::cout << (*cov) << std::endl;
					
			}
			catch (std::invalid_argument& ia) {
				std::cout << "That resulted in exception\n" 
					+ std::string(ia.what()) << std::endl;
					
				throw;
			}
			
			std::cout << "Try adding a site_seq : adding" << endl; 
			try {
				vector < site_seq > seqs;
				{
					int myints[] = {30}; 
					mct::site_seq seq (myints, myints + sizeof(myints) / sizeof(int) );
					copy(seq.begin(), seq.end(), ostream_iterator<int>(std::cout, "\t"));
					std::cout << std::endl;
					
					seqs.push_back(seq);
				}
				{
					int myints[] = {60}; 
					mct::site_seq seq (myints, myints + sizeof(myints) / sizeof(int) );
					copy(seq.begin(), seq.end(), ostream_iterator<int>(std::cout, "\t"));
					std::cout << std::endl;
					
					seqs.push_back(seq);
				}
				{
					int myints[] = {90}; 
					mct::site_seq seq (myints, myints + sizeof(myints) / sizeof(int) );
					copy(seq.begin(), seq.end(), ostream_iterator<int>(std::cout, "\t"));
					std::cout << std::endl;
					
					seqs.push_back(seq);
				}
				mlms.add(seqs);
				cout << "mlms has nsam " << mlms.nsam() << " and totalNsites " << mlms.totalNsites() << std::endl;
				cout << "mlms empty =  " << mlms.empty() << std::endl;
				
				cout << "mlms.nsites() is " << std::endl;
				vector < size_t > nsites = mlms.nsites();
				copy(nsites.begin(), nsites.end(), ostream_iterator<size_t>(std::cout, "\t"));
					std::cout << std::endl;
				
				cout << "toString is" << endl;
				cout << mlms.toString() << endl;
			
			}
			catch (std::invalid_argument& ia) {
				std::cout << "That resulted in exception\n" 
					+ std::string(ia.what()) << std::endl;
				
				throw;
			}		
			
			std::cout << "Try adding a site_seq : adding" << endl; 
			try {
				vector < site_seq > seqs;
				{
					int myints[] = {20}; 
					mct::site_seq seq (myints, myints + sizeof(myints) / sizeof(int) );
					copy(seq.begin(), seq.end(), ostream_iterator<int>(std::cout, "\t"));
					std::cout << std::endl;
					
					seqs.push_back(seq);
				}
				{
					int myints[] = {30}; 
					mct::site_seq seq (myints, myints + sizeof(myints) / sizeof(int) );
					copy(seq.begin(), seq.end(), ostream_iterator<int>(std::cout, "\t"));
					std::cout << std::endl;
					
					seqs.push_back(seq);
				}
				{
					int myints[] = {25}; 
					mct::site_seq seq (myints, myints + sizeof(myints) / sizeof(int) );
					copy(seq.begin(), seq.end(), ostream_iterator<int>(std::cout, "\t"));
					std::cout << std::endl;
					
					seqs.push_back(seq);
				}
				mlms.add(seqs);
				cout << "mlms has nsam " << mlms.nsam() << " and totalNsites " << mlms.totalNsites() << std::endl;
				cout << "mlms empty =  " << mlms.empty() << std::endl;
				
				cout << "mlms.nsites() is " << std::endl;
				vector < size_t > nsites = mlms.nsites();
				copy(nsites.begin(), nsites.end(), ostream_iterator<size_t>(std::cout, "\t"));
					std::cout << std::endl;
				
				cout << "toString is" << endl;
				cout << mlms.toString() << endl;
			
			}
			catch (std::invalid_argument& ia) {
				std::cout << "That resulted in exception\n" 
					+ std::string(ia.what()) << std::endl;
				throw;
			}	
			
			std::cout << "Try adding a site_seq : adding" << endl; 
			try {
				vector < site_seq > seqs;
				{
					int myints[] = {35}; 
					mct::site_seq seq (myints, myints + sizeof(myints) / sizeof(int) );
					copy(seq.begin(), seq.end(), ostream_iterator<int>(std::cout, "\t"));
					std::cout << std::endl;
					
					seqs.push_back(seq);
				}
				{
					int myints[] = {30}; 
					mct::site_seq seq (myints, myints + sizeof(myints) / sizeof(int) );
					copy(seq.begin(), seq.end(), ostream_iterator<int>(std::cout, "\t"));
					std::cout << std::endl;
					
					seqs.push_back(seq);
				}
				{
					int myints[] = {40}; 
					mct::site_seq seq (myints, myints + sizeof(myints) / sizeof(int) );
					copy(seq.begin(), seq.end(), ostream_iterator<int>(std::cout, "\t"));
					std::cout << std::endl;
					
					seqs.push_back(seq);
				}
				mlms.add(seqs);
				cout << "mlms has nsam " << mlms.nsam() << " and totalNsites " << mlms.totalNsites() << std::endl;
				cout << "mlms empty =  " << mlms.empty() << std::endl;
				
				cout << "mlms.nsites() is " << std::endl;
				vector < size_t > nsites = mlms.nsites();
				copy(nsites.begin(), nsites.end(), ostream_iterator<size_t>(std::cout, "\t"));
					std::cout << std::endl;
				
				cout << "toString is" << endl;
				cout << mlms.toString() << endl;
			
			}
			catch (std::invalid_argument& ia) {
				std::cout << "That resulted in exception\n" 
					+ std::string(ia.what()) << std::endl;
				throw;
			}	
			
			std::cout << "Try adding a site_seq : adding" << endl; 
			try {
				vector < site_seq > seqs;
				{
					int myints[] = {10}; 
					mct::site_seq seq (myints, myints + sizeof(myints) / sizeof(int) );
					copy(seq.begin(), seq.end(), ostream_iterator<int>(std::cout, "\t"));
					std::cout << std::endl;
					
					seqs.push_back(seq);
				}
				{
					int myints[] = {15}; 
					mct::site_seq seq (myints, myints + sizeof(myints) / sizeof(int) );
					copy(seq.begin(), seq.end(), ostream_iterator<int>(std::cout, "\t"));
					std::cout << std::endl;
					
					seqs.push_back(seq);
				}
				{
					int myints[] = {10}; 
					mct::site_seq seq (myints, myints + sizeof(myints) / sizeof(int) );
					copy(seq.begin(), seq.end(), ostream_iterator<int>(std::cout, "\t"));
					std::cout << std::endl;
					
					seqs.push_back(seq);
				}
				mlms.add(seqs);
				cout << "mlms has nsam " << mlms.nsam() << " and totalNsites " << mlms.totalNsites() << std::endl;
				cout << "mlms empty =  " << mlms.empty() << std::endl;
				
				cout << "mlms.nsites() is " << std::endl;
				vector < size_t > nsites = mlms.nsites();
				copy(nsites.begin(), nsites.end(), ostream_iterator<size_t>(std::cout, "\t"));
					std::cout << std::endl;
				
				cout << "toString is" << endl;
				cout << mlms.toString() << endl;
			
			}
			catch (std::invalid_argument& ia) {
				std::cout << "That resulted in exception\n" 
					+ std::string(ia.what()) << std::endl;
				throw;
			}	
			std::cout << "Try getting averages" << endl;
			try {
				boost::shared_ptr<const mct::SummaryStatistic> avs = mlms.getAverageRepeats();
				cout << "averages are" << endl;
				std::cout << (*avs) << std::endl;
					
			}
			catch (std::invalid_argument& ia) {
				std::cout << "That resulted in exception\n" 
					+ std::string(ia.what()) << std::endl;
				throw;
			}
			std::cout << "Try getting covariance " << endl;
			try {
				boost::shared_ptr<const mct::SummaryStatistic> cov = mlms.getCovariances();
				cout << "Covariances are" << endl;
				std::cout << (*cov) << std::endl;
					
			}
			catch (std::invalid_argument& ia) {
				std::cout << "That resulted in exception\n" 
					+ std::string(ia.what()) << std::endl;
				throw;
			}
			
			std::cout << "(should fail) Try getting stuff at locus ";
			try {
				
				size_t index = 4;
				cout << index << endl;
				boost::shared_ptr <MCTMicrosatLocus > m = mlms.at(index);
				
				throw std::logic_error("Should not be able to do that");
				
					
			}
			catch (std::exception& ee) {
				std::cout << "That resulted in exception\n" 
					+ std::string(ee.what()) << std::endl;
			}
			
			std::cout << "Try getting stuff at locus ";
			try {
				
				size_t index = 0;
				cout << index << endl;
				boost::shared_ptr <MCTMicrosatLocus > m = mlms.at(index);
				
				cout << "MCTMicrosatLocus toString is" << endl;
				cout << m->toString();
					
			}
			catch (std::invalid_argument& ia) {
				std::cout << "That resulted in exception\n" 
					+ std::string(ia.what()) << std::endl;
				throw;
			}
			
			std::cout << "Try getting stuff at locus ";
			try {
				
				size_t index = 1;
				cout << index << endl;
				boost::shared_ptr <MCTMicrosatLocus > m = mlms.at(index);
				
				cout << "MCTMicrosatLocus toString is" << endl;
				cout << m->toString();
					
			}
			catch (std::invalid_argument& ia) {
				std::cout << "That resulted in exception\n" 
					+ std::string(ia.what()) << std::endl;
				throw;
			}
			
			std::cout << "Try getting stuff at locus ";
			try {
				
				size_t index = 3;
				cout << index << endl;
				boost::shared_ptr <MCTMicrosatLocus > m = mlms.at(index);
				
				cout << "MCTMicrosatLocus toString is" << endl;
				cout << m->toString();
					
			}
			catch (std::invalid_argument& ia) {
				std::cout << "That resulted in exception\n" 
					+ std::string(ia.what()) << std::endl;
				throw;
			}
			
			
		}
		
		cout << "try a MultiLociMicrosat for 3 subpops but with no locus data" << endl;
			
		{
			size_t nsam = 2;
			// make an empty multiloci microsat table for nsam = 2
			boost::shared_ptr< hudson_ms::PopulationStructure > pop
				(new hudson_ms::PopulationStructure(nsam));
			pop->add(nsam);
			pop->add(nsam);
			MultiLociMicrosat mlms(pop);
			
			cout << "mlms has nsam " << mlms.nsam() << " and totalNsites " << mlms.totalNsites() << std::endl;
			cout << "mlms empty =  " << mlms.empty() << std::endl;
					
			try {
				cout << "toString is" << endl;
				cout << mlms.toString() << endl;
			}
			catch (std::invalid_argument& ia) {
				std::cout << "That resulted in exception\n" 
					+ std::string(ia.what()) << std::endl;
				throw;
			}
			std::cout << "Try getting averages (should give nans)" << endl;
			try {
				boost::shared_ptr<const mct::SummaryStatistic> avs = mlms.getAverageRepeats();
				cout << "Averages are" << endl;
				std::cout << (*avs) << std::endl;
					
			}
			catch (std::invalid_argument& ia) {
				std::cout << "That resulted in exception\n" 
					+ std::string(ia.what()) << std::endl;
				throw;
			}
			std::cout << "Try getting covariance " << endl;
			try {
				boost::shared_ptr<const mct::SummaryStatistic> cov = mlms.getCovariances();
				cout << "Covariances are" << endl;
				std::cout << (*cov) << std::endl;
					
			}
			catch (std::invalid_argument& ia) {
				std::cout << "That resulted in exception\n" 
					+ std::string(ia.what()) << std::endl;
				throw;
			}
			
			std::cout << "Try getting delta mus (this should give an error)" << endl;
			try {
				boost::shared_ptr<const mct::SummaryStatistic> 
					deltaMus = mlms.distanceGoldsteinDeltaMuSqBetweenPop();
				
				throw std::logic_error("Should not be able to do that");
					
			}
			catch (std::runtime_error& re) {
				std::cout << "That resulted in exception\n" 
					+ std::string(re.what()) << std::endl;
				
			}
		}
		std::cout << "\n\nEnd of testing\n" << std::endl;
	
	}
	catch (std::exception& e) {
	
		std::cout << "std::exception::\n" + std::string(e.what()) << std::endl;
	}
	catch (...) {
	
		std::cout << "Unknown error\n" << std::endl;
	}
	
}
