/* MCT - Markov Chains on Trees.


   Copyright (C) 2011, 2012 Jennifer Harlow

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
\brief Tools for doing various KDE routines

*/


#include "testing_toolsKDE.hpp"
#include "testing_tools.hpp"

//#include "log_memory.hpp"

//#include "test_routines.hpp"

#include "real_kde_mid_estimate.hpp"
#include "mcmc_kde.hpp"
#include "type_kde.hpp"
#include "functionestimator_kde.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <cmath>
//#include <numeric>
#include <iterator>
#include <cassert>


using namespace std;
using namespace subpavings;
using namespace subpavings::kde;	



std::vector < std::vector < cxsc::real > >& doGrid(size_t gps,
			const cxsc::ivector& box,
			std::vector < std::vector < cxsc::real > >& pts);

void outputGrid(const std::string& filename,
			const std::vector < std::vector < cxsc::real > >& pts);
			
void outputDensities(const std::string& filename,
					const std::vector < std::vector < cxsc::real > >& densities);

void getKDEGridDensitiesTwoD(const TypeKDE& tkde,
						const vector < vector < real > >& pts,
						std::vector < std::vector < real > >& kdeGridDensities);

void outputDensities(const std::string& filename,
					const std::vector < real >& densities);
					
void getKDEGridDensitiesOneD(const TypeKDE& tkde,
						const vector < vector < real > >& pts,
						std::vector < real > & kdeGridDensities);

						
// get an approximation to a kde
boost::shared_ptr < subpavings::PiecewiseConstantFunction > getKDEapprox(
				const std::vector < std::vector < double > >& container,
				const cxsc::ivector& box,
				const std::string& prefix,
				const std::string& spec,
				cxsc::real maxMeasure,
				const std::string& logFilename)
{
	try {
		
		string kdeLogFilename = prefix;
		{
			ostringstream oss;
			oss << "_log" << spec << ".txt";
			kdeLogFilename += oss.str();
			
		}
		string kdeResultsFilename = prefix;
		{
			ostringstream oss;
			oss << "_trace" << spec << ".txt";
			kdeResultsFilename += oss.str();
			
		}
		string gridFilename = prefix;
		{
			ostringstream oss;
			oss << "_grid" << spec << ".txt";
			gridFilename += oss.str();
			
		}
		string densitiesFilename = prefix;
		{
			ostringstream oss;
			oss << "_densities" << spec << ".txt";
			densitiesFilename += oss.str();
			
		}
		
		// for the mcmc kde
		size_t warm=5000;
		size_t step=10; 
		size_t num_batch = 50;
		size_t size_batch = 500; // 50 * 500 = 25,000 iterations
		unsigned int kdeSeed = 1234;
		
		// for the function estimation pq
		size_t maxLeaves = 1000000; 
		
		LOGGING_LEVEL logging = NOLOG;
				
		
		// log values used for kde approx to file
		{
			std::ostringstream stm;
			stm << "\n\nMCMC KDE using";
			stm << "\n\twarm\t" << warm;
			stm << "\n\tstep\t" << step;
			stm << "\tnum_batch\t" << num_batch;
			stm << "\tsize_batch\t" << size_batch;
			stm << "\tseed\t" << kdeSeed;
			stm << "\n\nand function estimation using";
			stm << "\n\tmaxMeasure\t" << maxMeasure;
			stm << "\n\tmaxLeaves\t" << maxLeaves;
			outputAHABCLogFileString(logFilename, stm.str());
		}
		
		clock_t start = clock();	
	
		MCMCKDE mkde(container);
		
		mkde.df(size_batch, num_batch, warm, step,
			kdeLogFilename,
			kdeResultsFilename,
			kdeSeed);
		
		clock_t end = clock();	
		
		double timing = (static_cast<double>(end-start)/CLOCKS_PER_SEC);	
		cout << "time to do KDE = " << timing << endl;
		
		/* use the mid-point estimator */
		RealKDEMidEstimator midest(mkde);
		FunctionEstimatorKDE estKDE(box, midest);
		
		
		start = clock();
		
		bool successfulPQ = estKDE.prioritySplit(maxMeasure, maxLeaves, 
											logging, kdeSeed);
		end = clock();
		timing = (static_cast<double>(end-start)/CLOCKS_PER_SEC);	
		cout << "time to do function approximation  = " << timing << endl;
		
		size_t l = estKDE.getRootLeaves();
		
		cout << "Number of leaves is " << l << endl;			
		if (!successfulPQ) {
			cout << "Note that PQ unsuccessful" << endl;
		}
					
		PiecewiseConstantFunction pcf 
				= estKDE.makePiecewiseConstantFunction();
		
		boost::shared_ptr < subpavings::PiecewiseConstantFunction > 
			pcfPtr(new subpavings::PiecewiseConstantFunction(pcf));
			
		
		int d = VecLen(box);
		if (d <= 2) {
			
			/*get a grid so that we can ouput stuff to make a visual of 
			* the kde itself */
			size_t gps = 100; // points to have along each dim of the box, for gridding

			/* grid the box */
			vector < vector < real > > pts;
		
			doGrid(gps, box, pts);
		
			outputGrid(gridFilename, pts);
		
		
			cout << "\nGet kde densities" << endl;
			
			if (d ==2 ) {
				 std::vector < std::vector < real > > kdeGridDensities;
			
							
				getKDEGridDensitiesTwoD(mkde,
								pts,
								kdeGridDensities);
				
				outputDensities(densitiesFilename, kdeGridDensities);
			}
			
			if (d == 1) {
				
				std::vector < real > kdeGridDensities;
				
				getKDEGridDensitiesOneD(mkde,
						pts,
						kdeGridDensities);
				
				outputDensities(densitiesFilename, kdeGridDensities);
			}
			
		}
		
					
		return pcfPtr;
		
	}
	catch (std::exception& e) {
	
		std::cout << "std::exception::\n" + std::string(e.what()) << std::endl;
		throw;
	}
	catch (...) {
	
		std::cout << "Unknown error\n" << std::endl;
		throw;
	}
}


// just do a kde
void doKDE(
				const std::vector < std::vector < double > >& container,
				const cxsc::ivector& box,
				const std::string& prefix,
				const std::string& spec,
				cxsc::real maxMeasure,
				size_t id,
				const std::string& logFilename)
{
	try {
		
		string kdeLogFilename = prefix;
		{
			ostringstream oss;
			oss << "_log" << spec << "_" << id <<".txt";
			kdeLogFilename += oss.str();
			
		}
		string kdeResultsFilename = prefix;
		{
			ostringstream oss;
			oss << "_trace" << spec << "_" << id <<".txt";
			kdeResultsFilename += oss.str();
			
		}
		string gridFilename = prefix;
		{
			ostringstream oss;
			oss << "_grid" << spec << "_" << id <<".txt";
			gridFilename += oss.str();
			
		}
		string densitiesFilename = prefix;
		{
			ostringstream oss;
			oss << "_densities" << spec << "_" << id <<".txt";
			densitiesFilename += oss.str();
			
		}
		
		// for the mcmc kde
		size_t warm=5000;
		size_t step=10; 
		size_t num_batch = 50;
		size_t size_batch = 500; // 50 * 500 = 25,000 iterations
		unsigned int kdeSeed = 1234;
		
		// for the function estimation pq
		size_t maxLeaves = 1000000; 
		
			
		
		// log values used for kde approx to file
		{
			std::ostringstream stm;
			stm << "\n\nMCMC KDE using";
			stm << "\n\twarm\t" << warm;
			stm << "\n\tstep\t" << step;
			stm << "\tnum_batch\t" << num_batch;
			stm << "\tsize_batch\t" << size_batch;
			stm << "\tseed\t" << kdeSeed;
			stm << "\n\nand function estimation using";
			stm << "\n\tmaxMeasure\t" << maxMeasure;
			stm << "\n\tmaxLeaves\t" << maxLeaves;
			outputAHABCLogFileString(logFilename, stm.str());
		}
		
		clock_t start = clock();	
	
		MCMCKDE mkde(container);
		
		mkde.df(size_batch, num_batch, warm, step,
			kdeLogFilename,
			kdeResultsFilename,
			kdeSeed);
		
		clock_t end = clock();	
		
		double timing = (static_cast<double>(end-start)/CLOCKS_PER_SEC);	
		cout << "time to do KDE = " << timing << endl;
		
		
		int d = VecLen(box);
		if (d <= 2) {
			
			/*get a grid so that we can ouput stuff to make a visual of 
			* the kde itself */
			size_t gps = 100; // points to have along each dim of the box, for gridding

			/* grid the box */
			vector < vector < real > > pts;
		
			doGrid(gps, box, pts);
		
			outputGrid(gridFilename, pts);
		
		
			cout << "\nGet kde densities" << endl;
			
			if (d == 2 ) {
				 std::vector < std::vector < real > > kdeGridDensities;
			
							
				getKDEGridDensitiesTwoD(mkde,
								pts,
								kdeGridDensities);
				
				outputDensities(densitiesFilename, kdeGridDensities);
			}
			
			if (d == 1) {
				
				std::vector < real > kdeGridDensities;
				
				getKDEGridDensitiesOneD(mkde,
						pts,
						kdeGridDensities);
				
				outputDensities(densitiesFilename, kdeGridDensities);
			}
			
		}
		
	}
	catch (std::exception& e) {
	
		std::cout << "std::exception::\n" + std::string(e.what()) << std::endl;
		throw;
	}
	catch (...) {
	
		std::cout << "Unknown error\n" << std::endl;
		throw;
	}
}



vector < vector < real > >& doGrid(size_t gps,
			const ivector& box,
			vector < vector < real > >& pts)
{
	/* grid the box */
	int lb = Lb(box);
	int dim = Ub(box) - lb + 1;
	vector < vector < real > > tmp(dim, vector < real >(gps, 0.0));
	
	for (int i = 0; i < dim; ++i) {
		real seg = diam(box[lb+i])/(gps-1.0);
		tmp[i][0] = Inf(box[lb+i]);
		tmp[i][gps-1] = Sup(box[lb+i]);
		
		for (size_t j = 1; j < gps-1; ++j) tmp[i][j] = tmp[i][j-1] + seg;			
	}
	tmp.swap(pts);
	return pts;

}

void outputGrid(const std::string& filename,
			const vector < vector < real > >& pts)
{
	subpavings::outputFileStart(filename);
	
	if (!pts.empty()) {
		ofstream os(filename.c_str(), ios::app);         // append
		
		if (os.is_open()) {
			size_t m = pts.size(); // cols
			size_t n = pts.front().size(); // rows
			for (size_t i = 0; i < n; ++i) {
				for (size_t j = 0; j < m-1; ++j) {
				
					os << pts[j][i] << "\t";
				}
				os << pts[m-1][i] << std::endl;
			}
					
			os.close();
		}
		else {
			std::cerr << "Error: could not open file named "
				<< filename << std::endl << std::endl;
		}
	}

}

// output densities for 2-d data (grid)
void outputDensities(const std::string& filename,
					const std::vector < std::vector < real > >& densities)
					
{
	subpavings::outputFileStart(filename);
	
	ofstream os(filename.c_str(), ios::app);         // append
	if (os.is_open()) {
		size_t m = densities.size(); // rows
		for (size_t i = 0; i < m; ++i) {
			ostream_iterator<real> out_it (os,"\t");
			copy ( densities[i].begin(), densities[i].end(), out_it );
			os << endl;
		}
				
		os.close();
	}
	else {
		std::cerr << "Error: could not open file named "
			<< filename << std::endl << std::endl;
	}
}

// output densities for 1-d data (one row)
void outputDensities(const std::string& filename,
					const std::vector < real >& densities)
					
{
	subpavings::outputFileStart(filename);
	
	ofstream os(filename.c_str(), ios::app);         // append
	if (os.is_open()) {
		ostream_iterator<real> out_it (os,"\n");
			copy ( densities.begin(), densities.end(), out_it );
		
		os.flush();		
		os.close();
	}
	else {
		std::cerr << "Error: could not open file named "
			<< filename << std::endl << std::endl;
	}
}

void getKDEGridDensitiesTwoD(const TypeKDE& tkde,
						const vector < vector < real > >& pts,
						std::vector < std::vector < real > >& kdeGridDensities)
{
	if (pts.size() != 2) {
		cerr << "getKDEGridDensitiesTwoD(...) : not 2 D" << endl;
		return;
	}
	
	clock_t starttime = clock();	
	
	
	size_t n = pts[0].size(); // columns
	size_t m = pts[1].size(); // rows
	
	std::vector < std::vector < real > > tmp(m, vector < real >(n));
		
	for (size_t i = 0; i < m; ++i) { // rows
	
		for (size_t j = 0; j < n; ++j) { // columns
		
			vector < real > r(1, pts[0][j]); // coord 1 is columns
			r.push_back(pts[1][i]); // coord 2 is rows
					
			tmp[i][j] = tkde.kde(r);
		}
	}
	
	kdeGridDensities.swap(tmp);
		
	clock_t endtime = clock();	
	
	double timing = (static_cast<double>(endtime-starttime)/CLOCKS_PER_SEC);	
	cout << "time to get kde densities = " << timing << endl;

}


void getKDEGridDensitiesOneD(const TypeKDE& tkde,
						const vector < vector < real > >& pts,
						std::vector < real > & kdeGridDensities)
{
	if (pts.size() != 1) {
		cerr << "getKDEGridDensitiesTwoD(...) : not 1 D" << endl;
		return;
	}
	
	clock_t starttime = clock();	
	
	
	size_t n = pts[0].size(); 
	
	std::vector < real > tmp(n);
		
	for (size_t j = 0; j < n; ++j) { 
		
		vector < real > r(1, pts[0][j]);
				
		tmp[j] = tkde.kde(r);
	}
	
	
	kdeGridDensities.swap(tmp);
		
	clock_t endtime = clock();	
	
	double timing = (static_cast<double>(endtime-starttime)/CLOCKS_PER_SEC);	
	cout << "time to get kde densities = " << timing << endl;

}


