/* MCT - Markov Chains on Trees.

   Copyright (C) 2011 Jenny Harlow

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
\brief ParameterParticleSet definitions.
*/

#include "parameter_particle_set.hpp"
#include "summary_statistic.hpp"

#include <boost/bind.hpp>

#include <boost/lambda/bind.hpp>
#include <boost/lambda/lambda.hpp>


#include <limits>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <iostream>
#include <stdexcept>

//#define MYDEBUG

#ifdef MYDEBUG
	#include <iostream>
#endif

using namespace aabc;

ParameterParticleSet::ParameterParticleSet() {}

ParameterParticleSet::ParameterParticleSet(size_t reserveSize)
{
		container.reserve(reserveSize);
}
			
ParameterParticleSet::ParameterParticleSet(
	const std::vector< boost::shared_ptr< const ParameterParticle > >& vec)
	: container(vec) {}
			

ParameterParticleSet::ParameterParticleSet(
				const ParameterParticleSet& other)
		: container(other.container) {}

ParameterParticleSet& ParameterParticleSet::operator=(ParameterParticleSet tmp)
{
	std::swap( container, tmp.container );

	return *this;
}
			
ParameterParticleSet::~ParameterParticleSet() 
{
	
	
}
			

const boost::shared_ptr< const ParameterParticle >& ParameterParticleSet::at(const size_t index) const
{
	return container.at(index);
}

boost::shared_ptr< const ParameterParticle >& ParameterParticleSet::at(const size_t index)
{
	return container.at(index);
}

void ParameterParticleSet::add(
					const boost::shared_ptr< const ParameterParticle >& toAdd)
{
	container.push_back(toAdd);
}

void ParameterParticleSet::add(
			const boost::shared_ptr < ParameterParticleSet >& toAdd)
{
	std::vector< boost::shared_ptr < const ParameterParticle > > otherContainer
		= toAdd->_getContainer();
	container.insert( container.end(),
					otherContainer.begin(), otherContainer.end() );
	
}

			
size_t ParameterParticleSet::size() const
{
	return container.size();
}

bool ParameterParticleSet::empty() const
{
	return container.empty();
}

bool ParameterParticleSet::sizeConsistencyCheck() const 
{
	bool result = true;
	if (!empty()) {
		ParameterParticleSet::const_iterator it = begin();
		size_t n = (*it)->size();
		++it;				
		while ( result && it < end()) {
			
			if ((*it)->size() != n) {
				result = false;
			}
			++it;
		}
	}
	return result;
}

#if(0)
boost::shared_ptr < ParameterParticleSet >
		ParameterParticleSet::makeParameterParticleSetStandardised() const
{
	//get all the standardised values as vec of vecs
	boost::shared_ptr < std::vector < std::vector < double > > > 
				stanvals = getAllValuesStandardised();
	
	size_t nss = size();
		
	// empty set with enough reserved space
	boost::shared_ptr < ParameterParticleSet > pset( new ParameterParticleSet(nss) );
	
	// construct and add standardised particles to the set
	if (!empty()) {
		
		size_t nstats = stanvals->size();
		
		for (size_t i = 0; i < nss; ++i) {
			
			std::vector < double > tmp;
			tmp.reserve(nstats);
			
			for (size_t j = 0; j < nstats; ++j) {
				tmp.push_back( (stanvals->at(j)).at(i) );
			}
			
			boost::shared_ptr < ParameterParticle > sptr( new ParameterParticle(tmp) );
			pset->add(sptr);
		}
	}
	
	return pset;
}
#endif

boost::shared_ptr < mct::SummaryStatisticSet >
				ParameterParticleSet::makeSummaryStatisticSet() const
{
	// empty set with enough reserved space
	boost::shared_ptr < mct::SummaryStatisticSet > 
						ss_set( new mct::SummaryStatisticSet(size()) );
	
	// construct and add summary statistics to the set
	ParameterParticleSet::const_iterator it = begin();
	while (it < end()) {
			
		boost::shared_ptr < mct::SummaryStatistic > 
			sptr( new mct::SummaryStatistic( (*it)->get_vector() ) );
		
		ss_set->add(sptr);
		
		++it;
	
	}
	
	return ss_set;
	
}

boost::shared_ptr < std::vector < std::vector < double > > >
				ParameterParticleSet::getAllValues() const
{
	boost::shared_ptr < std::vector < std::vector < double > > >
		results;
		
	if (!empty() ) {
		size_t n = (*begin())->size();
		
		std::vector < double > tmp;
		
		results = boost::shared_ptr < std::vector < std::vector < double > > >
					( new std::vector < std::vector < double > > (n, tmp));
		size_t set_size = size();
		for_each(results->begin(), results->end(), 
				boost::bind(&std::vector<double>::reserve, _1, set_size));
		
		#ifdef MYDEBUG
			std::cout << "n = " << n << std::endl;
				std::cout << "size of results is  " << results->size() << std::endl;
				std::cout << "results capacity is  " << results->capacity() << std::endl;
			for (std::vector < std::vector < double > >::iterator it = results->begin();
				it < results->end(); ++it) {
				std::cout << "Results capacity at " << distance(results->begin(), it) << " is " << (*it).capacity() << std::endl;
			}
		#endif
		
		for (ParameterParticleSet::const_iterator it = begin();
				it< end();
				++it) {
					
			std::vector < double > v = (*it)->get_vector();	
			size_t this_n = v.size();
			
			#ifdef MYDEBUG
				std::cout << "element " << distance(begin(), it) 
				<< "this_n = " << this_n << " and n = " << n <<std::endl;
			#endif	
			
			size_t extra = this_n > n ? this_n - n : 0;
			if (extra) n = this_n;
			
			// add in more containers if this_n > n, 
			// and stuff with nans for particles done so far
			while (extra > 0) {
				results->push_back( std::vector < double > 
						( distance( begin(), it ),
						std::numeric_limits<double>::quiet_NaN() ) );
				(results->back()).reserve( set_size );	
				extra--;
				#ifdef MYDEBUG
					std::cout << "n = " << n << " and this_n = " << this_n << std::endl;
					std::cout << "size of results now is  " << results->size() << std::endl;
				#endif
				
			}
					
			for (size_t i = 0; i < this_n; ++i) {
				(results->at(i)).push_back( (*it)->at(i) );
			}
			//stuff the rest with nans
			for (size_t i = this_n; i < n; ++i) {
				#ifdef MYDEBUG
					std::cout << "i = " << i << " is being stuffed" << std::endl;
				#endif
					
				(results->at(i)).push_back(std::numeric_limits<double>::quiet_NaN());
			}
					
		}
	}
	return results;
}


gsl_matrix* ParameterParticleSet::covariance() const
{
	// each inner vector is a vector of parameter values (one row)
	// there are N of these, for N observations (N rows, the ParameterParticles in the container)
	
	gsl_matrix* r = NULL;
	try {
		
		if (container.empty()) {
			throw std::runtime_error(
			std::string("ParameterParticleSet::covariance())")
			+ std::string(": No particles") );

		}
		
		size_t cols = (container.at(0))->size();
		size_t rows = container.size();
		
		r = gsl_matrix_alloc (cols, cols);
	
		std::vector < std::vector < double > > work(cols);
		for (size_t i = 0; i < cols; ++i) {
				
				std::vector< double > a(rows);
				
				for (size_t k = 0; k < rows; ++k) {
					a.at(k) = (container.at(k))->at(i);
				}
				// vectors of the mean, each element the same
				std::vector< double > a_bars( a.size(),
									std::accumulate(a.begin(),
									a.end(), 0.0)/rows ); 
				// subtract means from values
				std::transform(a.begin(), a.end(), a_bars.begin(),
							a.begin(), std::minus< double >());
				
				work[i] = a;
		}
	
		for (size_t i = 0; i < cols; ++i) {
			for (size_t j = 0; j < cols; ++j) {
				
				// summation of product over all k by innerproduct
				// div by N, not N-1, ie MLE version
				double v = inner_product((work[i]).begin(), (work[i]).end(),
				(work[j]).begin(), 0.0)/(rows); 
				
				gsl_matrix_set (r, i, j, v);;
				
			}
		}
		
		
		return r;
	}
	catch (...) {
		
		try {
			if ( r != NULL) gsl_matrix_free(r);
		}
		catch (...) {}
		throw;
	}
	
}

gsl_matrix* ParameterParticleSet::weightedCovariance() const
{
	// each inner vector is a vector of parameter values (one row)
	// there are N of these, for N observations (N rows, the ParameterParticles in the container)
	gsl_matrix* r = NULL;
	
	try {
		
		if (container.empty()) {
			throw std::runtime_error(
			std::string("ParameterParticleSet::weightedCovariance()")
			+ std::string(": No particles") );
		}

		size_t cols = (container.at(0))->size(); // n
		size_t rows = container.size(); // N
		
		r = gsl_matrix_alloc (cols, cols);

		size_t i = 0, j = 0;
		
		double totalWeight = accumulateWeights();
		
		if (0.0 == totalWeight) {
			throw std::runtime_error(
			std::string("ParameterParticleSet::weightedCovariance()")
			+ std::string(": Total weight == 0.0") );
		}
				
		for (i = 0; i < cols; ++i) {
			for (j = 0; j < cols; ++j) {
				
				std::vector< double > a(rows);
				std::vector< double > b(rows);
				std::vector< double > weights(rows);
				
				for (size_t k = 0; k < rows; ++k) {
					a.at(k) = (container.at(k))->at(i);
					b.at(k) = (container.at(k))->at(j);
					weights.at(k) = (container.at(k))->getWeight();
				}
				
				//weighted means x_bar_i and x_bar_j
				double x_bar_i = inner_product(a.begin(), a.end(),
				weights.begin(), 0.0)/(totalWeight); // normalise weights
				double x_bar_j = inner_product(b.begin(), b.end(),
				weights.begin(), 0.0)/(totalWeight); // normalise weights
				
				// make vector of means, all elements the same
				std::vector< double > a_bars( a.size(),
									x_bar_i );
				std::vector< double > b_bars( b.size(),
									x_bar_j );
				
				// subtract in place: a's - means x_bar_i		
				std::transform(a.begin(), a.end(), a_bars.begin(),
							a.begin(), std::minus< double >());
				// subtract in place: b's - means x_bar_j		
				std::transform(b.begin(), b.end(), b_bars.begin(),
							b.begin(), std::minus< double >());
							
				// mult in place (a's-x_bar_i) by weights
				std::transform(a.begin(), a.end(), weights.begin(),
							a.begin(), std::multiplies< double >());
				
				// summation of product over all k using inner product
				// if weights all 1/N this is equivlaent to 
				// div by N, not N-1, ie MLE version
				double v = inner_product(a.begin(), a.end(),
				b.begin(), 0.0)/(totalWeight); // normalise weights
				
				gsl_matrix_set (r, i, j, v);
				
			}
		}
		return r;
	}
	catch (...) {
		try {
			if ( r != NULL) gsl_matrix_free(r);
		}
		catch (...) {}
		throw;
	}
}

gsl_vector* ParameterParticleSet::weightedMean() const
{
	gsl_vector* m = NULL;
	
	try {
		// each inner vector is a vector of parameter values (one row)
		// there are N of these, for N observations (N rows, the ParameterParticles in the container)
		
		if (container.empty()) {
			throw std::runtime_error(
			std::string("ParameterParticleSet::weightedMean()")
			+ std::string(": No particles") );

		}

		size_t cols = (container.at(0))->size(); // n
		size_t rows = container.size(); // N
		
		m = gsl_vector_alloc (cols);

		double totalWeight = accumulateWeights();
		
		if (0.0 == totalWeight) {
			throw std::runtime_error(
			std::string("ParameterParticleSet::weightedMean()")
			+ std::string(": Total weight == 0.0") );
		}
				
		for (size_t i = 0; i < cols; ++i) {
							
			std::vector< double > a(rows);
			std::vector< double > weights(rows);
			
			for (size_t k = 0; k < rows; ++k) {
				a.at(k) = (container.at(k))->at(i);
				weights.at(k) = (container.at(k))->getWeight();
			}
			
			//weighted means x_bar_i
			double x_bar_i = inner_product(a.begin(), a.end(),
			weights.begin(), 0.0)/(totalWeight); // normalise weights
			
			gsl_vector_set (m, i, x_bar_i);
		}
		return m;
	}
	catch (...) {
		try {
			if (m != NULL) gsl_vector_free(m);
		}
		catch (...) {}
		throw;
	}
}

gsl_vector* ParameterParticleSet::mean() const
{
	gsl_vector* m = NULL;
	try {
		
		// each inner vector is a vector of parameter values (one row)
		// there are N of these, for N observations (N rows, the ParameterParticles in the container)
		if (container.empty()) {
			throw std::runtime_error(
			std::string("ParameterParticleSet::mean()")
			+ std::string(": No particles") );
		}

		size_t cols = (container.at(0))->size(); // n
		size_t rows = container.size(); // N
		
		m = gsl_vector_alloc (cols);

		for (size_t i = 0; i < cols; ++i) {
							
			double sum = 0.0;
			
			for (size_t k = 0; k < rows; ++k) {
				sum+= (container.at(k))->at(i);
			}
			
			//mean
			double x_bar_i = sum/rows; 
			
			gsl_vector_set (m, i, x_bar_i);
		}
		return m;
	}
	catch (...) {
		try {
			if (m != NULL) gsl_vector_free(m);
		}
		catch (...) {}
		throw;
	}
}


ParameterParticleSet::const_iterator ParameterParticleSet::begin() const
{
	return container.begin();
}

ParameterParticleSet::const_iterator ParameterParticleSet::end() const
{
	return container.end();
}

size_t ParameterParticleSet::pickOne(
					boost::function<double ()> rand_num_func) const
{
	double totalWeights = accumulateWeights();
	if (totalWeights > 0 && size() > 0) {
		
		double rand_num = rand_num_func();
		rand_num *= totalWeights; // rescale in case weights don't add to 1
		
		double stackedWeight = 0.0;
		size_t index = 0;
		
		for (std::vector< boost::shared_ptr< 
					const ParameterParticle > >::const_iterator 
					it = container.begin();
			it < container.end();
			++it, ++index) {
				
			stackedWeight += (*it)->getWeight();
			if (rand_num < stackedWeight) {
				break;
				
			}
		}
		
		return index;	
	}
	else {
		throw std::runtime_error(
			std::string("ParameterParticleSet::pickOne(boost::function<double ()>)")
			+ std::string(": No particle with weight > 0") );
	}
}


std::string ParameterParticleSet::toString() const
{
	std::ostringstream stm;
	
	for (std::vector< boost::shared_ptr< const ParameterParticle > >::const_iterator
				it = container.begin();
		it < container.end();
		++it) {
			
		stm << (*it)->toString() << "\n";
	}		
	return stm.str();
}

std::string ParameterParticleSet::plainString() const
{
	std::ostringstream stm;
	
	for (std::vector< boost::shared_ptr< const ParameterParticle > >::const_iterator
				it = container.begin();
		it < container.end();
		++it) {
			
		stm << (*it)->plainString() << "\n";
	}		
	return stm.str();
}

std::string ParameterParticleSet::plainString(
							const std::vector < double >& mults) const
{
	std::ostringstream stm;
	
	for (std::vector< boost::shared_ptr< const ParameterParticle > >::const_iterator
				it = container.begin();
		it < container.end();
		++it) {
			
		stm << (*it)->plainString(mults) << "\n";
	}		
	return stm.str();
	
}


void ParameterParticleSet::outputToFile( const std::string& s,
						bool append) const
{
	std::ofstream os;
        
	if (append) os.open(s.c_str(), std::ios::app);         // append
	else os.open(s.c_str()); // overwrite
        
    if (os.is_open()) {

        os << plainString();
        os.close();
    }
    else {
        std::cerr << "Error: could not open file named "
            << s << std::endl << std::endl;
    }
}

void ParameterParticleSet::outputToFile( 
					const std::vector < double >& mults,
					const std::string& s,
					bool append) const
{
	std::ofstream os;
        
	if (append) os.open(s.c_str(), std::ios::app);         // append
	else os.open(s.c_str()); // overwrite
        
    if (os.is_open()) {

        os << plainString(mults);
        os.close();
    }
    else {
        std::cerr << "Error: could not open file named "
            << s << std::endl << std::endl;
    }
}






//protected 


double ParameterParticleSet::accumulateWeights() const
{
	double result = 0.0;
	
	boost::function< double 
				(boost::shared_ptr < const ParameterParticle > ) >
					f1 = boost::bind(&ParameterParticle::getWeight, _1);
		
	//using namespace::boost::lambda;
	
	boost::lambda::placeholder1_type _p;
	
	std::for_each(container.begin(), container.end(), 
		result += boost::lambda::bind(f1, _p) );
		
	return result;
}

const std::vector< boost::shared_ptr < const ParameterParticle > >&
				ParameterParticleSet::_getContainer() const
{
	return container;
}		

// non member functions
std::ostream& aabc::operator<< (std::ostream& out, 
							const ParameterParticleSet& ps)
{
	out << ps.toString();
	return out;
}




