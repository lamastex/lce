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
\brief ParameterParticle definitions.
*/

#include "parameter_particle.hpp"

#include <algorithm> // for swap, transform
#include <functional> // for multiplies
#include <sstream>  // to be able to manipulate strings as streams
#include <iterator> // for ostream_iterator
#include <stdexcept>


using namespace aabc;

ParameterParticle::ParameterParticle(const std::vector < double >& values)
: n(values.size()), parameters(values), weight(0.0)
{}

ParameterParticle::ParameterParticle(const double a[], const size_t _n)
: n(_n), weight(0.0)
{
	if (NULL != a) {
		parameters.insert (parameters.begin(), a, a + n);
	}
	else {
		throw std::invalid_argument(
			"ParameterParticle::ParameterParticle(const double [], const size_t)");
	}
		
}

ParameterParticle::ParameterParticle(const size_t _n) 
		: n(_n), parameters(std::vector <double>(_n) ), weight(0.0) 
		{}

			
ParameterParticle::ParameterParticle(const gsl_vector * const v)
	: weight(0.0)
{
	if (NULL != v) {
		for (size_t i = 0; i < v->size; ++i) {
			parameters.push_back( gsl_vector_get(v, i) );
		}
		n = parameters.size();
	}
	else {
		throw std::invalid_argument(
			"ParameterParticle::ParameterParticle(const gsl_vector * const)");
	}
		
}

ParameterParticle::ParameterParticle(const ParameterParticle& other)
		: n(other.n), parameters(other.parameters), weight(other.weight) 
		{}


ParameterParticle& ParameterParticle::operator=(ParameterParticle tmp)
{
	std::swap( n, tmp.n );
	std::swap( weight, tmp.weight );
    std::swap( parameters, tmp.parameters );

	return *this;
}

ParameterParticle::~ParameterParticle()
{
		//std::cout << "ParameterParticle destructor called for " << toString() << std::endl; 
}


gsl_vector * ParameterParticle::get_gsl_vector() const
{
	gsl_vector * v = NULL;
	try {
		v = gsl_vector_calloc(n);
	
		// copy the current parameters into the new mean vector
		for (size_t i = 0; i < n; ++i) {
				
			gsl_vector_set(v, i, parameters.at(i));	
				
		}
		return v;
	}
	catch (...) {
		try { if (v != NULL) gsl_vector_free(v); } 
		catch (...) {}
		throw;
	}
}

std::vector < double > ParameterParticle::get_vector() const
{
	return parameters;	
}
		
double& ParameterParticle::at(const size_t index)
{
	return parameters.at(index);
}
			
const double& ParameterParticle::at(const size_t index) const
{
	return parameters.at(index);
}
			
double ParameterParticle::getWeight() const
{
	return weight;
}
			
void ParameterParticle::setWeight(const double w) const
{
	weight = w;
}
			
size_t ParameterParticle::size() const
{
	return n;
}	

ParameterParticle::const_iterator ParameterParticle::begin() const
{
	return parameters.begin();
}

ParameterParticle::const_iterator ParameterParticle::end() const
{
	return parameters.end();
}

ParameterParticle ParameterParticle::operator-(
					const ParameterParticle& b) const
{
	gsl_vector* a_vec = NULL;
	gsl_vector* b_vec = NULL;
	try {
		size_t n = size();
		if (n != b.size() ) {
			throw std::invalid_argument(
			"ParameterParticle ParameterParticle::operator-(const ParameterParticle&)");
		}
		a_vec = get_gsl_vector();
		b_vec = b.get_gsl_vector();
		
		gsl_vector_sub(a_vec, b_vec);
		
		ParameterParticle p(a_vec);
		
		gsl_vector_free(a_vec);
		gsl_vector_free(b_vec);
		
		return p;
	}
	catch (...) {
		try { if ( a_vec != NULL) gsl_vector_free(a_vec);}
		catch (...) {} // catch and swallow
		try {if ( b_vec != NULL) gsl_vector_free(b_vec);}
		catch (...) {} // catch and swallow
		throw;
	}
	
}


std::string ParameterParticle::toString() const
{
	std::ostringstream stm;
	
	stm << "(" << plainString() << ") weight: " << weight;
	return stm.str();
	
}

std::string ParameterParticle::plainString() const
{
	std::ostringstream stm;
	stm.precision(5);
	stm <<   std::fixed;
		
	std::ostream_iterator<double> out_it (stm, "\t");
	copy ( parameters.begin(), parameters.end(), out_it );
  
    std::string s = stm.str();
	size_t pos = s.find_last_of('\t');
	s.erase(pos, 1);


    return s;
}

std::string ParameterParticle::plainString(const std::vector < double >& mults) const
{
	if (n != mults.size() ) {
		throw std::invalid_argument( 
			"ParameterParticle::plainString(const std::vector < double >&)");
	}

	std::vector < double > temp_P = parameters;
	
	for (size_t i = 0; i < n; ++i) {
		temp_P[i] *= mults[i]; 
	}
	
	std::ostringstream stm;
	stm.precision(5);
	stm <<   std::fixed;
		
	std::ostream_iterator<double> out_it (stm, "\t");
	copy ( temp_P.begin(), temp_P.end(), out_it );
  
	std::string s = stm.str();
	size_t pos = s.find_last_of('\t');
	s.erase(pos, 1);

    return s;
}

// protected


// non member functions
std::ostream& aabc::operator<< (std::ostream &out, const ParameterParticle &p)
{
	out << p.toString();
	return out;
}


