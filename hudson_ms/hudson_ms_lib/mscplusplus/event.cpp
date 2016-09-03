/* mscplusplus - ms for c++.

   Copyright (C) 2011, 2012 Jenny Harlow
   
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
\brief HudsonMSDemogEvent definitions.
*/

#include "arg_factory.hpp"
#include "event.hpp"

#include <stdexcept>
#include <algorithm>
#include <iomanip>
#include <cmath>


using namespace hudson_ms;

HudsonMSEvent::HudsonMSEvent(double _time)
	: time(_time)
{
	if (time < 0.0) {
		throw std::invalid_argument(
			"HudsonMSEvent::HudsonMSEvent(double)" );
	}
}
	
double HudsonMSEvent::getTime() const
{
	return time;
}


std::string HudsonMSDemogEvent::msString() const
{
	std::ostringstream stm;
	stm << std::fixed;
	int prec = precEst(getTime());
	stm.precision(prec);
	outputDescription(stm);
	
	return stm.str();
}

bool HudsonMSDemogEvent::operator < (const hudson_ms::HudsonMSDemogEvent& rhs) const
{
	if (getTime() != rhs.getTime()) {
		return ( getTime() < rhs.getTime() );
	}
	else { // same times
		return ( getOrdering() < rhs.getOrdering() );
	}
}


int HudsonMSDemogEvent::getOrdering() const
{
	return default_ordering;
}

// rough and ready way to find an appropriate precision for printing values for ms
int HudsonMSDemogEvent::precEst(double param)
{
	int prec = 1;
	int maxprec = 10;
	double adjval = param*pow(10.0,prec);
	while ( (std::abs(adjval - floor(adjval)) > 0.0) && (prec < maxprec) ) {
		prec++;
		adjval *= 10;
	}
	return prec;
}

const int HudsonMSDemogEvent::default_ordering;



HudsonMSDemogEvent_N::HudsonMSDemogEvent_N(double _time, double _multNzero)
	: HudsonMSDemogEvent(_time), multNzero(_multNzero) 
{
	if (multNzero <= 0.0) {
		throw std::invalid_argument(
			"HudsonMSDemogEvent_N::HudsonMSDemogEvent_N(double, double)" );
	}
			
}

void HudsonMSDemogEvent_N::callEvent(const hudson_ms::HudsonMSArgFactory& builder) const
{
	builder.processEvent_N(multNzero);
}



int HudsonMSDemogEvent_N::getOrdering() const
{
	return ordering;
}

const int HudsonMSDemogEvent_N::ordering;

//-eN t x
std::ostringstream& HudsonMSDemogEvent_N::outputDescription(std::ostringstream& stm ) const
{
	stm << "-eN " << getTime();
	
	stm.precision(precEst(multNzero));
	
	stm << " " << multNzero; 
	
	return stm;
}



HudsonMSDemogEvent_n::HudsonMSDemogEvent_n(double _time, int _pop, double _multNzero)
	: HudsonMSDemogEvent(_time), multNzero(_multNzero)
{
	if (multNzero <= 0.0 || _pop < 0) {
		throw std::invalid_argument(
			"HudsonMSDemogEvent_n::HudsonMSDemogEvent_n(double, int, double)" );
	}
	pop = static_cast<size_t>(_pop);
}

void HudsonMSDemogEvent_n::callEvent(const hudson_ms::HudsonMSArgFactory& builder) const
{
	builder.processEvent_n(pop, multNzero);
}

int HudsonMSDemogEvent_n::getOrdering() const
{
	return ordering;
}

const int HudsonMSDemogEvent_n::ordering;

//-en t i x
// allow for pops being given as 1, 2, etc <- 0, 1, etc
std::ostringstream& HudsonMSDemogEvent_n::outputDescription(std::ostringstream& stm ) const
{
	stm << "-en " << getTime() << " " << (pop+1);
	
	stm.precision(precEst(multNzero));
	
	stm << " " << multNzero; 
	
	return stm;
}



HudsonMSDemogEvent_G::HudsonMSDemogEvent_G(double _time, double _newGrowthRate)
	: HudsonMSDemogEvent(_time), newGrowthRate(_newGrowthRate) {}
	
void HudsonMSDemogEvent_G::callEvent(const hudson_ms::HudsonMSArgFactory& builder) const
{
	builder.processEvent_G(getTime(), newGrowthRate);
}

//-eG t a
std::ostringstream& HudsonMSDemogEvent_G::outputDescription(std::ostringstream& stm ) const
{
	stm << "-eG " << getTime();
	
	stm.precision(precEst(newGrowthRate));
	
	stm << " " << newGrowthRate; 
	
	return stm;
}



HudsonMSDemogEvent_g::HudsonMSDemogEvent_g(double _time, int _pop, double _newGrowthRate)
	: HudsonMSDemogEvent(_time), newGrowthRate(_newGrowthRate)
{
	if (_pop < 0) {
		throw std::invalid_argument(
			"HudsonMSDemogEvent_g::HudsonMSDemogEvent_g(double, int, double)" );
	}
	pop = static_cast<size_t>(_pop);
}

void HudsonMSDemogEvent_g::callEvent(const hudson_ms::HudsonMSArgFactory& builder) const
{
	builder.processEvent_g(getTime(), pop, newGrowthRate);
}

//-eg t i a
// allow for pops being given as 1, 2, etc <- 0, 1, etc
std::ostringstream& HudsonMSDemogEvent_g::outputDescription(std::ostringstream& stm ) const
{
	stm << "-eg " << getTime() << " " << (pop+1);
	
	stm.precision(precEst(newGrowthRate));
	
	stm << " " << newGrowthRate; 
	
	return stm;
}


const int HudsonMSDemogEvent_j::ordering;

HudsonMSDemogEvent_j::HudsonMSDemogEvent_j(double _time, int _pop_i, int _pop_j)
	: HudsonMSDemogEvent(_time)
{
	if (_pop_i < 0 || _pop_j < 0) {
		throw std::invalid_argument(
			"HudsonMSDemogEvent_j::HudsonMSDemogEvent_j(double, int, int)" );
	}
	pop_i = static_cast<size_t>(_pop_i);
	pop_j = static_cast<size_t>(_pop_j);
}

void HudsonMSDemogEvent_j::callEvent(const hudson_ms::HudsonMSArgFactory& builder) const
{
	builder.processEvent_j(pop_i, pop_j);
}

int HudsonMSDemogEvent_j::getOrdering() const
{
	return ordering;
}

//-ej t i j
// allow for pops being given as 1, 2, etc -> 0, 1, etc
std::ostringstream& HudsonMSDemogEvent_j::outputDescription(std::ostringstream& stm ) const
{
	stm << "-ej " << getTime() << " " << (pop_i+1) << " " << (pop_j+1);
	
	return stm;
}


const int HudsonMSDemogEvent_s::ordering;

HudsonMSDemogEvent_s::HudsonMSDemogEvent_s(double _time, int _pop, double _probStaying)
	: HudsonMSDemogEvent(_time), probStaying(_probStaying)
{
	if (probStaying <= 0.0 || _pop < 0) {
		throw std::invalid_argument(
			"HudsonMSDemogEvent_s::HudsonMSDemogEvent_s(double, int, double)" );
	}
	pop = static_cast<size_t>(_pop);
}

void HudsonMSDemogEvent_s::callEvent(const hudson_ms::HudsonMSArgFactory& builder) const
{
	builder.processEvent_s(getTime(), pop, probStaying);
}

int HudsonMSDemogEvent_s::getOrdering() const
{
	return ordering;
}

//-es t i p
// allow for pops being given as 1, 2, etc <- 0, 1, etcstd::ostringstream& HudsonMSDemogEvent_s::outputDescription(std::ostringstream& stm ) const
std::ostringstream& HudsonMSDemogEvent_s::outputDescription(std::ostringstream& stm ) const
{
	stm << "-es " << getTime() << " " << (pop+1);
	
	stm.precision(precEst(probStaying));
	
	stm << " " << probStaying; 
	
	return stm;
}



HudsonMSDemogEvent_M::HudsonMSDemogEvent_M(double _time, double _newMigRate)
	: HudsonMSDemogEvent(_time), newMigRate(_newMigRate)
{
	if (newMigRate < 0.0) {
		throw std::invalid_argument(
			"HudsonMSDemogEvent_M::HudsonMSDemogEvent_M(double, double)" );
	}
}

void HudsonMSDemogEvent_M::callEvent(const hudson_ms::HudsonMSArgFactory& builder) const
{
	builder.processEvent_M(newMigRate);
}

//-eM t x
std::ostringstream& HudsonMSDemogEvent_M::outputDescription(std::ostringstream& stm ) const
{
	stm << "-eM " << getTime();
	
	stm.precision(precEst(newMigRate));
	
	stm << " " << newMigRate; 
	
	return stm;
}



HudsonMSDemogEvent_m::HudsonMSDemogEvent_m(double _time, int _pop_i, int _pop_j, double _newMigRate)
	: HudsonMSDemogEvent(_time), newMigRate(_newMigRate)
{
	if (newMigRate < 0.0 || _pop_i < 0 || _pop_j < 0) {
		throw std::invalid_argument(
			"HudsonMSDemogEvent_m::HudsonMSDemogEvent_m(double, int, int, double)" );
	}
	pop_i = static_cast<size_t>(_pop_i);
	pop_j = static_cast<size_t>(_pop_j);
}

void HudsonMSDemogEvent_m::callEvent(const hudson_ms::HudsonMSArgFactory& builder) const
{
	builder.processEvent_m(pop_i, pop_j, newMigRate);
}

//-em t i j x
// allow for pops being given as 1, 2, etc <- 0, 1, etc
std::ostringstream& HudsonMSDemogEvent_m::outputDescription(std::ostringstream& stm ) const
{
	stm << "-em " << getTime() << " " << (pop_i+1) << " " << (pop_j+1);
	
	stm.precision(precEst(newMigRate));
	
	stm << " " << newMigRate; 
	
	return stm;
}



HudsonMSDemogEvent_ma::HudsonMSDemogEvent_ma(double _time, const std::vector< std::vector < double > >& _newMigMatrix)
	: HudsonMSDemogEvent(_time), newMigMatrix(_newMigMatrix)
{
	for (std::vector < std::vector <double> >::const_iterator it = _newMigMatrix.begin();
		it < _newMigMatrix.end();
		++it) {
			
		if (*std::min_element(it->begin(), it->end()) < 0.0 ) {
			throw std::invalid_argument(
				"HudsonMSDemogEvent_ma::HudsonMSDemogEvent_ma(const std::vector < std::vector <double> >&)");
		}
	}
}

void HudsonMSDemogEvent_ma::callEvent(const hudson_ms::HudsonMSArgFactory& builder) const
{
	builder.processEvent_ma(newMigMatrix);
}

//-ema t npop M11 M12 .. M21, M22 ..
std::ostringstream& HudsonMSDemogEvent_ma::outputDescription(std::ostringstream& stm ) const
{
	stm << "-ema " << getTime() << " " << newMigMatrix.size();
	stm << " " << migrationMatrixString();
	return stm;
}

// it is annoying to duplicate this from pop structure - maybe later I can sort this out ...
std::string HudsonMSDemogEvent_ma::migrationMatrixString() const
{
	std::ostringstream stm;
		
	if (!newMigMatrix.empty()) {
		int prec = 1;
		int maxprec = 10;
		// rough and ready way to find an appropriate precision for printing mig value
		for (size_t i = 0; i < newMigMatrix.size(); ++i) {
			for (size_t j = 0; j < (newMigMatrix[i]).size(); ++j) {
			
				double adjval = (newMigMatrix[i][j])*pow(10.0,prec);
				
				while ((adjval - floor(adjval) > (1.0/maxprec)) && (prec < maxprec)) {
					prec++;
					adjval *= 10;
				}
			}
		}
		stm.precision(prec);
		stm << std::fixed;
		
		for (size_t i = 0; i < newMigMatrix.size(); ++i) {
			for (size_t j = 0; j < (newMigMatrix[i]).size(); ++j) {
			
				if (i != j) stm << " " << (newMigMatrix[i])[j];
				else stm << " x"; // x for diagonals
			}
		}
	}
	return stm.str();
}



HudsonMSDemogEventPtr::HudsonMSDemogEventPtr(
	boost::shared_ptr < const hudson_ms::HudsonMSDemogEvent > e)
	: eventPtr(e)
{
	if (!eventPtr) { // empty
		throw std::invalid_argument(
			"HudsonMSDemogEventPtr::HudsonMSDemogEventPtr(boost::shared_ptr < const HudsonMSDemogEvent >)");
	}
}

void HudsonMSDemogEventPtr::callEvent(const hudson_ms::HudsonMSArgFactory& builder) const
{
	eventPtr->callEvent(builder);
}

double HudsonMSDemogEventPtr::getTime() const 
{
	return eventPtr->getTime();
}

std::string HudsonMSDemogEventPtr::msString() const
{
	return eventPtr->msString();
}

bool HudsonMSDemogEventPtr::operator < (const hudson_ms::HudsonMSDemogEventPtr& rhs) const
{
	return  (*eventPtr) < (*(rhs.eventPtr) );
}
