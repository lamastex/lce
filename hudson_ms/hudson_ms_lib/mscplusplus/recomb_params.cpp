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



#include "recomb_params.hpp"

#include <cmath>
#include <stdexcept>
#include <sstream>
#include <iomanip>


using namespace hudson_ms;


// constructor
HudsonMSRecombParams::HudsonMSRecombParams()
	: 	rho(0.0), nsites(1), 
		conversion_parameter(0.0),
		conversion_track_len(0.0)
{}

HudsonMSRecombParams::HudsonMSRecombParams(	
									double _rho, 
									int _nsites,
									double _conversion_parameter,
									double _conversion_track_len)
	: 	rho(_rho), nsites(_nsites), 
		conversion_parameter(_conversion_parameter),
		conversion_track_len(_conversion_track_len)
{
	
	try {

		checkPars();
		
	}
	catch (std::invalid_argument& ia) {
		throw std::invalid_argument(
		std::string("HudsonMSRecombParams::HudsonMSRecombParams(...) : ")
				+ ia.what() );
	}
}

HudsonMSRecombParams::HudsonMSRecombParams(	
									double _rho, 
									int _nsites)
	: 	rho(_rho), nsites(_nsites), 
		conversion_parameter(0.0),
		conversion_track_len(0.0)
{
	
	try {

		checkPars();
		
	}
	catch (std::invalid_argument& ia) {
		throw std::invalid_argument(
		std::string("HudsonMSRecombParams::HudsonMSRecombParams(...) : ")
				+ ia.what() );
	}
}
HudsonMSRecombParams::~HudsonMSRecombParams() {}

double HudsonMSRecombParams::getRho() const
{
	return rho;
}
			
int HudsonMSRecombParams::getNsites() const
{
	return nsites;
}

double HudsonMSRecombParams::getConversionParameter() const
{
	return conversion_parameter;
}

double HudsonMSRecombParams::getConversionTrackLen() const
{
	return conversion_track_len;
}

std::string HudsonMSRecombParams::msString() const
{
	std::ostringstream stm;
	stm << std::fixed;
	
	if ( rho > 0.0 || conversion_parameter > 0.0) outputCrossover(stm);
	if ( conversion_parameter > 0.0 ) outputConversion(stm);
	
	return stm.str();
}

std::string HudsonMSRecombParams::toString() const
{
	std::ostringstream stm;
	stm << std::fixed << std::setprecision(precEst(rho));
	
	stm << "rho = " << rho;
	stm << ", nsites = " << nsites;
	stm << std::setprecision(precEst(conversion_parameter));
	stm << ", conversion_parameter = " << conversion_parameter;
	stm << std::setprecision(precEst(conversion_track_len));
	stm << ", conversion_track_len = " << conversion_track_len;
	
	return stm.str();
}

void HudsonMSRecombParams::checkPars() const
{
	if (rho < 0.0) {
		throw std::invalid_argument("_rho < 0.0");
	}
	if ((rho > 0.0 || conversion_parameter > 0.0) && nsites < 2) {
		throw std::invalid_argument("_nsites < 2");
	}
	if (conversion_parameter < 0.0) {
		throw std::invalid_argument("_conversion_parameter < 0.0");
	}
	if (conversion_track_len < 0.0) {
		throw std::invalid_argument("_conversion_track_len < 0.0");
	}
	if (conversion_parameter >  0.0 && conversion_track_len == 0.0) {
		throw std::invalid_argument(
		"_conversion_parameter > 0.0 and _conversion_track_len == 0.0");
	}
}

std::ostringstream& HudsonMSRecombParams::outputCrossover(
								std::ostringstream& stm)	const
{
	stm.precision(precEst(rho));
	stm << " -r " << rho << " " << nsites; 
	return stm;
}

std::ostringstream& HudsonMSRecombParams::outputConversion(
								std::ostringstream& stm)	const
{
	stm.precision(precEst(conversion_parameter));
	stm << " -c " << conversion_parameter;
	stm.precision(precEst(conversion_track_len));
	stm << " " << conversion_track_len;
	return stm;
}

int HudsonMSRecombParams::precEst(double val)	const
{
	int prec = 1;
	int maxprec = 10;
		
	// rough and ready way to find an appropriate precision for printing values
	double adjval = val*pow(10.0,prec);
	while ( (std::abs(adjval - floor(adjval)) > 0.0) && (prec < maxprec) ) {
		prec++;
		adjval *= 10;
	}
	
	return prec; 
}
