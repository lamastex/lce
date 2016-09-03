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


#include "segsite_params.hpp"

#include <cmath>
#include <stdexcept>
#include <sstream>
#include <iomanip>

using namespace hudson_ms;


// constructor
HudsonMSSegSiteParams::HudsonMSSegSiteParams()
	: 	theta(0.0), segsitesin(0), mfreq(1)
{}

HudsonMSSegSiteParams::HudsonMSSegSiteParams(double _theta, 
											int _segsitesin,
											int _mfreq)
	: 	theta(_theta), segsitesin(_segsitesin), mfreq(_mfreq)
{
	
	try {

		checkPars();
		
	}
	catch (std::invalid_argument& ia) {
		throw std::invalid_argument(
		std::string("HudsonMSSegSiteParams::HudsonMSSegSiteParams(...) : ")
				+ ia.what() );
	}
}

HudsonMSSegSiteParams::HudsonMSSegSiteParams(double _theta,
											int _mfreq)
	: 	theta(_theta), segsitesin(0), mfreq(_mfreq)
{
	
	try {

		checkPars();
		
	}
	catch (std::invalid_argument& ia) {
		throw std::invalid_argument(
		std::string("HudsonMSSegSiteParams::HudsonMSSegSiteParams(...) : ")
				+ ia.what() );
	}
}

HudsonMSSegSiteParams::HudsonMSSegSiteParams(double _theta)
	: 	theta(_theta), segsitesin(0), mfreq(1)
{
	
	try {

		checkPars();
		
	}
	catch (std::invalid_argument& ia) {
		throw std::invalid_argument(
		std::string("HudsonMSSegSiteParams::HudsonMSSegSiteParams(...) : ")
				+ ia.what() );
	}
}

HudsonMSSegSiteParams::HudsonMSSegSiteParams(int _segsitesin,
											int _mfreq)
	: 	theta(0.0), segsitesin(_segsitesin), mfreq(_mfreq)
{
	
	try {

		checkPars();
		
	}
	catch (std::invalid_argument& ia) {
		throw std::invalid_argument(
		std::string("HudsonMSSegSiteParams::HudsonMSSegSiteParams(...) : ")
				+ ia.what() );
	}
}

HudsonMSSegSiteParams::HudsonMSSegSiteParams(int _segsitesin)
	: 	theta(0.0), segsitesin(_segsitesin), mfreq(1)
{
	
	try {

		checkPars();
		
	}
	catch (std::invalid_argument& ia) {
		throw std::invalid_argument(
		std::string("HudsonMSSegSiteParams::HudsonMSSegSiteParams(...) : ")
				+ ia.what() );
	}
}

HudsonMSSegSiteParams::~HudsonMSSegSiteParams() {}

double HudsonMSSegSiteParams::getTheta() const
{
	return theta;
}
			
int HudsonMSSegSiteParams::getSegsitesin() const
{
	return segsitesin;
}

int HudsonMSSegSiteParams::getMfreq() const
{
	return mfreq;
}

std::string HudsonMSSegSiteParams::msString() const
{
	std::ostringstream stm;
	stm << std::fixed;
	
	if ( theta > 0.0 ) outputTheta(stm);
	if ( segsitesin > 0 ) outputSegsitesin(stm);
	if ( mfreq > 1) outputMfreq(stm);	
	
	return stm.str();
}

std::string HudsonMSSegSiteParams::toString() const
{
	std::ostringstream stm;
	stm << std::fixed << std::setprecision(precEst(theta));
	
	stm << "theta = " << theta;
	stm << ", segsitesin = " << segsitesin;
	stm << ", mfreq = " << mfreq;
	
	return stm.str();
}

void HudsonMSSegSiteParams::checkPars() const
{
	if (theta < 0.0) {
		throw std::invalid_argument("_theta < 0.0");
	}
	if (segsitesin < 0) {
		throw std::invalid_argument("_segsitesin < 0");
	}
	if (mfreq < 1) {
		throw std::invalid_argument("_mfreq < 1");
	}
}

std::ostringstream& HudsonMSSegSiteParams::outputTheta(
								std::ostringstream& stm)	const
{
	int prec = precEst(theta);
		
	stm.precision(prec);
	stm << " -t " << theta; 
	return stm;
}

std::ostringstream& HudsonMSSegSiteParams::outputSegsitesin(
								std::ostringstream& stm)	const
{
	stm << " -s " << segsitesin;
	return stm;
}

std::ostringstream& HudsonMSSegSiteParams::outputMfreq(
								std::ostringstream& stm)	const
{
	stm << " -F " << mfreq;
	return stm;
}

int HudsonMSSegSiteParams::precEst(double val)	const
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
