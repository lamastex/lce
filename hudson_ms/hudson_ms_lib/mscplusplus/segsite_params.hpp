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


#ifndef __HUDSON_MS_SEGSITE_PARAMS__
#define __HUDSON_MS_SEGSITE_PARAMS__

#include <string>
#include <iostream>


namespace hudson_ms {

	
	/*! \brief 
	A class for holding parameters relating to gene mutation. 
	* 
	This class is essentially a data structure plus some formatting
	for printing data values held, ie, to look like the parts
	of an ms command line relating to gene mutation parameters
	(see Hudson,R. 1983. Theor.Pop.Biol. 23: 183-201).  
	 
	*/
	class HudsonMSSegSiteParams {
		
		public:
		
			/*! \brief No-args constructor */
			HudsonMSSegSiteParams();
			
			/*! @name Constructors.
			
			\note There is no constructor for just parameters
			\a _theta and \a _segsitesin without specifying \a _mfreq.
			If you want to use just this combination then use the
			constructor with all three parameters and \a _mfreq = 1.

			See (Hudson,R. 1983. Generating samples under a Wright-Fisher neutral 
			model of genetic variation.  Theor.Pop.Biol. 23: 183-201)
			for more information about parameter values.
			* 
			Some parameter values are specified in terms of an 
			implicit \f$ N_0 \f$, 
			the effective subpopulation size. Parameter values
			are also given in the most general way, in terms of a
			a population of subpopulations.  If a single integral 
			population is being considered, this 
			is treated as a population consisting of a single
			subpopulation.  ie \a popnsams, \a popsizeratios and
			\a popgrowthrates will each have just one value,
			the value relating to the total
			population, and \f$ N_0 \f$ is the effective [total] 
			population size.  
			
			\param _theta is the mutation parameter specified 
			as \f$ \theta = 4N_0\mu \f$ where \f$ N_0 \f$ is the
			effective diploid population size and \f$ \mu \f$ is the neutral
			mutation rate for an <b>entire locus</b>.  
			\a _theta must be \f$ \geq 0\f$.
			\param _segsitesin can be used to specify a fixed number
			of segregating sites. If \a _segsitesin == 0, this parameter
			is irrelevant (an ms-like program using this with 
			\a _segsitesin = 0 should generate segsites randomly).
			If \a _segsitesin > 0, an ms-like program using this 
			should ensure that the number of segregating sites
			will be fixed to equal _segsitesin.
			\a _segsitesin must be \f$ \geq 0\f$.
			\param _mfreq is the minumum minor allele freq that
			there will be for any segregating site position
			in the segregating sites object made if
			an ms-like program uses this. \a _mfreq = 1
			should have no effect.  \a _mfreq must be \f$ \geq 1\f$.
			\pre \a _theta >= 0.0, 
			\a _segsitesin >= 0, _mfreq >= 1.*/
			//@{
			HudsonMSSegSiteParams(	double _theta, 
									int _segsitesin,
									int _mfreq);
			
			HudsonMSSegSiteParams(double _theta,
									int _mfreq);
			
			HudsonMSSegSiteParams(double _theta);
			
			HudsonMSSegSiteParams(int _segsitesin,
										int _mfreq);
			
			HudsonMSSegSiteParams(int _segsitesin);
			
			//@}
			
			~HudsonMSSegSiteParams();
			
			double getTheta() const;
			
			int getSegsitesin() const;
			
			int getMfreq() const;
			
			std::string msString() const;
			
			std::string toString() const;
			
						
		protected:
		
		private:
		
			void checkPars() const;

			std::ostringstream& outputTheta(
								std::ostringstream& stm)	const;
			
			std::ostringstream& outputSegsitesin(
								std::ostringstream& stm)	const;
			
			std::ostringstream& outputMfreq(
								std::ostringstream& stm)	const;
			
			int precEst(double val)	const;
								
			const double theta;
			
			const int segsitesin;
			
			const int mfreq;
			
	};
}

#endif
