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

#ifndef __HUDSON_MS_RECOMB_PARAMS__
#define __HUDSON_MS_RECOMB_PARAMS__

#include <string>
#include <iostream>


namespace hudson_ms {

	
	/*! \brief 
	A class for holding parameters relating to chromosone recombination,
	* ie cross-over and gene conversion. 
	* 
	This class is essentially a data structure plus some formatting
	for printing data values held, ie, to look like the parts
	of an ms command line relating to recombination parameters
	(see Hudson,R. 1983. Theor.Pop.Biol. 23: 183-201).  
	 
	*/
	class HudsonMSRecombParams {
		
		public:
		
			/*! \brief No-args constructor */
			HudsonMSRecombParams();
			
			/*! @name Constructors.
			
			\note There is no constructor for just parameters
			\a _conversion_parameter and \a _conversion_track_len
			(and \a _nsites) without specifying \a _rho.
			If you want to use just this combination then use the
			constructor with all four parameters and \a _rho = 0.

			See (Hudson,R. 1983. Generating samples under a Wright-Fisher neutral 
			model of genetic variation.  Theor.Pop.Biol. 23: 183-201)
			for more information about parameter values.
			
			Some parameter values are specified in terms of an 
			implicit \f$ N_0 \f$, 
			the effective subpopulation size. Parameter values
			are also given in the most general way, in terms of a
			a population of subpopulations.  If single integral population
			is being considered, this 
			is treated as a population consisting of a single
			subpopulation.  ie \a popnsams, \a popsizeratios and
			\a popgrowthrates will each have just one value,
			the value relating to the total
			population, and \f$ N_0 \f$ is the effective [total] 
			population size.  
			
			\param _rho is the crossover parameter.  _rho
			represents \f$ \rho = 4N_0r \f$ where \f$ r \f$ is the 
			probability of crossover per generation between the ends
			of the locus being simulated.  \a _rho must be \f$ \geq 0\f$.
			\param _nsites is the number of sites between which recombination
			can occur, which can be thought of as the number of
			base pairs (ie length) in the locus being simulated.
			If \a rho > 0.0, \a _nsites must be \f$ \geq 2\f$. 
			\param _conversion_parameter is the gene conversion parameter
			(see Hudson,R. 1983. Theor.Pop.Biol. 23: 183-201).
			If \a _rho > 0.0, \a _f should give the value for 
			\f$ \frac{g}{r} \f$ where \f$ g \f$
			is the probability of gene conversion initiating 
			between a pair
			of specified sites in a given chromosone and \f$ r \f$ is 
			the probability per generation of crossing over between
			adjacent \b sites (cf the definition for \a _rho).  If 
			\a _rho = 0.0, there can be gene conversion but no crossover 
			and \a _f should give the value of \f$ 4N_0g \f$ where
			\f$ g \f$ is the probability of gene conversion initiating 
			between a pair of specified sites in a given chromosone.
			\a _f must be \f$ \geq 0.0\f$.
			\param _conversion_track_len is the mean conversion track length (see 
			Hudson,R. 1983. Theor.Pop.Biol. 23: 183-201).
			\a _ftlen must be \f$ \geq 0.0\f$. 
			\pre \a _rho >= 0.0.
			\pre If \a _rho > 0.0 or \a _conversion_parameter > 0.0 then \a _nsites > 1.
			\pre \a _conversion_parameter >= 0.0. 
			\pre \a _conversion_track_len >= 0.0.
			\pre If \a _conversion_parameter > 0.0 then \a _conversion_track_len > 0.0.	*/
			//@{
			HudsonMSRecombParams(	double _rho, 
									int _nsites,
									double _conversion_parameter,
									double _conversion_track_len);
			
			HudsonMSRecombParams(	double _rho, 
									int _nsites);
			
			//@}
			
			~HudsonMSRecombParams();
			
			double getRho() const;
			
			int getNsites() const;
			
			double getConversionParameter() const;
			
			double getConversionTrackLen() const;
			
			std::string msString() const;
			
			std::string toString() const;
			
						
		protected:
		
		private:
		
			void checkPars() const;

			std::ostringstream& outputCrossover(
								std::ostringstream& stm)	const;
			
			std::ostringstream& outputConversion(
								std::ostringstream& stm)	const;
			
			int precEst(double val)	const;
								
			const double rho; 
									
			const int nsites;
			
			const double conversion_parameter;
			
			const double conversion_track_len;
			
	};
}

#endif
