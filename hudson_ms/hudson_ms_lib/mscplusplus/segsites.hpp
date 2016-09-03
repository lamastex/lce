/* mscplusplus - ms for c++.

   Copyright (C) 2011, 2012 Jenny Harlow
   
   This code was adapted for C++ from ms (Richard R. Hudson).
 
   ms - a program for generating samples under neutral models.
  
   Hudson, R. R. (2002).  Generating samples under a Wright-Fisher neutral 
   model of genetic variation.  Bioinformatics 18: 337-338.
  
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


#ifndef __HUDSON_MS_SEGSITES__
#define __HUDSON_MS_SEGSITES__

#include <vector>
#include <string>
#include <iostream>


namespace hudson_ms {

	
	/*! \brief 
	A class for holding seg site information for a 
	sample of \a nsam present-day individuals.  
	 
	*/
	class HudsonMSSegSites {
		
		public:
		
			/*! Constructor.
			
			\param _nsam is the number of individuals in the sample.
			\param _segsites is the fixed number of segsites specified
			by the user 
			(default -1 indicates no fixed number of segsites).*/
			HudsonMSSegSites(int _nsam, bool _showProbss);
			
			~HudsonMSSegSites();
			
			/*! \brief Set probss.*/
			void setProbss(double _probss);
			
			/*! \brief Add new positions to this.*/
			void addPositions(
					const std::vector < double >& newPos,
					const std::vector < std::vector < int > >& newSS);
			
			/*! \brief Output a description of this.
			 * 
			 * Output in same format as Hudson's seg sites.
			 * 
			 * If default precisions are used, output
			 * should be same as Hudson's.*/
			std::ostream& outputSegSites(std::ostream& os, int precPos = 4, int precProb = -1) const;
			
			/*! \brief Get a description of this as a string.
			 * 
			 * Same format as Hudson's seg sites.
			 * 
			 * If default precisions are used, precision in result
			 * should be same as Hudson's.*/
			std::string toString(int precPos = 4, int precProb = 6) const;
			
						
		protected:
		
		private:
		
			HudsonMSSegSites();

			//debug only
			//std::ostream& outputChrom(std::ostream& os, int ic = -1, int ic2 = -1) const;
			//std::ostream& outputChromo(std::ostream& os, const chromo& chr, int jseg = -1) const;
			//std::ostream& outputSegments(std::ostream& os, const std::vector < seg >& segments) const;

			size_t nsam;
			
			bool showProbss;
			
			double probss;
			
			std::vector < double > positions;
			
			// the segsites information, one inner vector for each individual
			// giving the seg sites at that individual
			// nsam inner vectors in total, each same length as positions
			std::vector < std::vector < int > > segsites;

	};
}

#endif
