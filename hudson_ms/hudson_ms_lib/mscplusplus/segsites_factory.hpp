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


#ifndef __HUDSON_MS_SEGSITES_FACTORY__
#define __HUDSON_MS_SEGSITES_FACTORY__

#include "segsites.hpp"
#include "arg.hpp"
#include "segsite_params.hpp"
#include "prng.hpp"

#include <boost/smart_ptr.hpp>

#include <vector>
#include <string>

namespace hudson_ms {
	
	/*! \brief 
	 A class for producing segregating sites for a sample using 
	 an equivalent method to that used in Hudson's program ms 
	(see Hudson,R. 1983. Generating samples under a Wright-Fisher neutral 
	model of genetic variation.  Theor.Pop.Biol. 23: 183-201).
	
	
	 */
	class HudsonMSSegSitesFactory {
		
		public:
		
			/*! Constructor.
			*/
			HudsonMSSegSitesFactory( 
				boost::shared_ptr < hudson_ms::HudsonMSPRNGen > _rr);
			
			~HudsonMSSegSitesFactory();
			
			/* should be same as ms*/
			/*! \brief Make a HudsonMSSegSites object.
			 
			Make and return a HudsonMSSegSites object using the given
			arg and other parameters.
			* 
			See (Hudson,R. 1983. Generating samples under a Wright-Fisher neutral 
			model of genetic variation.  Theor.Pop.Biol. 23: 183-201)
			for more information about parameter values.
			
			\param _arg_ptr is a pointer to an arg object representing
			the genealogical history of a sample of chromosones for a 
			locus.
			\param ssParam is a reference to a HudsonMSSegSiteParams
			object from which
			values for parameters controlling the gene mutation
			process (mutation parameter, fixed number of segsites
			if required, and minimum minor allele frequency if required
			can be obtained from this).
			\pre \a _arg_ptr must point to an arg containing at
			least one genealogical tree, and representing 
			a sample where there is at least one present-day individual
			(nsam > = 1), and where the number of sites in the locus
			is at least 2 (nsites >= 2). The mutation parameter given by 
			\a ssParams should be >= 0.0 and the fixed number of
			segregating sites required given by
			\a ssParams should be >= 0.  The minimum minor allele
			frequency given by \a ssParams must be >= 1 and <= nsam/2.*/
			boost::shared_ptr < hudson_ms::HudsonMSSegSites > 
				makeSegSites(
				const boost::shared_ptr < const hudson_ms::HudsonMSArg > _arg_ptr,
				const hudson_ms::HudsonMSSegSiteParams& ssParams) const;
						
		protected:
		
		private:
		
			void allocateSegSitesRandomly(
				std::vector < double >& treeTimes,
				std::vector < int >& segsites) const;
			
			void allocateFixedSegSites(
				std::vector < double >& treeTimes,
				std::vector < int >& segsites, double& probss) const;
				
			void addSegsites(
				const std::vector < double >& treeTimes,
				const std::vector < int >& segsites,
				int k, int segmnt) const;
			
			double genRand() const;
			
			double genRandGtZero() const;
			
			void initialise(
				const boost::shared_ptr < const hudson_ms::HudsonMSArg > _arg_ptr,
				const hudson_ms::HudsonMSSegSiteParams& ssParams) const;
			
			void clearWorkspace() const;
			
			int poisso(double u) const;
			
			double gasdev(double m, double v) const;
			
			std::vector < std::vector < int > > makeGametes(
						const std::vector < HudsonMSArg::node >& tree,
						double tt,
						int newsites) const; 
			
			int pickBranch(
					const std::vector < HudsonMSArg::node >& tree,
					double tt) const;

			bool isDescendant(const std::vector < HudsonMSArg::node >& tree,
								int tip, int node_picked ) const;

			
			std::vector < double > getNewPositions(int segsit,
									double beg, 
									double len) const;
			
			
			std::vector < int > allocatesSegSitesInOverTrees(
					const std::vector < double >& adjTreeTimes) const;
					
			// the initial configuration parameters
			boost::shared_ptr < hudson_ms::HudsonMSPRNGen > rr;
			
			// the workspace
			mutable boost::shared_ptr < const hudson_ms::HudsonMSArg > arg_ptr;
			mutable int nsegs;
			mutable int nsam;
			mutable int nsites;
			
			mutable double theta;
			mutable int segsitesin;
			mutable int mfreq;
			
			mutable int iset; // for gasdev
			mutable float gset; // for gasdev
	
			//mutable std::vector < double > positions;
			
			/* nsam inner vectors, 
			 * each containing segsite coding (0 or 1) for an individual.*/
			//mutable std::vector < std::vector < int > > list;
			
			mutable double segfac;
			
			// the segsites
			mutable boost::shared_ptr < hudson_ms::HudsonMSSegSites > ss_ptr;

	};
}

#endif
