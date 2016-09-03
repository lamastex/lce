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


/* *********  segtre_mig.c **********************************
*
*	This subroutine uses a Monte Carlo algorithm described in
*	Hudson,R. 1983. Theor.Pop.Biol. 23: 183-201, to produce
*	a history of a random sample of gametes under a neutral
*	Wright-Fisher model with recombination and geographic structure. 
*	Input parameters
*	are the sample size (nsam), the number of sites between
*	which recombination can occur (nsites), and the recombination
*	rate between the ends of the gametes (r). The function returns
*	nsegs, the number of segments the gametes were broken into
*	in tracing back the history of the gametes.  The histories of
*	these segments are passed back to the calling function in the
*	array of structures seglst[]. An element of this array,  seglst[i],
* 	consists of three parts: (1) beg, the starting point of
*	of segment i, (2) tree, which points to the first node of the
*	tree representing the history of the segment, (3) next, which
*	is the index number of the next segment.
*	     A tree is a contiguous set of 2*nsam nodes. The first nsam
*	nodes are the tips of the tree, the sampled gametes.  The other
*	nodes are the nodes ancestral to the sampled gametes. Each node
*	consists of an "abv" which is the number of the node ancestral to
*	that node, an "ndes", which is not used or assigned to in this routine,
*	and a "time", which is the time (in units of 4N generations) of the
*	node. For the tips, time equals zero.
*	Returns a pointer to an array of segments, seglst.

**************************************************************************/
#ifndef __HUDSON_MS_ARG_BUILDER__
#define __HUDSON_MS_ARG_BUILDER__

#include "recomb_params.hpp"
#include "arg.hpp"
#include "population_structure.hpp"
#include "event.hpp"
#include "prng.hpp"



#include <boost/smart_ptr.hpp>

#include <vector>
#include <stdexcept>

namespace hudson_ms {

	
	/*! \brief 
	 A class for building an ancestral recombination graph (ARG) given
	 information about the current-time population structure and 
	 past demographic events.  
	 
	This subroutine uses a Monte Carlo algorithm described in
	(Hudson,R. 1983. Generating samples under a Wright-Fisher neutral 
	model of genetic variation.  Theor.Pop.Biol. 23: 183-201), to produce
	a history of a random sample of gametes under a neutral
	Wright-Fisher model with recombination and geographic structure. 
	
	Input parameters
	The function returns
	nsegs, the number of segments the gametes were broken into
	in tracing back the history of the gametes.  The histories of
	these segments are passed back to the calling function in the
	array of structures seglst[]. An element of this array,  seglst[i],
	consists of three parts: (1) beg, the starting point of
	of segment i, (2) tree, which points to the first node of the
	tree representing the history of the segment, (3) next, which
	is the index number of the next segment.
	A tree is a contiguous set of 2*nsam nodes. The first nsam
	nodes are the tips of the tree, the sampled gametes.  The other
	nodes are the nodes ancestral to the sampled gametes. Each node
	consists of an "abv" which is the number of the node ancestral to
	that node, an "ndes", which is not used or assigned to in this routine,
	and a "time", which is the time (in units of 4N generations) of the
	node. For the tips, time equals zero.
	Returns a pointer to an array of segments, seglst.
	 */
	class HudsonMSArgFactory {
		
		public:
		
			/** \brief Constructor.
			
			The constructor parameters specify recombination
			(crossover and gene conversion), population structure, and
			the random number generation.
			
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
			population size.  The \a migration_matrix
			will similarly specify a one-row one-column matrix
			and the value actually given will be irrelevant.
			If a population with \f$ npop > 1 \f$
			subpopulations is being considered, \a popnsams, 
			\a popsizeratios and \a popgrowthrates will each have 
			\f$ npop \f$ values and the \a migration_matrix parameter
			will represent a matrix with \f$ npop \f$ rows and 
			\f$ npop \f$ columns. 
			
			Growth rates relate to exponential growth (or shrinkage) 
			of the subpopulations.
			Eg if the growth rate specified for a subpopulation is 
			\f$ \alpha \f$, then the population size at time \f$ t \f$, 
			before the present is \f$N(t) = N_0\exp^{\alpha t} \f$, 
			where time is measured in units of \f$ 4N_0 \f$ generations.
			\f$ \alpha < 0 \f$ indicates that the subpopulation was larger
			in the past than it is now (has been shrinking, in forward-time
			terms).
			 
			\param recombParam is a reference to a HudsonMSRecombParams
			object from which
			values for parameters controlling the chromosone recombinatino
			process (cross over and gene conversion, including the
			number of sites in the locus can be obtained from this).
			\param popnsams is an ordered collection of sub-population
			sample sizes (the number of chromosones sampled from
			each subpopulation) so that the number of elements in 
			\a popnsams is the number of subpopulations \f$ npop \f$
			and the element \a popnsams \f$ [i] \f$, \f$i = 0, \ldots npop-1 \f$
			is the number of samples (chromosones) 
			from the \f$ (i-i)^{th} \f$ subpopulation.
			\f$ \sum_{i=0}^{npop-1} popnams[i] = nsams \f$ where \f$ nsams \f$
			is the total sample size (total number of current-time 
			chromosones sampled).  
			\param popsizeratios is an ordered collection giving
			the ratios of the total size of 
			each subpopulation to the implicit \f$ N_0 \f$
			and the element \a popsizeratios \f$ [i] \f$, \f$i = 0, \ldots npop-1 \f$
			is the ratio of the effective population size of the
			the \f$ (i-i)^{th} \f$ subpopulation to the implicit \f$ N_0 \f$.
			\param popgrowthrates is an ordered collection giving
			the growth rates of each subpopulation
			and the element \a popgrowthrates \f$ [i] \f$, \f$i = 0, \ldots npop-1 \f$
			is the ratio of the effective population size of the
			the \f$ (i-i)^{th} \f$ subpopulation to the implicit \f$ N_0 \f$.
			\param migration_matrix is a specification of the 
			\b backwards migration matrix between subpopulations.  
			The rows of this matrix are given by the inner vectors in
			\a migration_matrix.  The value in the \f$ i^{th} \f$ row
			and \f$ j^{th} \f$ column of the migration matrix
			(the  \f$ j^{th} \f$ element in the \f$ i^{th} \f$ inner
			vector of \a migration_matrix) is \f$ 4N_0m_{ij} \f$ j^{th} \f$
			where \f$ m_{ij} \f$ is the fraction of the \f$ i^{th} \f$
			subpopulation made up of migrants from the \f$ j^{th} \f$
			subpopulation each generation (\f$i = 0, \ldots npop-1 \f$,
			\f$j = 0, \ldots npop-1 \f$). The diagonal elements
			of this matrix (the \f$ i^{th} \f$ element in the \f$ i^{th} \f$ inner
			vector of \a migration_matrix, \f$i = 0, \ldots npop-1 \f$)
			will be ignored.
			\param eList is a list of past demographic events, 
			represented as an ordered collection of 
			\link HudsonMSDemogEventPtr HudsonMSDemogEventPtrs\endlink.
			eList is expected to be ordered in ascending order
			of 'backwards time', ie if an event \f$ e_1 with associated
			time t_1 \f$ comes before an event \f$ e_2 with associated
			time t_2 \f$ in \a eList, then \f$ t_1 \leq t_2 \f$ and 
			\f$ e_1 \f$ took place closer to the present 
			than did \f$ e_2 \f$.*/
			
			HudsonMSArgFactory(
				const hudson_ms::HudsonMSRecombParams& recombParams,
				const boost::shared_ptr< const hudson_ms::PopulationStructure > pop,
				boost::shared_ptr < hudson_ms::HudsonMSPRNGen > _rr);
			
			
			~HudsonMSArgFactory();
			
			/* should be same as ms except when multiple demograhic events
			 * specified for the same time: this routine processes all of them
			 * before randomly scheduling more chomosone events, whereas ms
			 * will schedule another chromosone event even though this
			 * is guaranteed to never take place because there is a demographic
			 * event lined up for the exact current time.  Doing it my way
			 * ensures for example that the program will give the same
			 * results whether a whole migration matrix is reset or
			 * if individual elements are reset piecemeal.*/
			
			/* Also I have changed the way that cleft is calculated.
			 * This used to be set at the start and then 
			 * incremented/decremented in xover and
			 * commonAncestor according to changes in chromosones.  Now
			 * it is recalculated each time it is used.  On my tests
			 * this gives the same result as using the old methods, but
			 * the values for cleft are not totally identical (differences
			 * in smaller decimal places) and so 
			 * it is possible that there could be very small differences
			 * in the random numbers generated as cleft*genRand() and
			 * hence in the results */ 
			boost::shared_ptr < hudson_ms::HudsonMSArg > makeArg() const;
			
			
			// side effect = growths are set to 0
			void processEvent_N(double multNzero) const;

			void processEvent_n(size_t pop, double multNzero) const;
			
			void processEvent_G(double time, double newGrowthRate) const;
			
			void processEvent_g(double time, size_t pop, double newGrowthRate) const;
			
			// side effect - sets migration to other subpops from pop_i to zero
			void processEvent_j(size_t pop_i, size_t pop_j) const;
			
			// size of new subpop set to N0, no migration, no growth
			void processEvent_s(double time, size_t pop, double probStaying) const;

			void processEvent_M(double newMigRate) const;
			
			void processEvent_m(size_t pop_i, size_t pop_j, double newMigRate) const;
			
			void processEvent_ma(const std::vector< 
								std::vector < double > >& newMigMatrix) const;

			void processEvent_Coal(int cpop) const;
			
			void processEvent_Mig(double mig) const;

			void processEvent_Xover() const;
			
			
		protected:
		
		private:
		
			/*! A class for exceptions raised by population parameters in processing events.
			*/
			class Population_Index_Error : public std::runtime_error
			{
			   
			   public:
			   Population_Index_Error(std::string ss);
			   virtual ~Population_Index_Error() throw();
			   virtual const char* what() const throw();
			   
			};

			struct seg	{
				seg(int _beg, int _end, int _desc);
				int beg;
				int end;
				int desc;
			};

			struct chromo {
				chromo(int _pop, int _beg, int _end, int _desc); // one seg
				chromo(int _pop); // segs empty
				int links() const; //the number of links between beginning and end of chrom 
				int nseg() const;
				size_t pop;
				std::vector < seg > segments;
			};
			
			// cross over, migration, 
			class HudsonMSChromoEvent : public HudsonMSEvent
			{
				public:
				
					HudsonMSChromoEvent (double _time) : HudsonMSEvent(_time) {};

					/*! \brief Destructor */
					virtual ~HudsonMSChromoEvent () {};
			};

			// cross-over or gene conversion
			class HudsonMSChromoEvent_Xover : 
					public hudson_ms::HudsonMSArgFactory::HudsonMSChromoEvent
			{
				public :
				
					HudsonMSChromoEvent_Xover(double _time);
					
					void callEvent(const hudson_ms::HudsonMSArgFactory& builder) const;
				private :
					HudsonMSChromoEvent_Xover();
			};
			
			// migration
			class HudsonMSChromoEvent_Mig : 
					public hudson_ms::HudsonMSArgFactory::HudsonMSChromoEvent
			{
				public :
				
					HudsonMSChromoEvent_Mig(double _time, double _mig);
					
					void callEvent(const hudson_ms::HudsonMSArgFactory& builder) const;
				private :
					HudsonMSChromoEvent_Mig();
					double mig;
			};
			
			// coalescence
			class HudsonMSChromoEvent_Coal : 
					public hudson_ms::HudsonMSArgFactory::HudsonMSChromoEvent
			{
				public :
				
					HudsonMSChromoEvent_Coal(double _time, int _cpop);
					
					void callEvent(const hudson_ms::HudsonMSArgFactory& builder) const;
				private :
					HudsonMSChromoEvent_Coal();
					int cpop;
			};
			
			void setPars(const hudson_ms::HudsonMSRecombParams& recombParams,
				const boost::shared_ptr< const hudson_ms::PopulationStructure > pop);
			
			void initialise() const;
			
			void clearWorkspace() const;
			
			double getMig() const;
			
			bool checkMig() const;
			
			double genRand() const;
			
			double genRandGtZero() const;
			
			void acceptDemogEvent(const hudson_ms::HudsonMSDemogEventPtr& event) const;
			
			void acceptChromoEvent(boost::shared_ptr < 
					const hudson_ms::HudsonMSArgFactory::HudsonMSChromoEvent > event) const;
						
			int recombination() const;

			int cleftr() const;

			int cinr() const;

			void xover(int ic, int is) const;
			
			int commonAncestors(size_t c1, size_t c2) const;

			void pick2(int n, int& i,int& j) const;

			void pick2_chrom(size_t pop, int& pc1, int& pc2) const;

			bool isseg(int start, int c, int& psg) const;
			
			// total links in current chrom
			long int chrom_links() const;
			
			// cleft for current chrom
			double chrom_cleft() const;

			// the initial configuration parameters
			boost::shared_ptr < hudson_ms::HudsonMSPRNGen > rr;
			int nsites;
			double rho;
			size_t nPop;
			size_t nsam;
			double f;
			double track_len;
			std::vector < int > popNsams; // config
			std::vector < double > sizeRatios; //size
			std::vector < double > growthRates; // alphag
			
			/*
			The migration matrix is a backwards migration matrix where
			\f$ M_{ij} = 4N_0m_{ij} \f$ and \f$ m_{ij} \f$ is the fraction
			of subpopulation \f$ i \f$ made up of migrants <strong>from</strong>
			subpopulation \f$ j \f$ each generation.
			*/
			std::vector < std::vector < double > > migMat;

			std::vector< HudsonMSDemogEventPtr > demogEventList;

			
			// the workspace
			mutable size_t npop;
			mutable std::vector < int > config;
			mutable std::vector < double > size;
			mutable std::vector < double > alphag;

			mutable std::vector < std::vector < double > > migm; // migMat

			// the arg
			mutable boost::shared_ptr < hudson_ms::HudsonMSArg > arg_ptr;
			//mutable int nsegs;  // segments in the arg
						
			mutable double t;
			mutable std::vector < chromo > chrom; // set in initialise, changed in xover
			// starts with one chromo for each individual in each subpopulation, in order 
			// ie nsubsam_1 for subpop1, then nsubsam_2 for subpop2, etc

			//mutable long nlinks ; // set in initialise, changed in xover, commonAncestors
			
			mutable double r; // set in initialise
			
			mutable double prec; // all reset in each loop
			mutable double cin;
			mutable double clefta;
			mutable double prect;
	
			//mutable double cleft; // set in initialise, changed in xover, commonAncestors
			mutable double pc; // set in initialise
			mutable double lnpc ; // set in initialise
			mutable double rf; // set in initialise
			mutable double rft; // set in initialise


			mutable std::vector < double > tlast;

			mutable std::vector< HudsonMSDemogEventPtr >::const_iterator nextEventItr;

			//debug only
			std::ostream& outputChrom(std::ostream& os, int ic = -1, int ic2 = -1) const;
			std::ostream& outputChromo(std::ostream& os, const chromo& chr, int jseg = -1) const;
			std::ostream& outputSegments(std::ostream& os, const std::vector < seg >& segments) const;


	};
}

#endif
