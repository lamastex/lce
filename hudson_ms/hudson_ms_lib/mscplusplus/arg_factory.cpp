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


#include "arg_factory.hpp"

#include <string>
#include <cmath>
#include <numeric>
#include <algorithm>

#include <cassert>


//#define MYDEBUG

//#define MYDEBUG1

//#define MYDEBUG2

using namespace hudson_ms;


HudsonMSArgFactory::Population_Index_Error::Population_Index_Error(std::string ss) 
	: std::runtime_error(ss) {}

HudsonMSArgFactory::Population_Index_Error::~Population_Index_Error() throw() {}

const char* HudsonMSArgFactory::Population_Index_Error::what() const throw() 
{
	try {
		return std::runtime_error::what();
	}
	catch (std::exception& e) {
		std::cerr << "HudsonMSArgFactory::Population_Index_Error: Problem in what():\n" << (e.what()) << std::endl;
	}
	
}

HudsonMSArgFactory::seg::seg(int _beg, int _end, int _desc)
	: beg(_beg), end(_end), desc(_desc)
{}

HudsonMSArgFactory::chromo::chromo(int _pop,
								int _beg, int _end, int _desc)
	: pop(_pop)
{
	segments.push_back(seg(_beg, _end, _desc));
}

HudsonMSArgFactory::chromo::chromo(int _pop)
	: pop(_pop)
{}

int HudsonMSArgFactory::chromo::links() const
{
	int link = 0;
	if (!segments.empty()) {
		link = segments.back().end - segments.front().beg;
	}
	return link;
}

int HudsonMSArgFactory::chromo::nseg() const
{
	return segments.size();
}

// constructor not using converstion parameter and conversion track len
HudsonMSArgFactory::HudsonMSArgFactory(
			const hudson_ms::HudsonMSRecombParams& recombParams,
			const boost::shared_ptr< const hudson_ms::PopulationStructure > pop,
			boost::shared_ptr < hudson_ms::HudsonMSPRNGen > _rr	)
			: 	rr(_rr)
{
	try {
		setPars(recombParams, pop);
	}
	catch (std::invalid_argument& ia) {
		throw std::invalid_argument(
			std::string("HudsonMSArgFactory::HudsonMSArgFactory(...) : ")
			+ ia.what() );
	}
}


HudsonMSArgFactory::~HudsonMSArgFactory() {}


boost::shared_ptr < hudson_ms::HudsonMSArg > 
									HudsonMSArgFactory::makeArg() const
{
	try {
		initialise();

		#ifdef MYDEBUG
			std::cout << "\nAfter initialising" << std::endl;
		#endif
		
		while( chrom.size() > 1 ) {

			#ifdef MYDEBUG
				std::cout << "\nIn loop, nchrom = " << chrom.size() << std::endl;
			#endif
			
			long int nlinks = chrom_links();
			prec = nlinks*r;
			cin = nlinks*rf ;
			//assert ( abs(cleft - chrom_cleft()) < 0.0000000001 );
			//clefta = cleft*rft ;
			clefta = chrom_cleft()*rft ;
			prect = prec + cin + clefta ;

			#ifdef MYDEBUG1
				std::cout << "r = " << r << std::endl;
				std::cout << "f = " << f << std::endl;
				std::cout << "track_len = " << track_len << std::endl;
				std::cout << "rf = " << rf << std::endl;
				std::cout << "rtf = " << rft << std::endl;
				std::cout << std::endl;
				std::cout << "nlinks = " << nlinks << std::endl;
				std::cout << "cleft = " << chrom_cleft() << std::endl;
				std::cout << "prec = " << prec << std::endl;
				std::cout << "cin = " << cin << std::endl;
				std::cout << "clefta = " << clefta << std::endl;
				std::cout << "prect = " << prect << std::endl;
			#endif

			boost::shared_ptr < const hudson_ms::HudsonMSArgFactory::HudsonMSChromoEvent > chromEvent;
			
			//int cpop = -1; // record which pop has coalescent event
			//char event; // record what event is due to happen
			
			double mig = getMig(); 
			
			// check the mig makes sense
			if( (npop > 1) && ( mig == 0.0) && ( !(nextEventItr < demogEventList.end() ) ) 
														&& !checkMig() ) {
				throw std::logic_error(
				"HudsonMSArgFactory::makeArg() : Infinite coalescent time. No migration");
			}
			
			//int eflag = 0 ;  // indicator for whether an event has been scheduled

			if( prect > 0.0 ) {      // cross-over or gene conversion
				
				double rdum = genRandGtZero();
				double ttemp = -log( rdum)/prect ;
				
				if( (!chromEvent) || (ttemp < chromEvent->getTime() ) ) {
					chromEvent = boost::shared_ptr 
						< const hudson_ms::HudsonMSArgFactory::HudsonMSChromoEvent >
							( new hudson_ms::HudsonMSArgFactory::HudsonMSChromoEvent_Xover(ttemp) );

					#ifdef MYDEBUG1
						std::cout << "xover event created" << std::endl;
						std::cout << "prect = " << prect << "\trdum = " << rdum << std::endl;
						std::cout << "\ttime is = " << chromEvent->getTime() << std::endl;
					#endif
				}	
			}

			if  (mig > 0.0 ) {         // migration
				
				double rdum = genRandGtZero();
				double ttemp = -log( rdum)/mig ;
				
				if( (!chromEvent) || (ttemp < chromEvent->getTime() ) ) {
					chromEvent = boost::shared_ptr 
						< const hudson_ms::HudsonMSArgFactory::HudsonMSChromoEvent >
							( new hudson_ms::HudsonMSArgFactory::HudsonMSChromoEvent_Mig(ttemp, mig) );
					
					#ifdef MYDEBUG1
						std::cout << "migration event created" << std::endl;
						std::cout << "mig = " << mig << "\trdum = " << rdum << std::endl;
						std::cout << "\ttime is = " << chromEvent->getTime()  << std::endl;
					#endif

				}
			}
			
			for(size_t pop = 0; pop < npop ; ++pop) {     // coalescent
				
				double coal_prob = 1.0 * config[pop] * (config[pop]-1.0) ;
				
				if( coal_prob > 0.0 ) {
					
					double rdum = genRandGtZero();
					
					if( alphag[pop] == 0 ){
						double ttemp = -log( rdum )* size[pop] /coal_prob ;
						
						if( (!chromEvent) || (ttemp < chromEvent->getTime() ) ) {
							chromEvent = boost::shared_ptr 
								< const hudson_ms::HudsonMSArgFactory::HudsonMSChromoEvent >
									( new hudson_ms::HudsonMSArgFactory::HudsonMSChromoEvent_Coal(ttemp, pop) );
							#ifdef MYDEBUG1
								std::cout << "coalescent type 1 event created" << std::endl;
								std::cout << "coal_prob = " << coal_prob << "\trdum = " << rdum << std::endl;
								std::cout << "cpop = " << pop << std::endl;
								std::cout << "\ttime is = " << chromEvent->getTime()  << std::endl;
							#endif
						}
					}
					else {
						double arg_ptr  = 1.0 - 
							alphag[pop]*size[pop]*exp(-alphag[pop]
								*(t - tlast[pop] ) )* log(rdum) / coal_prob;
						
						if( arg_ptr > 0.0 ) {                          //if arg_ptr <= 0,  no coalescent within interval
							double ttemp = log( arg_ptr ) / alphag[pop]  ;
							
							if( (!chromEvent) || (ttemp < chromEvent->getTime() ) ) {
								chromEvent = boost::shared_ptr 
									< const hudson_ms::HudsonMSArgFactory::HudsonMSChromoEvent >
										( new hudson_ms::HudsonMSArgFactory::HudsonMSChromoEvent_Coal(ttemp, pop) );

								#ifdef MYDEBUG1
									std::cout << "coalescent type 2 event created" << std::endl;
									std::cout << "coal_prob = " << coal_prob << "\trdum = " << rdum << std::endl;
									std::cout << "cpop = " << pop << std::endl;
									std::cout << "\ttime is = " << chromEvent->getTime()  << std::endl;
								#endif

							}
						}
					}
				} //end if coal_prob > 0.0		
			} // end loop through pops

			if ( (!chromEvent) && !( nextEventItr < demogEventList.end() ) ) {
				std::string msg = std::string("HudsonMSArgFactory::makeArg() : ")
					+ std::string("infinite time to next event. ") + 
					std::string("Negative growth rate in last time interval ") +
					std::string("or non-communicating subpops.");
				
				throw std::logic_error(msg);
			}
			
			// if there are is at least one demographic event to go and either no chrom event   
			// or the first demographic event left is scheduled first
			
			#ifdef MYDEBUG1
				if (demogEventList.empty()) {
					std::cout << "nothing in demographic list" << std::endl;
				}
				else if (!(nextEventItr < demogEventList.end())) {
					std::cout << "nothing left in demographic list" << std::endl;
				}
				else {
					size_t still_to_go = distance(nextEventItr, demogEventList.end());
					std::cout << "size of demographic events list still to go is = " << still_to_go << std::endl;
				}
			#endif
			
			#ifdef MYDEBUG1
				std::cout << "\nAbout to check events" << std::endl;
				
			#endif
			
			if( nextEventItr < demogEventList.end() 
					&& ( !(chromEvent)  
						|| ((t + chromEvent->getTime()) >=  nextEventItr->getTime() ) ) ){
							
				t = nextEventItr->getTime();
				
				// consume all demographic events specified for this exact time together
				while ( (nextEventItr < demogEventList.end()) && 
					(nextEventItr->getTime() == t) ) {
						
					#ifdef MYDEBUG
						std::cout << "about to accept a demographic event, time " << t << std::endl;
					#endif
					
					acceptDemogEvent(*nextEventItr);
					
					nextEventItr++;
				}
			}
			else if (chromEvent) {  

				t += chromEvent->getTime();
				
				acceptChromoEvent(chromEvent);
				
			}// end else

		}  // end the main while loop

		#ifdef MYDEBUG
		{
			std::cout << "End buildArg, nsegs = " << arg_ptr->nsegs() << std::endl;
			arg_ptr->outputSegLst(std::cout);
			
		}
		#endif	
		
		//  finish off - fill in the end properties for each segl
		int segmnt=0; 
		int ns = arg_ptr->nsegs();
		for(int k=0; k < ns; ++k) { 
			
			int end = ( k < ns-1 ? arg_ptr->at(arg_ptr->at(segmnt).next).beg - 1 : nsites-1 );
			arg_ptr->at(segmnt).end = end;
			
			segmnt = arg_ptr->at(segmnt).next; 
		}
		
		boost::shared_ptr < hudson_ms::HudsonMSArg > tmp_arg_ptr = arg_ptr;

		// arg should be built
		clearWorkspace(); // clear the workspace
		
		return tmp_arg_ptr;
	}
	catch (Population_Index_Error& pie) {
		throw std::runtime_error(
		std::string("ArgFactory::makeArg() :\n") + std::string(pie.what()) );
	} 
}


// event processing


// side effect = growths are set to 0
void HudsonMSArgFactory::processEvent_N(double multNzero) const
{
	#ifdef MYDEBUG
		std::cout << "about process event N" << std::endl;
	#endif

	for(size_t pop = 0; pop < npop; ++pop) {
		size[pop]= multNzero ;
		alphag[pop] = 0.0 ;
	}
}

// side effect = growth is set to 0
void HudsonMSArgFactory::processEvent_n(size_t pop, double multNzero) const
{
	#ifdef MYDEBUG
		std::cout << "about process event n" << std::endl;
	#endif
	
	if (pop >= npop) {
		throw HudsonMSArgFactory::Population_Index_Error(
			"HudsonMSArgFactory::processEvent_n(size_t, double) : pop");
	}

	size.at(pop) = multNzero ;
	alphag[pop] = 0.0 ;
}

	
void HudsonMSArgFactory::processEvent_G(double time, 
										double newGrowthRate) const
{
	#ifdef MYDEBUG
		std::cout << "about process event G" << std::endl;
	#endif

	for(size_t pop =0; pop <npop; ++pop) {
		size[pop] = size[pop]*exp( -alphag[pop]*(time - tlast[pop]) ) ;
		alphag[pop]= newGrowthRate ;
		tlast[pop] = time ;
	}
}


void HudsonMSArgFactory::processEvent_g(double time, 
							size_t pop, double newGrowthRate) const
{
	#ifdef MYDEBUG
		std::cout << "about process event g" << std::endl;
	#endif
	
	if (pop >= npop) {
		throw HudsonMSArgFactory::Population_Index_Error(
			"HudsonMSArgFactory::processEvent_g(size_t, double) : pop");
	}

	size.at(pop) = size.at(pop)*exp( - alphag[pop]*(time - tlast[pop]) ) ;
	alphag.at(pop)= newGrowthRate ;
	tlast.at(pop) = time ;
}


// side effect - sets migration to other subpops from pop_i to zero
void HudsonMSArgFactory::processEvent_j(size_t pop_i, size_t pop_j) const
{
	#ifdef MYDEBUG
		std::cout << "about process event j" << std::endl;
	#endif

	if (pop_i >= npop && pop_j >= npop) {
		throw HudsonMSArgFactory::Population_Index_Error(
			"HudsonMSArgFactory::processEvent_j(size_t, size_t) : pop_i, pop_j");
	}
	if (pop_i >= npop) {
		throw HudsonMSArgFactory::Population_Index_Error(
			"HudsonMSArgFactory::processEvent_j(size_t, size_t) : pop_i");
	}
	if (pop_j >= npop) {
		throw HudsonMSArgFactory::Population_Index_Error(
			"HudsonMSArgFactory::processEvent_j(size_t, size_t) : pop_j");
	}

	config.at(pop_j) += config.at(pop_i) ;
	config[pop_i] = 0 ;
	for(size_t ic = 0; ic < chrom.size(); ic++) {
		if( chrom.at(ic).pop == pop_i ) chrom[ic].pop = pop_j ;
	}
	
	/* M_k,pop_i = 4N0m_k,pop_i where m_k,pop_i is fraction of subpop k
	 * made up of individuals from subpop pop_i each generation */
	for(size_t k = 0; k < npop; k++) {
		if( k != pop_i) {
			migm.at(k).at(pop_i) = 0.0 ;
		}
	}
}


void HudsonMSArgFactory::processEvent_s(double time, size_t pop, double probStaying) const
{
	#ifdef MYDEBUG
		std::cout << "about process event s" << std::endl;
	#endif

	if (pop >= npop) {
		throw HudsonMSArgFactory::Population_Index_Error(
			"HudsonMSArgFactory::processEvent_i(double, size_t, double) : pop");
	}
	
	/*split  pop i into two;p is the proportion from pop i, and 1-p from pop n+1  */
	
	// we need to make room for one sub more population
	npop ++;

	config.push_back(0);
	size.push_back(1.0);
	alphag.push_back (0.0);
	tlast.push_back(time);
	
	for (std::vector < std::vector < double > >::iterator it = migm.begin();
					it < migm.end(); ++it) {
		it->push_back(0.0);				
	}
	migm.push_back(std::vector < double > (npop, 0.0) );
	
	config.at(pop) = 0 ;
	size_t n = chrom.size();
	for(size_t ic = 0; ic < n; ++ic) {
		if( chrom[ic].pop == pop ) {
			if( genRand() < probStaying ) config[pop]++;
			else {
				chrom[ic].pop = npop-1 ;
				config[npop-1]++;
			}
		}
	}
}

void HudsonMSArgFactory::processEvent_M(double newMigRate) const
{
	#ifdef MYDEBUG
		std::cout << "about process event M" << std::endl;
	#endif

	for(size_t pop =0; pop < npop; ++pop) {
		for(size_t pop2 = 0; pop2 < npop; ++pop2) {
			if (pop == pop2) migm.at(pop).at(pop2) = 0.0;
			else migm.at(pop).at(pop2) = newMigRate/(npop-1.0);
		}
	}
}


void HudsonMSArgFactory::processEvent_m(size_t pop_i, size_t pop_j, double newMigRate) const
{
	#ifdef MYDEBUG
		std::cout << "about process event m" << std::endl;
	#endif

	if (pop_i >= npop && pop_j >= npop) {
		throw HudsonMSArgFactory::Population_Index_Error(
			"HudsonMSArgFactory::processEvent_m(size_t, size_t) : pop_i, pop_j");
	}
	if (pop_i >= npop) {
		throw HudsonMSArgFactory::Population_Index_Error(
			"HudsonMSArgFactory::processEvent_m(size_t, size_t) : pop_i");
	}
	if (pop_j >= npop) {
		throw HudsonMSArgFactory::Population_Index_Error(
			"HudsonMSArgFactory::processEvent_m(size_t, size_t) : pop_j");
	}
	
	migm.at(pop_i).at(pop_j)= newMigRate;
	
}


void HudsonMSArgFactory::processEvent_ma(const std::vector< 
									std::vector < double > >& newmigMatrix) const
{
	#ifdef MYDEBUG
		std::cout << "about process event ma" << std::endl;
	#endif

	migm = newmigMatrix;
	for(size_t pop = 0; pop < npop; ++pop) {
		// make sure diagonals of mig matrix are 0
		migm.at(pop).at(pop) = 0.0;
	}
}

void HudsonMSArgFactory::processEvent_Xover() const
{
	#ifdef MYDEBUG
		std::cout << "about process chromo event Xover : ";
	#endif

	double ran = genRand();
	if( ran < ( prec / prect ) ){ /*recombination*/
		#ifdef MYDEBUG
			std::cout << "recomb event" << std::endl;
		#endif
		
		int rchrom = recombination(); // the index of the chromosone recombd
		config[ chrom[rchrom].pop ] += 1 ;
	}
	else if( ran < (prec + clefta)/(prect) ){    /*  cleft event */
		#ifdef MYDEBUG
			std::cout << "cleft event" << std::endl;
		#endif
		int rchrom = cleftr();
		config[ chrom[rchrom].pop ] += 1 ;
	}
	else  {         /* cin event */
		#ifdef MYDEBUG
			std::cout << "cin event" << std::endl;
		#endif
		int rchrom = cinr();
		if( rchrom >= 0 ) config[ chrom[rchrom].pop ] += 1 ;
	}
}

void HudsonMSArgFactory::processEvent_Mig(double mig) const
{
	#ifdef MYDEBUG
		std::cout << "about process chromo event Mig" << std::endl;
	#endif

	double x = mig*genRand();
	double sum = 0.0 ;
	
	size_t i = 0;
	size_t n = chrom.size();
	for(i = 0; i < n; ++i) {
		
		sum = std::accumulate(migm[chrom[i].pop].begin(), migm[chrom[i].pop].end(), sum);
		
		if( x <sum ) break;
	}
	
	size_t migrant = i ;
	
	x = genRand() * std::accumulate(migm[chrom[i].pop].begin(), migm[chrom[i].pop].end(), 0.0);
	
	sum = 0.0;
	size_t j = 0;
	for(j=0; j < npop; ++j) {
	
		if( j != chrom[migrant].pop ) {
			sum += migm[chrom[migrant].pop][j];
			if( x < sum ) break;
		}
	}
	
	size_t source_pop = j;
	
	config[chrom[migrant].pop] -= 1;
	config[source_pop] += 1;
	chrom[migrant].pop = source_pop ;
}

void HudsonMSArgFactory::processEvent_Coal(int cpop) const
{
	#ifdef MYDEBUG
		std::cout << "about process chromo event Coal" << std::endl;
	#endif

	/* coalescent event */
	/* pick the two, c1, c2  */
	int c1 = 0;
	int c2 = 0;
	/* c1 and c2 are chrom's to coalesce */
	pick2_chrom( cpop, c1, c2);  /* c1 and c2 are updated */
	
	#ifdef MYDEBUG
		std::cout << "pop where coalescent event takes places is " << cpop << " and chromosones chosen are " << c1 << " and " << c2 << std::endl;
	#endif

	
	int dec = commonAncestors (c1, c2);
	
	config[cpop] -= dec ;
}


// private methods

void HudsonMSArgFactory::setPars(
				const hudson_ms::HudsonMSRecombParams& recombParams,
				const boost::shared_ptr< const hudson_ms::PopulationStructure > pop)
{
	
	if (!(rr)) {
		throw std::invalid_argument("_rr empty");
	}
	
	if (!(pop)) {
		throw std::invalid_argument("pop empty");
	}
	
	if (pop->totalNsam() < 1) {
		throw std::invalid_argument("nsam from pop < 1");
	}
	
	nsam = pop->totalNsam();
	
	popNsams = pop->makePopnsams();
	nPop = popNsams.size();
	
	sizeRatios = pop->makePopsizeratios();
	growthRates = pop->makePopgrowthrates();
	migMat = pop->makeMigration_matrix();
	demogEventList = pop->makeEventList();
		
	// sort the event list but keep relative ordering for equal events
	/* popstructures sort as well so this is double-doing */
	if (!demogEventList.empty()) {
		std::stable_sort( demogEventList.begin(), demogEventList.end() );
	}
	
	/* no checks on these - recombParams object should have checked them,
	except that we actually use a minimum of nsites = 2 always.  If rho
	or conversion parameter are > 0 then this should be the case anyway,
	so it just makes setting irrelevant parameters easier (avoids
	division by zero if nsites = 1) when 
	rho and conversion parameter are both 0, ie where there no 
	coalescence or crossover anyway.*/
	nsites = 2;
	int _nsites = recombParams.getNsites();
	if (_nsites > 2) nsites = _nsites;
	rho = recombParams.getRho();
	f = recombParams.getConversionParameter();
	track_len = recombParams.getConversionTrackLen();
	
	#ifdef MYDEBUG1
		std::cout << "In setPars, nsites = " << nsites << std::endl;
	#endif

	// make sure diagonals of mig matrix are 0
	for(size_t pop = 0; pop < nPop; ++pop) {
				
		migMat[pop][pop] = 0.0;
	}
	#ifdef MYDEBUG1
		std::cout << "end setPars\n " << std::endl;
	#endif
	
}

double HudsonMSArgFactory::getMig() const
{
	double mig = 0.0;
	
	/* each row i gives migration rate parameters into i from (col index) j
	 * each generation as 4N0m_ij, m_ij being the fraction of 
	 * subpop i made up of migrants from j each generation.
	 * So, this adds up total migrants into i each generation */
	
	#ifdef MYDEBUG1
			std::cout << "Checking mig:" << std::endl; 
	#endif
		
	for(size_t i = 0; i < npop; ++i) {
		double thisRowMigTotal = 0.0;
		for(size_t j = 0; j < npop; ++j) {
			if( j != i ) thisRowMigTotal += migm[i][j];
		}
		
		// accumulate total migration out of pop * subpop size
		mig += config[i] * thisRowMigTotal;
		#ifdef MYDEBUG1
			std::cout << std::endl;	
			std::cout << "i = " << i << "\tconfig[i] = " << config[i] 
				<< "\tthisRowMigTotal = " << thisRowMigTotal 
				<< "\tconfig[i] * thisRowMigTotal = " << config[i] * thisRowMigTotal << std::endl;
			std::cout << "mig = " << mig << std::endl;
		#endif

	}
	return mig;
}
	
bool HudsonMSArgFactory::checkMig() const
{
		
	int count = 0;
	for( size_t j = 0; j < npop; ++j) {
		if( config[j] > 0 ) count++;
	}
	
	return (count <= 1); // true if okay
}
	
	
// a rand > 0.0
double HudsonMSArgFactory::genRandGtZero() const
{
	return rr->genRandGtZero();
	
}

// a rand
double HudsonMSArgFactory::genRand() const
{
	return rr->genRand();
	
}

void HudsonMSArgFactory::acceptDemogEvent(const HudsonMSDemogEventPtr& event) const
{
	event.callEvent(*this);
}

void HudsonMSArgFactory::acceptChromoEvent(
		boost::shared_ptr < 
			const hudson_ms::HudsonMSArgFactory::HudsonMSChromoEvent > event) const
{
	if (event) {
		event->callEvent(*this);
	}
}

void HudsonMSArgFactory::initialise() const
{
	// copy these so we can change them
	// the workspace
	config = popNsams;
	size = sizeRatios;
	alphag = growthRates;
	migm = migMat;
	
	npop = nPop;
	
	r = rho/(nsites-1);
	
	t = 0;
	
	arg_ptr = boost::shared_ptr < HudsonMSArg >( new HudsonMSArg(nsam) );
	// makes a 'tree' with one segment, nsam nodes 	
	//nsegs = 1;
	
	pc = 1.0 ;
	if( f > 0.0 ) 	{
		pc = (track_len -1.0)/track_len ;
	}
	
	lnpc = log( pc ) ;
	
	rf = f/(nsites - 1);
	if( r > 0.0 ) rf = r*f ;
	
	rft = rf*track_len ;
	
	/*
	if( f > 0.0 ) 	pc = (track_len -1.0)/track_len ;
	else pc = 1.0 ;
	lnpc = log( pc ) ;
	cleft = nsam* ( 1.0 - pow( pc, (double)(nsites-1) ) ) ;
	if( r > 0.0 ) rf = r*f ;
	else rf = f /(nsites-1) ;
	rft = rf*track_len ;
	*/

	prec = 0.0;
	cin = 0.0;
	clefta = 0.0;
	prect = 0.0;

	// record last change for each pop in tlast
	std::vector < double > tmp1(npop, 0);
	tlast.swap(tmp1);

	
	nextEventItr = demogEventList.end();
	if (!demogEventList.empty()) {
		
		nextEventItr = demogEventList.begin();
	}

	std::vector < chromo > tmp3; 
	chrom.swap ( tmp3 );
	
	int ind = 0;
	for(size_t pop = 0; pop < npop; ++pop) {
		
		// make chromos, 1 seg each, each seg 0 - nsites-1
		for(int j =0 ; j < config[pop]; ++j) { //each subpop size is in config
			
			int seg_beg = 0;
			chrom.push_back( chromo(pop, seg_beg, nsites-1, ind) );
			++ind;
		}
	}
	
	assert(chrom.size() == nsam);
	//cleft = nsam* ( 1.0 - pow( pc, (double)(nsites-1) ) ) ;
	//cleft = chrom_cleft();
	assert ( chrom_cleft() == nsam* ( 1.0 - pow( pc, (double)(nsites-1) ) ) );
	
	//nlinks = static_cast<long>(nsam*(nsites-1));
	//assert(nlinks == chrom_links());
		
	#ifdef MYDEBUG1
		std::cout << "\nend of initialise" << std::endl;
		std::cout << "nsam = " << nsam << std::endl;
		std::cout << "npop = " << npop << std::endl;
		std::cout << "nsegs = " << arg_ptr->nsegs() << std::endl;
		std::cout << "nchrom = " << chrom.size() << std::endl;
		//std::cout << "nnodes[0] = " << nnodes[0] << std::endl;
		//std::cout << "chrom.size() = " << chrom.size() << std::endl;
		std::cout << "arg_ptr->nsegs() = " << arg_ptr->nsegs() << std::endl;
		
		std::cout << "r = " << r << std::endl;
		std::cout << "f = " << f << std::endl;
		std::cout << "track_len = " << track_len << std::endl;
		std::cout << "rf = " << rf << std::endl;
		std::cout << "rft = " << rft << std::endl;
		
		{
			std::cout << "\nmigration matrix is" << std::endl;
			for(size_t pop = 0; pop < npop; ++pop) {
				for(size_t pop2 = 0; pop2 < npop; ++pop2){
					std::cout << migm[pop][pop2] << "\t" << std::endl;
				}
				std::cout << std::endl;
			}
			std::cout << std::endl;	
		}
	#endif
	
}


void HudsonMSArgFactory::clearWorkspace() const
{
	arg_ptr = boost::shared_ptr < hudson_ms::HudsonMSArg >(); // empty
	
	std::vector < int >().swap(config);
	std::vector < double >().swap(size);
	std::vector < double >().swap(alphag);
	std::vector < std::vector < double > >().swap(migm);
	
	npop = 0.0;
	
	t = 0;
	
	//nsegs = 0;
	//std::vector < int >().swap(nnodes);
	
	r = 0.0;
	pc = 0.0 ;
	lnpc = 0.0;
	rf = 0.0;
	rft = 0.0;
	prec = 0.0;
	cin = 0.0;
	clefta = 0.0;
	prect = 0.0;

	std::vector < double >().swap(tlast);

	nextEventItr = demogEventList.end();
	
	std::vector < chromo >().swap(chrom);
	
	//nlinks = 0;
}


// chromosone manipulation routines


/*  recombination subroutine 
   Picks a chromosome and splits it in two parts. If the x-over point
   is in a new spot, a new segment is added to seglst and a tree set up
   for it.   */
int HudsonMSArgFactory::recombination() const
{
	#ifdef MYDEBUG
		std::cout << "starting recombination" << std::endl;
	#endif
	
	/* First generate a random x-over spot, then locate it as to chrom and seg. */
	long int nlinks = chrom_links();
	int spot = nlinks*genRand() + 1.;
	//int this_beg = 0;
    /* get chromosome # (ic)  */
	size_t ic = 0;
	for(ic = 0; ic < chrom.size() ; ic++) {
		//this_beg = chrom.at(ic).segments.front().beg;
		//int el = (chrom[ic].segments.back()).end  -  this_beg;
		int el = chrom[ic].links();
		if( spot <= el ) break;
		spot -= el ;
	}
	//int is = this_beg + spot - 1;
	assert(ic < chrom.size());
	int is = chrom.at(ic).segments.front().beg + spot - 1;

	#ifdef MYDEBUG
		std::cout << "ic = " << ic << "\t\tic is " << is << std::endl;
	#endif
	
	xover(ic, is);

	return(ic);	
}



int HudsonMSArgFactory::cleftr() const
{
	#ifdef MYDEBUG
		std::cout << "starting cleftr" << std::endl;
	#endif

	double x = 0.0;
	while( (x = chrom_cleft() * genRand() )== 0.0 )  ;
		
	double sum = 0.0;
	int ic = -1 ;
	
	while ( sum < x ) {
		++ic;
		sum +=  1.0 - pow( pc, chrom[ic].links() )  ;
		
	}
	
	double links = static_cast<double>(chrom[ic].links());
	
	int is = chrom[ic].segments.front().beg 
		+ floor( 1.0 + log( 1.0 - (1.0- pow( pc, links))*genRand() )/lnpc  )
		- 1  ;
	
	#ifdef MYDEBUG
		std::cout << "ic = " << ic << "\t\tic is " << is << std::endl;
	#endif

	xover( ic, is);
	
	return( ic) ;
}


int HudsonMSArgFactory::cinr() const
{
	#ifdef MYDEBUG
		std::cout << "starting cinr" << std::endl;
	#endif

	/* First generate a random x-over spot, then locate it as to chrom and seg. */
	long int nlinks = chrom_links();
	int spot = nlinks*genRand() + 1.;

	//int this_beg = 0.0;
    /* get chromosome # (ic)  */

	//int i = 0;
	int ic = 0;
	int n = chrom.size();
	for(ic = 0; ic < n ; ++ic) {
		//i = ic;

		//this_beg = chrom[ic].segments.front().beg;
		//int el = (chrom[ic].segments.back()).end - this_beg;
		int el = chrom[ic].links();
		if( spot <= el ) break;
		spot -= el ;
	}
	
	assert(ic < n);
	//assert(ic == i);
	//int is = this_beg + spot - 1;
	int is = chrom[ic].segments.front().beg + spot - 1;
	
	//int endic = (chrom[i].segments.back()).end ;
	int endic = (chrom[ic].segments.back()).end ;

	#ifdef MYDEBUG
		std::cout << "ic = " << ic << "\t\tic is " << is << std::endl;
	#endif
	
	xover(ic, is);

	int len = floor( 1.0 + log( genRand() )/lnpc ) ;
	
	if( is+len >= endic ) {
		return(ic) ;
	}  
	
	if( is+len < (chrom.back().segments.front()).beg ) {
	   commonAncestors(ic, chrom.size()-1);
	    return(-1) ;
	}
	
	#ifdef MYDEBUG
		std::cout << "back in cinr, is + len " << is+len << std::endl;
	#endif

	
	xover( chrom.size()-1, is+len) ;
	
	commonAncestors( ic,  chrom.size()-1);
	
	return(ic);	

}



void HudsonMSArgFactory::xover(int ic, int is) const
{
	#ifdef MYDEBUG
		std::cout << "\nIn xover, ic = " << ic << "\t\tis is "<< is << std::endl;
	#endif
	
	#ifdef MYDEBUG1
		long int nlinks = chrom_links();
		double cleft = chrom_cleft();
	#endif

	#ifdef MYDEBUG2
		std::cout << "\nTotal chrom at start is" << std::endl;
		outputChrom(std::cout, ic);
	#endif
	//chromosone ic, spot measured from chromosone begin = is
	
	//cleft -= 1 - pow(pc, links_ic) ;
	#ifdef MYDEBUG1
		double links_ic = static_cast<double>(chrom[ic].links());
		cleft -= 1 - pow(pc, links_ic) ;
    #endif
	
   /* get seg # (jseg)  */

	int jseg = 0;
	//while (is < chrom[ic].segments[jseg].end) jseg++
	for(jseg = 0; is >= chrom[ic].segments[jseg].end ; ++jseg) ;
	
	/* eg if spot was exactly at end of first segment, jseg will be 1
	 * and in that case in will be 0, indicating break was right at start
	 * in which case we have a neat break between existing segments
	 * otherwise break in the middle of a segment and we will end up 
	 * with part of the broken one in one place and part in another. */
	int in = 0;
	if( is >= chrom[ic].segments[jseg].beg ) in = 1;
	
	int lsg_old = chrom[ic].nseg() ;
	int newsg = lsg_old - jseg;

	#ifdef MYDEBUG
		std::cout << "jseg = " << jseg << " << in = " << in << "newsg = " << newsg << std::endl;
	#endif
	#ifdef MYDEBUG2
		std::cout << "\nChrom to be changed (with is = " << is << ") is" << std::endl;
		outputChromo(std::cout, chrom[ic], jseg);
	#endif
	
   /* copy last part of chrom to a new chrom  */

	/*segments for new chromo element is going to be newsg long,
	 * from the end of the current segment back to the one with index jseg */
	
	chrom.push_back( chromo(chrom[ic].pop ) );
	//assert(chrom.size() == nchrom+1);
	
	#ifdef MYDEBUG
		std::cout << "added new chromosone, nchrom = " << chrom.size() << std::endl;
	#endif

	std::vector < seg >::iterator start = chrom[ic].segments.begin();
	std::advance(start, jseg);
	
	//now put in the chromosone segments
	chrom.back().segments.insert(chrom.back().segments.end(),
			start, chrom[ic].segments.end());
	
	// if in, first seg in new ones has beg is + 1
	// else first seg in new ones has beg =  old seg at jseg beg
	if (in) {
		chrom.back().segments.front().beg = is + 1;
	}
	else {
		assert (chrom.back().segments.front().beg == chrom[ic].segments[jseg].beg);
	}
	
	assert(chrom.back().nseg() == newsg);
	assert(chrom.back().segments.front().end == chrom[ic].segments[jseg].end);
	assert(chrom.back().segments.front().desc == chrom[ic].segments[jseg].desc);

	#ifdef MYDEBUG1
		std::cout << "and copied in segments, size for these segments is " << chrom.back().segments.size() << std::endl;
		std::cout << "and links for these segments = " << chrom.back().links() << std::endl;
	#endif
	
	//double links_new = static_cast<double>(chrom.back().links());
	//cleft += 1.0 - pow( pc, links_new) ; // for the new section section
	#ifdef MYDEBUG1
		double links_new = static_cast<double>(chrom.back().links());
		cleft += 1.0 - pow(pc, links_new) ;
    #endif
	
	
	// and we take out either all the segments which went to the new 
	// one, or all less 1 if the break was in the middle of that segment
	
	int lsg_new = jseg + in; 
	
	std::vector < seg >::iterator er = chrom[ic].segments.begin();
	std::advance(er, lsg_new);
	
	chrom[ic].segments.erase (er,chrom[ic].segments.end());
	
	assert(chrom[ic].nseg() == lsg_new);

	#ifdef MYDEBUG1
		std::cout << "erased the copied segments, size for this chromosones segments = " << chrom[ic].nseg() << std::endl;
		std::cout << "and links for these segments = " << chrom[ic].links() << std::endl;

	#endif
	
	// if in, seg at jseg has end is, 
	if( in ) {
		chrom[ic].segments.back().end = is;
	}
		
	// surely this is just 1?
	#ifdef MYDEBUG1
		// not this is NOT the links for a single chrom!
		std::cout << "nlinks =  " << nlinks << std::endl;
		std::cout << "taking off links for change = " << chrom.back().segments.front().beg - chrom[ic].segments.back().end << " will leave"<< std::endl;
		std::cout << "nlinks = " << nlinks - (chrom.back().segments.front().beg - chrom[ic].segments.back().end) << std::endl;
	
		nlinks -= chrom.back().segments.front().beg - chrom[ic].segments.back().end;
		assert(nlinks == chrom_links());
	#endif
	
	//double links_ic_now = static_cast<double>(chrom[ic].links());
	//cleft += 1.0 - pow( pc, links_ic_now) ; // for what's left in the old ic section
	#ifdef MYDEBUG1
		double links_ic_now = static_cast<double>(chrom[ic].links());
		cleft += 1.0 - pow( pc, links_ic_now) ; // for what's left in the old ic section
    #endif
	
	#ifdef MYDEBUG1
		assert ( abs(cleft - chrom_cleft()) < 0.0000000001 );
	#endif
	
	
	#ifdef MYDEBUG2
		std::cout << "\nTotal chrom after change is" << std::endl;
		outputChrom(std::cout);
	#endif
	
	#ifdef MYDEBUG2
		std::cout << "Whole arg_ptr is" << std::endl;
		arg_ptr->outputSegLst(std::cout);
	#endif
	 
	if( in ) {
		int begs = chrom.back().segments.front().beg;
		
		int i = 0;
		int k = 0;
		int ns = arg_ptr->nsegs();
		for( i=0;
			(k < ns - 1) && (begs > arg_ptr->at(arg_ptr->at(i).next).beg - 1);
			i = arg_ptr->at(i).next, k++ ) ;
		 
		#ifdef MYDEBUG
			std::cout << "checking for new tree: begs = " << begs << "\ti = " << i << "\targ->at(i).beg = " << arg_ptr->at(i).beg << std::endl;
		#endif
	 
		if( begs != arg_ptr->at(i).beg ) { /* new tree  */

			#ifdef MYDEBUG
				std::cout << "adding new tree" << std::endl;
			#endif
			
			//nnodes.push_back(nnodes.at(i));
			
			arg_ptr->add( arg_ptr->at(i) ); // add a copy of this segl
			// so the tree, whatever state it is now in, is copied with it
			assert(arg_ptr->nsegs() == ns+1);
			
			arg_ptr->at(ns).beg = begs;
			arg_ptr->at(ns).next = arg_ptr->at(i).next;
			arg_ptr->at(i).next =  ns;
			
			//nsegs++; // increments nsegs

			#ifdef MYDEBUG1
				//std::cout << "nsegs = " << nsegs << " and nnodes.size() = " << nnodes.size() << std::endl;
				std::cout << "nsegs = " << arg_ptr->nsegs() << std::endl;
			#endif
			#ifdef MYDEBUG2
				std::cout << "Whole arg_ptr at end of xover after adding segl is" << std::endl;
				arg_ptr->outputSegLst(std::cout);
			#endif

			
		}
	}
}

/* common ancestor subroutine 
   Pick two chromosomes and merge them. Update trees if necessary. **/

// return the decrease in the number of chromosones (1 or 2)
int HudsonMSArgFactory::commonAncestors(size_t c1, size_t c2) const
{

	#ifdef MYDEBUG
		std::cout << "\nIn commonAncestor, c1 = " << c1 << "\t\tc2 = "<< c2 << std::endl;
	#endif
	#ifdef MYDEBUG1
		long int nlinks = chrom_links();
		double cleft = chrom_cleft();
	#endif

	#ifdef MYDEBUG2
		std::cout << "\nTotal arg_ptr at start is" << std::endl;
		arg_ptr->outputSegLst(std::cout);
		std::cout << "\nTotal chrom at start is" << std::endl;
		outputChrom(std::cout, c1, c2);
	#endif


	int seg1=0;
	int seg2=0;

	//std::vector < seg > segments(nsegs, seg(0,0,0));
	// this is bigger than we might need
	std::vector < seg > segments;
	
	int tseg = -1 ;

	int k = 0;
	int ns = arg_ptr->nsegs();
	for(int segmnt = 0; k < ns; segmnt = arg_ptr->at(segmnt).next, k++) {
		
		int start = arg_ptr->at(segmnt).beg;
		int yes1 = isseg(start, c1, seg1);  // seg1 updated by this
		int yes2 = isseg(start, c2, seg2);  // seg2 updated by this
		
		#ifdef MYDEBUG
			std::cout << "k = " << k << "\t\tsegmnt = "<< segmnt << std::endl;
			std::cout << "seg1 = " << seg1 << "\t\tyes1 = "<<yes1 << std::endl;
			std::cout << "seg2 = " << seg2 << "\t\tyes2 = "<<yes2 << std::endl;
		#endif

		
		if( yes1 || yes2 ) {
			tseg++;
			#ifdef MYDEBUG
				std::cout << "tseg is now = " << tseg << std::endl;
			#endif
			
			int end = ( k < ns-1 ? arg_ptr->at(arg_ptr->at(segmnt).next).beg-1 : nsites-1 ) ;
			
			segments.push_back( seg(arg_ptr->at(segmnt).beg, end, 0) );
			// set beg and end
			//segments.at(tseg).beg = arg_ptr->at(segmnt).beg;
			//segments.at(tseg).end = end ;

			#ifdef MYDEBUG1
				std::cout << "setting end for new seg at tseg to " << end << std::endl;
			#endif
			#ifdef MYDEBUG2
				std::cout << "The new segments now are" << std::endl;
				outputSegments(std::cout, segments);
			#endif


			if( yes1 && yes2 ) {
				
				// add a node to the tree here
				(arg_ptr->at(segmnt)).tree.push_back(HudsonMSArg::node());
				
				//nnodes[segmnt]++;

				#ifdef MYDEBUG
					//std::cout << "two yes's: adding a node: nnodes[" << segmnt << "] is now " << nnodes[segmnt] << std::endl;
					std::cout << "two yes's: adding a node: arg_ptr->at(" << segmnt << ").tree.size() is now " << arg_ptr->at(segmnt).tree.size() << std::endl;
				#endif

				//(arg_ptr->at(segmnt)).tree.at(nnodes.at(segmnt)).time = t;
					(arg_ptr->at(segmnt)).tree.back().time = t;

				#ifdef MYDEBUG1
					std::cout << "------------- time for new node is " << t << std::endl;
				#endif
				#ifdef MYDEBUG2
					std::cout << "Segl at " << segmnt << " is now" << std::endl;
					arg_ptr->outputSegl(std::cout, segmnt);
				#endif
				
				//assert(arg_ptr->at(segmnt).tree.size() == nnodes[segmnt] + 1);
				
				if( arg_ptr->at(segmnt).tree.size() >= (2*nsam-1) ) {
					// take the last segment out
					std::vector<seg>::iterator del = segments.begin();
					std::advance(del, segments.size() - 1);
					segments.erase(del);
					assert(static_cast<int>(segments.size()) == tseg);

					
					tseg--; //mrca
					
					#ifdef MYDEBUG
						std::cout << "mrca, tseg now " << tseg << std::endl;
					#endif
					
				}
				else { 
					//segments.at(tseg)).desc = nnodes[segmnt];
					segments.at(tseg).desc = (arg_ptr->at(segmnt)).tree.size() - 1;
					#ifdef MYDEBUG2
						std::cout << "new segments.at(" << tseg << ")).desc = (arg_ptr->at(" << segmnt << ").tree.size() - 1 = " <<  (arg_ptr->at(segmnt)).tree.size() - 1 << std::endl;
						outputSegments(std::cout, segments);
					#endif

				}
				/*
				(segments.at(tseg)).desc = nnodes[segmnt];
				#ifdef MYDEBUG2
					std::cout << "new segments.at(" << tseg << ")).desc = nnodes[" << segmnt << "] = " <<  nnodes[segmnt] << std::endl;
					outputSegments(std::cout, segments);
				#endif
				*/
				int desc1 = (chrom[c1].segments[seg1]).desc;
				(arg_ptr->at(segmnt).tree.at(desc1)).abv = (arg_ptr->at(segmnt)).tree.size() - 1;
												
				int desc2 = (chrom[c2].segments[seg2]).desc;
				(arg_ptr->at(segmnt).tree.at(desc2)).abv = (arg_ptr->at(segmnt)).tree.size() - 1;
				
				// and set descendants, lowest index on left
				assert(desc1 != desc2);
				(arg_ptr->at(segmnt).tree.back()).ldes = (desc1 < desc2 ? desc1 : desc2);
				(arg_ptr->at(segmnt).tree.back()).rdes = (desc1 < desc2 ? desc2 : desc1);
				// accumulate descendants
				(arg_ptr->at(segmnt)).tree.back().ndes = 
						arg_ptr->at(segmnt).tree.at(desc1).ndes
						+ arg_ptr->at(segmnt).tree.at(desc2).ndes;
		
				
				#ifdef MYDEBUG2
					std::cout << "identified a desc1 for chromo at c1 = chrom[c1].segments[seg1]).desc = chrom[" << c1 << "].segments[" << seg1 << "]).desc = " << desc1 << std::endl;
					std::cout << "identified a desc2 for chromo at c2 = chrom[c2].segments[seg2]).desc = chrom[" << c2 << "].segments[" << seg2 << "]).desc = " << desc2 << std::endl;
					std::cout << "altered arg_ptr segl at " << segmnt << " so that node at desc1 = " << desc1 << " has abv = (arg_ptr->at(segmnt)).tree.size() - 1 = (arg_ptr->at(" << segmnt << ") - 1 = " << (arg_ptr->at(segmnt)).tree.size() - 1 << std::endl;
					std::cout << "altered arg_ptr segl at " << segmnt << " so that node at desc2 = " << desc2 << " has abv = (arg_ptr->at(segmnt)).tree.size() - 1 = (arg_ptr->at(" << segmnt << ") - 1 = " << (arg_ptr->at(segmnt)).tree.size() - 1 << std::endl;
					std::cout << "arg_ptr segl at " << segmnt << " now is" << std::endl;
					arg_ptr->outputSegl(std::cout, segmnt);
					std::cout << std::endl;
				#endif


			}
			else { // only one yes
				int desc = ( yes1 ? chrom[c1].segments[seg1].desc :
									chrom[c2].segments[seg2].desc );
				segments.at(tseg).desc = desc;
				#ifdef MYDEBUG2
					if ( yes1 ) {
						std::cout << "Only yes1: identified a desc1 for chromo at c1 = chrom[c1].segments[seg1]).desc = chrom[" << c1 << "].segments[" << seg1 << "]).desc = " << desc << std::endl;
						std::cout << "new segments.at(" << tseg << ")).desc = desc1 = " << desc << std::endl;
					}
					else {
						std::cout << "Only yes2: identified a desc2 for chromo at c2 = chrom[c2].segments[seg2]).desc = chrom[" << c2 << "].segments[" << seg2 << "]).desc = " << desc << std::endl;
						std::cout << "new segments.at(" << tseg << ")).desc = desc2 = " << desc << std::endl;
					}
					outputSegments(std::cout, segments);
					
				#endif
			}
		}
	}
	
	
	#ifdef MYDEBUG1
		std::cout << "nlinks =  " << nlinks << std::endl;
		std::cout << "taking off links for c1 = " << chrom[c1].links() << " will leave"<< std::endl;
		std::cout << "nlinks = " << nlinks - chrom[c1].links() << std::endl;
		nlinks -= chrom[c1].links();
	#endif
		
	//cleft -= 1.0 - pow( pc, static_cast<double>(chrom[c1].links()) );
	#ifdef MYDEBUG1
		cleft -= 1.0 - pow( pc, static_cast<double>(chrom[c1].links()) ); 
    #endif
	
	//clear chrom[c1] segments
	std::vector < seg > tmp;
	chrom[c1].segments.swap( tmp );
	
	if( tseg < 0 ) {
		
		chrom[c1] = chrom.back();
		//chrom[c1].segments = chrom[nchrom-1].segments;
		//chrom[c1].nseg = chrom[nchrom-1].nseg;
		//chrom[c1].pop = chrom[nchrom-1].pop ;
		
		#ifdef MYDEBUG
			std::cout << "tseg < 0 " << tseg << "\t nchrom cut down to " << chrom.size()-1 << " by taking out the last one " << std::endl;
		#endif
		#ifdef MYDEBUG2
			outputChromo(std::cout, chrom.back());
		#endif
		
		if( c2 == chrom.size() - 1 ) {
			c2 = c1;
			#ifdef MYDEBUG2
				std::cout << "c1 == nchrom - 1, c2 = c1 = " << c1 << std::endl;
			#endif
		}

		//erase the last one
		std::vector <chromo>::iterator del = chrom.begin();
		std::advance(del, chrom.size() - 1);
		chrom.erase(del);
		//assert(chrom.size() == nchrom - 1);
		
		#ifdef MYDEBUG2
			std::cout << "\nTotal arg_ptr at end of ca is" << std::endl;
			arg_ptr->outputSegLst(std::cout);
			std::cout << "\nTotal chrom at end of ca is" << std::endl;
			outputChrom(std::cout, c1, c2);
		#endif

	}
	else {
		
		#ifdef MYDEBUG1
			std::cout << "tseg >= 0 : tseg + 1 = " << tseg + 1 << std::endl;
			std::cout << "          : nsegs = " << arg_ptr->nsegs() << std::endl;
			std::cout << "    segments.size() = " << segments.size() << std::endl;
		#endif

		
		// cut segments down to size
		assert(static_cast<size_t>(tseg + 1) <= segments.size());
		
		std::vector <seg>::iterator stop = segments.begin();
		std::advance(stop, tseg + 1);
		std::vector <seg> tmp(segments.begin(), stop);
		assert(tmp.size() == static_cast<size_t>(tseg + 1));
		
		segments.swap(tmp);
		assert(static_cast<size_t>(tseg + 1) == segments.size());
		
		#ifdef MYDEBUG2
			std::cout << "tseg >= 0 : new segments size tseg + 1 = " << tseg + 1 << std::endl;
		#endif
		
		// and c1 gets the new stuff
		chrom[c1].segments = segments;
		assert(chrom[c1].nseg() == tseg + 1) ;
		#ifdef MYDEBUG1
			std::cout << "tseg >= 0 : tseg + 1 = " << tseg + 1 << std::endl;
			std::cout << "    chrom[c1].nseg() = " << chrom[c1].nseg() << std::endl;
			std::cout << "links(c1) =  " << chrom[c1].links() << std::endl;
			
		#endif

		
		#ifdef MYDEBUG1
			std::cout << "nlinks =  " << nlinks << std::endl;
			std::cout << "adding on links for c1 = " << chrom[c1].links() << " will give"<< std::endl;
			std::cout << "nlinks = " << nlinks + chrom[c1].links() << std::endl;
			nlinks +=chrom[c1].links();
		#endif
		
		
	   	//cleft += 1.0 - pow( pc, static_cast<double>(chrom[c1].links()) );
		#ifdef MYDEBUG1
			cleft += 1.0 - pow( pc, static_cast<double>(chrom[c1].links()) ); 
		#endif
	   	#ifdef MYDEBUG2
			std::cout << "chrom[c1] = chrom[" << c1 << "] gets the new stuff and becomes" << std::endl;
			outputChromo(std::cout, chrom[c1]);
		#endif
	}
	
	#ifdef MYDEBUG1
		std::cout << "nlinks =  " << nlinks << std::endl;
		std::cout << "taking off links for c2 = " << chrom[c2].links() << " will give"<< std::endl;
		std::cout << "nlinks = " << nlinks - chrom[c2].links() << std::endl;
		nlinks -= chrom[c2].links();
	#endif
	
	//cleft -= 1.0 - pow( pc, static_cast<double>(chrom[c2].links()) );
	#ifdef MYDEBUG1
		cleft -= 1.0 - pow( pc, static_cast<double>(chrom[c2].links()) ); 
    #endif
	
	// c2 gets the stuff now at the end
	//chrom[c2].segments = chrom[nchrom-1].segments;
	//chrom[c2].nseg = chrom[nchrom-1].nseg;
	//chrom[c2].pop = chrom[nchrom-1].pop ;
	chrom[c2] = chrom.back();

	#ifdef MYDEBUG2
		std::cout << "Stuff now at the end moved to chrom[c2] = chrom[" << c2 << "] which is now" << std::endl;
		outputChromo(std::cout, chrom[c2]);
		std::cout << "and nchrom is now " << chrom.size() -1 << std::endl;
	#endif
	
	//erase the [now] last one again
	std::vector <chromo>::iterator del = chrom.begin();
	std::advance(del, chrom.size() -1);
	chrom.erase(del);
	//assert(chrom.size() == nchrom - 1);
	
	//cleft = chrom_cleft();
	#ifdef MYDEBUG1
		assert(nlinks == chrom_links());
		assert ( abs(cleft - chrom_cleft()) < 0.0000000001 );
	#endif
	
	if( tseg < 0) return( 2 );  /* decrease of nchrom is two */
	else return( 1 ) ;
}

/* Isseg: Does chromosome c contain the segment on seglst which starts at
	    start? psg is the segment of chrom[c] at which one is to begin 
	    looking, and this is updated.  **/
bool HudsonMSArgFactory::isseg(int start, int c, int& psg) const
{
	int ns = chrom[c].nseg();
	
	/*  ms note changed order of test conditions in following line on 6 Dec 2004 */
	for(  ; (psg < ns ) && ( (chrom[c].segments[psg]).beg <= start ) ; ++psg )
		if( (chrom[c].segments[psg]).end >= start ) return true;
	
	return false;
}

long int HudsonMSArgFactory::chrom_links() const
{
	long int clinks = 0;
	for(std::vector < chromo >::const_iterator it = chrom.begin();
		it < chrom.end();
		++it) {
			clinks += it->links();
	}
	return clinks;
}


double HudsonMSArgFactory::chrom_cleft() const
{
	double ccleft = 0;
	for(std::vector < chromo >::const_iterator it = chrom.begin();
		it < chrom.end();
		++it) {
			ccleft += 1.0 - pow( pc, (double)(it->links()) );
	}
	return ccleft;
}

void HudsonMSArgFactory::pick2(int n, int& i,int& j) const
{
	i = n * genRand() ;
	while( ( j = n * genRand() ) == i )
		;
}

void HudsonMSArgFactory::pick2_chrom(size_t pop, int& pc1, int& pc2) const
{
	int c1 = 0;
	int c2 = 0;
	
	pick2(config[pop], c1, c2);
	
	int cs = (c1>c2) ? c2 : c1;
	int cb = (c1>c2) ? c1 : c2 ;
	
	int i=0;
	int count = 0;

	for(;;){
		while( chrom[i].pop != pop ) i++;
		if( count == cs ) break;
		count++;
		i++;
	}
	
	pc1 = i;
	
	i++;
	count++;
	
	for(;;){
		while( chrom[i].pop != pop ) i++;
		if( count == cb ) break;
		count++;
		i++;
	}
	pc2 = i ;
}
	

// chromo event definitions

HudsonMSArgFactory::HudsonMSChromoEvent_Xover::HudsonMSChromoEvent_Xover(double _time)
	: HudsonMSChromoEvent(_time) {}


void HudsonMSArgFactory::HudsonMSChromoEvent_Xover::callEvent(
				const hudson_ms::HudsonMSArgFactory& builder) const
{
	builder.processEvent_Xover();
}

HudsonMSArgFactory::HudsonMSChromoEvent_Mig::HudsonMSChromoEvent_Mig(double _time, double _mig)
	: HudsonMSChromoEvent(_time), mig(_mig) {}


void HudsonMSArgFactory::HudsonMSChromoEvent_Mig::callEvent(
				const hudson_ms::HudsonMSArgFactory& builder) const
{
	builder.processEvent_Mig(mig);
}

HudsonMSArgFactory::HudsonMSChromoEvent_Coal::HudsonMSChromoEvent_Coal(double _time, int _cpop)
	: HudsonMSChromoEvent(_time), cpop(_cpop) {}


void HudsonMSArgFactory::HudsonMSChromoEvent_Coal::callEvent(
				const hudson_ms::HudsonMSArgFactory& builder) const
{
	builder.processEvent_Coal(cpop);
}

std::ostream& HudsonMSArgFactory::outputChrom(std::ostream& os, int ic, int ic2) const
{
	if (!chrom.empty()) {
		
		int n = chrom.size();
		
		for(int k = 0; k < n; ++k) {
			
			if ((k == ic) || (k == ic2)) os << "**";
			else os << "  ";
			os << "{" << chrom[k].pop << " : " << chrom[k].nseg() << "}" << std::endl;
			for (std::vector < seg >::const_iterator it = chrom[k].segments.begin(); 
						it < chrom[k].segments.end(); 
						++it) {
				os << "|" << it->desc << ": " << it->beg << " - " << it->end << "|  ";
			}
			os << std::endl;
		}
		os << std::endl;
	}
	return os;
}

std::ostream& HudsonMSArgFactory::outputChromo(std::ostream& os, const chromo& chr, int jseg) const
{
	os << "{" << chr.pop << " : " << chr.nseg() << "}" << std::endl;
	int k = 0;
	for (std::vector < seg >::const_iterator it = chr.segments.begin(); 
				it < chr.segments.end(); 
				++it) {
		if (k == jseg) os << " *";
		else os << "  ";
		os << "|" << it->desc << ": " << it->beg << " - " << it->end << "|  ";
		k++;
	}
	os << std::endl;
	return os;
}

std::ostream& HudsonMSArgFactory::outputSegments(std::ostream& os, const std::vector < seg >& segments) const
{
	for (std::vector < seg >::const_iterator it = segments.begin(); 
				it < segments.end(); 
				++it) {
		os << "  |" << it->desc << ": " << it->beg << " - " << it->end << "|  ";
	}
	os << std::endl;
	return os;
}
