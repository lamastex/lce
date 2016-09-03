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




#include "segsites_factory.hpp"


//#include <string>
#include <cmath>
//#include <numeric>
#include <algorithm>
//#include <iterator>
#include <stdexcept>
#include <cassert>


//#define MYDEBUG

//#define MYDEBUG1

//#define MYDEBUG2

//for making gametes
#define STATE1 1
#define STATE2 0


using namespace hudson_ms;



// constructor
HudsonMSSegSitesFactory::HudsonMSSegSitesFactory(
			boost::shared_ptr < hudson_ms::HudsonMSPRNGen > _rr)
			: rr(_rr)
{
	if (!(rr)) {
		throw std::invalid_argument(
		"HudsonMSSegSitesFactory::HudsonMSSegSitesFactory(...) : _rr empty");
	}
}


HudsonMSSegSitesFactory::~HudsonMSSegSitesFactory() {}


boost::shared_ptr < hudson_ms::HudsonMSSegSites > 
				HudsonMSSegSitesFactory::makeSegSites(
				const boost::shared_ptr < const hudson_ms::HudsonMSArg > _arg_ptr,
				const hudson_ms::HudsonMSSegSiteParams& ssParams) const
{
	try {
		initialise(_arg_ptr, ssParams);
	}
	catch (std::invalid_argument& ia) {
		throw std::invalid_argument(
		std::string("HudsonMSSegSitesFactory::makeSegSites(...) : ") + ia.what() );
	}
	
	#ifdef MYDEBUG
			std::cout << "\nin makeSS, nsam = " << nsam << ", segsites in = " << segsitesin;
			std::cout << ", nsites = " << nsites << std::endl;
			std::cout << "\nmfreq = " << mfreq << std::endl;
	#endif
	
	std::vector < double > treeTimes(nsegs);
	std::vector < int > segsites(nsegs, 0);
	
	double probss = 0.0;

	// now fill in treeTimes and segsites
	
	// if no segsites is specified, some mutation
	if ( (segsitesin == 0) && ( theta > 0.0) ) {
		
		allocateSegSitesRandomly(treeTimes, segsites);
		
	}
	
	// if some segsitesin specified
	// fill in probss as well
	else if( segsitesin > 0 ) {
		
		allocateFixedSegSites(treeTimes, segsites, probss);
	}
	
	// treeTimes and segsites, and probss if necessary, filled in
	// and actual positions and gamete states added to the *ss_ptr object
	
	// set the probss
	ss_ptr->setProbss(probss);
		
	// segsites should be built
	boost::shared_ptr < hudson_ms::HudsonMSSegSites > tmp_ss_ptr = ss_ptr;
	
	clearWorkspace(); // clear the workspace
		
	return tmp_ss_ptr;
}

// no segsites is specified, some mutation
void HudsonMSSegSitesFactory::allocateSegSitesRandomly(
				std::vector < double >& treeTimes,
				std::vector < int >& segsites) const
{	
	#ifdef MYDEBUG
		std::cout << "\nallocating sites randomly, segsitesin == 0 && theta > 0.0" << std::endl;
	#endif
	int segmnt = 0;
	for( int k = 0; k < nsegs; ++k) { 
		
		//if( mfreq > 1 ) ndes_setup( seglst[segmnt].ptree, nsam );
		
		//end = arg_ptr->at(segmnt).end;
		
		// tt should be tree time adjusted by len/nsites for this segment
		double treet = arg_ptr->treeTime(segmnt, mfreq);

		treeTimes[k] = treet;
					
		int len = arg_ptr->at(segmnt).len();
		double tt = (treet * len)/nsites;
		
		#ifdef MYDEBUG
			std::cout << "segmnt = " << segmnt << ", k = " << k << std::endl;
			std::cout << "tree time = "<< treet << std::endl;
		#endif
		
		int segsit = poisso( theta * tt );
		segsites[k] = segsit;
		
		#ifdef MYDEBUG
			std::cout << "back in allocateSegSitesRandomly, segsit = " << segsit << std::endl;
		#endif
		
		//add the segsites right away
		/* it would make more sense to do this after all the tree times calculations
		 * for both ways of allocating seg sites (randomly and fixed)
		 * but ms does it as k increments for randomly allocated segsites
		 * and so I do the same so that the random numbers
		 * generated for each part of the process are the same for me as ms*/
	
		addSegsites(treeTimes, segsites, k, segmnt);
					
		segmnt = arg_ptr->at(segmnt).next;
	}
}
	
// fixed segsites specified
void HudsonMSSegSitesFactory::allocateFixedSegSites(
				std::vector < double >& treeTimes,
				std::vector < int >& segsites, double& probss) const
{	
	std::vector < double > adjTreeTimes(nsegs);
		
	//pk = (double *)malloc((unsigned)(nsegs*sizeof(double)));
	//segsites = (int *)malloc((unsigned)(nsegs*sizeof(int)));
	//if( (pk==NULL) || (segsites==NULL) ) perror("malloc error. gensam.2");


	double tt = 0.0;
	
	int segmnt = 0;
	for(int k = 0; k < nsegs; k++) { 
	
		//if( mfreq > 1 ) ndes_setup( seglst[seg].ptree, nsam );
		
		double treet = arg_ptr->treeTime(segmnt, mfreq);
		
		int len = arg_ptr->at(segmnt).len();
		
		treeTimes[k] = treet;
		adjTreeTimes[k] = (treet * len)/nsites;
		
		tt += adjTreeTimes[k];
		
		segmnt = arg_ptr->at(segmnt).next;
	}

	
	if( theta > 0.0 ) { 
		double es = theta * tt;
		probss = exp( -es )*pow( es, (double) segsitesin) / segfac ;
		
	}

	if( tt > 0.0 ) {
		
		for (int k=0 ; k < nsegs ; k++) adjTreeTimes[k] /= tt ;
		
		// get some random allocation of segsitesin to segments
		segsites = allocatesSegSitesInOverTrees(adjTreeTimes);
	}
	// else segsites stays as a vector of nseg 0's 
	
	// treeTimes and segsites, and probss if necessary, filled in
	
	// now calculate actual positions and gamete states
	/* it would make more sense to do this after all the tree times calculations
	 * for both ways of allocating seg sites (randomly and fixed)
	 * but ms does it as k increments for randomly allocated segsites
	 * and so I do the same so that the random numbers
	 * generated for each part of the process are the same for me as ms*/
	segmnt = 0;
	for(int k = 0; k < nsegs; ++k) { 
		
		//add the segsites right away
		addSegsites(treeTimes, segsites, k, segmnt);
		
		segmnt = arg_ptr->at(segmnt).next; // move on
	}
  
}	
	
void HudsonMSSegSitesFactory::addSegsites(
				const std::vector < double >& treeTimes,
				const std::vector < int >& segsites,
				int k, int segmnt) const
{
	// now calculate actual positions and gamete states
	double nsinv = 1.0/nsites;
			
	if (segsites[k] > 0) {
		std::vector < HudsonMSArg::node > tree = arg_ptr->at(segmnt).tree;
	
		std::vector < std::vector < int > > gametes 
								= makeGametes(tree, treeTimes[k], segsites[k]);
	
		std::vector < double > newPos = 
			getNewPositions(segsites[k], 
							(arg_ptr->at(segmnt).beg)*nsinv, 
							(arg_ptr->at(segmnt).len())*nsinv);   

		//add to the segsites
		ss_ptr->addPositions(newPos, gametes);
	}
}


// a rand > 0.0
double HudsonMSSegSitesFactory::genRandGtZero() const
{
	return rr->genRandGtZero();
	
}

// a rand
double HudsonMSSegSitesFactory::genRand() const
{
	return rr->genRand();
	
}


void HudsonMSSegSitesFactory::initialise(
				const boost::shared_ptr < const hudson_ms::HudsonMSArg > _arg_ptr,
				const hudson_ms::HudsonMSSegSiteParams& ssParams) const
{
	
	if (!(_arg_ptr)) {
		throw std::invalid_argument("_arg_ptr empty");
	}
	
	nsegs = _arg_ptr->nsegs();
	nsam = _arg_ptr->nsams();
	nsites = _arg_ptr->nsites();
	
	#ifdef MYDEBUG
		std::cout << "nsites = " << nsites << std::endl;
	#endif
	
	if (nsegs == 0) {
		throw std::invalid_argument("nsegs in arg == 0");
	}
	if (nsam < 1) {
		throw std::invalid_argument("nsam in arg < 1");
	}
	if (nsites < 2) {
		throw std::invalid_argument("nsites in arg < 2");
	}
	
	mfreq = ssParams.getMfreq();
	//mfreq must be  <= nsam/2
	if (mfreq > nsam/2) {
		throw std::invalid_argument("ssParams.getMfreq() > nsam/2");
	}
	
	// no checks on the things ssParams is supposed to have checked
	theta = ssParams.getTheta();
	segsitesin = ssParams.getSegsitesin();
	
	arg_ptr = _arg_ptr;
	
	if(segsitesin > 0 && theta > 0.0 ) {
		segfac = 1.0 ; // calculate factorial
		for (int i = segsitesin; i > 1; --i) segfac *= i ;
	}
	
	iset = 0;
	gset = 0.0;
	
	/*
	std::vector < double > tmp1;
	tmp1.reserve(MAXSITES);
	tmp1.swap(positions);
	// positions is empty
	
	std::vector< int > tmp2;
	tmp2.reserve(MAXSITES);
	std::vector < std::vector < int > > tmp3(nsam, tmp2);
	tmp3.swap(list);
	// list contains nsam empty vectors
	*/
	
	bool showProbss= (segsitesin > 0 && theta > 0);
	
	ss_ptr = boost::shared_ptr < HudsonMSSegSites >( 
				new HudsonMSSegSites(nsam, showProbss) );
	
}

void HudsonMSSegSitesFactory::clearWorkspace() const
{
	arg_ptr = boost::shared_ptr < hudson_ms::HudsonMSArg >(); // empty
	ss_ptr = boost::shared_ptr < hudson_ms::HudsonMSSegSites >(); // empty
	
	nsegs = 0;
	nsam = 0;
	nsites = 0;
	
	theta = 0.0;
	segsitesin = 0;
	mfreq = 1;
		
	iset = 0;
	gset = 0.0;
	
	//std::vector < double >().swap(positions);
	//std::vector < std::vector< int > >().swap(list);
	segfac = 0;
}

int HudsonMSSegSitesFactory::poisso(double u) const
{
	int retvalue = 1;

	#ifdef MYDEBUG
		std::cout << "\n\tin poisso, u (time) = " << u << std::endl;
	#endif

	/*these horrible constants come from the ms code */
	if( u > 30. ) {
	    retvalue =  static_cast<int>( 0.5 + gasdev(u, u) ) ;
		
		#ifdef MYDEBUG
			std::cout << "\tu > 30, retvalue = (int)(0.5 + gasdev(u,u)) = " << retvalue << std::endl;
			if (retvalue < 0) std::cout << "\ti < 0, returning 0\n" << std::endl;
			else std::cout << "\ti >= 0, returning " << retvalue << std::endl;
		#endif
		
	    if( retvalue < 0 ) return( 0 ) ;
	    else return( retvalue ) ;
	}
	 
	double ru = genRand();
	double p = exp(-u);
	
	#ifdef MYDEBUG
		std::cout << "\tu <= 30, ru = " << ru << ", p = exp(-u) = " << p << std::endl;
		if (ru < p) std::cout << "\tru < p, returning 0" << std::endl;
	#endif
	
	if( ru < p) {
		return(0);
	}
	
	double cump = p;
	
	#ifdef MYDEBUG
		std::cout << "\tru >= p, retvalue = " << retvalue << ", cump = " << cump << std::endl;
	#endif
	
	while ( ru > ( cump += (p *= u/retvalue ) ) ) {
		retvalue++;
		#ifdef MYDEBUG
			std::cout << "\t\tretvalue = " << retvalue << ", cump = " << cump << std::endl;
		#endif
	}
	
	#ifdef MYDEBUG
		std::cout << "\treturning retvalue = " << retvalue << std::endl;
	#endif
	
	return(retvalue);
}

	
	
	

double HudsonMSSegSitesFactory::gasdev(double m, double v) const
{
	//int iset=0;
	//float gset;
	
	#ifdef MYDEBUG
		std::cout << "\n\t\tin gasdev, m = " << m << ", v = " << v << std::endl;
		std::cout << "\t\tiset = " << iset << ", gset = " << gset << std::endl;
	#endif
	
	if  (iset == 0) { // first time through, gset not set
		float r = 0.0;
		float v1 = 0.0;
		float v2 = 0.0;
			
		do {
			v1 = 2.0*genRand() - 1.0;
			v2 = 2.0*genRand() - 1.0;
			r = v1*v1 + v2*v2;
		} while (r >= 1.0);
		
		float fac = sqrt(-2.0*log(r)/r);
		
		gset = v1*fac;
		
		#ifdef MYDEBUG
			std::cout << "\t\troutine for iset = 0, v1 = " << v1 << ", v2 = " << v2 << std::endl;
			std::cout << "\t\tfac = " << fac << ", gset = " << gset << std::endl;
			std::cout << "\t\treturning m + sqrt(v)*v2*fac = " << m + sqrt(v)*v2*fac << std::endl;
		#endif
		
		iset = 1;
		
		return( m + sqrt(v)*v2*fac );
	} 
	else { // second time through, gset is set
		#ifdef MYDEBUG
			std::cout << "\t\troutine for iset = 1\n" << std::endl;
			std::cout << "\t\tgset = " << gset << std::endl;
			std::cout << "\t\treturning m + sqrt(v)*gset = " << m + sqrt(v)*gset << std::endl;
		#endif
		
		iset = 0; // allow gset to be reset next time
		
		return( m + sqrt(v)*gset ) ;
	}
}

//tt is tree time, unadjusted for len and nsites
//ns is at least as big as newsites
// post list has a new thingummmies in for each nsam
std::vector < std::vector < int > > HudsonMSSegSitesFactory::makeGametes(
			const std::vector < HudsonMSArg::node >& tree,
			double tt,
			int newsites) const 
{
	#ifdef MYDEBUG
		std::cout << "\nin makeGametes, newsites = " << newsites << ", tt = " << tt << std::endl;
	#endif

	std::vector < std::vector < int > > gametes(nsam);
	
	/* add newsites entries to each of the nsam segsite collections in list */
	for ( int j=0; j < newsites ;  ++j ) {
		
		#ifdef MYDEBUG
			std::cout << "\tj = " << j << ", about to pick node" << std::endl;
		#endif
		
		int node_picked = pickBranch(tree, tt);
		
		#ifdef MYDEBUG
			std::cout << "\n\tback in makeGametes with j = " << j << ", node = " << node_picked << std::endl;
		#endif
		
		
		for( int tip = 0; tip < nsam ; ++tip ) {
			int state = STATE2;
			if( isDescendant(tree, tip, node_picked) ) state = STATE1;
			gametes[tip].push_back(state); // add the new state
			
			#ifdef MYDEBUG
				std::cout << "\t\tgametes[" << tip << "][j] = " << gametes[tip][j] << std::endl;
			#endif 
			
		}
	}
	return gametes;
}

/*  pickBranch : returns a random branch from a tree. The probability of picking
              a particular branch is proportional to its duration. tt is total
	      time in tree.   */

int HudsonMSSegSitesFactory::pickBranch(
					const std::vector < HudsonMSArg::node >& tree,
					double tt) const
	//int nsam, mfreq;
	//struct node *ptree;
	//double tt;
{
	//double x, y, ran1();
	int lastbranch = 0 ;

	double x = genRand()*tt;
	
	#ifdef MYDEBUG
		std::cout << "\n\tin pickBranch, tt = " << tt << std::endl;
		std::cout << "\tx = " << x << std::endl;
	#endif

	
	double y = 0;
	
	for (int i=0; i < 2*nsam-2 ; i++ ) {
		int nd = tree.at(i).ndes;		
		if( ( nd >= mfreq ) && ( nd <= nsam - mfreq) ) {
			//y += (ptree + (ptree+i)->abv )->time - (ptree+i)->time ;
			y += tree.at(tree.at(i).abv).time - tree.at(i).time;
			lastbranch = i ;    // ms comment: changed 4 Feb 2010
			
			#ifdef MYDEBUG
				std::cout << "\t\ti = " << i << ", y = " << y << std::endl;
			#endif
		}
			
		if( y >= x ) {
			#ifdef MYDEBUG
				std::cout << "\ty >= x, returning i = " << i << std::endl;
			#endif
			return ( i );
		}
	}
	
	#ifdef MYDEBUG
		std::cout << "\treturning lastbranch = " << lastbranch << std::endl;
	#endif
	
	return( lastbranch );   // ms comment:  changed 4 Feb 2010
}

//  tdesn : returns true if tip is a descendant of node in tree, otherwise false.
bool HudsonMSSegSitesFactory::isDescendant(const std::vector < HudsonMSArg::node >& tree,
								int tip, int node_picked ) const
	
{
	/* work up from tip and either get to node_picked or overshoot*/
	int k = 0;
	for(k = tip ; k < node_picked ; k = tree.at(k).abv ) ;
	if( k == node_picked ) return true;
	else return false;
}

//locate(segsit,start*nsinv, len*nsinv,posit+ns);   
// ordran = ranvec and order
// fill in segsit new positions in the array with random numbers
//int ranvec(n,pbuf)
std::vector < double > HudsonMSSegSitesFactory::getNewPositions(
												int segsit, double beg, 
												double len) const
{
	#ifdef MYDEBUG
		std::cout << "\nin getNewPositions, segsit = " << segsit << ", beg = " << beg << ", len = " << len << std::endl;
	#endif
	
	// get the new positions as in interval (0, 1] and sort
	std::vector < double > newPos;
	
	#ifdef MYDEBUG
		std::cout << "\n\tallocating rands" << std::endl;
	#endif
	for(int i = 0; i < segsit; ++i) {
		
		newPos.push_back(genRand());
		
		#ifdef MYDEBUG
			std::cout << "\t(i = " << i << ", newPos[i] = " << newPos[i] << ")" << std::endl;
		#endif
	}
	#ifdef MYDEBUG
		std::cout << std::endl;
	#endif
	
	std::sort(newPos.begin(), newPos.end());
	
	#ifdef MYDEBUG
		std::cout << "\n\tafter ordering, new positions are:\n" << std::endl;
		int i;
		for(i=0; i<segsit; i++) {
			std::cout << "\t(i = " << i << ", newPos[i] = " << newPos[i] << ")" << std::endl;
			
		}
		std::cout << std::endl;
		int ii = 0;
	#endif
	//rescale
	for(std::vector < double >::iterator it = newPos.begin();
		it < newPos.end();
		++it) {
			
		#ifdef MYDEBUG
			std::cout << "\ti = " << ii << ", newPos[i] = " << newPos[ii] << ", newPos[i] = beg + ptr[i]*len = " << beg + (*it)*len << std::endl;
			ii++;
		#endif
		
		(*it) = ( beg + (*it)*len );
	}
	return newPos;
}

/*
	int
locate(n,beg,len,ptr)
	int n;
	double beg, len, *ptr;
{
	int i;

	ordran(n,ptr);
	for(i=0; i<n; i++)
		ptr[i] = beg + ptr[i]*len ;

}
*/

//mnmial(segsitesin,nsegs,pk,ss); ss = rv is an array, nsegs long
std::vector < int > HudsonMSSegSitesFactory::allocatesSegSitesInOverTrees(
					const std::vector < double >& adjTreeTimes) const
{
	
	std::vector < int > segsites(nsegs, 0);
	
	#ifdef MYDEBUG
		std::cout << "\n\tin allocatesSegSitesInOverTrees, segsites in = " << segsitesin << std::endl;
	#endif

	
	for(int i=0; i < segsitesin ; ++i) {
	   
		double x = genRand();
	   
		double tt = adjTreeTimes[0];
		
		int j=0; // incremented in loop
		
		#ifdef MYDEBUG
			std::cout << "\t\ti = " << i << std::endl;
			std::cout << "\t\tx = " << x << ", j = " << j << ", tt = " << tt << std::endl;
		#endif
	  
		while( (x > tt) && ( j < nsegs-1 ) ) {
			tt += adjTreeTimes[++j]; // accumulate times
			#ifdef MYDEBUG
				std::cout << "\t\t\tx > tt, tt now = " << tt << ", j = " << j << std::endl;
			#endif
		}
		
		segsites[j]++; // increment seg sites for this segment(tree)
		#ifdef MYDEBUG
			std::cout << "\t\trv[j] = segsites[" << j << "] = " << segsites[j] << std::endl;
		#endif
	}
	return segsites;
}



