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

#include "arg.hpp"

#include <iomanip>
#include <sstream>
#include <cmath>
#include <stdexcept>
#include <cassert>

//#define MYDEBUG

//#define MYDEBUG1

using namespace hudson_ms;

/*
bool HudsonMSArg::node::operator<(const HudsonMSArg::node& rhs) const
{
	return time < rhs.time;
}
*/

// set parent and descentants to -1, ndes to 1
HudsonMSArg::node::node()
	: abv(-1), rdes(-1), ldes(-1), ndes(1), time(0.0)
{}

// allocate space for nsam nodes in the tree
HudsonMSArg::segl::segl(int nsam)
	: beg(0), end(0), tree ( std::vector < HudsonMSArg::node > (nsam) ), next(-1)
{
	tree.reserve(2*nsam - 1);
}

int HudsonMSArg::segl::len() const
{
	return end - beg + 1;
}

HudsonMSArg::const_tree_iterator HudsonMSArg::segl::tree_begin() const
{
	return tree.begin();
}

HudsonMSArg::const_tree_iterator HudsonMSArg::segl::tree_end() const
{
	return tree.end();
}


HudsonMSArg::HudsonMSArg(int _nsam) : nsam(_nsam) 
{
	if (_nsam < 1) {
		throw std::invalid_argument(
			"HudsonMSArg::HudsonMSArg(...) : _nsam < 1");
	}
	seglst.push_back(segl(nsam));

}

const HudsonMSArg::segl& HudsonMSArg::at(const size_t index) const
{
	return seglst.at(index);
}

HudsonMSArg::segl& HudsonMSArg::at(const size_t index)
{
	return seglst.at(index);
}
		
void HudsonMSArg::add(const HudsonMSArg::segl& argSeg)
{
	seglst.push_back(argSeg);
}

int HudsonMSArg::nsegs() const
{
	return seglst.size();
}

int HudsonMSArg::nsams() const
{
	return nsam;
}

int HudsonMSArg::nsites() const
{
	int nsites = 0;
	if ( !seglst.empty()  ) {
		int i = 0;
		nsites = seglst.at(i).end + 1;
		while( seglst.at(i).next >= 0 ) { 
			i = seglst.at(i).next; 
			nsites = seglst.at(i).end + 1;
		}
	}
	return nsites;
}

std::ostream& HudsonMSArg::outputArg(std::ostream& os) const
{
	int mfreq = 1;
	int prec = 3;
	bool time = false;
	int precTime = 2*prec;
	os << toString(mfreq, prec, time, precTime);
	return os;
}

std::ostream& HudsonMSArg::outputArg(std::ostream& os, int prec, 
										bool time) const
{
	int mfreq = 1;
	int precTime = 2*prec;
	os << toString(mfreq, prec, time, precTime);
	return os;
}	

std::ostream& HudsonMSArg::outputArg(std::ostream& os, int prec, 
										bool time, int precTime) const
{
	int mfreq = 1;
	os << toString(mfreq, prec, time, precTime);
	return os;
}

std::ostream& HudsonMSArg::outputArg(std::ostream& os, bool time) const
{
	int mfreq = 1;
	int prec = 3;
	int precTime = 2*prec;
	os << toString(mfreq, prec, time, precTime);
	return os;
}

std::ostream& HudsonMSArg::outputArg(std::ostream& os, bool time, 
													int precTime) const
{
	int mfreq = 1;
	int prec = 3;
	os << toString(mfreq, prec, time, precTime);
	return os;
}

std::ostream& HudsonMSArg::outputArg(std::ostream& os, int prec) const
{
	int mfreq = 1;
	bool time = false;
	int precTime = 2*prec;
	os << toString(mfreq, prec, time, precTime);
	return os;
}

std::ostream& HudsonMSArg::outputArg(std::ostream& os, int mfreq, int prec) const
{
	bool time = true;
	int precTime = 2*prec;
	os << toString(mfreq, prec, time, precTime);
	return os;
}	

std::ostream& HudsonMSArg::outputArg(std::ostream& os, int mfreq, 
										int prec, int precTime) const
{
	bool time = true;
	os << toString(mfreq, prec, time, precTime);
	return os;
}

std::string HudsonMSArg::toString() const
{
	int mfreq = 1;
	int prec = 3;
	bool time = false;
	int precTime = 2*prec;
	return toString(mfreq, prec, time, precTime);
}

std::string HudsonMSArg::toString(int prec, bool time) const
{
	int mfreq = 1;
	int precTime = 2*prec;
	return toString(mfreq, prec, time, precTime);
}

std::string HudsonMSArg::toString(int prec, bool time, int precTime) const
{
	int mfreq = 1;
	return toString(mfreq, prec, time, precTime);
}

std::string HudsonMSArg::toString(bool time) const
{
	int mfreq = 1;
	int prec = 3;
	int precTime = 2*prec;
	return toString(mfreq, prec, time, precTime);
}

std::string HudsonMSArg::toString(bool time, int precTime) const
{
	int mfreq = 1;
	int prec = 3;
	return toString(mfreq, prec, time, precTime);
}

std::string HudsonMSArg::toString(int prec) const
{
	int mfreq = 1;
	bool time = false;
	int precTime = 2*prec;
	return toString(mfreq, prec, time, precTime);
}

std::string HudsonMSArg::toString(int mfreq, int prec) const
{
	bool time = true;
	int precTime = 2*prec;
	return toString(mfreq, prec, time, precTime);
}

std::string HudsonMSArg::toString(int mfreq, int prec, int precTime) const
{
	bool time = true;
	return toString(mfreq, prec, time, precTime);
}

std::string HudsonMSArg::toString(int mfreq, int prec, bool time, int precTime) const
{
	std::ostringstream oss;
	
	oss << std::fixed << std::setprecision(prec);
	
	if (!seglst.empty()) {
		
		int ns = nsegs();
		
		#ifdef MYDEBUG
			for(int k = 0; k < ns; ++k) {
				oss << "k = " << k << ", beg = " << at(k).beg 
				<< ", end = " << at(k).end << ", next = " 
				<< at(k).next << std::endl;
				
			}
		#endif

		double tmrca = 0.0;
		
		double ttot = argTimeAvgd(mfreq) ;
		
		int nsts = nsites();

		int segmnt = 0;
		for(int k = 0; k < ns; ++k) {
	
			#ifdef MYDEBUG
				oss << "\nNodes in tree are\n" << std::endl;
				for (std::vector < node >::const_iterator it = seglst[segmnt].tree.begin(); 
							it < seglst[segmnt].tree.end(); 
							++it) {
					oss << "node: abv = " << it->abv << "\ttime = " << it->time << std::endl;
				}
				oss << std::endl;
			#endif
			
			int end = seglst[segmnt].end;
			int start = seglst[segmnt].beg;
			
			if( (ns == 1) || ( ( start <= nsts/2) && ( end >= nsts/2 ) ) ) {
				tmrca = seglst[segmnt].tree[2*nsam-2].time ;
			}
			
			std::string stree = treeString( segmnt , prec) ;
			
			oss << stree << std::endl; 
			
			segmnt = seglst[segmnt].next; // go on to the next segment
		}
		
		oss << std::fixed << std::setprecision(precTime);
	
		if(time) oss << "time:\t" << tmrca << "\t" << ttot << std::endl;
	}
	
	return oss.str();
}


double HudsonMSArg::argTimeAvgd(int mfreq) const 
{
	double tt = 0.0;
	
	int ns = nsegs();
	int segmnt = 0;
	for( int k=0; k < ns; ++k) { 
		
		int len = seglst.at(segmnt).len();
		
		tt +=  (treeTime(segmnt, mfreq) * len);
		
		segmnt = seglst.at(segmnt).next;
	}
	return tt/nsites();
}

double HudsonMSArg::treeTime(int arg_index, int mfreq) const
{
	double t = 0;
	
	/* It should not matter whether we use mfreq method or not when 
	 * mfreq == 1 because the mfreq method with mfreq = 1 is effectively
	 * the same in that it should sum over all branches,
	 * but in practice there are differences in smaller decimal places in
	 * the results, so to ensure that the answers are the same 
	 * as used by ms, which does use the different methods, I
	 * have done the same here.*/
	
	if (mfreq < 2) {
		t = seglst.at(arg_index).tree.back().time ;
		for(int i = nsam; i < 2*nsam-1 ; ++i) {
			t += seglst.at(arg_index).tree.at(i).time ;
		}
	}
	
	else {
	
		for ( HudsonMSArg::const_tree_iterator it = seglst.at(arg_index).tree_begin();
			it < seglst.at(arg_index).tree_end();
			++it ) {
			// mfreq min is 1 so we won't try to use the mrca node (ndes = nsam)
			if( ( it->ndes >= mfreq )  && ( it->ndes <= nsam-mfreq) ) {
				//t += (ptree + (ptree+i)->abv )->time - (ptree+i)->time ;
				assert(it->abv >= 0);
				t += seglst.at(arg_index).tree.at(it->abv).time - it->time;
			}
		}
	}
	
	#ifdef MYDEBUG1
		if (mfreq == 1) {
			double tt = seglst.at(arg_index).tree.back().time ;
			for(int i = nsam; i < 2*nsam-1 ; ++i) {
				tt+= seglst.at(arg_index).tree.at(i).time ;
			}
			std::cout << std::setprecision(10);
			std::cout << "t using mfrq method = " << t << std::endl;
			std::cout << "t using no-mfreq method = " << tt << std::endl;
			assert(std::abs(tt - t) < 1.0E-8);
			cout.unsetf(ios::floatfield);
		}
	#endif
	
	return t;
}


HudsonMSArg::const_arg_iterator HudsonMSArg::begin() const
{
	return seglst.begin();
}

HudsonMSArg::const_arg_iterator HudsonMSArg::end() const
{
	return seglst.end();
}


std::string HudsonMSArg::treeString(size_t segmnt, int prec) const
{
	std::ostringstream stm;
	stm.precision(prec);
	stm << std::fixed;

	std::vector < node > tree = seglst[segmnt].tree;
	
	int len = seglst[segmnt].end -seglst[segmnt].beg + 1 ;

	stm << "[" << len << "]"; 
			
	parens( stm, tree, 2*nsam-2);
	
	return stm.str();
}

void HudsonMSArg::parens(std::ostringstream& stm,
	const std::vector < HudsonMSArg::node >& tree,
	int noden) const 
{
	
	if( tree[noden].ldes == -1 ) { // no descendents
		stm << noden + 1 << ":" << tree[tree[noden].abv].time;
	}
	else {
		stm << "(";
		parens( stm, tree, tree[noden].ldes ) ;
		stm << ",";
		parens(stm, tree, tree[noden].rdes ) ;
		
		if (tree[noden].abv == -1 ) { // root
			stm << ");";
		} 
		else {
		  double time = tree[tree[noden].abv].time - tree[noden].time ;
		  stm << "):" << time; 
		}
	}
}

std::ostream& HudsonMSArg::outputSegLst(std::ostream& os) const
{
	if (!seglst.empty()) {
		
		int ns = nsegs();
		
		for(int k = 0; k < ns; ++k) {
			
			os << "[" << seglst[k].beg << " - " << seglst[k].end << "] -> " << seglst[k].next << std::endl;
			for (std::vector < node >::const_iterator it = seglst[k].tree.begin(); 
						it < seglst[k].tree.end(); 
						++it) {
				os << "(" << it->time << ", -> " << it->abv << ")  ";
			}
			os << std::endl;
		}
		os << std::endl;
	}
	return os;
}

std::ostream& HudsonMSArg::outputSegl(std::ostream& os, int segmnt) const
{
	os << "[" << seglst[segmnt].beg << " - " << seglst[segmnt].end << "] -> " << seglst[segmnt].next << std::endl;
	for (std::vector < node >::const_iterator it = seglst[segmnt].tree.begin(); 
			it < seglst[segmnt].tree.end(); 
			++it) {
		os << "(" << it->time << ", -> " << it->abv << ")  ";
	}
	os << std::endl;
	return os;

}

