/* MCT - Markov Chains on Trees.

   Copyright (C) 2011 Jenny Harlow

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

/*! \file     
\brief LibSeq definitions.
*/


#include "libseq.hpp"

#include <Sequence/Coalescent/SimTypes.hpp>
#include <Sequence/Coalescent/Initialize.hpp>
#include <Sequence/Coalescent/DemographicModels.hpp>
#include <Sequence/Coalescent/Mutation.hpp>
#include <Sequence/stateCounter.hpp>

#include <gsl/gsl_randist.h>

#include <boost/bind.hpp>

#include <sstream>
//#include <exception>


using namespace mct;

			
Sequence::arg LibSeq::initialise_hist(
				const size_t nsam,
				const int nsites,
				const boost::shared_ptr < const mct::PRNGen >& r_ptr,
				const double g,
				const double t_begin,
				const double t_end,
				const double rho,
				const double size_at_end)
{
	/*
	single population, sample of $n_sample people, $nsites sites.
	n chromosomes, containing 1 segment.
	*/
	Sequence::marginal imarg = Sequence::init_marginal(nsam);
	
	/* define the sample 
	 * vector <int> (1, nsam) is a vector containing one element,
	 * nsam.
	 * In general, vector <int> (k, value) will give a vector of
	 * k repeats of value. Libseq is looking for a vector of
	 * population sizes, so for k populations, k elements and
	 * the i-th element (0 < i < k) is th i-th population size. 
	 */
	std::vector<Sequence::chromosome> sample =
			Sequence::init_sample(std::vector<int>(1,nsam),nsites);

	/*
	uni 	A binary function object (or equivalent) that returns a random deviate between a and b such that a <= x < b. a and b are the arguments to operator() of uni
	uni01 	A function object (or equivalent) whose operator() takes no arguments and returns a random deviate 0 <= x < 1.
	expo 	A unary function object whose operator() takes the mean of an exponential process as an argument and returns a deviate from an exponential distribution with that mean
	initialized_sample 	An initialized vector of chromosomes for a single population. 
	initialized_marginal 	An initialized marginal tree of the appropriate sample size for the simulation. 
	G 	The rate of exponential change in effective size. If G>0, the population grows exponentially (forwards in time). If G<0, it shrinks (again, forwards in time).
	t_begin 	The time in the past (in units of 4Ne generations) at which population size change begins (i.e., ends, moving forward in time)
	t_end 	The time in the past (in units of 4Ne generations) at which populations size change ends (begins forward in time)
	rho 	The population recombination rate 4N0r. 
	size_at_end 	At time t_end in the past, the population size is set to size_at_end. 
			If size_at_end = 1, the population is set to the same size that is was at t=0 (i.e. the beginning of the simulation).  
			If size_at_and < 0, the population size is not adjusted at t_end. In other words, it is left at whatever it grew or shrank to. 
	*/
	/* generate the ancestral recombination graph */
	return Sequence::exponential_change(
			boost::bind(&mct::PRNGen::uniform_double, r_ptr,_1,_2), // uniform_generator & uni,
			boost::bind(&mct::PRNGen::uniform01_double, r_ptr), // uniform01_generator & uni01,
			boost::bind(&mct::PRNGen::exponential_double, r_ptr,_1), //exponential_generator & expo,
			sample, //const std::vector<chromosome> & initialized_sample, 
			imarg, // const marginal & initialized_marginal,
			g, 
			t_begin, 
			t_end, 
			rho, 
			size_at_end); 
			
	
}


Sequence::SimData LibSeq::infinite_sites_sim_data(
				const int	nsites,
				const boost::shared_ptr < const mct::PRNGen >& r_ptr,
				const Sequence::arg& hist,
				const double& theta)
{
	//Apply mutations according to the infinitely-many sites scheme
	
	return Sequence::infinite_sites_sim_data(
				boost::bind(&mct::PRNGen::poisson_int,r_ptr,_1),
				boost::bind(&mct::PRNGen::uniform_double,r_ptr,_1,_2),
				nsites,
				hist,
				theta);
			
	
}
     
void LibSeq::startConversion (Sequence::PolyTable& poly,
			bool haveOutgroup, unsigned outgroup, bool strictInfSites,
			size_t this_nsites)
{
	unsigned segsites = poly.numsites();
	
	std::vector<double> newpositions;
	std::vector<unsigned> _pos_indexes;
	std::vector<std::string> newdata(poly.size());
	
	std::string _outgroup;
	if(haveOutgroup == true) {
		_outgroup = poly[outgroup];
	}
	else {
		_outgroup = poly[0]; // first one
	}
	
	for(unsigned j =0 ; j < segsites ; ++j) {
		Sequence::stateCounter counts;
		for(unsigned i = 0 ; i < poly.size() ; ++i) {
			counts (poly[i][j]);
		}
			
		if( (strictInfSites && counts.nStates()==2) || 
				(!strictInfSites && counts.nStates()>=2)) {
			// new position is going to be in interval [0,1]
			newpositions.push_back(poly.position(j)/this_nsites);
			_pos_indexes.push_back(j);
		}
	}
	
	for(unsigned i = 0 ; i < poly.size() ; ++i) {
		for(unsigned j = 0 ; j < _pos_indexes.size() ; ++j) {
			
			//skip sites where the ancestral information
			//is missing AND an outgroup sequence is present
			if(haveOutgroup==false ||
				(haveOutgroup==true && _outgroup[j] != 'N')) {
				
				if ( poly[i][_pos_indexes[j]] == 'N') {
					//missing data (leave untouched...)
					newdata[i] += 'N';
				}
				else if( poly[i][_pos_indexes[j]]
					!= _outgroup[_pos_indexes[j]] )	{
					//derived (maybe...)
					newdata[i] += '1';
				}
				else if ( poly[i][_pos_indexes[j]] 
					== _outgroup[_pos_indexes[j]] ) {
					//ancestral
					newdata[i] += '0';
				}
	
			}
		}
	}
	//assign takes care of setting non_const_access = true
	poly.assign(&newpositions[0],newpositions.size(),&newdata[0],newdata.size());
}




std::string LibSeq::nodeToString(
							const Sequence::node& snode)
{
	std::ostringstream stm;
	stm << "Time: " << snode.time;
	stm << ", abv: " << snode.abv;
	return stm.str();
	
}




std::string LibSeq::margToString(
							const Sequence::marginal& marg)
{
	std::ostringstream stm;
	stm << "beg: " << marg.beg;
	stm << ", nsam: " << marg.nsam;
	stm << ", nodes in tree: " << marg.tree.size();
	//note that nnodes is NOT a reliable source for nodes in the marg.tree
	//see libsequence documentation for Sequence::marginal::nnodes
	stm << ", nodes:\n";
	for (Sequence::marginal::const_iterator cit = marg.begin();
		cit < marg.end();
		++cit) {
		stm << nodeToString(*cit) << "\n";
	}	
	return stm.str();
	
}

std::string LibSeq::argToString(
							const Sequence::arg& arg)
{
	std::ostringstream stm;
	
	for (Sequence::arg::const_iterator cit = arg.begin();
		cit != arg.end();
		++cit) {
		stm << margToString(*cit) << "\n";
	}	
	return stm.str();
	
}

std::string LibSeq::polytableToString(const Sequence::PolyTable& polytable)
{
	std::ostringstream stm;
	polytable.print(stm);
	stm << "\n";
	return stm.str();
		
}

std::string LibSeq::argNewickString(
							const Sequence::arg& arg,
							const size_t prec)
{
	std::ostringstream stm;
	stm.precision(prec);
	stm <<   std::fixed;
		
	
	for (Sequence::arg::const_iterator cit = arg.begin();
		cit != arg.end();
		++cit) {
		makeNewick(cit, stm);	
		stm << "\n";
	}	
	return stm.str();
	
}

std::ostream & LibSeq::makeNewick(Sequence::arg::const_iterator m,
						std::ostream & o)
{
	Sequence::marginal::const_iterator mi = m->begin();
	int nsam = m->nsam;
	std::vector < int > left(2*nsam-1,-1);
	std::vector < int > right(left);
	std::vector< Sequence::node > tree;

	for( int i = 0; i < 2*nsam-2 ;++i ) {
		if( left[ (mi+i)->abv ] == -1 ) {
			left[(mi+i)->abv] = i;
		}
		else {
			right[ (mi+i)->abv] = i;
		}
	}
	
	parens( mi, left, right, 2*nsam-2,o);
	return o;
	
}


std::ostream & LibSeq::parens( 
							Sequence::marginal::const_iterator mi,
							const std::vector < int >& left,
							const std::vector < int >& right,
							const int & noden,
							std::ostream & o)
{
	double time;
 
	if( left[noden] == -1 ) {
		assert( (mi+((mi+noden)->abv))->time >= 0. );
		o << noden+1 << ':' <<(mi+((mi+noden)->abv))->time;
	}
	else {
		o << '(';
		parens( mi, left, right, left[noden], o ) ;
		o << ',';
		parens( mi, left, right, right[noden], o ) ;
		if( (mi+noden)->abv == -1 ) {
			o << ");";
		}
		else {
			time = (mi + (mi+noden)->abv )->time - (mi+noden)->time ;
			assert(time >= 0.);
			o << "):"<<time;
		}
	}
	return o;
}

// non_member function
bool mct::seq_node_comp_height(const Sequence::node& n1, const Sequence::node& n2)
{

	return n1.time < n2.time;
	
}

