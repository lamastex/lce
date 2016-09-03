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
*	of segment i, (2) ptree, which points to the first node of the
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
#ifndef __HUDSON_MS_ARG__
#define __HUDSON_MS_ARG__


#include <iostream>

#include <vector>
#include <string>


namespace hudson_ms {

	
	class HudsonMSArg {
		
		public:
		
			
		// unlike original hudson nodes, mine keep track of descendents
		struct node {
			node(); // no-args constructor
			int abv; // parent
			int rdes; // right descendant
			int ldes; // left descendant
			int ndes; // total number of descendants
			float time;
			//bool operator<(const HudsonMSArg::node& rhs) const;
		} ;
		
		typedef std::vector < node >::const_iterator const_tree_iterator;


		struct segl {
			segl(int nsam);
			int beg;
			int end;
			std::vector < HudsonMSArg::node > tree;
			int next;
			int len() const;
			HudsonMSArg::const_tree_iterator tree_begin() const;
			HudsonMSArg::const_tree_iterator tree_end() const;
			
		};
		
		typedef std::vector < segl >::const_iterator const_arg_iterator;
		
		
		/*! \brief Constructor. */
		HudsonMSArg(int _nsam);
		
		/*! \brief Destructor. */
		~HudsonMSArg() {};
		
		/*! \brief Get a const ref to an indexed element.*/
		const HudsonMSArg::segl& at(const size_t index) const;
		
		/*! \brief Get a ref to an indexed element.*/
		HudsonMSArg::segl& at(const size_t index);
		
		/*! \brief Add an element.*/
		void add(const HudsonMSArg::segl& argSeg);

		/*! \brief Get the number of elements.*/
		int nsegs() const;
		
		/*! \brief Get the number of individuals.*/
		int nsams() const;
		
		/*! \brief Get the number of sites.*/
		int nsites() const;
		
		/*! \brief Get the arg as a string.*/
		std::string toString() const;
		
		std::string toString(int prec, bool time) const;
		
		std::string toString(int prec, bool time, int precTime) const;
		
		std::string toString(bool time) const;
		
		std::string toString(bool time, int precTime) const;
		
		std::string toString(int prec) const;
		
		std::string toString(int mfreq, int prec) const;
		
		std::string toString(int mfreq, int prec, int precTime) const;
		
		/*! \brief Output the arg.*/
		std::ostream& outputArg(std::ostream& os) const;
								
		std::ostream& outputArg(std::ostream& os, 
								int prec, bool time) const;
		
		std::ostream& outputArg(std::ostream& os, 
								int prec, bool time, int precTime) const;
								
		std::ostream& outputArg(std::ostream& os, bool time) const;
		
		std::ostream& outputArg(std::ostream& os, bool time, 
												int precTime) const;
		
		std::ostream& outputArg(std::ostream& os, int prec) const;
		
		std::ostream& outputArg(std::ostream& os, int mfreq,
								int prec) const;
		
		std::ostream& outputArg(std::ostream& os, int mfreq,
								int prec, int precTime) const;
		
		const_arg_iterator begin() const;
		
		const_arg_iterator end() const;
		
		double argTimeAvgd(int mfreq) const;
		
		double treeTime(int arg_index, int mfreq) const;
		
		//debug only
		std::ostream& outputSegLst(std::ostream& os) const;
		//debug only
		std::ostream& outputSegl(std::ostream& os, int segmnt) const;
		
		protected:
		
		private:

		HudsonMSArg();
		
		std::string treeString(size_t segmnt, int prec) const;
		
		std::string toString(int mfreq, int prec, bool time) const;
		
		std::string toString(int mfreq, int prec, 
							bool time, int precTime) const;
		
		void parens(std::ostringstream& stm,
						const std::vector < node >& tree,
						//const std::vector < int >& descl,
						//const std::vector < int >& descr,
						int noden) const;
		
		int nsam;
		
		std::vector < HudsonMSArg::segl > seglst;
		/*
			struct segl {
			segl(int nsam);
			int beg;
			std::vector < node > tree;
			int next;
		*/
	};

	
}

#endif
