/*
LCE: A C++ class library for lumped coalescent experiments

Copyright (C) 2005--2011 Raazesh Sainudiin and Kevin Thornton
Copyright (C) 2009--2011 Jennifer Harlow and Brendan Bycroft

Copying Permission Statement:

LCE is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or (at
your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
See http://www.gnu.org/copyleft/gpl.html
*/

/*!
\mainpage

\image html "nCoalescentExperimetGraph.png" "A Kingman's n-coalescent experiments graph.  An observed multiple sequence alignment of the mother experiment and its offspring are shown on the left.  The corresponding formalities are shown on the right."

- \ref mainpage_sec_overview
- \ref mainpage_sec_credits

\section mainpage_sec_overview Overview

LCE is a C++ class library for lumped coalescent experiments that builds on libsequence (C++ class library for population genetics http://molpopgen.org/software/libsequence.html , boost (Boost provides free peer-reviewed portable C++ source libraries http://www.boost.org/) and GSL (Gnu Scientific Library http://www.gnu.org/software/gsl/).  
LCE, libsequence and GSL are distributed under the terms of the GNU General Public License (GPL) and boost is distributed under the boost license.

lce0.1 is under construction but full source docs are here. 

The goal of lce is to provide a code-base for inference in population geentics that takes advantage of a set of diverse tools, including extended arithmetics, computational algebraic algorithms and graph algorithms.

LCE library provides the following features:
- MCT: Markov Chains on Trees -- for simulating the evolution of repetitive and non-repetitive DNA on demographically and spatially structured ancestral recombination graphs and conducting Approximate Bayesian Computations.
- EBC: Approximate Bayesian Computation Done Exactly -- Generating realisations for lumped coalescent processes, Integrations over irreducible graphs via Markov bases, etc. 
- Generating structured ancestral recombination graphs (via libsequence and ms).

To learn how to use libsequence, boost and the GSL go to the appropriate help documents.

\section mainpage_sec_credits History and Credits


\version 0.1
\author 1. Raazesh Sainudiin (r.sainudiin ATADDRESSSYMBOL math.canterbury.ac.nz), 
Biomathematics Research Centre, 
Dept. of Maths and Stats, 
University of Canterbury, 
Private Bag 4800, 
Christchurch, New Zealand, 

\author 2. Kevin Thornton (krthornt ATADDRESSSYMBOL uci.edu) 
Department of Ecology and Evolution
UC Irvine
321 Steinhaus Hall
Irvine, CA 92697

\subsection mainpage_sec_subcredits Module Authors
- MCT: Brendan Bycroft, Jennifer Harlow and Raazesh Sainudiin
- EBC: Jennifer Harlow, Raazesh Sainudiin and Kevin Thornton 

The work on LCE started during Raazesh Sainudiin's and Kevin Thornton's post-doctoral work at Cornell.
Raazesh Sainudiin continued working on LCE during his subsequent post-doctoral work at The University of Oxford.
Many colleagues have contributed to the realization of LCE. Special thanks go to: 
- Ruriko Yoshida, Michael Stillman, James Booth, Peter Donnelly, Robert Griffiths and Gil McVean

For the latest news and up to date software contact http://www.math.canterbury.ac.nz/~r.sainudiin/codes/lce/index.shtml.

\date April 2011
*/
