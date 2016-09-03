# LCE: A C++ Class Library for Lumped Coalescent Experiments

### Authors: Jennifer Harlow, Raazesh Sainudiin and Kevin Thornton 
### CC 2005-2014

This library was written for the research needs of the following works:

* Experiments with the Site Frequency Spectrum, Raazesh Sainudiin, Kevin Thornton, Jennifer Harlow, James Booth, Michael Stillman, Ruriko Yoshida, Robert Griffiths, Gilean McVean and Peter Donnelly, [Bulletin of Mathematical Biology, Volume 73, Number 4, 829-872, 2011](http://www.springerlink.com/content/0748966716753484/).
* Finding the best resolution for the Kingman-Tajima coalescent: theory and applications, Raazesh Sainudiin, Tanja Stadler and Amandine Véber, Journal of Mathematical Biology, Volume 70, Issue 6, pp 1207-1247, 2015 ([preprint PDF 384KB](http://lamastex.org/preprints/SixCoalv4.pdf)). The final publication is available at Springer via [10.1007/s00285-014-0796-5](http://dx.doi.org/10.1007/s00285-014-0796-5).
* Jenny Harlow's MSc Thesis

## Overview
LCE is a C++ class library for lumped coalescent experiments that builds on libsequence (C++ class library for population genetics http://molpopgen.org/software/libsequence.html , boost (Boost provides free peer-reviewed portable C++ source libraries http://www.boost.org/) and GSL (Gnu Scientific Library http://www.gnu.org/software/gsl/). LCE, libsequence and GSL are distributed under the terms of the GNU General Public License (GPL) and boost is distributed under the boost license.

lce0.1 works for our research needs but is not autoconfiscated yet.

The goal of lce is to provide a code-base for inference in population geentics that takes advantage of a set of diverse tools, including extended arithmetics, computational algebraic algorithms and graph algorithms.

### LCE library provides the following features:

* MCT: Markov Chains on Trees -- for simulating the evolution of repetitive and non-repetitive DNA on demographically and spatially structured ancestral recombination graphs and conducting Approximate Bayesian Computations.
* EBC: Approximate Bayesian Computation Done Exactly -- Generating realisations for lumped coalescent processes, Integrations over irreducible graphs via Markov bases, etc.
* Generating structured ancestral recombination graphs (via libsequence and ms).

To learn how to use libsequence, boost and the GSL go to the appropriate help documents.

## History and Credits

Version: 0.1 


## Module Authors

* MCT: Brendan Bycroft, Jennifer Harlow and Raazesh Sainudiin
* EBC: Jennifer Harlow, Raazesh Sainudiin and Kevin Thornton

The work on LCE started during Raazesh Sainudiin's and Kevin Thornton's post-doctoral work at Cornell. Raazesh Sainudiin continued working on LCE during his subsequent post-doctoral work at The University of Oxford. It was significantly improved and extended by Jenny Harlow. Many other colleagues have contributed to the realization of LCE. Special thanks go to:

Ruriko Yoshida, Michael Stillman, James Booth, Peter Donnelly, Robert Griffiths and Gil McVean

This software is copy-lefted and distributed under the terms of the 
GNU General Public License (GPL).
see http://www.gnu.org/copyleft/gpl.html


## TO COMPILE

You are assumed to have the following GNU General Public Licensed 
or boost licensed libraries already installed
1. libsequence http://molpopgen.org/software/libsequence.html

2. GSL http://www.gnu.org/software/gsl/

3. boost http://www.boost.org/

Then you need to be brave, stubbrn and read the Makefiles in each directory to figure out the process.
