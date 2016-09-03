/* MCT - Markov Chains on Trees.

   Copyright (C) 2009-11 Brendan Bycroft <brb44@student.canterbury.ac.nz>
   Copyright (C) 2010-11 Jenny Harlow

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
\brief Doxygen documentation.
*/

//mct namespace documentation
/*! \namespace mct

\brief The namespace for MCT classes, structs, and subroutines.

*/

//mct namespace documentation
/*! \namespace aabc

\brief The namespace for classes, structs, and subroutines for
population Monte Carlo approximate Bayesian computation 
('adaptive abc').

*/


/** @mainpage MCT: Markov Chains on Trees


<ul>
<li>\ref main_sec_overview</li>
<li>\ref main_sec_key_terms</li>
<li>\ref main_sec_invocation</li>
<li>\subpage models</li>
<li>\subpage args</li>
<li>\subpage libsequence</li>
</ul>

<HR>

@section main_sec_overview Overview

This is a standalone program which can run a number of time-continuous 
Markov chain models over an ancestral recombination graph (ARG) and 
simulate a probabalistic genetic data at one locus for each 'present-day' 
individual in the ARG.

@section main_sec_key_terms Key terms

\b MCT simulates genetic data for a single \b locus on the chromosone 
for a sample of present-day individuals.  

The \b locus covers covers a number \b n_sites.  A \b site can be one
individual nucleotide site or the site for another unit 
of genetic data, such as a microsatellite.

The evolution of a 
the genetic data at a locus for a sample of present day individuals 
is described with an 
<b>ancestral recombination graph</b>.

An <b>ancestral recombination graph</b> (ARG) is a DAG 
(directed acyclic graph) whose nodes correspond to genotypes and 
edges correspond to either recombination or coalescence events 
<a href="http://biowiki.org/AncestralRecombinationGraph">(Biowiki)</a>.

The nodes in the ARG at time 0 are the 'present day' individuals or the
sample to which the ARG relates.  In this documentation the number of 
present day individuals is referred to as the sample size or \b nsam.

An ARG can be seen as a sequence of marginal trees, each tree relating 
to a portion of the locus on which no recombination has occurred 
(i.e. each marginal tree describes only coalescent events). 
See @subpage args for more information on ARGs and MCT. 

The points where recombination takes place divide the total 
locus into non-overlapping \b partitions,
each of which can be associated with a marginal tree.  


\b MCT can simulate an nsam-large sample of n_site-long genetic sequence
data
by applying a user-specified Markov chain \b model to each site. The 
models are independent of the partition boundaries,  
i.e. different models can be applied to different sites within a 
section of the locus covered by one marginal tree. See @subpage models 
for more information on Markov chain models and MCT.

A graphical description for a sample size of nsam=3, with 4 different 
partitions over 25 sites is shown below.  

Each partition of the locus is characterised by a marginal
tree, the site number beginning of subsection of DNA to which 
that tree applies, and the number of sites in that subsection.  

The nsam=3 sample nodes (present-day
individual) are denoted by s<sub>1</sub>, s<sub>2</sub>, s<sub>3</sub> 
and the nsam-12 internal nodes are denoted n<sub>a</sub> and n<sub>b</sub>.
(In the case of nsam=3 there are only two general tree shapes, depending on
which two sample nodes coalesce first, but each tree can 
have different branch lengths.)

\image html sites_and_trees_and_partions.png "" width=7cm
\image latex sites_and_trees_and_partions.png "" width=7cm

Each site x (i.e. a column of `-'s in the graphical depiction) is
associated with a Markov chain model m<sub>x</sub> that defines how the 
in a child node is simulated given the parent node and 
the time between the nodes (branch lengths in the tree).
In the diagram above, the models are applied
as follows:
<ul>
<li>sites 0-3: model m<sub>1</sub></li>
<li>sites 4-6: model m<sub>2</sub></li>
<li>sites 7-9: model m<sub>3</sub></li>
<li>sites 10-13: model m<sub>2</sub></li>
<li>sites 14-16: model m<sub>3</sub></li>
<li>sites 17-21: model m<sub>4</sub></li>
<li>sites 22-24: model m<sub>2</sub></li>
</ul>


The marginal tree and the model combine to give the simulated data for 
the site for each node (the present day individual sample nodes
and their ancestors represented by internal nodes) in the tree. 
A one-site
simulation for a nucleotide model is depicted below.

\image html one_site_simulation.png "" height=5cm
\image latex one_site_simulation.png "" height=5cm

For any node, the simulations for each site are combined in site order
to give the simulated sequence of genetic data for that node.  
The simulations of
interest are for the sample nodes or present day individuals (not the 
internal nodes).  

A possible simulated sample, using nucleotide models, for
the 3 individuals and 25 sites described above is:

\code
CTAGGGTCCCTAGGAGATTCCAAATCCGTC
ACAGGGTGTTATCGAGTTCGTGCAAGTCGC
ACTGGGCGTACTCGGGTTACTGCATCTCGC
\endcode 







 

@section main_sec_invocation Invocation

After compilation to the executable \b mct the program can be invoked as

./mct -m \e model_sites_file | -p [-a \e newick_file] [-n \e sample_size] 
[--seed \e prng_seed] [-q]
 
@subsection invocation_options Invocation options

<ul>
<li>-m \e model_sites_file  The file specifying the models to use and
how models map to sites.  Required option unless -p switch for 
model parameters only is used.  See \ref invocation_model_sites_file </li>
<li>-p    Switch to specify that only details of available models and
their default parameters should be printed to standard output.  
-m and -p should not be used in the same invocation.  All other options
will be ignored if this switch is used.</li>
<li>-n \e sample_size    The number of individuals in the sample.
Optional (defaults to 10).</li>
<li>-a \e newick_file    A file giving an ARG in newick format.  
Optional: If 
this option is used the ARG used by the program is read in from the
specified file. If not, the ARG is generated within the program
using \ref libsequence "libsequence".  See \ref invocation_newick_file </li> 
<li>--seed \e prng_see    Pseudo-random number generator seed.
Optional.</li>
<li>-q    Switch to supress detailed program output.</li>
</ul>

@subsection invocation_model_sites_file The model sites file

A file, for example a txt file, specifying the models to use and how
models map to sites. 

See @subpage models 
for more information on Markov chain models and MCT.

@subsection invocation_newick_file The newick tree file

A file specifying an ARG in newick format to be used in the program.
Each line specifies the marginal tree for n sites where n is given
by [n] at the beginning of the line.

See @subpage args for more information on ARGs and MCT. 

@subsection invocation_examples Invocation examples

<tt>./mct -m modelsites2.txt -n 15</tt>\n
Invoke \b mct to use model data in the file modelsites2.txt and
a sample of 15 individuals.  


<tt>./mct -m modelsites2.txt</tt>\n
Invoke \b mct to use model data in the file modelsites2.txt and
a sample of the default number of individuals (10).  


<tt>./mct -m modelsites2.txt -a newick.tree -n 15</tt>
\anchor invoke_newick_file \n
Invoke \b mct to use model data in the file modelsites2.txt, the
ARG specified in newick.tree, and a sample of 15 individuals.  


<tt>./mct -m modelsites2.txt -n 15 -q</tt>\n
Invoke \b mct to use model site in the file modelsites2.txt and
a sample of 15 individuals, and suppress detailed output.  


<tt>./mct -m modelsites2.txt -n 15 -seed 1234</tt>\n
Invoke \b mct to use model data in the file modelsites2.txt and
a sample of 15 individuals, and prng seed 1234.  


<tt>./mct -p</tt>\anchor invoke_params_only \n
Invoke \b mct to output only details of available models and 
their parameters.  

*/

//------------------------------------------------------------------------

/*!
@page models MCT, Markov chain models and model site files

A Markov chain model m<sub>x</sub> defines how the site-data in a child
node in a coalescence tree is simulated given the site-data in its parent 
node and the time between the nodes (branch lengths in the tree).

\anchor models_sites A site is defined with reference to the model 
applying to that site. In general, "site" refers to the unit relevant 
to the model. If a model is a nucleotide model a site refers to a 
single nucleotide in the DNA and the model is described as a nucleotide 
model.  A model could also be a microsatellite model.  

The state space of the model is related to the type of model:  A
nucleotide model has a state space of 4, being the 4 nucleotides A, G,
C and T.  A microsatellite model would have a larger state space.  
See the \ref models_statespace "model state space parameter" below.  

The models currently available can be printed to standard output
by invoking \b mct with the \b -p switch 
(see \ref invoke_params_only "invocation for model parameters").

An example of this output is shown below:

\code
Model parameters:

Model AQB: State Space (ss) = 40
	avgmr = 1
	d0 = 0.5
	d1 = 0.005
	d2 = 0.03
	threskappa = 0
	u0 = 0.5
	u1 = 0.005
	u2 = 0.03

Model EXP: State Space (ss) = 40
	alphaD = 0.302
	alphaU = 0.2
	avgmr = 1
	gammaD = 4e-07
	gammaU = 3.1e-06
	lambdaD = 1.06
	lambdaU = 1.06

Model JC69: State Space (ss) = 4

Model PL1: State Space (ss) = 40
	avgmr = 1
	bdryPHi = 0.5
	bdryPLo = 0.5
	m = 1
	p = 0
	s = 0.8752
	u = 0.6246
	v = 0.01542

Model PL2: State Space (ss) = 40
	avgmr = 0.761
	bdryPHi = 0.5
	bdryPLo = 0.5
	m = 0.5475
	p = 0
	s = 0.7638
	u = 0.8158
	v = 0.03947

Model REV: State Space (ss) = 40

Model SMM: State Space (ss) = 40

\endcode

The models to be used for each site of interest must be specified
in a modelsites file.  

A modelsites file is file, for example a txt file, 
specifying the models to use and how
models map to sites. Each line of the file gives a recognised model 
name and, in square brackets, the integer number of contiguous sites 
that will use that model, followed by other model-specific parameters. 
The following line gives the model to use for the following specified 
number of contiguous sites.

The "recognised model names" are those used in the ./mct -p output, 
for example "JC69".  If \a x contiguous sites use the JC69 model then one
line of the model sites file will be "JC69 [\a x]".  

Each model is related to a different set of parameters.  The parameters
for each model are shown in the ./mct -p output.  Default values
for each parameter are available, as shown, and will be used unless
an overriding value is given in the modelsites file.  For example, to 
override the \p avgmr parameter for the PL1 model and set it to, say,
1.2, include "avgmr=1.2" on the relevant line of the model sites file.

A special parameter, which applies to all models, is the 
\anchor models_statespace <b>state space</b> 
(\p ss). This is the number of possible different (mutually exclusive)
states for the site to which the model applies.  For example, for the
Jukes and Cantor 1969 (JC69) nucleotide model the possible states are the 
neucleotides A, G, C or T.  For the JC69 model the default of ss=4 cannot
in fact be overridden.  To override the state space for other models, 
say to use 30 rather than 40, include "ss=30" on the relevant line of 
the modelsites file (see also \ref models_sites "model sites" above).

The \anchor number_of_sites number of sites to to be used for the 
simulation is derived from
the modelsites file as the sum of the number of sites specified
on each line of the file. 

For example the file below, which might be called modelsites2.txt,
specifies models for 30 sites (9 + 3 + 5 + 3 + 10 = 30).  
   
\e modelsites2.txt
\code
JC69 [9]
SMM [3] ss=30
JC69 [5]
PL1 [3] ss=35 avgmr=1.2
JC69 [10]
\endcode

Malformed modelsites files, including unrecognised model names and/or 
parameters in the model sites file will cause
an exception to be thrown.  


*/

//------------------------------------------------------------------------

/*!
@page args MCT and ARGs

There are two ways of providing the  ancestral recombination graph
(ARG) used by MCT:

<ul>
<li>The graph to be used can be specified in an input file.  
See \ref args_sec_input_file </li>
<li>The graph can be simulated within the program.  
See \ref args_sec_simulated </li>
</ul>

Unless the -q switch is used on invocation, \b mct will print the ARG used
to standard output once it has been read from the file/simulated within
the program.
 
@section args_sec_input_file Specifiying the ARG in an input file

An input file specifying the ARG in newick format can be given using 
the \b -a switch
in the command line invocation of the program 
(see \ref invoke_newick_file "invoking with an arg file).  

(No formats other than newick are currently supported.)

By default each line of the file specifies the marginal tree for a 
single site unless the number of sites is specified as \a n sites 
where \a n is given by [\a n] at the beginning of the line.

For example, the line
\code
[1](2:1.059,((1:0.042,4:0.042):0.311,(3:0.058,5:0.058):0.296):0.706);
\endcode
explicitly states that this marginal tree applies to one site, and is
interpreted in the same way as the line
\code
(2:1.059,((1:0.042,4:0.042):0.311,(3:0.058,5:0.058):0.296):0.706);
\endcode
where the number of sites = 1 by default.

The total number of sites covered by the ARG is the sum of the number
 of sites specified on each line of the file.  This should match the 
 \ref number_of_sites "total number of sites derived from the model sites file".
 A mismatch will result in an exception being thrown.

For example the file below, which might be called tree.newick,
specifies 11 marginal trees covering 30 sites 
(3+ 4 + 2 + 2 + 1 + 1 + 1 + 5 + 2 + 7 + 2 = 30).  
   
\e tree.newick
\code
[3]((4:0.042,(1:0.033,2:0.033):0.009):0.713,(3:0.058,5:0.058):0.698);
[4]((4:0.042,(1:0.033,2:0.033):0.009):0.923,(3:0.058,5:0.058):0.907);
[2]((4:0.042,(1:0.033,2:0.033):0.009):0.169,(3:0.058,5:0.058):0.153);
[2](2:0.965,((1:0.042,4:0.042):0.169,(3:0.058,5:0.058):0.153):0.754);
(2:1.059,((1:0.042,4:0.042):0.169,(3:0.058,5:0.058):0.153):0.848);
(2:1.059,((1:0.042,4:0.042):0.313,(3:0.058,5:0.058):0.297):0.704);
[1](2:1.059,((1:0.042,4:0.042):0.311,(3:0.058,5:0.058):0.296):0.706);
[5]((1:0.042,4:0.042):1.017,(2:0.258,(3:0.058,5:0.058):0.200):0.802);
[2]((1:0.042,4:0.042):1.017,(2:0.258,(3:0.058,5:0.058):0.200):0.802);
[7]((1:0.042,4:0.042):0.727,(2:0.258,(3:0.058,5:0.058):0.200):0.511);
[2](1:0.769,(2:0.258,(3:0.058,(4:0.034,5:0.034):0.024):0.200):0.511);
\endcode

@section args_sec_simulated Simulating the ARG

If a newick-format ARG file is not specified, then the ARG will be
simulated using the 
<a href="http://molpopgen.org/software/libsequence.html">libsequence</a> 
library (see \ref libsequence).  

*/

//------------------------------------------------------------------------

/*!
@page libsequence The libsequence library

See <a href="http://molpopgen.org/software/libsequence.html">libsequence</a>
and <a href="http://www.molpopgen.org/software/libsequence/doc/html">libsequence documentation</a>

Libsequence is a C++ library for population genetics developed by Kevin Thornton.

Libsequence members used in MCT are denoted with the namespace Sequence:: 
and include:

<ul>
<li>\anchor libsequence_arg <b>Sequence::arg</b>\n
An ARG represented as an ordered list of Sequence::marginal marginal trees.</li>
<li>\anchor libsequence_marginal <b>Sequence::marginal</b>\n
The genealogy of a portion of a chromosome on which no recombination has occurred, 
represented as a coalescent tree.</li>
</ul>

*/




