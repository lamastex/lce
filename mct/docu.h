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

\brief The namespace for \b MCT classes, structs, and subroutines.

*/

//mct namespace documentation
/*! \namespace aabc

\brief The namespace for classes, structs, and subroutines for
population Monte Carlo approximate Bayesian computation 
('adaptive abc').

*/

//reg namespace documentation
/*! \namespace reg

\brief The namespace for classes, structs, and subroutines related
to local-linear regression of simulated parameter values, following
Thornton, 'Automating approximate Bayesian computation by local linear
regression (BMC Genetics (2009), 10:35) (see also Beaumont, Zhang and 
Balding, 'Approximating Bayesian computation in population genetics',
Genetics (2002), 162, pp.2025-2035).

*/

//mct_utilities namespace documentation
/*! \namespace mct_utilities

\brief The namespace for some general utilities to help to use
the library, particularly outputting and manipulating collections
of results.

*/


/** @mainpage MCT: Markov Chains on Trees


<ul>
<li>\ref main_sec_overview</li>
<li>\ref main_sec_key_terms</li>
<li>\subpage page_mct</li>
<li>\subpage page_aabc</li>
<li>\subpage page_reg</li>
</ul>

<HR>

@section main_sec_overview Overview

This is a library of classes and methods related to the simulation and
analysis of genetic data including for statistical inference

The \b mct namespace includes types concerned with the representation
simulation and analysis of genetic data in nucleotide or microsatellite
form:
<ul>
<li>Types representing <b>summaries of genetic data</b> for multiple loci
such as mct::MultiLociPolyTable and mct::MultiLociMicrosat.</li>
<li>Types holding \b collections of genetic data types, such as 
mct::MultiLociPolyTableSet and mct::MultiLociMicrosatSet.</li>
<li>Types capable of \b simulating the genetic data types, such as 
mct::MultiLociPolyTableSampler and mct::MultiLociMicrosatSampler.</li>
<li>mct::MCTSampleSequence, representing <b>nucleotide and/or 
microsatellite-repeat encoded genetic data</b>.</li>
<li>SummaryStatistic, holding <b>summary statistics</b>
for objects of interest, including the summaries of genetic data.</li>
</ul> 

See \ref page_mct for more information.

The \b aabc namespace contains types for population 
Monte Carlo approximate Bayesian computation ('adaptive abc'), 
a technique for statistical inference about the 
values of parameters controlling population genetic processes.  The namespace
includes:
<ul>
<li>aabc::ParameterParticle, representing a set of parameter values.</li>
<li>aabc::ParameterPrior, representing a Bayesian prior for
a parameter.</li>
<li>The aabc::PopMCabcSimulatorAdaptive type for performing 
population Monte Carlo approximate Bayesian computation.</li>
</ul> 

See \ref page_aabc for more information.

The \b reg namespace contains types for local-linear regression 
of simulated parameter values, following
Thornton, 'Automating approximate Bayesian computation by local linear
regression (BMC Genetics (2009), 10:35) (see also Beaumont, Zhang and 
Balding, 'Approximating Bayesian computation in population genetics',
Genetics (2002), 162, pp.2025-2035).

See \ref page_reg for more information.

@section main_sec_key_terms Key terms


A \b locus is the location of a gene or DNA sequence on a chromosone.  

A \b locus covers covers a number of \b sites.  In this documentation 
the number of sites in a locus is referred to as \b nsites.  
A site can be one individual nucleotide site or the site for another unit 
of genetic data, such as a microsatellite.

Population genetic inference typically uses genetic data from a \b sample of 
present-day individuals to make inferences about the processes controlling
the evolution of the populations to which those individuals belong.
In this documentation the number of 
present day individuals is referred to as the sample size or \b nsam.

The evolution of a 
the genetic data at a locus for a sample of present day individuals 
is described with an 
<b>ancestral recombination graph</b>.

An <b>ancestral recombination graph</b> (ARG) is a DAG 
(directed acyclic graph) whose nodes correspond to genotypes and 
edges correspond to either recombination or coalescence events.

The nodes in the ARG at time 0 are the 'present day' individuals or the
sample to which the ARG relates.

See @subpage page_args for more details.

*/

//------------------------------------------------------------------------

/*!
@page page_mct Representing and sampling genetic data

<ul>
<li>\ref mct_sec_representing</li>
<li>\ref mct_sec_simulating</li>
<li>\ref mct_sec_analysing</li>
<li>\subpage page_mct_our_mct</li>
<li>\subpage libsequence</li>
</ul>

<HR>

@section mct_sec_representing Representing genetic data

@subsection mct_representing_subsection_encoding Polymorphic nucleotide site and microsatellite repeat encodings

The genetic data in a section of DNA can be represented (encoded)
as a sequence of \b nucleotides (eg "... AGCTAAGCCCTA ...").  For
some purposes a \b microsatelite encoding may be used.

\todo More about microsatellite encoding.

For population genetics, where we have a sample of nsam present-day
individuals, the focus of interest is usually in \e differences in
the sequences between those individuals.  Differences can be 
encoded as <b>nucleotide polymorphism</b> data.  In the case where there
can only be maximum of two different states for a particular site
in the sample ('single nucleotide polymorphisms, or SNPs), 
we can also encode this as binary (0/1) data.

\todo more on polymorphic site data

The \b MCT type mct::MultiLociPolyTable represents genetic data, encoded
as polymorphic site data, for some sample of present day individuals 
over multiple loci.
Each locus in the MultiLociPolyTable covers some number of sites, which
can differ between the loci represented in the MultiLociPolyTable.

Different sub-types of %mct::MultiLociPolyTable may exist, differing
in the source (underlying model) and form of the polymorphic data that 
they hold.  The \b MCT
types mct::MultiLociSimData
and mct::MultiLociPolySites are both sub-types of %mct::MultiLociPolyTable.



The \b MCT type mct::MultiLociMicrosat represents genetic data, encoded
as microsatellite repeats, for some sample of present day individuals 
over multiple loci. Each locus in the MultiLociMicrosat covers some 
number of sites, which can differ between the loci represented
in the MultiLociMicrosat.


Where the individuals in a sample come from different subpopulations,
the \b MCT types representing sample sequences (of polymorphic nucleotide
data or microsatellite data) maintain a relationship between the ordering of
the individual sequences for the sample and an ordering of the subpopulations.
i.e for \f$ m \f$ subpopulations,
with indices \f$ 0,\,1,\, \ldots \, m-1 \f$ and \f$ nsam^{(j)} \f$ samples from
subpopulation \f$ j \f$ , \f$ j = 0,\,1,\, \ldots \, m-1 \f$ , the first \f$ nsam^{(0)} \f$
individual sequences in the sample are those for the the \f$ nsam^{(0)} \f$ individuals from
the first subpopulation, next \f$ nsam^{(1)} \f$ sequences are those 
for the \f$ nsam^{(1)} \f$ individuals from the second subpopulation, etc. 


@subsection mct_representing_subsection_arg Ancestral recombination graphs

An ancestral recombination graph (ARG) represents the genealogical history
of a sample.  See \ref page_args for more details.  ARGs are relevant in 
\b MCT because methods for simulating sample sequences usually depend on
on information ('real', or also simulated) in this genealogical history  

The type mct::PopulationLabeledARG represents an ARG with subpopulation
labels, ie a representation of an ARG and information relating
to each present-day individual in the ARG to that individual's
subpopulation.  


@section mct_sec_simulating Simulating genetic data

\b MCT aims to be flexible in the way that it simulates genetic data. 
To avoid ending up with an rag-bag of unrelated types and processes, this
flexibility is built around a limited number of key types  common to 
the different methods:

<ul>
<li>A single type to describe a population structure;</li>
<li>A single mct::PopulationLabeledARG type representing an
ARG;</li>
<li>A limited number of generalised types for holding
and manipulating the results of different simulation methods.</li>
</ul>

@subsection mct_representing_subsection_population_structure Population structure

\b MCT simulates genetic data given a population structure.  The 
population structure describes the present 'shape' of the sample and the
relevant demographic history of the subpopulations from which that sample derives:

<ul>
<li>The number of subpopulations in the overall population;</li>
<li>The number of present-day individuals in the sample coming from each subpopulation;</li>
<li>The growth rate(s) assumed to apply to each subpopulation;</li>
<li>Migration between the subpopulations;</li>
<li>Other demographic events including changes in growth rates, migration etc, and 
subpopulation splits or merges;</li>
</ul>

At its simplest a population structure consists of at least one present-day individual
from a single population.   

\todo Add a link to mscplusplus

@subsection mct_simulating_subsection_args Simulating ancestral recombination graphs

In all but the very simplest cases, there are many possible genealogical 
histories that could have resulted in the given present-day sample structure.  The first
step in \b MCT's simulation is to identify an ARG to form the basis for further
simulation of mutations and hence obtain nucleotide or microsatellite sequence data.
 
The type mct::ArgFactory is a collection of different methods for making 
\link mct::PopulationLabeledARGs mct::PopulationLabeledARG\endlink.  These
include:
<ul>
<li>Simulating the %mct::PopulationLabeledARG using mscplusplus.  
This is the most flexible method of simulating a random ARG because 
it can deal with population structures containing more than one subpopulation
and can include past demograhic events in the simulation process.</li>
<li>Simulating the %mct::PopulationLabeledARG using ARG simulation
routines in \link libsequence libsequence\endlink library</li>.  This
is only suitable for simple population structures with only one subpopulation
and similarly can take into account only a limited past demographic history.
<li>Parsing an ARG description from string (see \ref page_args_section_tree_file_format).  
This method can be used to fix the ARG used for further simulations.</li>
</ul>

@subsection mct_simulating_subsection_seqs Simulating sample sequences

Given some (fixed or randomly simulated) %mct::PopulationLabeledARG
\b MCT uses various sampler types to simulate sets of the types representing
genetic data: 
<ul>
<li>mct::MultiLociSimDataSampler and mct::MultiLociPolySitesSampler
both create multiple replications of simulated polymorphic site data 
as a mct::MultiLociPolyTableSet.  The %mct::MultiLociPolyTableSet type 
is a container for a collection of mct::MultiLociPolyTable types.</li>
<li>mct::MultiLociMicrosatSampler creates multiple replications of 
simulated microsatellite repeats data as a mct::MultiLociMicrosatSet.
The %mct::MultiLociMicrosatSet</li> type is a container for a 
collection of \link mct::MultiLociMicrosats mct::MultiLociMicrosat\endlink.</li>
</ul>

\b MCT's approach to flexibility is illustrated in the ways in which 
polymorphic nucleotide site data sample sequences can be simulated.  There are 
two types of sampler (%mct::MultiLociSimDataSampler 
and %mct::MultiLociPolySitesSampler), both creating
mct::MultiLociPolyTableSets, the generalised container for polymorphism table data.

mct::MultiLociSimDataSampler uses an <b>infinite sites model</b> to simulate
SNP-encoded polymorphic nucleotide site data (the infinite site model used is 
that implemented in the \link libsequence libsequence\endlink
library).

The %mct::MultiLociPolySitesSampler (for nucleotide polymorphic data) and the 
%mct::MultiLociMicrosatSampler both use a <b>finite sites model</b> to simulate
polymorphic site data and microsatelite repeats data respectively.  
The finite site models implemented in  \b MCT are Markov
chain models. The Markov chain-based modelling methods within \b MCT can 
simulate an nsam-large sample of nsite-long 
genetic sequence
data (in nucleotide or microsatellite format) by applying a 
user-specified Markov chain \b model to each site in a locus.  
See \ref page_mct_our_mct for more information about \b MCT's Markov
chain models.  



@section mct_sec_analysing Summarising and analysing genetic data

Polymorphic nucleotide data, SNP data, and microsatellite repeats data are compressions 
of the information in DNA data.  Often it is necessary to further compress
this information.  

\todo. More on mct::SummaryStatistics and mct::DescriptiveStatistics  


*/

//------------------------------------------------------------------------

/*!

@page page_mct_our_mct MCT, Markov chain models and model site files

<li>\ref mct_our_mct_simulation</li>
<li>\ref mct_our_mct_MCTSampleSequence</li>
<li>\ref mct_our_mct_samplers</li>
<li>\subpage page_mct_our_mct_models</li>
<li>\subpage page_mct_our_mct_model_site_file</li>

<HR>

@section mct_our_mct_simulation Simulation with MCT

Given an ARG containing one or more marginal trees 
(see \ref page_args) representing the genealogical history
of a sample of size nsam for a locus covering nsite sites, 
MCT can simulate an nsam-large sample of nsite-long genetic sequence
data by applying a user-specified Markov chain \b model to each site. 

These Markov-chain models are all <b>finite site models.</b>

Different models can be applied to different sites within a 
section of the locus covered by one marginal tree. See \ref page_mct_our_mct_models
for more information on Markov chain models and MCT.

The mct::Model type represents a single \b MCT Markov chain model.  The
set of models for an entire locus is represented by a mct::ModelSet.

A ModelSet is created using an ordered collection of model data 
(see mct::ModelData) describing the models to be applied to the sites
on the locus.  Model data includes the specifc type of model 
(eg JC69), model parameters, and
the number of continguous sites within the locus covered by the model. 
The models are applied to the locus in the order in which they are
specified in the ordered collection of model data.  
The number of sites, say \f$ nsites_1 \f$, specified in the 
first ModelData in the model
collection \f$ m_1 \f$ tells the ModelSet that model \f$ m_1 \f$
should be applied to the first \f$ nsites_1 \f$ sites in the locus. 
The number of sites, say \f$ nsites_2 \f$, specified in the 
second ModelData in the model
collection \f$ m_2 \f$ tells the ModelSet that model \f$ m_2 \f$
should be applied to the next \f$ nsites_2 \f$ sites in the locus, etc.

Because specifying model configurations can be somewhat complicated
\b MCT includes a mct::ModelSetConfigBuilder type to build up the 
configuration for a ModelSet  within a program.  Using 
ModelSetConfigBuilder means that direct interaction with ModelData
can be avoided.  
	
The ModelSetConfigBuilder can either parse a string giving the model
configurations 
(see mct::ModelSetConfigBuilder::parseModelConfigsFromString and
\ref page_mct_our_mct_model_site_file) or
the configuration can be built up by progressively adding 
individual model specifications to the configuration (see
mct::ModelSetConfigBuilder::addModelConfig()).

The configuration can be retrieved from the ModelSetConfigBuilder
and passed direction to the ModelSet constructor (see
mct::ModelSetConfigBuilder::getModelSetConfigs()).

The code extract
below shows the use of a ModelSetConfigBuilder to specify
a very simple one-model configuration using only the JC69 model 
and the creation of the 
ModelSet using this configuration. 

See the documentation for 
mct::ModelSetConfigBuilder and mct::ModelSet for more details. 
	 

\code
	...
	// create a random number generator
	boost::shared_ptr < mct::PRNGen > r ( new mct::PRNGenGSL(seed) );
	
	// use a ModelSetConfigBuilder to build the model data
	mct::ModelSetConfigBuilder configBuilder;
	configBuilder.addModelConfig("JC69", nsites);
	
	// create our models
	mct::ModelSet models(configBuilder.getModelSetConfigs(), r);

\endcode


@section mct_our_mct_MCTSampleSequence Sample sequences from MCT Markov chain models

Sample sequences are simulated using the \b MCT Markov chain models in 
a mct::ModelSet using a mct::PartitionedARGSequenceGenerator 
(the 'partitioned arg' refers to partitioning of the locus
by recombination events and hence the need to take account of different
marginal trees in these different sections as well
as the different models that can be applied to each site as specified
in the %mct::ModelSet.  A %mct::PartitionedARGSequenceGenerator is 
configured with information about the ARG in the form of a 
mct::PopulationLabelledArg and can repeatedly create sample sequences
using this ARG and a given mct::ModelSet - 
see mct::PartitionedARGSequenceGenerator::createSampleSeqs().  Thus sequences
can be created using the same ARG but different configurations of models.

The sequences created by an %mct::PartitionedARGSequenceGenerator are in 
the form of the mct::MCTSampleSequence type.  This type holds the underlying
'raw' sequence data and - crucially - can decode it according to the
type of model that was used to generate it.  This is important because
the mct::PartitionedARGSequenceGenerator can deal with a model configuration
that includes either nucleotide or microsatellite models or mixture of
both nucleotide and microsatellite models.  An %mct::MCTSampleSequence 
provides a number of methods for accessing the sequence data it holds 
in different formats. 

The extract below illustrates using a model configuration that alternates
between microsatellite models (SMM, REV) and nucleotide models (JC69).  
The %mct::PartitionedARGSequenceGenerator partition_arg_seq_gen has 
already been created (configured with an ARG for a sample of nsam = 5
present-day individuals) and is used here to 
generate the sequences.  The sequences are then output in 'raw' format
(encoded by the models) and as a string of characters decoded from the raw
sequences according to the model type.    

\code

	...
	mct::ModelSetConfigBuilder configBuilder;
		
	configBuilder.addModelConfig("SMM", 2); // model_factory.cpp for options
	configBuilder.addModelConfig("JC69", 3);
	configBuilder.addModelConfig("REV", 2); // model_factory.cpp for options
	configBuilder.addModelConfig("JC69", 3);
	
	boost::shared_ptr < mct::ModelSet > models_ptr	( 
		new mct::ModelSet(configBuilder.getModelSetConfigs(), r_models) );
		
	boost::shared_ptr < mct::MCTSampleSequence > sample_seqs_ptr
					= partition_arg_seq_gen.createSampleSeqs(models_ptr);

	/* print the sequences of the sample in full form */
	std::cout << "\nsequence in raw format (with labels) is" << std::endl;
	std::cout << sample_seqs_ptr->toString << std::endl;

	/* print the sequences of the sample population in the phylip format */
	std::cout << "\nsequence in stringPhylipFormat() is" << std::endl;
	std::cout << sample_seqs_ptr->stringPhylipFormat();
	 ...

\endcode

The result might be something like

\code

sequence in raw format (with labels) is
0( 47   20      3       3       3       28      39      1       1       3 )
0( 47   19      0       2       3       28      11      2       1       0 )
0( 46   19      0       2       3       28      11      2       1       0 )
1( 46   19      0       2       3       28      43      1       1       0 )
1( 46   19      0       2       3       28      43      1       1       0 )


sequence in stringPhylipFormat() is
.47.20.TTT.28.39.CCT
.47.19.AGT.28.11.GCA
.46.19.AGT.28.11.GCA
.46.19.AGT.28.43.CCA
.46.19.AGT.28.43.CCA

\endcode

If only one type of model is used other sequence formats can be used. If
only nucleotide models are used then the sequence is all nucleotide data 
and we can ask for it to be converted into a 
libsequence::Sequence::PolySites object (see mct::snpTableFormat()).
Sample output showing three different
formats from one mct::MCTSampleSequence object created using just the JC69
nucleotide model is shown below (nsam = 5, nsites = 10)

\code

sequence in raw format (with labels) is
0( 3    3       0       0       0       3       2       1       2       2 )
0( 0    2       2       0       3       3       1       1       1       0 )
0( 0    1       2       0       3       3       1       1       1       0 )
1( 0    2       2       0       3       3       2       1       1       0 )
1( 0    2       2       0       3       3       0       1       1       0 )


sequence in stringPhylipFormat() is
TTAAATGCGG                         
AGGATTCCCA                         
ACGATTCCCA                         
AGGATTGCCA                         
AGGATTACCA                         

getting sequence in Sequence::PolySites format
polysites format polysites.print(...) gives                    
1       2       3       5       7       9       10             
T       T       A       A       G       G       G              
A       G       G       T       C       C       A              
A       C       G       T       C       C       A              
A       G       G       T       G       C       A              
A       G       G       T       A       C       A              

\endcode

Similarly, using just microsatellite models we can convert the sequence into
collections of microsatellite repeats 
(see mct::MCTSampleSequence::microsatRepeatsFormat()).

\code
sequence in raw format (with labels) is
0( 34   15      29      41      15      44      31      32      39      14 )
0( 34   15      28      42      15      44      16      32      39      10 )
0( 34   15      28      42      15      44      16      32      39      10 )
1( 34   15      29      41      16      44      28      32      39      10 )
1( 35   15      29      41      16      44      28      32      39      10 )


sequence in stringPhylipFormat() is
.34.15.29.41.15.44.31.32.39.14.    
.34.15.28.42.15.44.16.32.39.10.    
.34.15.28.42.15.44.16.32.39.10.    
.34.15.29.41.16.44.28.32.39.10.    
.35.15.29.41.16.44.28.32.39.10.    

getting sequence in microsat repeats form and printing the repeats
repeats format is                          
( 34    15      29      41      15      44      31      32      39      14 )

( 34    15      28      42      15      44      16      32      39      10 )

( 34    15      28      42      15      44      16      32      39      10 )

( 34    15      29      41      16      44      28      32      39      10 )

( 35    15      29      41      16      44      28      32      39      10 )

\endcode



@section mct_our_mct_samplers Samplers using MCT Markov chain models

The mct::MultiLociPolySitesSampler (for polymorphic nucleotide data) and the 
mct::MultiLociMicrosatSampler (for microsatellite repeats data) both 
use the \b MCT Markov
chain models for the sample sequence simulation, creating
\link mct::MultiLociPolyTableSets mct::MultiLociPolyTableSet\endlink
and \link mct::MultiLociMicrosatSets mct::MultiLociMicrosatSet\endlink
respectively.  The samplers deal with much of the complexity described above.
The intermediary step of using a %mct::PartitionedARGSequenceGenerator
to create a %mct::MCTSampleSequence and then converting the sequence data 
into the correct format and turning it into mct::MultiLociMicrosat or
mct::MultiLociPolyTable ojects for the set is all transparent to the user.

*/

//------------------------------------------------------------------------

//*** This section not used any more - excluded from markup recognised by doxygen ***
 

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



//------------------------------------------------------------------------

/*!
@page page_mct_our_mct_models MCT's Markov chain models

@section page_mct_our_mct_models_section_mc Markov chain models

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

@section page_mct_our_mct_models_section_models_available MCT models

The models currently available can be printed to standard output
by using mct::ModelFactory::printAvailableModels().

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

@section page_mct_our_mct_models_section_nuc_ex Example using nucleotide models

A graphical description for a sample size of nsam=3, with 4 different 
partitions over 25 sites is shown below.  

Each partition of the locus is characterised by a marginal
tree, the site number beginning the subsection of the locus to which 
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



\todo Microsatellite example


*/


//------------------------------------------------------------------------

/*!
@page page_mct_our_mct_model_site_file Model site files

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

Each model is related to a different set of parameters.  Default values
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
@page page_args Ancestral recombination graphs

ul>
<li>\ref page_args_section_newick_format</li>
<li>\ref page_args_section_population_structure</li>
<li>\ref page_args_section_tree_file_format</li>
</ul>

<HR>

The evolution of a 
the genetic data at a locus for a sample of present day individuals 
is described with an 
<b>ancestral recombination graph</b>.

An <b>ancestral recombination graph</b> (ARG) is a DAG 
(directed acyclic graph) whose nodes correspond to genotypes and 
edges correspond to either recombination or coalescence events 
<a href="http://biowiki.org/AncestralRecombinationGraph">(Biowiki)</a>.

The terminal nodes in the ARG at time 0 are the 'present day' 
individuals or the
sample to which the ARG relates.  For a sample of nsam present-day 
individuals there will be nsam nodes in the ARG at time 0.

If recombination has taken place within the locus, there will be different
genealogical histories for each section of the locus separated 
by the points at which recombination took place: The ARG can be seen 
as a sequence of marginal trees, each tree relating 
to a portion of the locus on which no recombination has occurred 
(i.e. each marginal tree describes only coalescent events). The points 
where recombination takes place divide the total 
locus into non-overlapping \b partitions,
each of which is associated with a marginal tree.  


@section page_args_section_newick_format Newick format

ARGs are often described or encoded in an adaption of newick tree format.
This shows each marginal tree in newick format, in the order in which
those trees apply to the locus.  Because it is also useful to know
the length of the \segment of the locus (number of sites within the locus)
to which each tree applies this
information may be prepended to the newick description of each tree.  In
the coding used in \b MCT the length of the segment is shown in square
brackets and the description for each marginal tree ends with a semicolon.
The following encoding would describe 
an ARG for nsam = 5 present-day individuals for locus covering nsites = 30
(3 + 4 + 2 + 2 + 1 + 1 + 1 + 5 + 2 + 7 + 2 = 30) 
in total in which 10 recombination events have taken place (making
11 segments each with its own marginal tree):

\anchor complex_tree
[3]((4:0.042,(1:0.033,2:0.033):0.009):0.713,(3:0.058,5:0.058):0.698);

[4]((4:0.042,(1:0.033,2:0.033):0.009):0.923,(3:0.058,5:0.058):0.907);

[2]((4:0.042,(1:0.033,2:0.033):0.009):0.169,(3:0.058,5:0.058):0.153);

[2](2:0.965,((1:0.042,4:0.042):0.169,(3:0.058,5:0.058):0.153):0.754);

[1](2:1.059,((1:0.042,4:0.042):0.169,(3:0.058,5:0.058):0.153):0.848);

[1](2:1.059,((1:0.042,4:0.042):0.313,(3:0.058,5:0.058):0.297):0.704);

[1](2:1.059,((1:0.042,4:0.042):0.311,(3:0.058,5:0.058):0.296):0.706);

[5]((1:0.042,4:0.042):1.017,(2:0.258,(3:0.058,5:0.058):0.200):0.802);

[2]((1:0.042,4:0.042):1.017,(2:0.258,(3:0.058,5:0.058):0.200):0.802);

[7]((1:0.042,4:0.042):0.727,(2:0.258,(3:0.058,5:0.058):0.200):0.511);

[2](1:0.769,(2:0.258,(3:0.058,(4:0.034,5:0.034):0.024):0.200):0.511);

In this encoding only the nodes corresponding to the root nodes - the
five present-day individuals are labelled.  The labels are the numbers
(1, 2, 3, 4, 5) which are followed by a colon and then the time in that
tree branch.  The internal nodes of the tree are not labelled.

@section page_args_section_population_structure Spatial population structure

The present-day individuals represented by the nodes of an ARG may come 
from different subpopulations where there has been sufficient inter-subpopulation
'communication' to allow there to be a single common ancestor (root node) in
in each marginal tree).  Inter-subpopulation communication can include
migration between subpopulations, or subpopulations merging or splitting at
some point in their genealogical history.

To related individuals/nodes in the tree to subpopulations, it is useful
to consider the subpopulations as an ordered collection, each one with
a unique index in the collection and to relate the nodes at time 0 to
the subpopulations in the same order.  i.e for \f$ m \f$ subpopulations,
with indices \f$ 0,\,1,\, \ldots \, m-1 \f$ and \f$ nsam^{(j)} \f$ samples from
subpopulation \f$ j \f$ , \f$ j = 0,\,1,\, \ldots \, m-1 \f$ , the nodes
numbered \f$ 1,\ldots\, nsam^{(1)} \f$ are the \f$ nsam^{(0)} \f$ individuals from
subpopulation index 0, the nodes numbered 
\f$ nsam^{(0)} + 1,\ldots\, nsam^{(0)} + nsam^{(1)} \f$ 
are the \f$ nsam^{(1)} \f$ individuals from subpopulation index 1 etc. 

For example, if the ARG above related to two subpopulations with three present-day
individuals from the first subpopulation and two from the second
subpopulation, the nodes numbered 1, 2, 3 would represent the 
first subpopulation's individuals in the ARG and the nodes numbered 4, 5 
would represent the second subpopulation's individuals in the ARG.


@section page_args_section_tree_file_format Newick tree file format

\b MCT can make an ARG by reading a file specifying the ARG in 
the adapted newick format described above (no formats other than 
newick are currently supported).  
See mct::ARGFactory::makeHistFromNewickString()

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
 of sites specified on each line of the file.  

For example the file below, which might be called tree.newick,
specifies the same ARG of 11 marginal trees covering 30 sites shown 
\link complex_tree above\endlink. 
   
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

*/

//------------------------------------------------------------------------

//*** This section not used any more - excluded from markup recognised by doxygen ***

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


//------------------------------------------------------------------------

/*!
@page page_aabc Adaptive approximate Bayesian computation

The \c aabc namespace contains types for population 
Monte Carlo approximate Bayesian computation ('adaptive abc' or AABC), 
a technique for statistical inference about the 
values of parameters controlling population genetic processes.  

\todo More on aabc

In AABC, the focus of interest is not the simulation of DNA data itself
but on approximating a posterior distribution for parameters relevant
to the evolution of the population given the observed data.

\b MCT provides an aabc::ParameterParticle type to hold  
particular values of for a parameter vector.  The %aabc::ParameterParticle
can hold values for one or multiple values and is designed to 
be as general as possible:  subject to the restriction that the 
values are held as doubles, the %aabc::ParameterParticle can represent
values for parameters for many different kinds of processes. 

aabc::ParameterPrior is abstract type for specifying a Bayesian
prior for a single element of a vector of parameter. 
Two concrete implementations
are provided:  aabc::ParameterPriorUniform (a uniform prior, configured
with the upper and lower bounds for the required uniform distribution)
and aabc::ParameterPriorNormal (a normal prior, configured with the required
mean and standard deviation).  A prior for a parameter vector is specified
with a aabc::ParameterPriorSet, an ordered collection of aabc::ParameterPriors.

\b MCT provides a very generalised population 
Monte Carlo approximate Bayesian computation simulator, 
aabc::PopMCabcSimulatorAdaptive.  This is designed to be as
independent as possible of the type of process of interest, the 
parameterisation of that process, and and the type of observed data 
on which the posterior distribution of interest is conditioned.  

This independence is achieved by configuring the aabc::PopMCabcSimulatorAdaptive
with a aabc::ParameterPrior and a aabc::ParameterParticleValidator.  
The prior is responsible for generating possible parameter particles
and the %aabc::ParameterParticleValidator is responsible for quantifying how
'good' they are, with reference to the observed data and returning
the results of the validation in non-process-specific general format 
(aabc::ParameterParticleValidationResult).  The validator 
'decouples' the %aabc::PopMCabcSimulatorAdaptive from the details of 
the process of interest, what roles the elements of a %aabc::ParameterParticle
 fulfil in that process,
 and how the 'goodness' of a %aabc::ParameterParticle
is assessed.  

The aabc::ParameterParticleValidator is an abstract type which must simply
be capable of providing validation results for a 
%aabc::ParameterParticle or a set of aabc::ParameterParticles
(aabc::ParameterParticleSet).  The aabc::ParameterParticleValidationResult type
is a simple data structure that the aabc::PopMCabcSimulatorAdaptive can use
without any knowledge of the validation process. 

Specific aabc::ParameterParticleValidator subtypes provide the link between
the %aabc::PopMCabcSimulatorAdaptive and specific processes and types of observed data.
In the case of inferences about parameters for a population genetics process, for
example, an aabc::ParameterParticleValidatorMultiLociPolyTable is configured with 
a mct::MultiLociPolyTable reference object and a mct::MultiLociPolyTableSample.
It uses the sampler to simulate data with each particle it is given to validate
and calculates a distance between the simulated data and the reference object
using the summary statistics from each.  The distances it returns to 
the %aabc::PopMCabcSimulatorAdaptive as part of the 
aabc::ParameterParticleValidationResult are all that that object needs to
decide whether to 
retain a %aabc::ParameterParticle in the population and what weight to give it. 


\internal
(In fact, the validation does not necessarily have to have any reference to an
observed data object, ie the approximation sought need not be to
to a posterior distribution.)


*/

//------------------------------------------------------------------------


/*!
@page page_reg Local linear regression

Types for performing local linear regression on sets of
aabc::ParameterParticles are in the \c reg namespace.  
See Thornton, 'Automating approximate Bayesian computation by local linear
regression', BMC Genetics (2009), 10:35 and also see also Beaumont, Zhang and 
Balding, 'Approximating Bayesian computation in population genetics',
Genetics (2002), 162, pp.2025-2035.  

These types are closely based on code supplied by Kevin Thornton.

These types have not been extensively tested or used.  Caveat emptor!

\section reg_section_comparison Comparison of results to Thornton

Thornton's code takes input from files whereas the \b MCT types can
take input directly from other \b MCT objects.  There will therefore
be some small differences in results even if the local linear
regression methods are equivalent, because of the rounding of values 
in the process of being written to file which will be experienced 
when using Thornton's code.  

\internal
\subsection reg_subsection_tan The tan transformation reg::TanTransformer

Testing so far has given the same results as Kevin's, subject to
the comment about rounding, when no transformation or the log transformation
are used.  The implementation of the tan transformation  (see
Hamilton, G., Stoneking, M., Excoffier., L. (2005), "Molecular
 analysis reveals tighter social regulation of immigration in
 patrilocal populations than in matrilocal populations."
 PNAS, 102, pp. 746-7480) is however a little different to Kevin's.  He
 uses mins and maxs from \b data referred to as 'prior', ie
 the parameter values to be transformed before being regressed (which means that
 the transformation of an individual value depends on the range of data 
 in the set to which that value belongs).  This is not what 
 Hamilton et al use - they use (or say they do) the 
 actual min and max bounds <b>on the prior as the distribution
 from which possible parameter values are generated </b>.  Since
 the aim of the whole tan transformation is to ensure that the 
 posterior values remain within the bounds of the \b prior I have followed
 what I understand to be Hamilton et al's approach.

*/

//------------------------------------------------------------------------


/*!
@page libsequence The libsequence library

See <a href="http://molpopgen.org/software/libsequence.html">libsequence</a>
and <a href="http://www.molpopgen.org/software/libsequence/doc/html">libsequence documentation</a>

Libsequence is a C++ library for population genetics developed by Kevin Thornton.
\b MCT uses Libsequence for some representations and routines.

Libsequence members used in MCT are denoted with the namespace Sequence:: 
and may include:

<ul>
<li>\anchor libsequence_arg <b>Sequence::arg</b>\n
An ARG represented as an ordered list of Sequence::marginal marginal trees.</li>
<li>\anchor libsequence_marginal <b>Sequence::marginal</b>\n
The genealogy of a portion of a chromosome on which no recombination has occurred, 
represented as a coalescent tree.</li>
</ul>

*/

//------------------------------------------------------------------------

/*!


