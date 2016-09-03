/*! \file FseqCompatibles.cpp
    \brief main function for finding which unconditional F sequences are compatible with a SFS
*/
//Copyright(c) 2005, Raazesh Sainudiin and Kevin Thornton --  Wed Aug 17 23:44:58 EDT 2005
#include <ebc_sfstypes.hpp>
#include <ebc_graphtypes.hpp>
#include <ebc_moves.hpp>
#include <ebc_summaries.hpp>
#include <ebc_output.hpp>
#include <ebc_sfsinteg.hpp>
#include <ebc_sfsdata.hpp>
#include <ebc_params.hpp>
#include <ebc_precomp_TsPs.hpp>
#include <ebc_process_options.hpp>
#include <ebc_fseq.hpp>
#include <ebc_fxzseq_named_graph.hpp>
#include <ebc_fxzseq_efficient_graph.hpp>
#include <ebc_fseq_efficient_graph.hpp>
#include <ebc_fseq_general_graph.hpp>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>

#include <time.h> // for timing

//#define MOVE_TYPE 1

using namespace std;
using namespace boost;

// quick hack helper function
bool checkStateSequence(const StateSequence& sequence, const state_type& xstar)
{

    bool retValue = true;

    size_t length = (*(*sequence.begin())).size() - 1;

    state_type summary(length, 0);

    size_t i = 2;

    while (i < length && retValue) {
       if (xstar[i] > 0) {
            size_t summ = 0;
            StateSequenceConstIt it;
            for (it = sequence.begin(); it < sequence.end(); it++) {
                //the actual state is *(*it)
                summ += (**it)[i];
                            }
            if (summ == 0)
                retValue = false;
        }
        i++;
    }

    return retValue;
}


/*
Main generates a sample (importance sampling or straight sampling as specfied in shell script)
of f sequences for unconditional xstar for nsam from shell
Then takes each SFS in data file specified in shell script, generates 1 sfs just as
a quick hack to get the iterator to the xstar for that SFS, then takes each unconditional
f sequence and checks if it is compatible with this SFS's xstar.  For each compatible F sequence accumulates number of importance samples (only relevant if importance sampling), KC probability, and number of unique fsequences.  Output unconditional and conditional totals and proportions.

// You have to do a hack to compile this, or change makefile

*/

int main(int argc, char **argv)
{
    //set up a results file
    string filename_Res = "UncInCondResults_n10_all.txt";

    // get and print out the process id - helps manage long-running programs
    int myPid = getpid();
    cout << "This process id is " << myPid << endl;



    CheckUserDefinedDoubleTruncation();

    // useful statements for pretty printing of numbers
    std::cout << std::scientific << std::showpoint << std::setprecision(6);
    std::cout << std::fixed << std::showpoint << std::setprecision(6);

    // Getting options from command line arguments
    params p = process_options(argc,argv);

    const int nsam = p.SfsSz+1;    //sample size to simulate

    // vector of sfs data (these are typically simulated sfs data -- can also be real data)
    vector<sfs_array_type> DataSfspectra;

    //Read the Sfs data from the DataFile
    ReadSfsData(DataSfspectra, p);

    //gsl setup
    gsl_rng * rgsl = gsl_rng_alloc(gsl_rng_mt19937);
    if(!p.USER_SEED) p.seed = time( NULL ); //(unsigned short)time( NULL ) ;
    gsl_rng_set (rgsl, p.seed);
    //std::cout << "SEED = " << p.seed << std::endl;

    p.NumOfSfsData = DataSfspectra.size();

    // set up structures we may need to create graph
    /*
    two of the possible graphs are here, we only need one of them
    mainly depending on whether we are importance sampling (graphVeryEfficient)
    or not (graphEfficient)
    */

    f_x_z_efficient_graph graphEfficient; // a graph on product space, unlabelled

    f_efficient_graph graphVeryEfficient; // a graph on f-space, unlabelled

    state_typeSet extFStateSet; // set of unique f states
    state_typeSet extZStateSet; // set of unique z states

    ExtSetItPairVertexMap extStatesPairVertexMap; // map from set pairs of iterators to vertices
    ExtSetItVertexMap extStateVertexMap; // map from set iterators to vertices

    // parameters to control importance sampling or not and importance sample size
    bool importanceSampling = p.boolImportanceSampling;
    size_t importanceSampleSize = p.numImportanceSamples;

    cout << "Output for n = " << (nsam) << ", data file " << p.DataFile;
    if (importanceSampling) {
        cout  << " with Importance Sampling, sample size " << importanceSampleSize;
    }
    else {
        cout << " and sampling to make sum(P(fsequence)) = 1";
    }
    cout << " and rgsl seed is " << p.seed << "\n" << endl;


    // set up the structures we need to record fsequences with
    XStarSequencesMap xstarToSequencesInfoMap; //map from xstar states to fsequence data
    /* Each unique xstar can be associated with a number of fsequences.  To achieve this,
    each unique xstar maps to a another map (ie, key is xstar, value is another map).
    This other map contains representations of all the fsequences associated with that xstar.
    Each of these fsequences is mapped to information about that fsequence (keys are fsequences,
    values are xstar-specific information). The representation of each f-sequence used as the
    key in this map is just an iterator into a set of all the unique
    f-sequences but the information is particular to the xstar: it gives the counts (number of
    times that fsequence was generated for that xstar - all generation is stochastic - and
    the probabilities (proposal probability and Kingman probability) - see FseqInfo class
    */

    StateSequenceSet fsequenceSet; // set to store fsequences in
    /* To save space, each fsequence is held only once, in a set.  Each unique fsequence can be
    associated with many different x-stars by the xstars.  Each xstar only needs to know
    the iterators or pointers to the fsequences it is associated with.
    In fact, to save even more space, an fsequence itself is representated as a sequence of
    iterators to unique f-states.  Eg the same state (5,0,0,0,0) will be the final state for
    every f-sequence when n=5, and similarly other states will appear in more than one
    fsequence.  So, the unique possible states are held in a set and an fsequence is represented
    by a sequence of iterators into the fstate set.

    Note also that the fsequence is the forward fsequence and does not include the state
    in epoch 1 when one common ancestor subtends all nsam leaves.  So, the element
    at index [i], i = 0,..., nsam-2 in the fsequence represents the state at epoch i+2.
    */


    size_t totalFSImpSample = 0; // only used if importance sampling
    XStarSequencesMapIt xstarMapItUnc; // define iterator variable
    //gsl_ran_discrete_t *gslpdfstructUnc;

    sfs_array_type obsSfs0(nsam-1);// sfs of all zeros

    //vector<StateSequenceSetIt> fseqsUnc;
    //vector<double> fseqProbsUnc;

    double accAllUncKCProb = 0;
    size_t accAllUncImpSampleSize = 0;


    //uncontrolled f-sequence generation
    if (importanceSampling) {
        boost::tie(xstarMapItUnc, totalFSImpSample)
                    = generateImportanceSample(graphVeryEfficient,
                                    extStateVertexMap,
                                    extFStateSet,
                                    extZStateSet,
                                    rgsl,
                                    obsSfs0,
                                    xstarToSequencesInfoMap,
                                    fsequenceSet,
                                    nsam,
                                    importanceSampleSize);
       if (totalFSImpSample == 0) {
            cout << "Error trying to do importance sampling - aborting for sfs "
                    << endl;
        }
    }


    else { // not importance sampling
        bool successfulGeneration = false;
        boost::tie(xstarMapItUnc, successfulGeneration)
                        = generateAllFSeq(graphEfficient,
                                    extStatesPairVertexMap,
                                    extFStateSet,
                                    extZStateSet,
                                    rgsl,
                                    obsSfs0,
                                    xstarToSequencesInfoMap,
                                    fsequenceSet,
                                    nsam);

        if (!successfulGeneration) {
            cout << "Error trying to generate all fsequences for sfs "
                    << endl;
        }
    }

    size_t uniqueUncondFseq = xstarMapItUnc->second.size();

    cout << "Done uncontrolled Fs, size " << xstarMapItUnc->second.size() << endl;

    // we now have an iterator into a map of xstars->collections of fsequences with info
    //cout << "about to loop over " << xstarMapItUnc->second.size()
    //                          << " fsequences for this xstar " << endl;
    for (StateSequenceSetItSummaryMapIt fsMapIt = xstarMapItUnc->second.begin();
                fsMapIt != xstarMapItUnc->second.end(); fsMapIt++) {
        // StateSequence is a vector of iterators to a set of f-states
        // what we have here is an iterator to that state sequence in a set


        StateSequence fs = *(fsMapIt->first);
        accAllUncKCProb += fsMapIt->second.kcProb;
        accAllUncImpSampleSize += fsMapIt->second.countImpSample;

        //size_t impCount = fsMapIt->second.countImpSample;
        //cout << "impCount is " << impCount << endl;
        //cout << "kc prob is " << kcProbFSeq << "\n" << endl;
        //cout << "proposal prob is " << proposalProbFSeq << "\n" << endl;

    }


    ofstream fileRes;

    fileRes.open (filename_Res.c_str(), ios::trunc);
    if (fileRes.is_open()) {
        fileRes << "Process id " << myPid << endl;
        fileRes << "Output for n = " << (nsam) << ", Sfs data set size " << p.NumOfSfsData
        << ", data file " << p.DataFile << endl;
        if (importanceSampling) {
            fileRes << " with Importance Sampling, sample size "
            << importanceSampleSize << "\n";
        }
        else {
            fileRes << " and sampling to make sum(P(fsequence)) = 1" << "\n";
        }
        fileRes << "rgls seed is " << p.seed << "\n\n";
        fileRes << std::fixed << std::showpoint << std::setprecision(5) << std::endl;


        fileRes << "xStar" << "\t"
            << "uncFs" << "\t" << "uncKCPr" << "\t" << "uncISsize" << "\t"
            << "uInCFs" << "\t" << "uIncKCPr" << "\t" << "uIncISsize" << "\t"
            << "p Fs" << "\t"
            << "p KCpr" << "\t"
            << "p ISsize" << endl;

         fileRes.close();


    }
    cout << "Starting controlled SFS" << endl;

    // go through each of the sfs in the data set in turn
    for (unsigned CNTR=0; CNTR < p.NumOfSfsData; CNTR++)
    {

        //Load the Observed (Simulated) Sfs and manipulate it
        sfs_array_type obsSfs(p.SfsSz);
        obsSfs = DataSfspectra[CNTR];

        //SfsInfo SFSI(DataSfspectra[CNTR].size());
        //SFSI.Inform(DataSfspectra[CNTR]);
        //Print the Observed Statistics and sfs and Pis
        //std::cout << SFSI;



        bool canContinue = true; // controls whether loop continues after fsequences
        XStarSequencesMapIt xstarMapIt; // define iterator variable
        // this gets its value in the following code

        size_t tempimportanceSampleSize = 1;

        // get an iterator to a map from xstars to fsequences
        if (importanceSampling) {

            boost::tie(xstarMapIt, totalFSImpSample)
                        = generateImportanceSample(graphVeryEfficient,
                                        extStateVertexMap,
                                        extFStateSet,
                                        extZStateSet,
                                        rgsl,
                                        obsSfs,
                                        xstarToSequencesInfoMap,
                                        fsequenceSet,
                                        nsam,
                                        tempimportanceSampleSize);



            if (totalFSImpSample == 0) {
                cout << "Error trying to do importance sampling - aborting for sfs " << endl;
                canContinue = false; // can't carry on for this sfs
            }

        }

        else { // not importance sampling

            bool successfulGeneration = false;

            boost::tie(xstarMapIt, successfulGeneration)
                            = generateAllFSeq(graphEfficient,
                                        extStatesPairVertexMap,
                                        extFStateSet,
                                        extZStateSet,
                                        rgsl,
                                        obsSfs,
                                        xstarToSequencesInfoMap,
                                        fsequenceSet,
                                        nsam);

            if (!successfulGeneration) {
                cout << "Error trying to generate all fsequences for sfs " << endl;
                canContinue = false; // can't carry on for this sfs
            }
        }

        if (!canContinue) {
            continue;  // stop processing for this sfs, go to next one
        }


        // go through the unconditional f-seqs and find who is compatible
        size_t uniqueUncInCondFseq = 0;
        double accUncInCondKCProb = 0;
        double accUncInCondImpSampleSize = 0;

        state_type thisXstar = xstarMapIt->first;
        string xStarStr = stateToString(thisXstar);


        size_t jjj = 1;
        for (StateSequenceSetItSummaryMapIt fsMapIt = xstarMapItUnc->second.begin();
                    fsMapIt != xstarMapItUnc->second.end(); fsMapIt++) {
            // StateSequence is a vector of iterators to a set of f-states
            // what we have here is an iterator to that state sequence in a set

            StateSequence thisFseq = *(fsMapIt->first);


            bool isOkay = checkStateSequence(thisFseq, thisXstar);
            if (isOkay) {
                uniqueUncInCondFseq++;
                accUncInCondKCProb += (fsMapIt->second.kcProb);
                accUncInCondImpSampleSize += fsMapIt->second.countImpSample;
            }
            jjj++;
        }


        cout << xStarStr << "\t" << scientific << endl;
        cout << uniqueUncondFseq << "\t"
             << accAllUncKCProb << "\t" << accAllUncImpSampleSize << endl;
        cout << uniqueUncInCondFseq << "\t"
             << accUncInCondKCProb << "\t" << accUncInCondImpSampleSize << endl;
        cout << (1.0*uniqueUncInCondFseq)/(1.0*uniqueUncondFseq) << "\t"
             << accUncInCondKCProb/accAllUncKCProb << "\t"
             << (1.0*accUncInCondImpSampleSize)/(1.0*accAllUncImpSampleSize) << endl;

        fileRes.open (filename_Res.c_str(), ios::app);
        if (fileRes.is_open()) {

            fileRes << xStarStr << "\t" << scientific
             << uniqueUncondFseq << "\t"
             << accAllUncKCProb << "\t" << accAllUncImpSampleSize << "\t"
             << uniqueUncInCondFseq << "\t"
             << accUncInCondKCProb << "\t" << accUncInCondImpSampleSize << "\t"
             << (1.0*uniqueUncInCondFseq)/(1.0*uniqueUncondFseq) << "\t"
             << accUncInCondKCProb/accAllUncKCProb << "\t"
             << (1.0*accUncInCondImpSampleSize)/(1.0*accAllUncImpSampleSize) << endl;
            fileRes.close();
                    }
        else {
            cout << "Could not open " << filename_Res << endl;
        }


    }



    gsl_rng_free(rgsl); //------------------------------- freeing the RNG
}
