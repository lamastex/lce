/*! \file TajDPvalsOutputDs.cpp
    \brief main function for Tajima's D statistics with some Ds sent to file
*/

// You have to do a hack to compile this, or change makefile

// SFS to do D's for is hard coded


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

/*
Main takes one hard-coded SFS and generates MAXNTRIES D's using conditional and unconditional null distribution, and for first MAXOUTPUT Pvalues outputs individual Ds as well as p-values.

You have to do a hack to compile this, or change makefile

SFS to do D's for is hard coded

*/

int main(int argc, char **argv)
{
    // get and print out the process id - helps manage long-running programs
    int myPid = getpid();
    cout << "This process id is " << myPid << endl;

    const int nsam = 10;    //sample size to simulate
    int MAXOUTPUT = 1000; // number of tries to give the Ds for recording to file
    int MAXNTRIES = 100000; // number of tries to give the Ds for p-values

    // make up an SFS to do the Ds on
    string sfsIdStr = "t100_sfs899"; // a name for the SFS (used in output file name


    // ***************************************************************
    string thisSfsStr = "137 7 1 1 3 0 1 16 11"; // hard coded hack
    // and also have to do this - bad bad bad
    sfs_array_type thisSfs(nsam-1);// sfs of all zeros sfsToXstar
    thisSfs[2] = 6; thisSfs[3] = 3; thisSfs[4] = 12; thisSfs[5] = 6;
    thisSfs[6] = 0; thisSfs[7] = 0; thisSfs[8] = 84;
    // ***************************************************************


    CheckUserDefinedDoubleTruncation();

    // useful statements for pretty printing of numbers
    std::cout << std::scientific << std::showpoint << std::setprecision(6);
    std::cout << std::fixed << std::showpoint << std::setprecision(6);

    // Getting options from command line arguments
    params p = process_options(argc,argv);


    // vector of sfs data (these are typically simulated sfs data -- can also be real data)
    vector<sfs_array_type> DataSfspectra;

    //Read the Sfs data from the DataFile
    ReadSfsData(DataSfspectra, p);

    //gsl setup
    gsl_rng * rgsl = gsl_rng_alloc(gsl_rng_mt19937);
    if(!p.USER_SEED) p.seed = time( NULL ); //(unsigned short)time( NULL ) ;
    gsl_rng_set (rgsl, p.seed);
    //std::cout << "SEED = " << p.seed << std::endl;

    p.NumOfSfsData = (p.NumOfSfsData <= DataSfspectra.size() ) ? p.NumOfSfsData : DataSfspectra.size();


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

    cout << "Output for n = " << (nsam) << ", making Ds for " << thisSfsStr;
    if (importanceSampling) {
        cout  << " with Importance Sampling, sample size " << importanceSampleSize;
    }
    else {
        cout << " and sampling to make sum(P(fsequence)) = 1";
    }
    cout << " and rgsl seed is " << p.seed << "\n" << endl;

    ofstream outfile;



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
    gsl_ran_discrete_t *gslpdfstructUnc;

    sfs_array_type obsSfs0(nsam-1);// sfs of all zeros
    vector<StateSequenceSetIt> fseqsUnc;
    vector<double> fseqProbsUnc;

    {//uncontrolled f-sequence generation
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

        // we now have an iterator into a map of xstars->collections of fsequences with info
        //cout << "about to loop over " << xstarMapItUnc->second.size()
        //                          << " fsequences for this xstar " << endl;
        double accKCProb = 0;
        double accPropProb = 0;


        for (StateSequenceSetItSummaryMapIt fsMapIt = xstarMapItUnc->second.begin();
                    fsMapIt != xstarMapItUnc->second.end(); fsMapIt++) {
            // StateSequence is a vector of iterators to a set of f-states
            // what we have here is an iterator to that state sequence in a set
            StateSequence fs = *(fsMapIt->first);
            fseqsUnc.push_back(fsMapIt->first);
            double kcProbFSeq = (fsMapIt->second.kcProb);
            fseqProbsUnc.push_back(kcProbFSeq);

            accKCProb += kcProbFSeq;
            double proposalProbFSeq = (fsMapIt->second.proposalProb);
            accPropProb += proposalProbFSeq;
            /*
            size_t impCount = fsMapIt->second.countImpSample;
            cout << "impCount is " << impCount << endl;
            cout << "kc prob is " << kcProbFSeq << "\n" << endl;
            cout << "proposal prob is " << proposalProbFSeq << "\n" << endl;
            */
        }
        assert(fseqsUnc.size()==fseqProbsUnc.size());
        size_t K=fseqsUnc.size();
        //cout << "Acc kc prob for uncontrolled xstar is " << accKCProb << "\n\n" << endl;
        //cout << "Acc proposal prob for uncontrolled xstar is " << accPropProb << "\n\n" << endl;
        //cout << "There are " << K << " many unique f-seqs in Sampler Set\n\n" << endl;
        gslpdfstructUnc = gsl_ran_discrete_preproc (K, &fseqProbsUnc[0]);
    }//end of uncontrolled f-sequence generation

    //std::cout << "about to go through SFS ..." << std::endl;


    cout << "MAXNTRIES is " << MAXNTRIES << endl;
    cout << "MAXOUTPUT is " << MAXOUTPUT << endl;

    //Load the test Sfs and manipulate it
    sfs_array_type obsSfs(nsam-1);
    obsSfs = thisSfs;
    SfsInfo SFSI(obsSfs.size());
    SFSI.Inform(obsSfs);
    //Print the Observed Statistics and sfs and Pis
    //std::cout << SFSI;

    ofstream fileDs;
    string filename_Ds = "FileOfDs";
    filename_Ds = filename_Ds + sfsIdStr + ".txt";

    fileDs.open (filename_Ds.c_str(), ios::trunc);
    if (fileDs.is_open()) {
        fileDs << "Process id " << myPid << endl;
        fileDs << "Output for n = " << nsam << ", Unconditional Ds for "
            << thisSfsStr << " with MAXNTRIES = " << MAXNTRIES
            << " and MAXOUTPUT = " << MAXOUTPUT << ",\n";
        if (importanceSampling) {
            fileDs << " with Importance Sampling, sample size "
            << importanceSampleSize << "\n";
        }
        else {
            fileDs << " and sampling to make sum(P(fsequence)) = 1" << "\n";
        }
        fileDs << "rgls seed is " << p.seed << "\n";
        fileDs << std::fixed << std::showpoint << std::setprecision(5) << std::endl;
        fileDs << "\nD for this SFS is " << SFSI.D << '\n' << std::endl;

        fileDs.close();
    }

    bool canContinue = true; // controls whether loop continues after fsequences
    XStarSequencesMapIt xstarMapIt; // define iterator variable
    // this gets its value in the following code

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
                                    importanceSampleSize);



        if (totalFSImpSample == 0) {
            cout << "Error trying to do importance sampling - aborting for sfs "
                    << thisSfsStr << endl;
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
            cout << "Error trying to generate all fsequences for sfs "
                    << thisSfsStr << endl;
            canContinue = false; // can't carry on for this sfs
        }
    }

    if (canContinue) {

        // assuming we can continue ...
        // we now have an iterator into a map of xstars->collections of fsequences with info
        //cout << "about to loop over " << xstarMapIt->second.size()
        //                            << " fsequences for this xstar " << endl;
        double accKCProb = 0;
        double accPropProb = 0;

        vector<StateSequenceSetIt> fseqs;
        vector<double> fseqProbs;

        for (StateSequenceSetItSummaryMapIt fsMapIt = xstarMapIt->second.begin();
                    fsMapIt != xstarMapIt->second.end(); fsMapIt++) {

            // StateSequence is a vector of iterators to a set of f-states
            // what we have here is an iterator to that state sequence in a set


            StateSequence fs = *(fsMapIt->first);
            fseqs.push_back(fsMapIt->first);

            /*
            size_t cherryStat = getCherryCoalescentStat(fs);
            size_t unevenStat = getUnevenCoalescentStat(fs);
            cout << "f sequence\n";
            printStateSequence(fs);
            cout << "has cherryStat " << cherryStat << endl;
            cout << "and unevenStat " << unevenStat << endl;
            */
            double kcProbFSeq = (fsMapIt->second.kcProb);
            fseqProbs.push_back(kcProbFSeq);

            accKCProb += kcProbFSeq;
            double proposalProbFSeq = (fsMapIt->second.proposalProb);
            accPropProb += proposalProbFSeq;
            /*
            size_t impCount = fsMapIt->second.countImpSample;
            cout << "impCount is " << impCount << endl;
            cout << "kc prob is " << kcProbFSeq << "\n" << endl;
            cout << "proposal prob is " << proposalProbFSeq << "\n" << endl;
            */
        }

        assert(fseqs.size()==fseqProbs.size());
        size_t K=fseqs.size();
        //cout << "Acc kc prob for this xstar is " << accKCProb << "\n\n" << endl;
        //cout << "Acc proposal prob for this xstar is " << accPropProb << "\n\n" << endl;
        //cout << "There are " << K << " many unique f-seqs in Sampler Set\n\n" << endl;
        gsl_ran_discrete_t *gslpdfstruct;
        gslpdfstruct = gsl_ran_discrete_preproc (K, &fseqProbs[0]);

        int TajimasD_2pvalue = 0;
        int TajimasDUnc_2pvalue = 0;
        int TajimasD_pvalue = 0;
        int TajimasDUnc_pvalue = 0;

        vector<double> unconditionalDcoll;
        vector<double> conditionalDcoll;


        sfs_array_type vaSFS(nsam-1);
        //Get the total tree size in TEES[i] and the epoch-times in PEES[i]
    //        make_set_of_pees_tees(rgsl, PEES, TEES, growth_rate, p);
        //cout << "about to try"<< endl;
        for (int tries = 1; tries <= MAXNTRIES; tries++)
          {                 // try a bunch of times
            //cout << "about to try"<< endl;
            valarray<double> EpochTimes;
            EpochTimes = make_StdNeutral_EpochTimes(rgsl, p, EpochTimes);
            //cout << "EpochTimes is "<< endl;
            //stdoutArray(EpochTimes);

            // uncontrolled p-values
            int proposed_index_Unc = static_cast<int>(gsl_ran_discrete (rgsl, gslpdfstructUnc));
            //cout << proposed_index_Unc << endl;
            StateSequence fsUnc = *(fseqsUnc[proposed_index_Unc]);
            //cout << "f sequence Unc trial " << tries << "\n" << endl;
            //printStateSequence(fsUnc);
            valarray<double> lineageLenItonsUnc;
            lineageLenItonsUnc = epochTimesProdFseq(lineageLenItonsUnc,EpochTimes,fsUnc,nsam);
            //cout << "\n Lineage Lengths are "<< endl;
            //printVA1D(lineageLenItonsUnc);

            gsl_ran_multinomial (rgsl, nsam-1, ((unsigned) SFSI.S), &lineageLenItonsUnc[0],
                                 (unsigned *)(const_cast<int *> (&vaSFS[0])));

            SfsInfo SFSIsimUnc(vaSFS.size());
            SFSIsimUnc.Inform(vaSFS);
            //printVA1D(vaSFS);

            if(abs(SFSIsimUnc.D) >= abs(SFSI.D))
                TajimasDUnc_2pvalue++;
            if(SFSI.D >= 0) {
                if(SFSIsimUnc.D >= SFSI.D)
                    TajimasDUnc_pvalue++;
            }
            else {
                if(SFSIsimUnc.D <= SFSI.D)
                    TajimasDUnc_pvalue++;
            }
            if (tries < MAXOUTPUT)
                unconditionalDcoll.push_back(SFSIsimUnc.D);

            // controlled p-values
            int proposed_index = static_cast<int>(gsl_ran_discrete (rgsl, gslpdfstruct));
            //cout << proposed_index << endl;
            StateSequence fs = *(fseqs[proposed_index]);
            //cout << "f sequence trial " << tries << "\n" << endl;
            //printStateSequence(fs);
            valarray<double> lineageLenItons;
            lineageLenItons = epochTimesProdFseq(lineageLenItons,EpochTimes,fs,nsam);
            //cout << "\n Lineage Lengths are "<< endl;
            //printVA1D(lineageLenItons);
            gsl_ran_multinomial (rgsl, nsam-1, ((unsigned) SFSI.S), &lineageLenItons[0],
                                 (unsigned *)(const_cast<int *> (&vaSFS[0])));
            SfsInfo SFSIsim(vaSFS.size()); SFSIsim.Inform(vaSFS);
            //printVA1D(vaSFS);
            if (abs(SFSIsim.D) >= abs(SFSI.D))
                TajimasD_2pvalue++;
            if (SFSI.D >= 0) {
                if(SFSIsim.D >= SFSI.D)
                    TajimasD_pvalue++;
            }
            else {
                if (SFSIsim.D <= SFSI.D)
                    TajimasD_pvalue++;
            }

            if (tries < MAXOUTPUT)
                conditionalDcoll.push_back(SFSIsim.D);

            //cout << SFSIsim.D << endl;

            //size_t cherryStat = getCherryCoalescentStat(fs);
            //size_t unevenStat = getUnevenCoalescentStat(fs);

        }

        double TajimasD_2p = ((double)(TajimasD_2pvalue))/(double(MAXNTRIES));
        double TajimasDUnc_2p = ((double)(TajimasDUnc_2pvalue))/(double(MAXNTRIES));
        double TajimasD_p = ((double)(TajimasD_pvalue))/(double(MAXNTRIES));
        double TajimasDUnc_p = ((double)(TajimasDUnc_pvalue))/(double(MAXNTRIES));
        //cout << "TajimasD_pvalues Un/controlled = \n";

        fileDs.open (filename_Ds.c_str(), ios::app);
        if (fileDs.is_open()) {
            fileDs << "uncondD\tcondD" << endl;
            vector<double>::iterator uit = unconditionalDcoll.begin();
            vector<double>::iterator cit = conditionalDcoll.begin();
            for (uit, cit;
                    (uit < unconditionalDcoll.end() || cit < conditionalDcoll.end())
                    ; uit++, cit++) {

                    fileDs << *uit << '\t' << *cit << '\n';
            }
            fileDs << '\n' << "SFSI.D" << "\t"
             << "TajDUnc_2p" << "\t" << "TajD_2p" << "\t"
             << "TajDUnc_p" << "\t" << "TajD_p" << endl;
            fileDs << SFSI.D << "\t"
             << TajimasDUnc_2p << "\t" << TajimasD_2p << "\t"
             << TajimasDUnc_p << "\t" << TajimasD_p << endl;
            fileDs.close();
        }
        else {
            cout << "Could not open " << fileDs << endl;
        }


        gsl_ran_discrete_free (gslpdfstruct);

        gsl_ran_discrete_free (gslpdfstructUnc);
    }

    gsl_rng_free(rgsl); //------------------------------- freeing the RNG
}
