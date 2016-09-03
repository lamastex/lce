/*! \file TajDPvals1.cpp
    \brief main function for p-values from Tajima's D statistics 1 Locus at a time (less RAM)
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

int main(int argc, char **argv)
{
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

    cout << "Output for n = " << (nsam) << ", Sfs data set size " << p.NumOfSfsData;
    if (importanceSampling) {
        cout  << " with Importance Sampling, sample size " << importanceSampleSize;
    }
    else {
        cout << " and sampling to make sum(P(fsequence)) = 1";
    }
    cout << " and rgsl seed is " << p.seed << "\n" << endl;

    ofstream outfile;

    // make a file for output and put general info into it
    string filename = "fsequencestats";
    outfile.open (filename.c_str(), ios::out);
    if (outfile.is_open()) {
        outfile << "Process id " << myPid << endl;
        outfile << "Output for n = " << (p.SfsSz+1)
                << ", Sfs data set size " << p.NumOfSfsData;
        if (importanceSampling) {
            outfile << " with Importance Sampling, sample size "
            << importanceSampleSize << "\n";
        }
        else {
            outfile << " and sampling to make sum(P(fsequence)) = 1" << "\n";
        }
        outfile << "rgls seed is " << p.seed;
        outfile.close();
    }
    else {
        cout << "Error with likelihood output file " << filename << endl;
    }

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
   sfs_array_type obsSfs0(p.SfsSz);// sfs of all zoros
   gsl_ran_discrete_t *gslpdfstructUnc;
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


    /*
    // Could set up containers to summarise results in here

    */

    // go through each of the sfs in the data set in turn
    for (unsigned CNTR=0; CNTR < p.NumOfSfsData; CNTR++)
    {

        //std::cout << "SFS index "<< CNTR << std::endl;

        //Load the Observed (Simulated) Sfs and manipulate it
        sfs_array_type obsSfs(p.SfsSz);
        obsSfs = DataSfspectra[CNTR];
        SfsInfo SFSI(DataSfspectra[CNTR].size());
        SFSI.Inform(DataSfspectra[CNTR]);
        //Print the Observed Statistics and sfs and Pis
        //std::cout << SFSI;


        f_x_z_efficient_graph graphEfficientLoc; // a graph on product space, unlabelled
        f_efficient_graph graphVeryEfficientLoc; // a graph on f-space, unlabelled
        state_typeSet extFStateSetLoc; // set of unique f states
        state_typeSet extZStateSetLoc; // set of unique z states
        // map from set pairs of iterators to vertices
        ExtSetItPairVertexMap extStatesPairVertexMapLoc; 
        ExtSetItVertexMap extStateVertexMapLoc; // map from set iterators to vertices
        XStarSequencesMap xstarToSequencesInfoMapLoc; //map from xstar states to fsequence data
        StateSequenceSet fsequenceSetLoc; // set to store fsequences in for each locus

        bool canContinue = true; // controls whether loop continues after fsequences
        size_t totalFSImpSample = 0; // only used if importance sampling
        XStarSequencesMapIt xstarMapItLoc; // define iterator variable
        // this gets its value in the following code

        /* what are we doing here?  Many different sfs can have the same xstar or the
        same 'equivalent' xstar (eg x star (0,1,0,0,0) is equivalent to  (0,0,0,0,0)
        because it would be compatible with exactly the same complete collection fsequences).
        To save regenerating fsequences unnecessarily, we make have map of unique xstars
        where each unique xstar is mapped to its own collection of fsequences (mapped to
        their xstar-specific counts, probabilities etc as described above).

        When we 'generate' a sample of fsequences for a given sfs, we want to end up with
        an iterator to the map of unique xstars to fsequences, ie an iterator to a (key, value)
        pair where the key will be the xstar for that sfs, or an equivalent xstar and the value
        will be a representation of a sample of fsequences associated with that xstar.
        What the 'generate' functions do is first look whether the sfs's xstar,
        or an equivalent one, is already in the xstar map.
        If so, they just return an interator to that xstar's (key,value) pair in the map.
        If not, they generate a sample of fsequences for that xstar, add the new (key, value)
        pair to the map and return the iterator to this new pair in the map.
        The 'sample' can be all the fsequences associated with the xstar.  If
        importance sampling is used, a specified number of fsequences is generated which may
        or may not contain all the fsequences associated with that xstar.
        */


        // get an iterator to a map from xstars to fsequences
        if (importanceSampling) {

            boost::tie(xstarMapItLoc, totalFSImpSample)
                        = generateImportanceSample(graphVeryEfficientLoc,
                                        extStateVertexMapLoc,
                                        extFStateSetLoc,
                                        extZStateSetLoc,
                                        rgsl,
                                        obsSfs,
                                        xstarToSequencesInfoMapLoc,
                                        fsequenceSetLoc,
                                        nsam,
                                        importanceSampleSize);



            if (totalFSImpSample == 0) {
                cout << "Error trying to do importance sampling - aborting for sfs "
                        << CNTR << endl;
                canContinue = false; // can't carry on for this sfs sample
            }

        }

        else { // not importance sampling

            bool successfulGeneration = false;

            boost::tie(xstarMapItLoc, successfulGeneration)
                            = generateAllFSeq(graphEfficientLoc,
                                        extStatesPairVertexMapLoc,
                                        extFStateSetLoc,
                                        extZStateSetLoc,
                                        rgsl,
                                        obsSfs,
                                        xstarToSequencesInfoMapLoc,
                                        fsequenceSetLoc,
                                        nsam);

            if (!successfulGeneration) {
                cout << "Error trying to generate all fsequences for sfs "
                        << CNTR << endl;
                canContinue = false; // can't carry on for this sfs sample
            }
        }

        if (!canContinue) {
            continue;  // stop processing for this sfs, go to next one
        }
        // assuming we can continue ...
        // we now have an iterator into a map of xstars->collections of fsequences with info
        //cout << "about to loop over " << xstarMapItLoc->second.size()
        //                            << " fsequences for this xstar " << endl;
        double accKCProb = 0;
        double accPropProb = 0;

        vector<StateSequenceSetIt> fseqs;
        vector<double> fseqProbs;

        for (StateSequenceSetItSummaryMapIt fsMapIt = xstarMapItLoc->second.begin();
                    fsMapIt != xstarMapItLoc->second.end(); fsMapIt++) {

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
        int MAXNTRIES = 100000;
        sfs_array_type vaSFS(p.SfsSz);
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
            gsl_ran_multinomial (rgsl, p.SfsSz, ((unsigned) SFSI.S), &lineageLenItonsUnc[0],
                                 (unsigned *)(const_cast<int *> (&vaSFS[0])));
            SfsInfo SFSIsimUnc(vaSFS.size()); SFSIsimUnc.Inform(vaSFS);
            //printVA1D(vaSFS);
            if(abs(SFSIsimUnc.D) >= abs(SFSI.D)) TajimasDUnc_2pvalue++;
            if(SFSI.D >= 0) {if(SFSIsimUnc.D >= SFSI.D) TajimasDUnc_pvalue++;}
            else {if(SFSIsimUnc.D <= SFSI.D) TajimasDUnc_pvalue++;}
            
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
            //std::valarray<double> MPs = lineageLenItons/(lineageLenItons.sum());
            //printVA1D(MPs);
            gsl_ran_multinomial (rgsl, p.SfsSz, ((unsigned) SFSI.S), &lineageLenItons[0],
                                 (unsigned *)(const_cast<int *> (&vaSFS[0])));
            SfsInfo SFSIsim(vaSFS.size()); SFSIsim.Inform(vaSFS);
            //printVA1D(vaSFS);

            if(abs(SFSIsim.D) >= abs(SFSI.D)) TajimasD_2pvalue++;
            if(SFSI.D >= 0) {if(SFSIsim.D >= SFSI.D) TajimasD_pvalue++;}
            else {if(SFSIsim.D <= SFSI.D) TajimasD_pvalue++;}

            //cout << SFSIsim.S << "\t" << SFSIsim.Pi << "\t" << SFSIsim.D << endl;

            //size_t cherryStat = getCherryCoalescentStat(fs);
            //size_t unevenStat = getUnevenCoalescentStat(fs);

          }
        double TajimasD_2p = ((double)(TajimasD_2pvalue))/(double(MAXNTRIES));
        double TajimasDUnc_2p = ((double)(TajimasDUnc_2pvalue))/(double(MAXNTRIES));
        double TajimasD_p = ((double)(TajimasD_pvalue))/(double(MAXNTRIES));
        double TajimasDUnc_p = ((double)(TajimasDUnc_pvalue))/(double(MAXNTRIES));
        //cout << "TajimasD_pvalues Un/controlled = \n";
        cout << SFSI.S << "\t" << SFSI.Pi << "\t" << SFSI.D << "\t"
             << TajimasDUnc_2p << "\t" << TajimasD_2p << "\t"
             << TajimasDUnc_p << "\t" << TajimasD_p << endl;
        gsl_ran_discrete_free (gslpdfstruct);
    } // end loop over SFS data sets

    // would summarise results over SFS here

    gsl_ran_discrete_free (gslpdfstructUnc);
    gsl_rng_free(rgsl); //------------------------------- freeing the RNG
}
