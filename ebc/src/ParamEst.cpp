/*! \file ParamEst.cpp
    \brief main function for parameter estimation
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

#include <gsl/gsl_sf_gamma.h>

#include <time.h> // for timing

#define MOVE_TYPE 1

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

    OneMove OneMoveTemp(p.SfsSz);// temp variable of the struct OneMove
    // vector of moves
    vector<OneMove> MovesVector;
    //Read the moves from the MovesFile *****************************
    ReadMoves(MovesVector, p);

    //gsl setup
    gsl_rng * rgsl = gsl_rng_alloc(gsl_rng_mt19937);
    if(!p.USER_SEED) p.seed = time( NULL ); //(unsigned short)time( NULL ) ;
    gsl_rng_set (rgsl, p.seed);
    //std::cout << "SEED = " << p.seed << std::endl;

    double prob_min_sieve = 1e300;

    p.NumOfSfsData = (p.NumOfSfsData <= DataSfspectra.size() ) ? p.NumOfSfsData : DataSfspectra.size();

    //Load and print the Parameter Grid
    std::vector<output_params> output_GridParams;
    PrintLoadGridParams(p, output_GridParams);
    // after this call, the grid contains (numRho x numGrowth x numTheta) output_param types
    // each holding the relevant value for rho, growth, theta parameters
    // numRho, numGrowth, numTheta depend on the parameters specified for the run

    double volume = p.theta_per_locus_incr * p.growth_rate_incr;
    unsigned trueIndexX = 10; // with a 101 X 101 theta X growth_rate grid,  true parameter 10,0
    double coverageCriteria = 0.95; // coverage criteria to use

    // ie, TrueIndex says where in the parameter grid the output_params
    // with the true parameter values is
    // recombination rate per locus
    double true_rho_per_locus = output_GridParams[trueIndexX].rho_per_locus;
    // exponential growth rate
    double true_growth_rate =output_GridParams[trueIndexX].growth_rate;
    // scaled mutation rate theta per locus
    double true_theta_per_locus = output_GridParams[trueIndexX].theta_per_locus;

    // grid size
    unsigned gridSize = (((p.rho_max - p.rho_min)/p.rho_incr) + 1)
                            * p.growth_number_of_points
                            * p.theta_number_of_points;
    double lnGridSize = log(gridSize); // natural log

    // set up structures we need to create graph
    f_x_z_named_graph graphNamed; // a graph on product space using inefficient labelled implementation

    f_x_z_efficient_graph graphEfficient; // a graph on product space, unlabelled

    f_efficient_graph graphVeryEfficient; // a graph on f-space, unlabelled

    state_typeSet extFStateSet; // set of unique f states
    state_typeSet extZStateSet; // set of unique z states

    ExtSetItPairVertexMap extStatesPairVertexMap; // map from set pairs of iterators to vertices
    ExtSetItVertexMap extStateVertexMap; // map from set iterators to vertices

    // parameters to control importance sampling or not and importance sample size
    bool importanceSampling = p.boolImportanceSampling;
    size_t importanceSampleSize = p.numImportanceSamples;

    cout << "Output for n = " << (nsam) << ", Sfs data set size " << p.NumOfSfsData
            << ", number of replicates per fsequence (pees) " << p.Num_Of_Replicates;
    if (importanceSampling) {
        cout  << " with Importance Sampling, sample size " << importanceSampleSize;
    }
    else {
        cout << " and sampling to make sum(P(fsequence)) = 1";
    }
    cout << "\nGrid size is " << gridSize << " and rgsl seed is " << p.seed << "\n" << endl;

    ofstream outfile;

    // make a file for output and give it header lines
    string startFileName = p.LklOutFile;
    string filename = startLikelihoodOutputFile(outfile, myPid,
                                importanceSampling, importanceSampleSize,
                                p, coverageCriteria,
                                true_rho_per_locus, true_growth_rate, true_theta_per_locus);


    // set up the structures we need to record fsequences with
    XStarSequencesMap xstarToSequencesInfoMap; //map from xstar states to fsequence data
    StateSequenceSet fsequenceSet; // set to store fsequences in

    std::cout << "about to go through SFS, including generating Fseqns ..." << std::endl;

    // for timing
    clock_t startTime = 0, endTime = 0;
    clock_t accSeqTime = 0;
    clock_t accNonSeqTime = 0;
    clock_t accPeesTime = 0;
    clock_t accLineageTime = 0;
    clock_t accIntegrateTime = 0;
    clock_t accRGenTime = 0;
    clock_t accTime = 0;

    // Containers to summarise results in
    std::vector<output_params> resultsMaxLikelihoods; // maximum likelihood results
    std::vector<double> resultsKLDivergence; // Kullback-Leibler divergence
    size_t resultsGrowthCoverageAccumulation = 0; // growth parameter coverage indicator
    size_t resultsThetaCoverageAccumulation =0; // theta parameter coverage indicator
    size_t resultsJointGrowthThetaCoverageAccumulation =0; // joint coverage indicator

    for (unsigned CNTR=0; CNTR < p.NumOfSfsData; CNTR++) //loop over the sfs data sets
    {

        // only needed for timing
        clock_t thisNonSeqTime = 0;
        clock_t thisSeqTime = 0;
        clock_t thisPeesTime = 0;
        clock_t thisLineageTime = 0;
        clock_t thisIntegrateTime = 0;
        clock_t thisRGenTime = 0;

        if (p.NOISY == 10) { // timing
            std::cout << "SFS index "<< CNTR << std::endl;
            // timing
            startTime = clock();

        }

        //std::vector< output_results2 > output;
        //std::vector< output_LklShi > output;
        // vectors to store various likelihood values
        std::vector< double > LKL_MTr;
        std::vector< double > LKL_MDT;
        std::vector< double > LKL_MPT;
        std::vector< double > LKL_MPD;
        std::vector< double > LKL_Trs;

        //Load the Observed (Simulated) Sfs and manipulate it
        sfs_array_type obsSfs(p.SfsSz);
        obsSfs = DataSfspectra[CNTR];

        sfs_sieve sfsSieve;// the basic "sieve" or set of sfs that satisfy a given set of summaries
        sfs_sieve::iterator SfsSieve_Iter; // the corresponding iterator

        // used only for initial sieve construction
        // -- basically heating the observed (MLE) multinomail probs
        //std::vector<double> ObsMultinomPs(p.SfsSz);
        //SfsToMultinomialPs(obsSfs, ObsMultinomPs, 1 );
        std::vector<double> UnfMultinomPs(p.SfsSz);

        for(unsigned i=0; i < p.SfsSz; i++)
            UnfMultinomPs[i] = 1.0 / ((double) p.SfsSz ); // uniform distribution

        // prior alphas in the Dirichlet conjugate to the multinomial pi's
        vector<double> MeanPostDirich(p.SfsSz, p.prior_alpha);

        MeanPost_UnifDirichletPrior(obsSfs, MeanPostDirich, 1.0, 0);
        // used only for initial sieve construction
        // -- basically heating the observed (MLE) multinomail probs
        // MeanPost_UnifDirichletPrior(obsSfs, MeanPostDirich, p.heat, 1);
        // posterior (with prior alphas) in the Dirichlet conjugate
        // to the multinomial pi's giving Law of SFS

        // prior alphas in the Dirichlet conjugate to the multinomial pi's
        vector<double> AlphasPostDirich(p.SfsSz, p.prior_alpha);

        SfsInfo SFSI(DataSfspectra[CNTR].size());
        SFSI.Inform(DataSfspectra[CNTR]);
        //Print the Observed Statistics and sfs and Pis
        if(p.NOISY!=0 && p.NOISY!=10 && p.NOISY!=56){
            std::cout << SFSI;
        }
        //cout << SFSI.S << "\t" << SFSI.Pi << "\t" << SFSI.D << endl;
        ////////////--BEGIN SFS SIEVE-BUILDING (DETERMINISTIC/STOCHASTIC)--////////////

        // building the sfsSieve by starting the chain at the observed sfs

        double Observed_Sieve_Prob=0.;

        switch(p.sieve_Construct_Type)
        {
            case 2 :
                Observed_Sieve_Prob = build_sieve(rgsl, sfsSieve, obsSfs, MovesVector,
                                        MovesVector.size(), MeanPostDirich,
                                        &prob_min_sieve, p,
                                        p.sieve_NeededSamples, p.sieve_MaxTrials );
                //ObsMultinomPs, &prob_min_sieve, p, p.sieve_NeededSamples, p.sieve_MaxTrials );
            break;

            case 1 :
                Observed_Sieve_Prob = expand_uniformly_sieve(rgsl, sfsSieve, obsSfs,
                                        MovesVector, MovesVector.size(),
                                        UnfMultinomPs,
                                        &prob_min_sieve, p,
                                        p.sieve_NeededSamples, p.sieve_MaxTrials );
            break;

            case 0 :
                Observed_Sieve_Prob = graph_sieve(sfsSieve, obsSfs, MovesVector,
                                        MovesVector.size(), UnfMultinomPs,
                                        &prob_min_sieve, p,
                                        p.sieve_NeededSamples, p.sieve_MaxTrials );
            break;

            case -1 :
            //not making sieve at all
            break;

            default :
                std::cerr << "ERROR: p.sieve_Construct_Type = 0,1,2 or -1  BUT it is "
                        << p.sieve_Construct_Type << std::endl;
            exit(0);
            break; // !!!
        } // end switch

        if(p.NOISY>99) {
            cout << "# built sieve sfsSieve by method " << p.sieve_Construct_Type
                << std::endl
                << " with average Multinomial Prob = " << Observed_Sieve_Prob
                << "\t and Poisson Multinomial Prob = "
                << (Observed_Sieve_Prob + log(gsl_ran_poisson_pdf (SFSI.uS, double(SFSI.uS) )))
                << std::endl;
        }

        unsigned sfsSieveSize = sfsSieve.size();

        ///////////--END OF SFS SIEVE-BUILDING (DETERMINISTIC/STOCHASTIC)--//////////

        ////////////////--BEGIN COALECENT SIMULATINS OF GENEALOGIES--///////////////

        // timing
        if (p.NOISY == 10) {
            endTime = clock();
            thisNonSeqTime = endTime - startTime;
            startTime = clock();
        }

        bool canContinue = true; // controls whether loop continues after fsequences
        size_t totalFSImpSample = 0; // only used if importance sampling
        XStarSequencesMapIt xstarMapIt;

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
                        << CNTR << endl;
                canContinue = false; // can't carry on for this sfs sample
            }

        }

        else {

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
                        << CNTR << endl;
                canContinue = false; // can't carry on for this sfs sample
            }
        }

        if (p.NOISY == 10) { // timing
            endTime = clock();
            thisSeqTime = endTime - startTime;
            accSeqTime += thisSeqTime;
        }

        if (p.NOISY == 10) {    // timing
            cout << "end generated sequences for index " << CNTR << endl;
        }

        //timing
        if (importanceSampling && (p.NOISY == 10)) {
            cout << "generated " << totalFSImpSample << " fsequences" << endl;
        }

        if (!canContinue) {
            continue;  // stop processing for this sfs, go to next one
        }

        // assuming we can continue ...

        // timing
        startTime = clock();

        double rho = 0.0;   //rho=4Nr = recombination rate across region
        unsigned num_Of_EpochTimes = p.Num_Of_Replicates;

        //std::cout << "about to enter loops over rho"<< std::endl;

        //for(rho = p.rho_min; rho < p.rho_min + p.rho_incr; rho += p.rho_incr) {
        for(rho = p.rho_min; rho <= p.rho_max; rho += p.rho_incr) {

            //std::cout << "inside rho-loop with rho = " << rho << std::endl;

            //std::cout << "about to enter loops over growth"<< std::endl;

            for( int growth_i = 0 ; growth_i < p.growth_number_of_points ; ++growth_i) { // growth loop

                double growth_rate = (double(growth_i) * p.growth_rate_incr) + p.growth_rate_min;
                //std::cerr << "inside growth-loops with growth = " << growth_rate << std::endl;

                Tees TEES;
                Pees PEES;

                //debug
                //std::cout << "growth rate "<< growth_rate << "growth_i = "
                   // << growth_i << ")" << std::endl;

                // debug - temporary rgsl
                //cout << "making pees and tees with temporary rgsl seeded with 1111" << endl;

                //Get the total tree size in TEES[i] and the epoch-times in PEES[i]
                make_set_of_pees_tees(rgsl, PEES, TEES, growth_rate, p);

                if(TEES.size() != PEES.size()) {
                    std::cerr << std::endl << "Tees and Pees not same size!" << std::endl;
                    exit(1);
                }

                //std::cout << "about to enter loops over theta"<< std::endl;

                // theta loop
                for( int theta_i = 0 ; theta_i < p.theta_number_of_points ; ++theta_i) {

                    double theta_per_locus = (double(theta_i) * p.theta_per_locus_incr)
                                                    + p.theta_per_locus_min;
                    double theta_per_site = theta_per_locus/(double(p.length));

                    double sieve_integral_on_curr_arg=0.0;
                    double sieve_integral_on_all_args=0.0;
                    double sieve_ShanInfo_on_curr_arg=0.0;

                    // loop over each compatible f_sequence for this sfs
                    // xstarBool.first is an iterator pointing to a pair (xstar, fsequenceMap)
                    // fsequenceMap is a map from an iterator pointing to an fsequence to a
                    // struct FseqInfo
                    // so we want to iterate over this map

                    //cout << "about to loop over " << xstarMapIt->second.size()
                      //              << "fsequences for this xstar " << endl;
                    for (StateSequenceSetItSummaryMapIt fsMapIt = xstarMapIt->second.begin();
                                fsMapIt != xstarMapIt->second.end(); fsMapIt++) {

                        // if we are using importance sampling, we will have to remember to
                        // take into account the count of each fsequence as well as
                        // the probabilities

                        double kcProbFSeq = (fsMapIt->second.kcProb);
                        double proposalProbFSeq = (fsMapIt->second.proposalProb);

                        // StateSequence is a vector of iterators to a set of f-states
                        // what we have here is an iterator to that state sequence in a set

                        StateSequence fs = *(fsMapIt->first);

                        // variables incremented in the for loop over PEES
                        double sieve_integral_on_curr_fseq = 0.0;

                        // timing
                        clock_t timePees = 0;
                        clock_t lineageCalcTime = 0;
                        clock_t integrateCalcTime = 0;
                        clock_t rGenTime = 0;

                        for(size_t sampled_TeesPees = 0;
                            sampled_TeesPees < PEES.size(); sampled_TeesPees++) {

                            // only needed for timing
                            clock_t thisFSeqPeesTime = 0;
                            clock_t thisFSeqRGenTime = 0;
                            clock_t thisFSeqIntegrateCalcTime = 0;
                            clock_t thisFSeqLineageCalcTime = 0;
                            clock_t startThisPees = clock();

                            // set up an array for lineageLenItons
                            valarray<double> lineageLenItons;

                            // timing
                            clock_t startLineageCalcTime = clock();

                            // calculate lineage lengths
                            lineageLenItons = epochTimesProdFseq(lineageLenItons,
                                                                PEES[sampled_TeesPees],
                                                                fs,
                                                                nsam);

                            if (p.NOISY == 10) {
                                // timing
                                clock_t endLineageCalcTime = clock();
                                thisFSeqPeesTime += endLineageCalcTime - startThisPees;
                                thisFSeqLineageCalcTime
                                        = endLineageCalcTime - startLineageCalcTime;
                                lineageCalcTime += thisFSeqLineageCalcTime;
                                startThisPees = clock();
                            }
                            clock_t startIntegrateCalcTime = clock();

                            // integration
                            integrate_over_sieve(sfsSieve, lineageLenItons, p,
                                        sieve_integral_on_curr_arg, sieve_ShanInfo_on_curr_arg);

                            if (p.NOISY == 10) { //timing
                                clock_t endIntegrateCalcTime = clock();
                                thisFSeqPeesTime += endIntegrateCalcTime - startThisPees;
                                thisFSeqIntegrateCalcTime
                                    = endIntegrateCalcTime - startIntegrateCalcTime;
                                integrateCalcTime += thisFSeqIntegrateCalcTime;
                                startThisPees = clock();
                            }
                            clock_t startRGenTime = clock();

                            // poisson probability
                            double Theta_dependent_Poisson_Prob
                                = gsl_ran_poisson_pdf
                                    (SFSI.uS,(theta_per_site*TEES[sampled_TeesPees]));

                            if (p.NOISY == 10) { // timing

                                clock_t endRGenTime = clock();
                                thisFSeqPeesTime += endRGenTime - startThisPees;
                                thisFSeqRGenTime = endRGenTime - startRGenTime;
                                rGenTime += thisFSeqRGenTime;
                                startThisPees = clock();
                            }

                            // bringing it all together
                            double ln_Theta_dependent_Poisson_Prob
                                    = std::log(Theta_dependent_Poisson_Prob);

                            double ln_sieve_integral = log(sieve_integral_on_curr_arg);

                            double ln_sieve_integral_Poisson_Prob = ln_sieve_integral
                                                    + ln_Theta_dependent_Poisson_Prob;

                            if( (ln_Theta_dependent_Poisson_Prob < MOST_NEGATIVE_DOUBLE_ALLOWED)
                                || (ln_sieve_integral < MOST_NEGATIVE_DOUBLE_ALLOWED)      ) {

                                if(p.NOISY>1000) {
                                    std::cerr << std::endl
                                        << "WARNING log(prob) too small for screen !!!"
                                        << " resetting to "
                                            << MOST_NEGATIVE_DOUBLE_ALLOWED << std::endl;
                                }

                                ln_sieve_integral_Poisson_Prob = MOST_NEGATIVE_DOUBLE_ALLOWED;
                            }

                            sieve_integral_on_curr_fseq
                                += std::exp(ln_sieve_integral_Poisson_Prob);

                            if (p.NOISY == 10) { // timing

                                clock_t endThisPees = clock();

                                thisFSeqPeesTime += endThisPees - startThisPees;
                                timePees += thisFSeqPeesTime;
                            }

                        } // end for sampled TeesPees over PEES

                        if (p.NOISY == 10) { // timing
                            thisNonSeqTime += timePees;
                            thisPeesTime += timePees;
                            thisLineageTime += lineageCalcTime;
                            thisIntegrateTime += integrateCalcTime;
                            thisRGenTime += rGenTime;
                            startTime = clock();
                        }

                        // MC integration over num_Of_EpochTimes many EpochTimes
                        // but this is constant for all sequences

                        // if we are doing importance sampling, we need to take account of
                        // fsequence counts and proposal probabilities
                        if (importanceSampling) {

                            sieve_integral_on_curr_fseq *=
                                ((static_cast<double>(fsMapIt->second.countImpSample)
                                            /totalFSImpSample) / proposalProbFSeq);
                        }

                        // whether importance sampling or not [adjusted] sieve integral*P(f)
                        // is added to the integral over all f sequences
                        sieve_integral_on_all_args += (sieve_integral_on_curr_fseq * kcProbFSeq);

                    } // end replicates_per_theta (f-sequence) loop

                    LKL_Trs.push_back(sieve_integral_on_all_args);

                    if(p.NOISY==56 ) {
                        std::cerr << "\rgrowth = " << growth_rate << "  theta = " << theta_per_locus
                              << "\t" << num_Of_EpochTimes
                              << "  Replicate " << (CNTR+p.DATA_id) << std::flush;
                    }

                }// end theta loop
            }//end growth loop
        }//end rho loop

        // we did each theta num_Of_EpochTimes times
        unsigned actual_Num_Of_Replicates_per_theta = num_Of_EpochTimes;
        // for each fsequence that we used in the integral
        if (importanceSampling) {
            actual_Num_Of_Replicates_per_theta *= totalFSImpSample;
        }
        else {
            actual_Num_Of_Replicates_per_theta *= xstarMapIt->second.size();
        }

        if(actual_Num_Of_Replicates_per_theta > 0 ) {

            // the following containers are filled by SummarizeDensityGrid call
            std::vector<double> shanInfoX;
            std::vector<double> coverageProb;
            std::vector<double> maxX;
            std::vector<unsigned> maxIndexX;

            // summarise the density grid - fills in MaxIndexX, MaxX, CoverageProb, ShanInfoX
            SummarizeDensityGrid(LKL_Trs, volume, maxIndexX, maxX, coverageProb, trueIndexX, shanInfoX, p);

            // calculate the Kullback-Leibler divergence measure
            // reference measure is the uniform prior over the grid
            // which 1/gridSize
            // so log(reference measure) = log(1/gridSize) = -lnGridSize
            // Kullback-Leibler divergence simplifies to
            // (-Shannons Entropy + lnGridSize)/lnGridSize
            double divergenceKL = (-shanInfoX[2] + lnGridSize)/lnGridSize;

            // send to the output file we started right at the top

            unsigned indexSFS = CNTR+p.DATA_id;

            // output results for this SFS
            outputResultsForOneSFS(filename,
                            indexSFS, actual_Num_Of_Replicates_per_theta,
                            maxIndexX, output_GridParams,
                            maxX, shanInfoX,
                            coverageProb, coverageCriteria,
                            sfsSieveSize,
                            true_rho_per_locus, true_growth_rate, true_theta_per_locus,
                            divergenceKL);

            // store results we want to summarise
            resultsMaxLikelihoods.push_back(output_GridParams[maxIndexX[2]]);
            resultsKLDivergence.push_back(divergenceKL);
            if (coverageProb[0] <= coverageCriteria)
                resultsGrowthCoverageAccumulation ++;
            if (coverageProb[1] <= coverageCriteria)
                resultsThetaCoverageAccumulation ++;
            if (coverageProb[2] <= coverageCriteria)
                resultsJointGrowthThetaCoverageAccumulation++;

        }

        if (p.NOISY == 10) { // timing
            endTime = clock();
            thisNonSeqTime += (endTime - startTime);
            accNonSeqTime += thisNonSeqTime;
            accPeesTime += thisPeesTime;
            accLineageTime += thisLineageTime;
            accIntegrateTime += thisIntegrateTime;
            accRGenTime += thisRGenTime;

            std::cout << "finished generating index " << CNTR << ", took "
            << ((static_cast<double>(thisSeqTime)) / CLOCKS_PER_SEC)
                << "s sequence time, "
                << ((static_cast<double>(thisNonSeqTime)) / CLOCKS_PER_SEC)
                << "s non sequence time" << std::endl;
            cout << ((static_cast<double>(thisSeqTime + thisNonSeqTime)) / CLOCKS_PER_SEC)
                << "s in total" << std::endl;
            std::cout << "Pees took "
                    << ((static_cast<double>(thisPeesTime)) / CLOCKS_PER_SEC)
                    << "s" << std::endl;
                std::cout << "of which  "
                    << ((static_cast<double>(thisLineageTime)) / CLOCKS_PER_SEC)
                    << "s was over lineage";
                std::cout << " and "
                << ((static_cast<double>(thisIntegrateTime)) / CLOCKS_PER_SEC)
                    << "s was over Integrate" << std::endl;
                std::cout << " and "
                << ((static_cast<double>(thisRGenTime)) / CLOCKS_PER_SEC)
                    << "s was over Poisson" << std::endl;
        }

    } // end loop over SFS data sets

    // summarise results
    outputSummaryResults(filename,
                        resultsMaxLikelihoods, resultsKLDivergence,
                        resultsGrowthCoverageAccumulation,
                        resultsThetaCoverageAccumulation,
                        resultsJointGrowthThetaCoverageAccumulation,
                        true_rho_per_locus, true_growth_rate, true_theta_per_locus,
                        coverageCriteria);

    if (p.NOISY == 10) {    // timing
        accTime = accSeqTime + accNonSeqTime;
        std::cout << "finished all, sequence time took "
            << ((static_cast<double>(accSeqTime)) / CLOCKS_PER_SEC)
            << "s" << std::endl;
        std::cout << "non sequence time took "
            << ((static_cast<double>(accNonSeqTime)) / CLOCKS_PER_SEC)
            << "s" << std::endl;
        std::cout << "total time took "
            << ((static_cast<double>(accTime)) / CLOCKS_PER_SEC)
            << "s" << std::endl;
        std::cout << "of which  "
            << ((static_cast<double>(accPeesTime)) / CLOCKS_PER_SEC)
            << "s was over Pees" << std::endl;
        std::cout << "of which  "
            << ((static_cast<double>(accLineageTime)) / CLOCKS_PER_SEC)
            << "s was over Lineage" << std::endl;
        std::cout << "of which  "
            << ((static_cast<double>(accIntegrateTime)) / CLOCKS_PER_SEC)
            << "s was over Integrate" << std::endl;
        std::cout << "of which  "
            << ((static_cast<double>(accRGenTime)) / CLOCKS_PER_SEC)
            << "s was over poisson random number generator" << std::endl;
    }



    gsl_rng_free(rgsl); //------------------------------- freeing the RNG
}
