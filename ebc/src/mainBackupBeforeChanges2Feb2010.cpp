/*! \file main.cpp
	\brief main function of the program
*/
//Copyright(c) 2005, Raazesh Sainudiin and Kevin Thornton --  Wed Aug 17 23:44:58 EDT 2005
#include<ebc_sfstypes.hpp>
#include<ebc_moves.hpp>
#include<ebc_summaries.hpp>
#include<ebc_output.hpp>
#include<ebc_sfsinteg.hpp>
#include<ebc_sfsdata.hpp>
#include<ebc_params.hpp>
#include<ebc_precomp_TsPs.hpp>
#include<ebc_process_options.hpp>
#include<ebc_fseq.hpp>
#include<main.hpp>

#include <gsl/gsl_sf_gamma.h>

#define MOVE_TYPE 1

using namespace std;


int main(int argc, char **argv)
{

  CheckUserDefinedDoubleTruncation();

  // pretty printing of numbers
  std::cout << std::scientific << std::showpoint << std::setprecision(6);
  std::cout << std::fixed << std::showpoint << std::setprecision(6);

  // Getting options from command line arguments
  params p = process_options(argc,argv);

  ofstream LklOutFileStream;
  /*
  LklOutFileStream.open(p.LklOutFile.c_str(), ios::out); // opening OutFile in append mode -- clean start!!!
  if(!LklOutFileStream)
  {
    std::cerr << "Unable to open OutFile " << p.LklOutFile << " for writing :(" << std::endl;
    exit(1001);
  }
  LklOutFileStream.close(); // closing for now
  */

  //std::vector< output_results2 > output;
  //std::vector< output_LklShi > output;
  // vectors to store various likelihood values
  std::vector< double > LKL_MTr;
  std::vector< double > LKL_MDT;
  std::vector< double > LKL_MPT;
  std::vector< double > LKL_MPD;
  std::vector< double > LKL_Trs;


  OneMove OneMoveTemp(p.SfsSz);// temp variable of the struct OneMove
  // various site frequency spectrum or sfs arrays as valarrays
  sfs_array_type va(p.SfsSz), vaProp(p.SfsSz), vaBest(p.SfsSz), VA(p.SfsSz);
  //vector< vector<double> > MultiProbVector;// stores the multinomial probabilities
  //vector< vector<double> >::iterator MultiProb_Iter;
  sfs_sieve SfsSieve;// the basic "sieve" or set of sfs that satisfy a given set of summaries
  sfs_sieve::iterator SfsSieve_Iter; // the corresponding iterator

  // vector of sfs data (these are typically simulated sfs data -- can also be real data)
  vector<sfs_array_type> DataSfspectra;

  //Read the Sfs data from the DataFile
  ReadSfsData(DataSfspectra, p);

  // vector of moves
  vector<OneMove> MovesVector;
  //Read the moves from the MovesFile
  ReadMoves(MovesVector, p);

  //gsl setup
  gsl_rng * rgsl = gsl_rng_alloc(gsl_rng_mt19937);
  if(!p.USER_SEED) p.seed = time( NULL ); //(unsigned short)time( NULL ) ;
  gsl_rng_set (rgsl, p.seed);
  //std::cout << "SEED = " << p.seed << std::endl;

  double prob_min_sieve = 1e300;
  vector<unsigned> vaUS(p.SfsSz);
  vector<unsigned> vaPropUS(p.SfsSz);
  sfs_array_type ObsSfs(p.SfsSz);
  std::vector<output_params> output_GridParams;

  p.NumOfSfsData = (p.NumOfSfsData <= DataSfspectra.size() ) ? p.NumOfSfsData : DataSfspectra.size();

    const int nsam = p.SfsSz+1;    //sample size to simulate

    //std::cerr << "about to Generate Fseqns"<< std::endl;
    typedef std::map<fseq_type, Count_PRw_ISw, LexicoSorting<fseq_type> > fseq_map;
    fseq_map an_fseq_map;
    fseq_map::iterator fseq_map_Iter;
    std::pair<fseq_map::iterator, bool>  fseq_map_Bool;
    fseq_type fs((nsam-1)*(nsam-1));
    double Probfs = 0.0;
    double PropProbfs = 0.0;
    double TotalProbfs = 0.0;
    Count_PRw_ISw CnPwIw(1, Probfs, PropProbfs);
    int iii=0;
    for(iii=0; (iii< 10000) && (TotalProbfs < 0.999999999); iii++) {
      //std::cout << iii << '\n';
      Probfs = 0.0;
      PropProbfs = 0.0;
      GenerateFSequence(rgsl, fs, Probfs, PropProbfs, nsam);
      fseq_map_Bool = an_fseq_map.insert(make_pair(fs, CnPwIw));
      if(fseq_map_Bool.second) {
	//have successfully inserted a new f-seqn and update its KC_Prob
        (fseq_map_Bool.first)->second.PRw = Probfs;
        (fseq_map_Bool.first)->second.ISw = PropProbfs;
	TotalProbfs += Probfs;//increment the total prob of all f-seqns visited
      }
      //increment the number of visits to this particular f-seqn
      else (fseq_map_Bool.first)->second.Count += 1;
    }

    //std::cout << "NumOfSfsData = " << p.NumOfSfsData << '\n';
  for (unsigned CNTR=0; CNTR < p.NumOfSfsData; CNTR++) //loop over the sfs data sets
  {
    std::cout << std::scientific << std::showpoint << std::setprecision(6);
    std::cerr << std::scientific << std::showpoint << std::setprecision(6);
    //Load and Print the Parameter Grid during first run only:
    if(CNTR==0)  PrintLoadGridParams(p, output_GridParams);
    SfsSieve.clear(); // emptying out the Sieve
    //Load the Observed Sfs and Manipulate it
    ObsSfs = DataSfspectra[CNTR];
    // used only for initial sieve construction -- basically heating the observed (MLE) multinomail probs
    //std::vector<double> ObsMultinomPs(p.SfsSz);
    //SfsToMultinomialPs(ObsSfs, ObsMultinomPs, 1 );
    std::vector<double> UnfMultinomPs(p.SfsSz);
    for(unsigned i=0; i < p.SfsSz; i++) UnfMultinomPs[i] = 1.0 / ((double) p.SfsSz ); // uniform distribution
    vector<double> MeanPostDirich(p.SfsSz, p.prior_alpha); // prior alphas in the Dirichlet conjugate to the multinomial pi's
    MeanPost_UnifDirichletPrior(ObsSfs, MeanPostDirich, 1.0, 0);
    // used only for initial sieve construction -- basically heating the observed (MLE) multinomail probs
    //MeanPost_UnifDirichletPrior(ObsSfs, MeanPostDirich, p.heat, 1);
    // posterior (with prior alphas) in the Dirichlet conjugate to the multinomial pi's giving Law of SFS
    vector<double> AlphasPostDirich(p.SfsSz, p.prior_alpha); // prior alphas in the Dirichlet conjugate to the multinomial pi's
    SfsInfo SFSI(DataSfspectra[CNTR].size());
    SFSI.Inform(DataSfspectra[CNTR]);
    //Print the Observed Statistics and sfs and Pis
    if(p.NOISY!=0 && p.NOISY!=56){
	    std::cout << SFSI;
    }
    //cout << SFSI.S << "\t" << SFSI.Pi << "\t" << SFSI.D << endl;
    /////////////////////////////////////////////////////////////--BEGIN SFS SIEVE-BUILDING (DETERMINISTIC/STOCHASTIC)--///////////////////////

    // building the SfsSieve by starting the chain at the observed sfs

    double Observed_Sieve_Prob=0.;

    switch(p.sieve_Construct_Type)
    {
      case 2 :
	Observed_Sieve_Prob = build_sieve(rgsl, SfsSieve, ObsSfs, MovesVector, MovesVector.size(),
		    				MeanPostDirich, &prob_min_sieve, p, p.sieve_NeededSamples, p.sieve_MaxTrials );
		    				//ObsMultinomPs, &prob_min_sieve, p, p.sieve_NeededSamples, p.sieve_MaxTrials );
	break;

      case 1 :
        Observed_Sieve_Prob = expand_uniformly_sieve(rgsl, SfsSieve, ObsSfs, MovesVector, MovesVector.size(),
		    				UnfMultinomPs, &prob_min_sieve, p, p.sieve_NeededSamples, p.sieve_MaxTrials );
	break;

      case 0 :
        Observed_Sieve_Prob = graph_sieve(SfsSieve, ObsSfs, MovesVector, MovesVector.size(),
		    				UnfMultinomPs, &prob_min_sieve, p, p.sieve_NeededSamples, p.sieve_MaxTrials );
	break;
      case -1 :
	//not making sieve at all
	break;

      default :
	std::cerr << "ERROR: p.sieve_Construct_Type = 0,1,2 or -1  BUT it is " << p.sieve_Construct_Type << std::endl;
	exit(0);
	break;
    }

    if(p.NOISY>99) {
	    cout << "# built sieve SfsSieve by method " << p.sieve_Construct_Type << std::endl
		 << " with average Multinomial Prob = " << Observed_Sieve_Prob  << "\t and Poisson Multinomial Prob = "
	         << (Observed_Sieve_Prob + log(gsl_ran_poisson_pdf (SFSI.uS, double(SFSI.uS) ))) << std::endl;
    }

    unsigned SfsSieveSize = SfsSieve.size();
    //getchar();
    //continue;
    /////////////////////////////////////////////////////////////--END OF SFS SIEVE-BUILDING (DETERMINISTIC/STOCHASTIC)--///////////////////////

    /////////////////////////////////////////////////////////////--BEGIN COALECENT SIMULATINS OF GENEALOGIES--///////////////////////


    //std::cout << SFSI;
    std::vector<unsigned> DefiniteItons = SFStoDefiniteItons(ObsSfs);
    std::valarray<double> DefiniteItonsDbl(p.SfsSz);
    SFStoDefiniteItonsDbl(DefiniteItonsDbl, ObsSfs);
    //printVA1D(DefiniteItonsDbl);

    const unsigned Num_Of_Replicates_per_theta = p.MAX_RS_NUM;
    unsigned Num_Of_EpochTimes=0;
    unsigned Actual_Num_Of_Replicates_per_theta=0;

    double volume = p.theta_per_locus_incr * p.growth_rate_incr;
    std::vector<double> ShanInfoX;
    unsigned TrueIndexX=41;
    std::vector<double> CoverageProb;
    std::vector<double> MaxX;
    std::vector<unsigned> MaxIndexX;

    double growth_rate;
    double theta_per_site;
    double theta_per_locus;
    double rho;   //rho=4Nr = recombination rate across region
    double littler;

    Tees TEES;
    Pees PEES;
    std::vector<double> ImpWeightsD;
    double ImpWtC;
    double Theta_dependent_Poisson_Prob = 0.0;
    double ln_Theta_dependent_Poisson_Prob = 0.0;
    double ln_sieve_integral = 0.0;
    double ln_sieve_integral_Poisson_Prob = 0.0;
    double avgprob = 0.0;


    typedef std::vector< fseq_type > fseq_vec;
    typedef std::vector< Count_PRw_ISw > fseq_CPI;
    fseq_vec::iterator fseq_vec_Iter;
    fseq_vec an_fseq_vec;
    fseq_CPI an_fseq_CPI;
    std::valarray<size_t> fsidx(nsam-1);
    std::valarray<double> fsx(nsam-1);
    double DfISum=DefiniteItonsDbl.sum();
    for(size_t i=0; i < size_t(nsam-1); i++) fsidx[i]=(i+((nsam-2)*(nsam-1)));

    //std::cout << "about to print data-compatible fseqn's from an fseq_map of size " << an_fseq_map.size() << "  effort " << iii << "  TotProb " << TotalProbfs <<"\n"; getchar();
    for(fseq_map_Iter=an_fseq_map.begin( ); fseq_map_Iter != an_fseq_map.end( ); ++fseq_map_Iter )
    {
      fs = fseq_map_Iter->first;
      fsx = fs[fsidx];//std::cout << "fsx \n";printVA1D(fsx);
      //std::cout << "DfI \n"; printVA1D(DefiniteItonsDbl);
      fsx *= DefiniteItonsDbl;//std::cout << "fsxDfI \n";printVA1D(fsx);
      if(fsx.sum() == DfISum)
      {
	an_fseq_vec.push_back(fs);
	an_fseq_CPI.push_back(fseq_map_Iter->second);
	/*
        std::cout << "fs  \n";printVA2D( fs, (nsam-1));
        std::cout << "\t\t"<< fseq_map_Iter->second.PRw << "\t" << fseq_map_Iter->second.ISw
		  << "\t" << double(fseq_map_Iter->second.Count)/double(iii) << "\n";
		  */
      }
      //else {std::cout << "no match : fsxDfI \n"; printVA1D(fsx); std::cout << "fs  \n";printVA2D( fs, (nsam-1)); }
    }
    /*
    for(size_t vi=0; vi < an_fseq_vec.size(); vi++)
    {
      std::cout << "fs  \n";printVA2D(an_fseq_vec[vi], (nsam-1));
      std::cout << "\t\t"<< an_fseq_CPI[vi].PRw << "\t" << an_fseq_CPI[vi].ISw
	                        << "\t" << double(an_fseq_CPI[vi].Count)/double(iii) << "\n";
    }*/

    //std::cerr << "about to enter loops "<< std::endl;
    for(rho = p.rho_min; rho <= p.rho_max; rho += p.rho_incr)
    {
      //std::cerr << "inside rho-loop with rho = " << rho << std::endl;
      littler = (rho > 0. && p.length > 1) ?  rho/double(p.length-1) : 0.;
      for( int growth_i = 0 ; growth_i < p.growth_number_of_points ; ++growth_i)
      {
	growth_rate = (double(growth_i) * p.growth_rate_incr) + p.growth_rate_min;
        //std::cerr << "inside growth-loops with growth = " << growth_rate << std::endl;
        unsigned TeesPeesSz;
        if(!TEES.empty()) TEES.clear();
        if(!PEES.empty()) PEES.clear();
	if(!ImpWeightsD.empty()) ImpWeightsD.clear();
	ImpWtC = 0.0;

	//Get the total tree size in TEES[i] and the epoch-times in PEES[i]
	make_set_of_pees_tees(rgsl, PEES, TEES, growth_rate, p);
        if(TEES.size()==PEES.size()) TeesPeesSz=PEES.size();
        else {std::cerr << std::endl << "Tees and Pees not same size!" << std::endl; exit(1);}
	//cout << "Size of TEES" << std::endl << TEES.size() << std::endl;
        //cout << "TEES we have" << std::endl; for(unsigned i=0; i<TeesPeesSz; i++) { cout << (TEES[i]) << "\t";  } cout << std::endl;
	//double integral_on_Average_arg=0.0;
	double sieve_integral_on_Average_arg=0.0;
	double sieve_ShanInfo_on_Average_arg=0.0;
  	valarray<double> MH_Average_pis(0.0,p.SfsSz);
  	double MH_Average_TotTree = 0.0;

	//cout << "theta_number_of_points "<< p.theta_number_of_points<<"\n";
        for( int theta_i = 0 ; theta_i < p.theta_number_of_points ; ++theta_i)
        {
	  theta_per_locus = (double(theta_i) * p.theta_per_locus_incr) + p.theta_per_locus_min;
	  theta_per_site = theta_per_locus/(double(p.length));
	  //thet_per_site specific Poisson(theta_per_site*TEES[i])-weighted averaging of TEES and PEES
	    double sieve_integral_on_curr_arg=0.0;
	    double sieve_integral_on_all_args=0.0;
	    double sieve_ShanInfo_on_curr_arg=0.0;
	    double sieve_integral_on_curr_fseq=0.0;
	    double ProbfSeq = 0.0;
	    Actual_Num_Of_Replicates_per_theta=0;

	    for(size_t vi=0; vi < an_fseq_vec.size(); vi++)
	    {
	      ProbfSeq = an_fseq_CPI[vi].PRw;
	      sieve_integral_on_curr_fseq=0.0;
	      Num_Of_EpochTimes=0;
	      fseq_type fSeqCurr = an_fseq_vec[vi];
	      valarray<double> LineageLenItons(nsam-1);
	      for(size_t sampled_TeesPees=0; sampled_TeesPees< PEES.size(); sampled_TeesPees++)
	      {
	        LineageLenItons = EpochTimesProdFseq(PEES[sampled_TeesPees], fSeqCurr, nsam);
		//printVA1D(LineageLenItons);
	        integrate_over_sieve(SfsSieve, LineageLenItons, p, sieve_integral_on_curr_arg, sieve_ShanInfo_on_curr_arg);
	        Theta_dependent_Poisson_Prob = gsl_ran_poisson_pdf(SFSI.uS,(theta_per_site*TEES[sampled_TeesPees]));
	        ln_Theta_dependent_Poisson_Prob = std::log(Theta_dependent_Poisson_Prob);
	        ln_sieve_integral = log(sieve_integral_on_curr_arg);
	        ln_sieve_integral_Poisson_Prob = ln_sieve_integral+ln_Theta_dependent_Poisson_Prob;
	        if( (ln_Theta_dependent_Poisson_Prob < MOST_NEGATIVE_DOUBLE_ALLOWED)
	  	  || (ln_sieve_integral < MOST_NEGATIVE_DOUBLE_ALLOWED)      )
	        {
	          if(p.NOISY>1000) std::cerr << std::endl << "WARNING log(prob) too small for screen !!! resetting to "
			            << MOST_NEGATIVE_DOUBLE_ALLOWED << std::endl;
	          ln_sieve_integral_Poisson_Prob = MOST_NEGATIVE_DOUBLE_ALLOWED;
	        }
	        sieve_integral_on_curr_fseq += std::exp(ln_sieve_integral_Poisson_Prob);
		Num_Of_EpochTimes++;
		Actual_Num_Of_Replicates_per_theta++;
	        //std::cout << std::endl << "sieve_integral_on_curr_fseq    = " << sieve_integral_on_curr_fseq << std::endl;
	      }
	      sieve_integral_on_curr_fseq /= double(Num_Of_EpochTimes);//MC integration over Num_Of_EpochTimes many EpochTimes
	      sieve_integral_on_all_args += (sieve_integral_on_curr_fseq * ProbfSeq);
	    }// Replicates_per_theta loop
	    sieve_integral_on_all_args /= double(Num_Of_EpochTimes);
            LKL_Trs.push_back(sieve_integral_on_all_args);
	  if(p.NOISY==56 )
	  {
            std::cerr << "\rgrowth = " << growth_rate << "  theta = " << theta_per_locus
	       	      << "\t" << TeesPeesSz << "  Replicate " << (CNTR+p.DATA_id) << std::flush;
	  }
	}// theta loop
      }//growth loop
    }//rho loop

    //print out the remaining crap
    /*
     if (!output.empty())
      {
	LklOutFileStream.open(p.LklOutFile.c_str(),ios::app);
  	LklOutFileStream << std::scientific << std::showpoint << std::setprecision(6);
	for(unsigned output_index=0; output_index < output.size(); output_index++)
	  {
	    LklOutFileStream << output[output_index] << std::endl;
	  }
	LklOutFileStream.close();
      }
      */
    TrueIndexX=11;// a 101 X 101 theta X growth_rate grid with true parameter 10,0
    /*
    SummarizeDensityGrid(LKL_MTr, volume, MaxIndexX, MaxX,CoverageProb, TrueIndexX, ShanInfoX, p);
    //std::cout << "argmax_i X " << '\t' << "max_i X " << '\t' << "Coverage " <<
	//      << '\t' << "ShanInformation " << '\n';
    //std::cout << " : aVgs : " << '\t'
    std::cerr << std::endl;
    //std::cout << std::endl;
    std::cout << CNTR+p.DATA_id << '\t' << "MTr" << '\t' << p.Num_Of_Replicates
	      << '\t' << MaxIndexX[2]
	      << '\t' << output_GridParams[MaxIndexX[2]] << '\t' << MaxX[2] << '\t' << ShanInfoX[2] << '\t'
              << std::fixed << std::showpoint << std::setprecision(2)
	      << CoverageProb[0] << ' ' << CoverageProb[1] << ' ' << CoverageProb[2] << ' '
	      << ( (CoverageProb[0] < 0.95) ? " Y" : " N")
              << ( (CoverageProb[1] < 0.95) ? " Y" : " N") << ( (CoverageProb[2] < 0.95) ? " Y" : " N");
    std::cout << std::scientific << std::showpoint << std::setprecision(6);
    std::cout << '\t' << SfsSieveSize << std::endl;
    */
    if(Actual_Num_Of_Replicates_per_theta > 0 )
    {
      SummarizeDensityGrid(LKL_Trs, volume, MaxIndexX, MaxX,CoverageProb, TrueIndexX, ShanInfoX, p);
      std::cout << CNTR+p.DATA_id
	      //<< '\t' << "Trs"
	      << '\t' << Actual_Num_Of_Replicates_per_theta
	      << '\t' << MaxIndexX[2]
	      << '\t' << output_GridParams[MaxIndexX[2]] << '\t' << MaxX[2] << '\t' << ShanInfoX[2] << '\t'
              << std::fixed << std::showpoint << std::setprecision(2)
	      << CoverageProb[0] << ' ' << CoverageProb[1] << ' ' << CoverageProb[2] << ' '
	      << ( (CoverageProb[0] < 0.95) ? " 1" : " 0")
              << ( (CoverageProb[1] < 0.95) ? " 1" : " 0") << ( (CoverageProb[2] < 0.95) ? " 1" : " 0");
    std::cout << std::scientific << std::showpoint << std::setprecision(6);
    std::cout << '\t' << SfsSieveSize << std::endl;
    }
    LKL_MTr.clear(); LKL_Trs.clear();
  }
  gsl_rng_free(rgsl); //------------------------------- freeing the RNG
}
