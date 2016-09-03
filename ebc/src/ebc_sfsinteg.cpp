/*
 * Copyright (C) 2005--2009 Raazesh Sainudiin and Kevin Thornton
 *
 * This file is part of lce, a C++ class library for lumped coalescent experiments.
 *
 * lce is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */
/*! \file ebc_sfsinteg.cpp
    \brief methods for constructing an sfs sieve and integrating over it

    We construct the sfs sieve and operate over it.
*/
#include<ebc_sfsinteg.hpp>

/*! checking if the bounds on the number screen are system-specific

    want to make sure that our MOST_NEGA(POSI)TIVE_DOUBLE_ALLOWED larg(small)er than system's number screen limits otherwise reset them
*/
void CheckUserDefinedDoubleTruncation()
{// This function will be called to check the lower/upper bound on doubles MOST_NEGA/POSITIVE_DOUBLE_ALLOWED
  if( !(MOST_NEGATIVE_DOUBLE_ALLOWED > MOST_NEGATIVE_DOUBLE) && (MOST_POSITIVE_DOUBLE_ALLOWED < MOST_POSITIVE_DOUBLE) )
    //std::cerr << "MOST_NEGA(POSI)TIVE_DOUBLE_ALLOWED = " << MOST_NEGATIVE_DOUBLE_ALLOWED << " > "
    //      << MOST_NEGATIVE_DOUBLE << " = " << MOST_NEGATIVE_DOUBLE << std::endl;
  {
    std::cerr << std::endl << "ERROR: unable to set MOST_NEGA(POSI)TIVE_DOUBLE_ALLOWED within " << std::endl
        << "[ " << MOST_NEGATIVE_DOUBLE << " , " << MOST_POSITIVE_DOUBLE << " ] " << std::endl
        << " FIX :  make MOST_NEGA(POSI)TIVE_DOUBLE_ALLOWED larg(small)er !!! " << std::endl;
    exit(1);
  }
}


void PisToMultinomialPs(std::vector<double>& Multinom_Ps, const std::valarray<double>& Pis, double heat, int PRINT )
{//Normalized
  std::valarray<double> Pis_tmp = Pis/(Pis.sum());
  Pis_tmp = Pis_tmp + heat/double(Pis.size()); // heat many units of uniform (1/N,..., 1/N) is added
  double accum = 0.0;
  if( heat<0.) { std::cerr << "heat < 0 encountered in PisToMultinomialPs() " << std::endl; exit(12); }//error check
  for(unsigned i=0; i<Multinom_Ps.size(); i++) {
    Multinom_Ps[i] = Pis_tmp[i];
    if( (Multinom_Ps[i]<=0.00000000000000000001) ) Multinom_Ps[i] = 0.00000000000000000001;
    accum += Multinom_Ps[i];
  }
  for(unsigned i=0; i<Multinom_Ps.size(); i++) Multinom_Ps[i] /= accum;
  if(PRINT) printMultiProb(Multinom_Ps, double(PRINT)); // PRINT is a multiplying factor; if PRINT =1.0 it prints Ps
}

void MultinomialPsToSfs(sfs_array_type& va, int S, const std::valarray<double>& Pis, int PRINT )
{//sfs Need not sum to S -- OK since we fit inside Sieve later
  double accum = Pis.sum();
  double Sd = double(S);
  for(unsigned i=0; i<va.size(); i++) va[i] = (int) ( rint( (Pis[i]*Sd)/accum ) );
  if(PRINT) printValarray1D(va);
}

void SfsToMultinomialPs(const sfs_array_type& va, std::vector<double>& Multinom_Ps, int PRINT )
{//Not Normalized
  for(unsigned i=0; i<va.size(); i++) {
    Multinom_Ps[i] = ( (double) va[i] );
    if( Multinom_Ps[i]<=0.) Multinom_Ps[i]=0.00000000000000000001;
    if( Multinom_Ps[i]<0.)
    {
      std::cerr << "p_i < 0 encountered for index " << i << std::endl;
      exit(10);
    }
  }
  if(PRINT) printMultiProb(Multinom_Ps, double(PRINT)); // PRINT is a multiplying factor; if PRINT =1.0 it prints Ps
}

void MeanPost_UnifDirichletPrior(const sfs_array_type& va, std::vector<double>& Multinom_Ps, double alpha_i, int PRINT )
{//Normalized
  double accum=0.0;
  for(unsigned i=0; i<va.size(); i++) {
    Multinom_Ps[i] = double(va[i]) + alpha_i;
    accum = accum+Multinom_Ps[i];
  }
  for(unsigned i=0; i<va.size(); i++) Multinom_Ps[i] /= accum;
  if(PRINT) printMultiProb(Multinom_Ps, double(PRINT)); // PRINT is a multiplying factor; if PRINT =1.0 it prints Ps
}

/*! build the sfs sieve using a Dirichlet Sampler, i.e. Type_Of_Sieve_Construction=2

 We run a Metropolis-Hastings Markov chain with local sfs proposals from adding random moves
 from Markov Basis with the target stationary distribution on sfs sieve given by the multinomial probabilities in Ps
*/
double build_sieve(gsl_rng* rgsl,
           sfs_sieve & SfsSieve,
           sfs_array_type& vaInit,
           std::vector<OneMove> & MovesVector,
           unsigned NumberOfMoves,
           const std::vector<double> & Ps,
           double * prob_min_sieve,
           const params& p, unsigned NeededSamples, unsigned MaxTrials)
{
  sfs_sieve::iterator SfsSieve_Iter;
#ifndef NDEBUG
    int SsSfs = Ss(vaInit);
    int PISfs = PI(vaInit);
    int Obs_First_Plus_Last = First_Plus_Last(vaInit);
#endif
  unsigned SfsSz = vaInit.size();
  double prob_va=0.;
  double prob_va_Prop, prob_min = *prob_min_sieve;
  sfs_array_type va(SfsSz), vaProp(SfsSz);
  std::vector<double> MPs(SfsSz);
  for(unsigned i=0; i<SfsSz; i++) MPs[i] = Ps[i];


  if (p.NOISY>99) std::cerr << "# NeededSamples = " << NeededSamples << "\tMaxTrials = " << MaxTrials << std::endl;
  unsigned Trials=0;
  unsigned Samples;
  double sfs_integral = 0.0;
  for (Samples = 0; ( (SfsSieve.size() < NeededSamples) && (Trials < MaxTrials) ); Trials++)
  {
    if(Trials==0)
    { //initialize the chain
      va = vaInit;
      // the cast hackery is for GSL compatibility only
      prob_va = gsl_ran_multinomial_lnpdf (SfsSz, &MPs[0],  (unsigned *)(const_cast<int *> (&va[0])));
      SfsSieve.insert(make_pair(va,prob_va));
      Samples++;
    }
    else
    { // propose an sfs
      //std::cerr << "trying to make a move" << std::endl;
      Make_A_Valid_Move(rgsl, va, vaProp, MovesVector, NumberOfMoves);
      //std::cerr << "move proposed" << std::endl; printValarray1D(vaProp);
      prob_va_Prop = gsl_ran_multinomial_lnpdf (SfsSz, &MPs[0], (unsigned *)(const_cast<int *> (&vaProp[0])));
      //Metropolis-Hastings Step
      if ( log(gsl_rng_uniform_pos (rgsl)) <= prob_va_Prop - prob_va )
      { // move accepted
        va = vaProp;
        prob_va = prob_va_Prop;
        SfsSieve.insert(make_pair(va,prob_va));
        Samples++;
        //std::cerr << "move accepted" << std::endl; printValarray1D(va);
      }
    }
    prob_min = (prob_va < prob_min) ? prob_va : prob_min;
    //std::cerr << PI(va) << " " << Ss(va) << " :  "; printValarray1D(va);
    assert( PI(va) == PISfs);
    assert( Ss(va) == SsSfs);
    assert( (p.MoveType == 3 && First_Plus_Last(va) == Obs_First_Plus_Last) || (p.MoveType != 3 ));
    sfs_integral += prob_va;
  }
  *prob_min_sieve = prob_min;

  if(p.NOISY == 1) {
    std::cout << "# sfs Sieve:" << std::endl;
    for ( SfsSieve_Iter = SfsSieve.begin( ) ; SfsSieve_Iter != SfsSieve.end( ) ; SfsSieve_Iter++ ) {
      std::cout << "# ";
      printValarray1Dcout(SfsSieve_Iter->first);
    }
  }
  if(p.NOISY == 2) {
    std::cerr << "# sfs Sieve First and Last Elements:" << std::endl;
    SfsSieve_Iter = SfsSieve.begin( ) ;
    printValarray1D(SfsSieve_Iter->first);
    SfsSieve_Iter = SfsSieve.end( ) ;
    SfsSieve_Iter--;
    printValarray1D(SfsSieve_Iter->first);
  }
  if (p.NOISY>99) std::cerr << "# Trials = " << Trials << std::endl << "# Samples = " << Samples
      << std::endl << "# Sieve Size = " << SfsSieve.size() << std::endl;
  //return (sfs_integral/(double(Trials)));
  return (sfs_integral);
}

/*! build the sfs sieve stochastically with uniform target distribution, i.e. Type_Of_Sieve_Construction=1
*/
double expand_uniformly_sieve(gsl_rng* rgsl,
           sfs_sieve & SfsSieve,
           sfs_array_type& vaInit,
           std::vector<OneMove> & MovesVector,
           unsigned NumberOfMoves,
           const std::vector<double> & Ps,
           double * prob_min_sieve,
           const params& p, unsigned NeededSamples, unsigned MaxTrials)
{
  sfs_sieve::iterator SfsSieve_Iter;
#ifndef NDEBUG
    int SsSfs = Ss(vaInit);
    int PISfs = PI(vaInit);
    int Obs_First_Plus_Last = First_Plus_Last(vaInit);
#endif
  unsigned SfsSz = vaInit.size();
  double prob_va=0.;
  double prob_va_Prop, prob_min = *prob_min_sieve;
  sfs_array_type va(SfsSz), vaProp(SfsSz);
  std::vector<double> MPs(SfsSz);
  for(unsigned i=0; i<SfsSz; i++) MPs[i] = Ps[i];


  if (p.NOISY>99) std::cerr << "# NeededSamples = " << NeededSamples << "\tMaxTrials = " << MaxTrials << std::endl;
  unsigned Trials=0;
  unsigned Samples;
  double sfs_integral = 0.0;
  for (Samples = 0; ( (SfsSieve.size() < NeededSamples) && (Trials < MaxTrials) ); Trials++)
  {
    if(Trials==0)
    { //initialize the chain
      va = vaInit;
      // the cast hackery is for GSL compatibility only
      prob_va = gsl_ran_multinomial_lnpdf (SfsSz, &MPs[0],  (unsigned *)(const_cast<int *> (&va[0])));
      SfsSieve.insert(make_pair(va,prob_va));
      Samples++;
    }
    else
    { // propose an sfs
      //std::cerr << "trying to make a move" << std::endl;
      Make_A_Valid_Move(rgsl, va, vaProp, MovesVector, NumberOfMoves);
      //std::cerr << "move proposed" << std::endl; printValarray1D(vaProp);
      prob_va_Prop = gsl_ran_multinomial_lnpdf (SfsSz, &MPs[0], (unsigned *)(const_cast<int *> (&vaProp[0])));
      // move recorded
      va = vaProp;
      prob_va = prob_va_Prop;
      SfsSieve.insert(make_pair(va,prob_va));
      Samples++;
      //std::cerr << "move accepted" << std::endl; printValarray1D(va);
      }
    }
    prob_min = (prob_va < prob_min) ? prob_va : prob_min;
    //std::cerr << PI(va) << " " << Ss(va) << " :  "; printValarray1D(va);
    assert( PI(va) == PISfs);
    assert( Ss(va) == SsSfs);
    assert( (p.MoveType == 3 && First_Plus_Last(va) == Obs_First_Plus_Last) || (p.MoveType != 3 ));
    sfs_integral += prob_va;

  *prob_min_sieve = prob_min;

  if(p.NOISY == 1) {
    std::cerr << "# sfs Sieve:" << std::endl;
    for ( SfsSieve_Iter = SfsSieve.begin( ) ; SfsSieve_Iter != SfsSieve.end( ) ; SfsSieve_Iter++ ) {
      printValarray1Dcout(SfsSieve_Iter->first);
    }
  }
  if(p.NOISY == 2) {
    std::cerr << "# sfs Sieve First and Last Elements:" << std::endl;
    SfsSieve_Iter = SfsSieve.begin( ) ;
    printValarray1D(SfsSieve_Iter->first);
    SfsSieve_Iter = SfsSieve.end( ) ;
    SfsSieve_Iter--;
    printValarray1D(SfsSieve_Iter->first);
  }
  if (p.NOISY>99) std::cerr << "# Trials = " << Trials << std::endl << "# Samples = "
          << Samples << std::endl << "# Sieve Size = " << SfsSieve.size() << std::endl;
  std::cerr << "\t" << SfsSieve.size() << std::endl;
  //return (sfs_integral/(double(Trials)));
  return (sfs_integral);
}

/*! sieve is built deterministically with uniform target distribution, i.e. Type_Of_Sieve_Construction=0

here prob_va is used to flag whether the sfs has already had its neighbors swept and stored
*/
double graph_sieve(sfs_sieve & SfsSieve,
           sfs_array_type& vaInit,
           std::vector<OneMove> & MovesVector,
           unsigned NumberOfMoves,
           const std::vector<double> & Ps,
           double * prob_min_sieve,
           const params& p, unsigned NeededSamples, unsigned MaxTrials)
{
  sfs_sieve::iterator SfsSieve_Iter;
  //this object encodes the return type of an insert operation done to a map, ie our sfs_filter SFS_Filter
  std::pair<sfs_sieve::iterator, bool>  SieveIterBool;
  //this vector stores the valid neighborhood of an sfs (the neighborhood of x includes x here)
  //std::vector<sfs_array_type> vaNbhd;
  std::vector<sfs_array_type> SfssToProcess;
  std::vector<sfs_array_type> SfssToProcessNext;
#ifndef NDEBUG
    int SsSfs = Ss(vaInit);
    int PISfs = PI(vaInit);
    int Obs_First_Plus_Last = First_Plus_Last(vaInit);
#endif
  unsigned SfsSz = vaInit.size();
  double prob_va=0.;
  double prob_va_Prop, prob_min = *prob_min_sieve;
  sfs_array_type va(SfsSz), vaProp(SfsSz);
  std::vector<double> MPs(SfsSz);
  for(unsigned i=0; i<SfsSz; i++) MPs[i] = Ps[i];

  if (p.NOISY>99) std::cerr << "# NeededSamples = " << NeededSamples << "\tMaxTrials = " << MaxTrials << std::endl;
  unsigned Trials=0;
  unsigned Samples;
  int Nbhd_Sweeps=0;
  bool Nbhd_Exists=true;
  double sfs_integral = 0.0;
  //for (Samples = 0; ( (SfsSieve.size() < NeededSamples) && (Trials < MaxTrials) ); Trials++)
  for (Samples = 0; ( (SfsSieve.size() < NeededSamples) && (Nbhd_Sweeps <= int(MaxTrials)) && Nbhd_Exists ); Trials++)
  {
    if(Trials==0)
    { //initialize the chain
      va = vaInit;
      //make sure Sieve is empty to begin with
      if(!SfsSieve.empty()) {std::cerr << "Sfs Sieve is expected to be empty here ... exiting " << std::endl; exit(9);}
      //store va in SfsSieve for processing -- neighborhood sweeping
      SfssToProcess.push_back(va);
      Nbhd_Exists=true;//va is a neighbor of the empty set -- initial Sfs Sieve
      //Samples++;
    }
    else
    {
      Nbhd_Exists=false;//First assume that SfssToProcess is contained in SfsSieve
      for(unsigned sfs_i=0; sfs_i < SfssToProcess.size(); sfs_i++)
      {
    va = SfssToProcess[sfs_i];
        assert( PI(va) == PISfs);
        assert( Ss(va) == SsSfs);
        assert( (p.MoveType == 3 && First_Plus_Last(va) == Obs_First_Plus_Last) || (p.MoveType != 3 ));
        SieveIterBool = SfsSieve.insert(make_pair(va,prob_va));
    if(SieveIterBool.second)
    { //the insertion into Sieve was successful ==> an sfs whose neighborhood has not been sweept yet
          //The Sieve only contains sfs whose neighborhoods had been swept already
      // find the valid neighboring moves of sfs va and store in SfssToProcessNext
          //std::cerr << "trying to make moves around " << std::endl;
      //printValarray1D(va);
      //vaNbhd.clear();
          //Make_Nbhd_Valid_Moves(va, vaNbhd, MovesVector, NumberOfMoves);
          //std::cerr << "the neighboring moves are " << std::endl;
          //for(unsigned Nbhd_i=0; Nbhd_i < vaNbhd.size(); Nbhd_i++) printValarray1D(vaNbhd[Nbhd_i]);
          //std::cerr << "moves for future processing" << std::endl;
          Make_Nbhd_Valid_Moves(va, SfssToProcessNext, MovesVector, NumberOfMoves);
      Nbhd_Sweeps ++;
          Nbhd_Exists=true;
          //for(unsigned Nbhd_i=0; Nbhd_i < SfssToProcessNext.size(); Nbhd_i++) printValarray1D(SfssToProcessNext[Nbhd_i]);
          //prob_va_Prop = gsl_ran_multinomial_lnpdf (SfsSz, &MPs[0], (unsigned *)(const_cast<int *> (&vaProp[0])));
          //move recorded
    }
        if( p.NOISY>99 && ((Nbhd_Sweeps % 1) == 0)  )
        {
          std::cerr << "\rExpanding Graph : " << Nbhd_Sweeps << "/" << NeededSamples
            << " .  |  " << Samples << std::flush;
    }
      }
      //if(!Nbhd_Exists) std::cerr << "\nIsolated Component!!!" << std::endl;
      Samples += SfssToProcessNext.size();//Samples here is the number of valid neighborhood sweeps up to now
      SfssToProcess.clear(); // This has been processed
      SfssToProcess = SfssToProcessNext; // This has to be processed next
      SfssToProcessNext.clear(); // This has been processed
    }
    prob_min = (prob_va < prob_min) ? prob_va : prob_min;
    sfs_integral += prob_va;
  }
  *prob_min_sieve = prob_min;

  if(p.NOISY == 1) {
    std::cerr << "# sfs Sieve:" << std::endl;
    for ( SfsSieve_Iter = SfsSieve.begin( ) ; SfsSieve_Iter != SfsSieve.end( ) ; SfsSieve_Iter++ ) {
      printValarray1Dcout(SfsSieve_Iter->first);
    }
  }
  if(p.NOISY == 2) {
    std::cerr << "# sfs Sieve First and Last Elements:" << std::endl;
    SfsSieve_Iter = SfsSieve.begin( ) ;
    printValarray1D(SfsSieve_Iter->first);
    SfsSieve_Iter = SfsSieve.end( ) ;
    SfsSieve_Iter--;
    printValarray1D(SfsSieve_Iter->first);
  }
  if (p.NOISY>99) std::cerr << "# Trials = " << Trials << std::endl << "# Samples = "
          << Samples << std::endl << "# Sieve Size = " << SfsSieve.size() << std::endl;


  //std::cerr << "\t" << SfsSieve.size() << std::endl;

  if( p.NOISY>99 )
    std::cerr << "\rExpanding Graph : " << SfsSieve.size() << "/" << NeededSamples << " .  |  " << Samples << std::endl;
  //return (sfs_integral/(double(Trials)));
  return (sfs_integral);
}

bool integrate_over_sieve(sfs_sieve & SfsSieve, const std::valarray<double>& Pis,
        const params& p, double & sfs_sieve_integral, double & sfs_sieve_ShanInfo)
{
  sfs_sieve_integral = 0.;
  sfs_sieve_ShanInfo = 0.;
  sfs_sieve::iterator SfsSieve_Iter;
  for(SfsSieve_Iter=SfsSieve.begin(); SfsSieve_Iter!=SfsSieve.end(); ++SfsSieve_Iter )
  {
    //double MultinomProb = gsl_ran_multinomial_pdf (p.SfsSz, &Pis[0], (unsigned *)(const_cast<int *> (&(SfsSieve_Iter->first)[0])));
    double MultinomProb = exp(Multinomial_Prob_Ln(p.SfsSz, &Pis[0], (unsigned *)(const_cast<int *> (&(SfsSieve_Iter->first)[0]))));
    sfs_sieve_integral += MultinomProb;
    sfs_sieve_ShanInfo += log(MultinomProb);

    //std::cout << MultinomProb << '\n'; getchar();
    //MultinomPROBS_PopMC.push_back( MultinomProb_PopMC );
   }
  //std::cout << std::endl << sfs_sieve_integral << '\t' << sfs_sieve_integral/double(SfsSieve.size()) << '\t' << sfs_sieve_ShanInfo << '\n';
  return true;
}

double Multinomial_Prob_Ln(const size_t K, const double p[], const unsigned int n[])
{
  size_t k;
  unsigned int N = 0;
  double log_pdf = 0.0;
  double norm = 0.0;
  for (k = 0; k < K; k++) N += n[k];
  log_pdf = gsl_sf_lngamma (N + 1);/* Note: n! == gamma(n+1) */
  for (k = 0; k < K; k++)
  {
    log_pdf -= gsl_sf_lngamma (n[k] + 1);
    norm += p[k];
  }
  //the following is different from GSL -- as we only care about the prob of entries with >0 count
  for (k = 0; k < K; k++) if(n[k]>0) log_pdf += log (p[k] / norm) * n[k];
  return log_pdf;
}

double Multinomial_Unordering_Factor_Ln(const size_t K, const unsigned int n[])
{
  size_t k;
  unsigned int N = 0;
  double log_pdf = 0.0;

  for (k = 0; k < K; k++) N += n[k];
  /* Note: n! == gamma(n+1) */
  log_pdf = gsl_sf_lngamma (N + 1);
  for (k = 0; k < K; k++) log_pdf -= gsl_sf_lngamma (n[k] + 1);
  return log_pdf;
}

double Multinomial_Ordered_Lnpdf(const size_t K, const double p[], const unsigned int n[])
{
  size_t k;
  double log_pdf = 0.0;
  double norm = 0.0;
  for (k = 0; k < K; k++) norm += p[k];
  for (k = 0; k < K; k++) if(n[k]>0) log_pdf += log (p[k] / norm) * n[k];
  return log_pdf;
}



double integrate_over_conditional_sfs(gsl_rng* rgsl, const sfs_array_type& vaInit,
        std::vector<OneMove> & MovesVector, unsigned NumberOfMoves, const std::valarray<double>& Pis,
        const params& p)
{
  sfs_sieve SfsSieve_temp;
  sfs_sieve::iterator SfsSieve_Iter;
  int SsSfs = Ss(vaInit);
  unsigned SsSfsUS = ((unsigned) SsSfs);
#ifndef NDEBUG
    int PISfs = PI(vaInit);
    int Obs_First_Plus_Last = First_Plus_Last(vaInit);
#endif
  unsigned SfsSz = vaInit.size();
  double prob_va, prob_va_Prop, prob_min_sieve;
  sfs_array_type va(SfsSz), vaProp(SfsSz), vaTmp(SfsSz), vaPisToSfs(SfsSz);
  std::vector<double> Multinom_Ps(SfsSz);
  std::vector<double> Multinom_Ps_hot(SfsSz);
  std::vector<double> LogMultinomProbs;
  unsigned MULTINOMIAL_JUMPS=0;
  //double average_tajd=0.0;

  if (p.NOISY>99)
  {
    std::cerr << "Integrating Over conditional Sfs for the simulated ARG and the corresponding Pis " << "lkl_sim" << std::endl;
    std::cerr << "Pis To Multinomial Pis -- COLD weight :" << std::endl;
    //PisToMultinomialPs(Multinom_Ps, Pis, 0.0001, SsSfs );//actually heat is LOW at 0.0001
    PisToMultinomialPs(Multinom_Ps, Pis, 0., SsSfs );//actually heat is LOW at 0.0001
    std::cerr << std::endl << "Pis To Multinomial Pis -- HOT weight :" << std::endl;
    PisToMultinomialPs(Multinom_Ps_hot, Pis, p.heat, SsSfs );
  }
  else
  {
    PisToMultinomialPs(Multinom_Ps, Pis, 0., 0 );//actually heat is LOW at 0.0001 -- or the cold chain drags along -- I should really do simulated tempering here!!!!!!
    PisToMultinomialPs(Multinom_Ps_hot, Pis, p.heat, 0 );
  }
  double sfs_integral=0.0;
  int integration_samples=0;
  for (int Num_Starts=0; Num_Starts <= p.Num_Multinomial_Starts; Num_Starts++ )
  {
    if(Num_Starts== 0)
    { //build Sfs_Sieve_temp based on Multinom_Ps_hot weights to throw multinomial darts for integration when Num_Starts>0
      if(p.NOISY>99)
      {
        std::cerr << "Building the ARG-specific Sfs Sieve of size  = " << p.sieve_temp_NeededSamples << " with heat = " << p.heat << std::endl;
    std::cerr << "and initialized by the following sfs in 'observed-sfs-initialized Sfs Sieve' that is closest to the current 'Pis To Sfs' :" << std::endl;
    printValarray1D(vaInit);
    MultinomialPsToSfs(vaPisToSfs, SsSfs, Pis, 1 );
        std::cerr << std::endl;
      }
      va = vaInit;// in the case of 0 multinomial jumps the initial sfs = vaInit
      build_sieve(rgsl, SfsSieve_temp, va, MovesVector, NumberOfMoves, Multinom_Ps_hot, &prob_min_sieve, p, p.sieve_temp_NeededSamples, p.sieve_temp_MaxTrials);
     //cout << "Initial Hot Sfs Sieve Prob for the arg(s) = " << Initial_Hot_Sieve_Prob << std::endl;
    }
    else
    { //va is the multinomial realization with weights Multinom_Ps
      if(MULTINOMIAL_JUMPS)
      {
    gsl_ran_multinomial (rgsl, SfsSz, SsSfsUS, &Multinom_Ps[0], (unsigned *)(const_cast<int *> (&va[0])));
      }
      else va = vaPisToSfs;
    }

    if(p.NOISY>99) vaTmp=va;
    //initialize va to the closest sfs to the multinomial realization in current Sieve
    if ((SfsSieve_Iter = SfsSieve_temp.lower_bound(va)) == SfsSieve_temp.end()) SfsSieve_Iter--;
    va = SfsSieve_Iter->first;
    if(p.NOISY>99)
    {
      std::cerr << "Multinomial Jump Start Num. = " << Num_Starts;
      std::cerr << " with reinitializing sample sfs, currently-closest sfs, and current 'Pis To Sfs' are :- w/ Hamming Dist = ";
      std::cerr << (abs(va-vaTmp)).sum() << "\t" << (abs(va-vaPisToSfs)).sum() << std::endl;
      printValarray1D(vaTmp);
      printValarray1D(va);
      printValarray1D(vaPisToSfs);
    }
    //double prob_va_closest_MObs = SfsSieve_Iter->second;
    prob_va = gsl_ran_multinomial_lnpdf (SfsSz, &Multinom_Ps[0],  (unsigned *)(const_cast<int *> (&va[0])));
    //printValarray1D(va);

    //for each multinomail starting va integrate Num_Sfs_Integrations times
  int Trials, int_Trials;
    for (Trials=0, int_Trials=0; integration_samples <= p.Num_Sfs_Integrations; Trials++)
    {
      // propose an sfs
      Make_A_Valid_Move(rgsl, va, vaProp, MovesVector, NumberOfMoves);
      //printValarray1D(vaProp);

      //Metropolis-Hastings Step
      prob_va = gsl_ran_multinomial_lnpdf (SfsSz, &Multinom_Ps[0],  (unsigned *)(const_cast<int *> (&va[0])));
      prob_va_Prop = gsl_ran_multinomial_lnpdf (SfsSz, &Multinom_Ps[0], (unsigned *)(const_cast<int *> (&vaProp[0])));
      if ( log(gsl_rng_uniform_pos (rgsl)) <= prob_va_Prop - prob_va )
      {
        //std::cerr << "move accepted" << std::endl; printValarray1D(vaProp);
        va = vaProp;
        prob_va = prob_va_Prop;
        SfsSieve_temp.insert(make_pair(va,prob_va));
      }
      assert( PI(va) == PISfs);
      assert( Ss(va) == SsSfs);
      assert( (p.MoveType == 3 && First_Plus_Last(va) == Obs_First_Plus_Last) || (p.MoveType != 3 ));
      if( (Trials >= p.BURNIN_SFS) && ( (Trials % p.THIN_OUT_SFS) ==0 ) )
      {
        int_Trials++;
        integration_samples++;
        LogMultinomProbs.push_back( exp(prob_va) );
    //std::cout << prob_va << std::endl;
        //average_tajd += tajd_SFS(SsSfsUS, va);
      }
    }//end of for loop in Num_Sfs_Integrations
  }//end of for loop in Num_Multinomial_Starts
  sfs_integral = gsl_stats_mean (const_cast<double *>(&LogMultinomProbs[0]), 1, unsigned(LogMultinomProbs.size()));
  double simstderr = gsl_stats_sd_m(const_cast<double *>(&LogMultinomProbs[0]), 1, unsigned(LogMultinomProbs.size()), sfs_integral);
  if(p.NOISY == 1) {
    std::cerr << "# sfs Sieve temp:" << std::endl;
    for ( SfsSieve_Iter = SfsSieve_temp.begin( ) ; SfsSieve_Iter != SfsSieve_temp.end( ) ; SfsSieve_Iter++ ) {
      printValarray1D(SfsSieve_Iter->first);
    }
  }
  if(p.NOISY == 2) {
    std::cerr << "# sfs Sieve temp's First and Last Elements:" << std::endl;
    SfsSieve_Iter = SfsSieve_temp.begin( ) ;
    printValarray1D(SfsSieve_Iter->first);
    SfsSieve_Iter = SfsSieve_temp.end( ) ;
    SfsSieve_Iter--;
    printValarray1D(SfsSieve_Iter->first);
  }
  if (p.NOISY>99)
  {
    std::cerr << "# Integration Samples = " << integration_samples << "\t";
    std::cerr << "# Sieve temp Size = " << SfsSieve_temp.size() << std::endl;
  }
  SfsSieve_temp.clear();

  //std::cout << sfs_integral << "\t" << log(sfs_integral);
  //std::cout << "\t" << simstderr << std::endl;
  //cout << "average Tajima's D = " << average_tajd/(double(integration_samples)) << std::endl;
  return sfs_integral;
}

//the grid is assumed uniform for the density X, otherwise volume should be a vector of doubles of length = X.size()
void SummarizeDensityGrid(std::vector<double> & X, const double volume, std::vector<unsigned> & MaxIndexX, std::vector<double> & MaxX,
        std::vector<double> & CoverageProb, const unsigned TrueIndexX, std::vector<double> & ShanInfoX, params& p)
{
 if(!CoverageProb.empty()) CoverageProb.clear();
 if(!MaxX.empty()) MaxX.clear();
 if(!MaxIndexX.empty()) MaxIndexX.clear();
 //std::cout << '\n' << "X :" << '\n';
 //std::transform(X.begin(), X.end(), std::ostream_iterator<double>(std::cout, " "), std::bind2nd(std::multiplies<double>(), 1.0));
 unsigned NumThetas = unsigned(p.theta_number_of_points);
 unsigned NumGrowths = unsigned(p.growth_number_of_points);

 std::transform(X.begin(), X.end(), X.begin(), std::bind2nd(std::multiplies<double>(),volume));
 //std::cout << '\n' << "X * volume :" << '\n';
 //std::transform(X.begin(), X.end(), std::ostream_iterator<double>(std::cout, " "), std::bind2nd(std::multiplies<double>(), 1.0));

 double Nc = std::accumulate(X.begin(),X.end(),0.,kahan_sum<double>());
 //std::cout << std::endl << "Kahan Ns = " << Nc << std::endl;

//Alternatively we can use the recursion in GSL's mean algorithm and multiply by X.size()
 //Nc = gsl_stats_mean (const_cast<double *>(&X[0]), 1, unsigned(X.size()));
 //Nc *= double(X.size());
// std::cout << std::endl << "GSL   Ns = " << Ns << std::endl;

//std::cout << '\n' << "X.size() = " << X.size() << '\n';
 std::transform(X.begin(), X.end(), X.begin(), std::bind2nd(std::divides<double>(),Nc));
 //std::transform(X.begin(), X.end(), std::ostream_iterator<double>(std::cout, " "), std::bind2nd(std::multiplies<double>(), 1.0));
 if(p.NOISY==66)
 {
   std::cout << '\n' << "X * volume / Nc :" << '\n';
   for(size_t i=0; i < NumGrowths; i++)
   {
     std::transform(X.begin()+(i*NumThetas), X.begin()+(i*NumThetas + NumThetas),
           std::ostream_iterator<double>(std::cout, " "), std::bind2nd(std::multiplies<double>(), 1.0));
     std::cout << std::endl;
   }
 }
 Nc = std::accumulate(X.begin(),X.end(),0.);
 //std::cout << '\n' << "Sum = " << Nc << " |Sum-1.0|  = " << fabs(Nc-1.) <<'\n';
 assert( fabs(Nc-1.) < 1e-10);

 std::vector<double> Xtheta(NumThetas,0.0);
 for(unsigned G = 0; G < NumGrowths; G++) for(unsigned T = 0; T < NumThetas; T++) Xtheta[T] += X[(G*NumThetas)+T];
 //Nc = std::accumulate(Xtheta.begin(),Xtheta.end(),0.,kahan_sum<double>());
 //std::cout << '\n' << "Sum(Xtheta) = " << Nc << " |Sum-1.0|  = " << fabs(Nc-1.) <<'\n';
 //std::transform(Xtheta.begin(), Xtheta.end(), Xtheta.begin(), std::bind2nd(std::divides<double>(),Nc));
 //CoverageProbShannonInfo(Xtheta, MaxIndexX, MaxX, CoverageProb, TrueIndexX%NumGrowths, ShanInfoX, p);
 CoverageProbShannonInfo(Xtheta, MaxIndexX, MaxX, CoverageProb, TrueIndexX%NumThetas, ShanInfoX, p);

 if(p.NOISY==66)
 {
   std::cout << '\n' << "Xtheta * volume / Nc :" << '\n';
   std::transform(Xtheta.begin(), Xtheta.begin()+NumThetas,
           std::ostream_iterator<double>(std::cout, " "), std::bind2nd(std::multiplies<double>(), 1.0));
   std::cout << std::endl;
   std::cout << "TrueIndex = " << TrueIndexX << '\t' << " TrueIndexX%NumThetas= "
             << TrueIndexX%NumThetas  << " TrueIndexX / NumThetas = " << TrueIndexX / NumThetas << '\n';
 }

 std::vector<double> Xgrowth(NumGrowths,0.0);
 for(unsigned G = 0; G < NumGrowths; G++) for(unsigned T = 0; T < NumThetas; T++) Xgrowth[G] += X[(G*NumThetas)+T];
 CoverageProbShannonInfo(Xgrowth, MaxIndexX, MaxX, CoverageProb, TrueIndexX / NumThetas, ShanInfoX, p);

 if(p.NOISY==66)
 {
   std::cout << '\n' << "  Xgrowth * volume / Nc :" << '\n';
   std::transform(Xgrowth.begin(), Xgrowth.begin()+NumGrowths,
           std::ostream_iterator<double>(std::cout, " "), std::bind2nd(std::multiplies<double>(), 1.0));
   std::cout << std::endl;
 }

 CoverageProbShannonInfo(X, MaxIndexX, MaxX, CoverageProb, TrueIndexX, ShanInfoX, p);

}

void CoverageProbShannonInfo(const std::vector<double> & X, std::vector<unsigned> & MaxIndexX, std::vector<double> & MaxX,
        std::vector<double> & CoverageProb, const unsigned TrueIndexX, std::vector<double> & ShanInfoX, params& p)
{
  double Nc=0.;
  std::multimap<double, std::pair<double,size_t> > Xmulmap; // pdf, <cdf, OriginalIndex>


  if(!Xmulmap.empty()) Xmulmap.clear();
  for(size_t i = 0; i < X.size(); i++) Xmulmap.insert(std::make_pair(X[i], ( std::make_pair(0.0,i) ) ));

  std::multimap<double, std::pair<double,size_t> >::iterator mpos;
  for(Nc=0., mpos = Xmulmap.begin(); mpos != Xmulmap.end(); ++mpos)
  {
    Nc += mpos->first; (mpos->second).first = Nc;
    unsigned index = (mpos->second).second;
    if(index == TrueIndexX)
    {
      CoverageProb.push_back(std::max(0., (1.0 - (mpos->second).first) ));
      //std::cout << mpos->first << ' ' << index << ' ' << CoverageProb << '\n';
    }
  }
  mpos--;
  //std::cout << "argmax_i X = " <<
  (MaxIndexX.push_back((mpos->second).second));
  //std::cout << "max_i X = " <<
  MaxX.push_back((mpos)->first);
  //std::cout << std::endl;

  //Shannon's Entropy Computation
  double logF=0.;
  for(Nc=0., mpos = Xmulmap.begin(); mpos != Xmulmap.end(); ++mpos)
  {
    if( (logF = log(mpos->first) ) < MOST_NEGATIVE_DOUBLE_ALLOWED)
    {
      if(p.NOISY>1000) std::cerr << std::endl << "WARNING log(prob) too small for screen !!! resetting to "
                       << MOST_NEGATIVE_DOUBLE_ALLOWED << std::endl;
      logF=MOST_NEGATIVE_DOUBLE_ALLOWED;
    }
    Nc += (mpos->first * -logF);
  }
  ShanInfoX.push_back(Nc);
}
