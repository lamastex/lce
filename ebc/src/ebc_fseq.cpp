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
/*! \file ebc_fseq.cpp
    \brief methods for unvintaged (sized) n-coalescent and associated shape statistics
*/


#include <ebc_fseq.hpp>



#define SPLITS_DEBUG 0
#define SPLITS_DEBUG1 0
#define TREEOUT_DEBUG 0
#define TIMESOUT_DEBUG 1
#define FSEQOUT_DEBUG 1


/*! generate fsequences of Kingman's unlabeled n-coalescent (i.e. the unvintaged and sized n-coalescent) forward in time

*/



/*! generate fsequences of Kingman's unlabeled n-coalescent (i.e. the unvintaged and sized n-coalescent) forward in time
*/
void GenerateFSequence(gsl_rng* rgsl, fseq_type& fs, double& Probfs, double& PropProbfs, size_t nsam)
{
  if(nsam < 2) {std::cerr << "\nnsam must be > 1 in GenerateFSequence\n"; exit(0);}
  //Row indices run from 0 to n-2 --
  //0 --> 2 lineages at EPOCH 2, ..., n-2 --> n lineages at EPOCH n
  size_t NROWS=(nsam-1);
  // Column represent the frequency of the leaves subtended by current lineages --
  // 0 --> 1 leaf ,..., n-2 --> n-1 leaves
  size_t NCOLS=(nsam-1);
  //Total length of an fseqn as an (nsam-1)X(nsam-1) matrix in vector form
  size_t NELMS=(NROWS*NCOLS);
  //tracks a valarray of NROWS many indices 0,1,...,NROWS-1
  //for easy access of elements of fseq fs via valarray vecortized index calls
  std::valarray<size_t> fsIndexVecCurr(NROWS);
  for(size_t i=0; i < NROWS; i++) fsIndexVecCurr[i]=i;
  //initialize the fseq fs as 0.0's
  fs = 0.0;
  if(fs.size() != NELMS) fs.resize(NELMS,0.);//ensure fs is the right size

  //EPOCH 2 from EPOCH 1 with f_1 = (0 0 0 0 ... 0 nsam) -- {we don't track f_1}
  //initial split of the ancestral lineage subtending nsam many leaves
  //into two lineages subtending split1 and split2 many leaves
  size_t split0 = nsam;// # leaves subtended by the lineage that's about to be split
  size_t split1 = 1 + gsl_rng_uniform_int (rgsl, split0-1); // size of the left split
  size_t split2 = split0 - split1; // size of the right split
  //proposal probability of the Fseqn -- this will depend on the observed SFS x
  double ProposalFseq = (split1 == split2) ? 1.0/(double(nsam-1)) : 2.0/(double(nsam-1));
  //Kingman Coalescent probability of the Fseqn
  double ProbabKCFseq = (split1 == split2) ? 1.0/(double(nsam-1)) : 2.0/(double(nsam-1));
  //std::cout << "rnd int ::::" << split1 << "\t" << split2 << std::endl;
  //Update Fseq of EPOCH 2 according to the realized split1 and split2
  fs[ fsIndexVecCurr[split1-1] ] += 1.0;
  fs[ fsIndexVecCurr[split2-1] ] += 1.0;
  fs[(nsam-2)*NROWS + split1-1] = 1.0;
  fs[(nsam-2)*NROWS + split2-1] = 1.0;
  //printVA2D(fs,NROWS);
  //---------End of EPOCH 2

  //Begin EPOCH 3, EPOCH 4, ..., EPOCH nsam-1
  //-- EPOCH nsam with prob. 1 is (nsam,0,0,..,0)
  for(size_t c_EPOCH = 2; c_EPOCH < nsam-1; c_EPOCH++)
  {
    //std::cout << "EPOCH : " << c_EPOCH << "  BEGIN " << std::endl;
    //printVA2D(fs,NROWS);
    //preprocess the probabilities in current f-vector
    //to draw the lineage subtending >= 2 leaves to split based on its frequency
    size_t beginIndx = (c_EPOCH-2)*NROWS;
    gsl_ran_discrete_t * fvec_pdfstruct = gsl_ran_discrete_preproc (NROWS-1, &fs[beginIndx+1]);
    size_t draw = gsl_ran_discrete (rgsl, fvec_pdfstruct);
    // # leaves subtended by the lineage that's about to be split
    split0 = 2 + draw;
    ProposalFseq *= (fs[beginIndx+1+draw]/(double(c_EPOCH) - fs[beginIndx])); //normalize the f's by the total number of lineages
    //ProposalFseq *= gsl_ran_discrete_pdf(draw, fvec_pdfstruct);
    //std::cout << "draw = " << draw << "\n";
    //for(size_t iii=0; iii < NROWS-1; iii++)
    //      std::cout << iii << " : " << gsl_ran_discrete_pdf(iii, fvec_pdfstruct) << "\t"; std::cout << "\n";
    //std::cout << "gsl_ran_discrete_pdf = " <<  gsl_ran_discrete_pdf(draw, fvec_pdfstruct) << "\n";
    //std::cout << " ProposalFseq = " << (fs[beginIndx+1+draw]/(double(c_EPOCH) - fs[beginIndx])) << "\n";
    //getchar();
    //assert ( fabs(gsl_ran_discrete_pdf(draw, fvec_pdfstruct) - (fs[beginIndx+1+draw]/(double(c_EPOCH) - fs[beginIndx])) ) < 0.000001 );
    ProbabKCFseq *= fs[beginIndx+1+draw];
    gsl_ran_discrete_free(fvec_pdfstruct);//freeing the pdf struct
    //end of picking the lineage to split
    split1 = 1 + gsl_rng_uniform_int (rgsl, split0-1); // size of the left split
    split2 = split0 - split1; // size of the right split
    //proposal probability of the Fseqn -- this will LATER depend on the observed SFS x
    ProposalFseq *= (split1 == split2) ? 1.0/(double(split0 - 1)) : 2.0/(double(split0 - 1));
    //Kingman Coalescent probability of the Fseqn
    ProbabKCFseq *= (split1 == split2) ? 1.0/(double(nsam - c_EPOCH)) : 2.0/(double(nsam - c_EPOCH));
    //std::cout << "rnd int ::::" << split1 << "\t" << split2 << std::endl;
    //Update Fseq of c_EPOCH + 1 according to the realized split1 and split2
    fs[beginIndx+NROWS+fsIndexVecCurr] = fs[beginIndx+fsIndexVecCurr];//copy f-seq from c_EPOCH to f-seq of c_EPOCH+1
    fs[beginIndx + NROWS + 1 + draw] -= 1.0; // take away one count of the lineage that was just split
    fs[beginIndx + NROWS + split1 - 1] += 1.0;// add one count to account for the left split
    fs[beginIndx + NROWS + split2 - 1] += 1.0;// add one count to account for the right split
    fs[(nsam-2)*NROWS + split1-1] = 1.0;
    fs[(nsam-2)*NROWS + split2-1] = 1.0;
    //std::cout << "EPOCH : " << c_EPOCH << "    END " << std::endl;
    //printVA2D(fs,NROWS);
    //std::cout << "next fv" << std::endl;
    //getchar();
  }
  Probfs = ProbabKCFseq;
  PropProbfs = ProposalFseq;
  //printVA2D(fs,NROWS);
  //std::cout << ProbabKCFseq << "\t" << ProposalFseq << "\n";
  //std::cout << "*******" << std::endl;
}


/*! return (EpochTimes^Transpose * fseq), i.e. the lengths of lineages leading to singleton, doubleton,...,(n-1)-ton mutations

    Input 1: the EpochTimes, where EpochTimes[nsam-2] = t_n = Time to 1st Coal event, ..., EpochTimes[0] = t_2 = Time to last Coal event.
    Input 2: the f-sequence fs under unvintaged and sized (Kingman's unlabeled) n-coalescent.
        Input 3: nsam is the sample size n in the n-coalescent.
        Output: LineageLenItons = (EpochTimes^Transpose * fseq)
*/
std::valarray<double> EpochTimesProdFseq(std::valarray<double> & EpochTimes, fseq_type & fs, size_t nsam)
{
  size_t NROWSCOLS=(nsam-1);
  //tracks a valarray of NROWS many indices 0,1,...,NROWS-1
  //for easy access of elements of fseq fs via valarray vecortized index calls
  std::valarray<size_t> fsIndexVecCurr(NROWSCOLS);
  std::valarray<double> LineageLenItons(NROWSCOLS);
  std::valarray<double> Tmp(NROWSCOLS);
  for(size_t i=0; i < NROWSCOLS; i++) fsIndexVecCurr[i]=i;
  for(size_t c_EPOCH = 2; c_EPOCH < nsam; c_EPOCH++)
  {
    size_t IndexCurr = (c_EPOCH-2)*NROWSCOLS;
    //std::cout << "EPOCH : " << c_EPOCH << "   Times = " << EpochTimes[c_EPOCH-2] << std::endl;
    Tmp = fs[IndexCurr+fsIndexVecCurr];
    Tmp *= (EpochTimes[c_EPOCH-2]);
    LineageLenItons += Tmp;
    //std::cout << "LineageLenItons :\n";printVA1D(LineageLenItons);
  }
  //std::cout << "fseq :\n";printVA2D(fs,NROWSCOLS);
  LineageLenItons[0] += (EpochTimes[nsam-2] * double(nsam));//taking care of the epoch with n lineages with t_n epoch time
  //std::cout << "Epoch Times : \n"; printVA1D(EpochTimes);
  //std::cout << "Final LineageLenItons :\n";printVA1D(LineageLenItons);getchar();
  return LineageLenItons;
}


/*! This returns the shape statistic DefiniteItons (as a vector of unsigned ints) that gives the topology information in the observed sfs

    the shape statistic DefiniteItons gives all information available in the observed sfs regarding the
    equivalence classes of binary coalescent trees upto f-sequences.
    When the i-th element of the sfs X is positive, i.e. X[i]>0, then we know that there had to have been
    a lineage in the underlying coalescent tree (or equivalently in the corresponding f-sequence) that lead to the
    i-ton mutation(s) observed in X[i].  Thus, DefiniteItons is the vector such definite i-ton mutation subtending lineages.
    Of course all of this is only valid in the current formulation when recombination rate is 0 and we have the same realization
    of the binary n-coalescent at all sites.
*/
std::vector<unsigned> SFStoDefiniteItons(const sfs_array_type & X)
{
  std::vector<unsigned> DefiniteItons;
  DefiniteItons.clear();
  unsigned START = 2; // since singletons and doubletons are on every tree with 3 or more leaves
  for (unsigned i=START; i<X.size(); i++)
  {
    //std::cout << " X[" << i << "] = " << X[i] << ' ' ;
    if(X[i] > 0) DefiniteItons.push_back(i+1);
  }
  /*
  std::cout << std::endl;
  for (unsigned i=0; i<DefiniteItons.size(); i++)
      //std::cout <<  DefiniteItons[i] << ' ' ;
      std::cout << " I[" << i << "] = " << DefiniteItons[i] << ' ' ;
  std::cout << std::endl;
  */
  return DefiniteItons;
}

/*! makes the shape statistic DefiniteItons of the observed sfs X as an Indicator vector of doubles of length equal to the length of the sfs X

    Once again this shape statistic DefiniteItons as an Indicator vector of doubles of length equal to the length of the sfs X
    gives the topology information in the observed sfs as explained in SFStoDefiniteItons in a form that is
    more conducive to matrix operations in the sequel.  So, in this form of DefiniteItons, when the observed
    sfs X[i]>0 then the shape statistic DefiniteItons[i]=1.0 else DefiniteItons[i]=0.0;
*/
void SFStoDefiniteItonsDbl(std::valarray<double> & DefiniteItons, const sfs_array_type & X)
{
  DefiniteItons = 0.0;
  if(DefiniteItons.size() != X.size()) DefiniteItons.resize(X.size(),0.);
  for (unsigned i=0; i<X.size(); i++)
  {
    if(X[i] > 0) DefiniteItons[i]=1.0;
    else DefiniteItons[i]=0.0;
  }
}

/*!  The simpleminded O(n) loop-based method to obtain a sample from n bins, based on probability P[]
*/
unsigned SampleBins(std::valarray<double> & p, double r)
{
   unsigned i;
   unsigned Z=p.size();

   for (i=1; i<Z; i++) p[i]+=p[i-1];//make cdf out of pdf
   if(fabs(p[Z-1]-1.0) > 1e-5 ) for (i=0; i<Z; i++) p[i] /= p[Z-1];//Normalize if necessary
   if (fabs(p[Z-1]-1.0) > 1e-5) {std::cerr << "Sum P != 1.. in SampleBins"; exit (-1);}//ensure CDF OK
   for (i=0; i<Z; i++) if (r<p[i]) break;
   return (i);
}


std::valarray<double>& make_StdNeutral_EpochTimes( gsl_rng* rgsl, params& p, std::valarray<double>& EpochTimes)
{
  double rcoal=0.0;//coalescent rate
  double tcoal=0.0;//coalescent time
  unsigned NSAM=0;
  unsigned nsam=p.SfsSz+1;
  unsigned SfsSz=p.SfsSz;

  EpochTimes.resize(SfsSz, 0.0);

  NSAM=nsam;
  while( NSAM > 1 )
    {
      rcoal = double(NSAM*(NSAM-1));
      tcoal = - std::log( gsl_rng_uniform(rgsl) )/rcoal;
      //EpochTimes[nsam-2] = t_n = Time to 1st Coal event, ...,
      //EpochTimes[0] = t_2 = Time to last Coal event
      EpochTimes[NSAM-2]=tcoal;
      NSAM--;
    }
   EpochTimes = EpochTimes * double(p.length); //records the n-1 epoch-times for each replicate

   return EpochTimes; //return by reference
}


/*! make a set of Tees := total tree size = sum of all branch-lengths and Pees := the n-1 epoch-times
    for each replicate simulation of the n-coalescent
*/
void make_set_of_pees_tees( gsl_rng* rgsl,
               Pees & PEES,
               Tees & TEES,
               const double & G,//growth_rate
               params& p)
{
  double t=0.0;//tre height
  double T=0.0;//tree size
  double rcoal=0.0;
  double tcoal=0.0;
  unsigned NSAM=0;
  unsigned nsam=p.SfsSz+1;
  unsigned SfsSz=p.SfsSz;
  std::valarray<double> EpochTimes(0.0, SfsSz); // this keeps track of the pi's for the best arg so far in the M-H sampler

  unsigned MC_ARGTRYNUM=1;
  while(MC_ARGTRYNUM <= p.Num_Of_Replicates)
  {
    t=0.0;
    T=0.0;
    NSAM=nsam;
    while( NSAM > 1 )
      {
    rcoal = double(NSAM*(NSAM-1));
    //tcoal = expo(1./rcoal);
    tcoal = - std::log( gsl_rng_uniform(rgsl) )/rcoal;
    if(G > 0.0) tcoal = std::log(1. + G * std::exp(-G * (t)) * tcoal)/G;
    t += tcoal;
    T += tcoal*(NSAM);
    EpochTimes[NSAM-2]=tcoal;//EpochTimes[nsam-2] = t_n = Time to 1st Coal event, ..., EpochTimes[0] = t_2 = Time to last Coal event
    NSAM--;
      }
    TEES.push_back( (T * double(p.length)) );//records the total tree size = sum of all branch-lengths for each replicate
    PEES.push_back( (EpochTimes * double(p.length)) );//records the n-1 epoch-times for each replicate
    //std::cout << "Epoch times: \n"; printVA1D( std::valarray<double>(EpochTimes * double(p.length)) ); getchar();
    MC_ARGTRYNUM++;
  }
}
