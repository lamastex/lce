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
#include<ebc_sfsparticles.hpp>

//This function takes in a vector of SFS (SFS_Particles) and makes each SFS "walk/wander about" according to a M-H chain with
//temperarute-specific hot target (Temperature-heated Pis, our proposal distribution) while keeping track of all the current 
//and visited states.  Then it draws NumberOfParticles many SFS according to the proposal distribution.  Then it reweights them
//according to the importance weights given by the cold target (the actual target we want to sample from) and finally resamples
//NumberOfParticles many SFS from the normalized reweighted importance weights.  Thus, at the end SFS_Particles is distributed
//according to the actual cold target Pis.
void iterate_sfs_particles_Repel(gsl_rng* rgsl, sfs_particlesV& SFS_Particles, unsigned ParticleWalkLength, 
		unsigned NumberOfParticles, double Temperature,
		std::vector<OneMove> & MovesVector, unsigned NumberOfMoves, const std::valarray<double>& Pis,
		const params& p, const SfsInfo& SFSInfo)
{

  //some summaries we might need
  int SsSfs = SFSInfo.S;
#ifndef NDEBUG 
    int PISfs = SFSInfo.Pi;
    int Obs_First_Plus_Last = First_Plus_Last(SFSInfo.SFS);
#endif

  //some variables we will need
  double prob_va, ln_prob_va, ln_prob_va_Prop;
  double prob_va_cold; 
  sfs_array_type va(p.SfsSz), vaProp(p.SfsSz), vaTmp(p.SfsSz), vaPisToSfs(p.SfsSz);
  std::vector<double> Multinom_Ps(p.SfsSz);
  std::vector<double> Multinom_Ps_hot(p.SfsSz);
  unsigned Trials;

  //here we get the hot and cold targets
  if (!p.NOISY) 
  {
    //  This is our actual target distribution -- cold target
    PisToMultinomialPs(Multinom_Ps, Pis, 0., 0 );
    //  This is our proposal target distribution -- hot target
    PisToMultinomialPs(Multinom_Ps_hot, Pis, Temperature, 0 );
    //  We typically use a M-H chain targetting a hotter target
    //  to anneal over the cold target -- seems to help get over 
    //  local minima/minimize wandering into th tais of the actual target
  }
  else 
  {
    std::cerr << "Integrating Over conditional Sfs for the simulated ARG and the corresponding Pis " << "lkl_sim" << std::endl;
    std::cerr << "Pis To Multinomial Pis -- COLD weight -- ACTUAL target:" << std::endl;
    PisToMultinomialPs(Multinom_Ps, Pis, 0., SsSfs );
    std::cerr << std::endl << "Pis To Multinomial Pis -- HOT weight -- PROPOSAL target :" << std::endl;
    PisToMultinomialPs(Multinom_Ps_hot, Pis, Temperature, SsSfs );
  }

  // this object is the instance of a container that stores the lexicographically 
  // sorted SFS, and its (Visit_Counts, importance_weight, proposal_weight)
  sfs_filter SfsFilter;
  sfs_filter::iterator SfsFilter_Iter;
  //NOTE for later-- ideally we should be passing the filter across iterations of this function to minimize Filter
  //reconstruction time from scratch -- but for this to work we need to master smart_pointers so that we can store the same SFS in two containers
  //-- one ordered lexicographically as now and another ordered by the importance_weight.  This way, we can periodically prune the SFS with the lowest
  //importance_weights and prevent the container from reaching too high a size
  
  //We will clear away the filter in this scheme 
  SfsFilter.clear();
  //this object encodes the return type of an insert operation done to a map, ie our sfs_filter SFS_Filter
  std::pair<sfs_filter::iterator, bool>  FilterIterBool;
  
  //for each SFS in SFS_Particles we will do a M-H walk of length ParticleWalkLength and record all the unique SFS states we visit through the map SfsFilter
  for (unsigned particle_index=0; particle_index < SFS_Particles.size(); particle_index++)// could make SFS_Particles.size() smaller in this loop
  {// for each of the particles -- do the following :
    // 1. make sure it is inserted into the filter right away
    //store SFS in va
    va = SFS_Particles[particle_index];
    //compute the proposal target probability at the SFS va
    prob_va = gsl_ran_multinomial_pdf (p.SfsSz, &Multinom_Ps_hot[0],  (unsigned *)(const_cast<int *> (&va[0])));	
    //construct the object CnPwIw of class Count_PRw_ISw WITH (1,prob_va,0.), ie, visit count = 1, prposal weight = prob_va, and importance weight = 0. 
    Count_PRw_ISw CnPwIw(1,prob_va,0.);
    //insert va and its CnPwIw into the SfsFilter
    FilterIterBool = SfsFilter.insert(make_pair(va,CnPwIw));
    if(FilterIterBool.second)
    { // if the insertion were successful then update the importance weight, CnPwIw.ISw to the cold/hot ratio (prob_va_cold/prob_va)
      prob_va_cold = gsl_ran_multinomial_pdf (p.SfsSz, &Multinom_Ps[0],  (unsigned *)(const_cast<int *> (&va[0])));
      (FilterIterBool.first)->second.ISw = (prob_va_cold/prob_va); // update the importance sampling weight when SFS is inserted into map (filter)
    }
    else
    {// otherwise save the computation of (prob_va_cold/prob_va), as it would have already been done during first successful insertion
     // BUT remember to increment the number of visits to this particular SFS -- this information is mostly kept track of for fun/debug/insight for now
      (FilterIterBool.first)->second.Count += 1; // increment the number of visits to this particukar SFS
    }
    //begin a random M-H walk for particle va=SFS_Particles[particle_index] of walk-length ParticleWalkLength with target Multinom_Ps_hot  
    for (Trials=0; Trials <= ParticleWalkLength; Trials++)
    {
      // propose an sfs
      Make_A_Valid_Move(rgsl, va, vaProp, MovesVector, NumberOfMoves); //printValarray1D(vaProp);
      //Metropolis-Hastings Step 
      ln_prob_va = gsl_ran_multinomial_lnpdf (p.SfsSz, &Multinom_Ps_hot[0],  (unsigned *)(const_cast<int *> (&va[0])));
      ln_prob_va_Prop = gsl_ran_multinomial_lnpdf (p.SfsSz, &Multinom_Ps_hot[0], (unsigned *)(const_cast<int *> (&vaProp[0])));
      if ( log(gsl_rng_uniform_pos (rgsl)) <= ln_prob_va_Prop - ln_prob_va ) 
      {
        //std::cerr << "move accepted" << std::endl; printValarray1D(vaProp);
       	va = vaProp;
        ln_prob_va = ln_prob_va_Prop;
        prob_va = exp(ln_prob_va);	
	//we insert the newly visited SFS into our SfsFilter -- same as the initial insertions described earlier
	Count_PRw_ISw CnPwIw(1,prob_va,0.);
        FilterIterBool = SfsFilter.insert(make_pair(va,CnPwIw));
	if(FilterIterBool.second)
	{
          prob_va_cold = gsl_ran_multinomial_pdf (p.SfsSz, &Multinom_Ps[0],  (unsigned *)(const_cast<int *> (&va[0])));
	  (FilterIterBool.first)->second.ISw = (prob_va_cold/prob_va); // update the importance sampling weight when SFS is inserted into map (filter)
	}
	else
	{
	  //(FilterIterBool.first)->second.ISw += 1; // increment the number of visits to this particukar SFS
	  (FilterIterBool.first)->second.Count += 1; // increment the number of visits to this particukar SFS
	}
      } 
      assert( PI(va) == PISfs);
      assert( Ss(va) == SsSfs);
      assert( (p.MoveType == 3 && First_Plus_Last(va) == Obs_First_Plus_Last) || (p.MoveType != 3 ));
    }//end of loop for particlewalklength
  }//end of for loop in particle_index and hence end of "Building" of SfsFilter
  std::cout << "Completed SfsFilter with " << SfsFilter.size() << "  unique states" << std::endl;
    
  std::cout << "Next we extract the SfsFilter::Iterators into a vector and make a vector of Proposal weights for sampling" << std::endl;
  std::vector<double> proposal;// the hot prposal target we have been targeting with the M-H walk before
  std::vector<sfs_filter::iterator> SfsFilterTracker;
  for ( SfsFilter_Iter = SfsFilter.begin( ) ; 
		  SfsFilter_Iter != SfsFilter.end( ) ; ++SfsFilter_Iter ) 
  {
    //std::cout << SfsFilter_Iter->second.Count << "\t"; 
    //std::cout << SfsFilter_Iter->second.PRw << "\t"; 
    //std::cout << SfsFilter_Iter->second.ISw << "\n"; 
    //printValarray1D(SfsFilter_Iter->first);
    SfsFilterTracker.push_back(SfsFilter_Iter);
    proposal.push_back(SfsFilter_Iter->second.PRw);
  }

  std::cout << "Now we sample  " << NumberOfParticles << " many samples from proposal of size " << (size_t)(proposal.size()) << std::endl;
  gsl_ran_discrete_t * proposal_pdfstruct = gsl_ran_discrete_preproc ((size_t)(proposal.size()), &proposal[0]);
  // this vector stores the indices of the proposal which also indexes the vector<SfsFilter::Iterator> SfsFilterTracker
  std::vector<unsigned> ProposedIndices; 
  // the importance weight: the cold (our true and final target) OVER hot target (temperature -specifc proposal)
  std::vector<double> importance_weights; 
  
  //Finally, we resample with importance weight corrections and push them into SFS_Particles if Temperature 
  //is high else we use the samples from the proposal = target directly.
  SFS_Particles.clear();//clear out the previous particles
  while(ProposedIndices.size() < NumberOfParticles)
  {
    unsigned sample_id = gsl_ran_discrete(rgsl, proposal_pdfstruct);
    //std::cout << sample_id << "\t";
    //push back the index into the vector ProposedIndices
    ProposedIndices.push_back(sample_id);
    //When Temperature is large the target is different from the proposal, so we have to save the importance weights for resampling
    if(Temperature > 1.0e-10)
    {
      //find the corresponding iterator of the SfsFilter stored in the vector<SfsFilter::Iterator> SfsFilterTracker
      SfsFilter_Iter = SfsFilterTracker[sample_id];
      //push back the corresponding pre-computed importance weight into the vector importance_weights
      importance_weights.push_back(SfsFilter_Iter->second.ISw);
    }
    //when the proposal is this cold it is basically the target and we can use the samples from the proposal directly
    else SFS_Particles.push_back(  SfsFilterTracker[sample_id]->first );

  }
  std::cout << std::endl << "Sampled : " << ProposedIndices.size() << "\t" << importance_weights.size() << std::endl;
  gsl_ran_discrete_free (proposal_pdfstruct);

  //when Temperature is NOT cold we need to resample from the normalized importance weights as the prposal != target
  if(Temperature > 1.0e-10)
  {
   std::cerr << "making gsl_discrete PDF importance_weights_pdfstruct with " << (size_t)(importance_weights.size()) << std::endl;
    gsl_ran_discrete_t * importance_weights_pdfstruct = gsl_ran_discrete_preproc ((size_t)(importance_weights.size()), &importance_weights[0]);
    while(SFS_Particles.size() < NumberOfParticles)
    {
      unsigned sample_id = gsl_ran_discrete(rgsl, importance_weights_pdfstruct);
      SFS_Particles.push_back(  SfsFilterTracker[ProposedIndices[sample_id]]->first );
    }
    //sfs_integral = gsl_stats_mean (const_cast<double *>(&LogMultinomProbs[0]), 1, unsigned(LogMultinomProbs.size()));
    //double simstderr = gsl_stats_sd_m(const_cast<double *>(&LogMultinomProbs[0]), 1, unsigned(LogMultinomProbs.size()), sfs_integral);
    gsl_ran_discrete_free (importance_weights_pdfstruct);
  }

  // we finally clear out the filter 
  SfsFilter.clear();
  //-- BUT there is a lot more to be gained in recycling the filter itself across iterations of this function
  //But to accomplish this in a manner that will only keep the currently relevant SFS we need smart pointers 
  //that allow simultaneous sorting of SFS lexicographically as well as only keeping the ones with the the 
  //higher weights when some threshold size for the filter is reached... This is for future work
}

