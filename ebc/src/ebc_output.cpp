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
/*! \file ebc_output.cpp
        \brief methods for output management
*/
#include <ebc_output.hpp>
#include <ebc_fseq_general_graph.hpp> // to be able to use general graphs utilities

#include <sstream>  // to be able to manipulate strings as streams

using namespace std;

output_params::output_params( double r, double gr, double t ) : rho_per_locus(r), growth_rate(gr), theta_per_locus(t)
{
}

std::ostream & output_params::output(std::ostream& o) const
{
  o << rho_per_locus << "\t" << growth_rate << "\t" << theta_per_locus ;
  return o;
}

output_Lkl::output_Lkl( double LNL ) : lnl(LNL)
{
}

std::ostream & output_Lkl::output(std::ostream& o)
{
  o << lnl;
  return o;
}

output_LklShi::output_LklShi( double LNL, double SHI ) : lnl(LNL), shi(SHI)
{
}

std::ostream & output_LklShi::output(std::ostream& o)
{
  o << lnl << "\t" << shi;
  return o;
}

output_results::output_results( double r, double gr, double t, double LNL ) : rho_per_locus(r), growth_rate(gr), theta_per_locus(t), lnl(LNL)
{
}

std::ostream & output_results::output(std::ostream& o)
{
  o << rho_per_locus << "\t" << growth_rate << "\t" << theta_per_locus << "\t" << lnl ;
  return o;
}

output_results2::output_results2( double r, double gr, double t, double SHI, double LKL ) : rho_per_locus(r), growth_rate(gr), theta_per_locus(t), shi(SHI), lkl(LKL)
{
}

std::ostream & output_results2::output(std::ostream& o)
{
  o << rho_per_locus << "\t" << growth_rate << "\t" << theta_per_locus << "\t" << shi << "\t" << lkl ;
  return o;
}

/*! print one dimensional sfs_array_type to std::cout*/
void printValarray1Dcout (const sfs_array_type& va)
{
  for(unsigned i=0; i<va.size()-1; i++) {
    std::cout << ( (i>0) ? std::setw(3) : std::setw(1) ) << va[i] << " ";
  }
  std::cout << std::setw(3) << va[va.size()-1] << std::endl;
}

/*! print one dimensional sfs_array_type to some user-specified std::ostream o and return o*/
std::ostream & printValarray1Dostream (std::ostream & o, const sfs_array_type& va)
{
  for(unsigned i=0; i<va.size()-1; i++) {
    o << ( (i>0) ? std::setw(3) : std::setw(1) ) << va[i] << " ";
  }
  o << std::setw(3) << va[va.size()-1] << std::endl;
  return o;
}

/*! print one dimensional sfs_array_type to std::cerr*/
void printValarray1D (const sfs_array_type& va)
{
  for(unsigned i=0; i<va.size()-1; i++) {
    std::cerr << ( (i>0) ? std::setw(3) : std::setw(1) ) << va[i] << " ";
  }
  std::cerr << std::setw(3) << va[va.size()-1] << std::endl;
}

/*! Printing Grid of Parameters from Set Parameters*/
void PrintLoadGridParams(params& p, std::vector<output_params> & output_GridParamsPTP)
{
  std::ofstream PrmOutFileStream;
  std::ifstream PrmInpFileStream;
  bool APPEND=false;
  PrmInpFileStream.open(p.PrmOutFile.c_str(), std::ios::in); // opening OutFile in read mode -- file must exist
  if(!PrmInpFileStream)
  {//file doesn't exist already
    PrmInpFileStream.close(); // closing for now
    std::cerr << "Unable to open PrmOutFile " << p.PrmOutFile << " for reading :( -- may not exist --> creating one now ... " << std::endl;
    PrmOutFileStream.open(p.PrmOutFile.c_str(), std::ios::out); // opening OutFile in append mode -- clean start!!!
    if(!PrmOutFileStream)
    {
      std::cerr << "Unable to open OutFile " << p.PrmOutFile << " for writing :(" << std::endl;
      exit(1001);
    }
    APPEND=true;
  }
  PrmOutFileStream.close(); // closing for now
  if(APPEND) std::cout << std::endl << "Printing Grid of Parameters from Set Parameters" << std::endl;

  for(double rho = p.rho_min; rho <= p.rho_max; rho += p.rho_incr)
  {
   double littler = (rho > 0. && p.length > 1) ?  rho/double(p.length-1) : 0.;


   for( int growth_i = 0 ; growth_i < p.growth_number_of_points ; ++growth_i)
   {

    double growth_rate = (double(growth_i) * p.growth_rate_incr) + p.growth_rate_min;

    //Now we loop over theta
    for( int theta_i = 0 ; theta_i < p.theta_number_of_points ; ++theta_i)
    {
      double theta_per_locus = ( double(theta_i) * p.theta_per_locus_incr) + p.theta_per_locus_min;
      //double theta_per_site = theta_per_locus/(double(p.length));

      output_GridParamsPTP.push_back
      (
        output_params( littler, growth_rate, theta_per_locus)
      );

      if ( APPEND && (output_GridParamsPTP.size() == p.BUFFER_SIZE) )
      {
        //print out contents of output in bursts of BUFFER_SIZE
        PrmOutFileStream.open(p.PrmOutFile.c_str(),std::ios::out|std::ios::app);
        PrmOutFileStream << std::scientific << std::showpoint << std::setprecision(6);
        for(unsigned output_index=0; output_index < output_GridParamsPTP.size(); output_index++)
        {
          PrmOutFileStream << output_GridParamsPTP[output_index] << std::endl;
        }
        PrmOutFileStream.close();
      }
    }// theta loop
   }//growth loop
  }//rho loop
  //print out the remaining crap
  if (APPEND && !output_GridParamsPTP.empty())
  {
    PrmOutFileStream.open(p.PrmOutFile.c_str(),std::ios::out|std::ios::app);
    PrmOutFileStream << std::scientific << std::showpoint << std::setprecision(6);
    for(unsigned output_index=0; output_index < output_GridParamsPTP.size(); output_index++)
    {
      PrmOutFileStream << output_GridParamsPTP[output_index] << std::endl;
    }
    PrmOutFileStream.close();
  }
}

/*! Start a file for likelihoods output.*/
std::string startLikelihoodOutputFile(ofstream& outfile,
                                const int myPid,
                                const bool importanceSampling,
                                const size_t importanceSampleSize,
                                const params& p,
                                double coverageCriteria,
                                const double& true_rho_per_locus,
                                const double& true_growth_rate,
                                const double& true_theta_per_locus)
{

    // make a file for output and give it a header line
    std::ostringstream stm;
    stm << p.LklOutFile;    // start of the filename
    if (importanceSampling) {
        stm << "ImpSamp_" << importanceSampleSize;
    }
    stm << "_n_" << (p.SfsSz+1) << "_d_" << p.NumOfSfsData
                << "_t_" << p.Num_Of_Replicates << ".txt";
    std::string filename = stm.str();

    outfile.open (filename.c_str(), ios::out);
    if (outfile.is_open()) {
        outfile << "Process id " << myPid << endl;
        outfile << "Output for n = " << (p.SfsSz+1) << ", Sfs data set size "
            << p.NumOfSfsData << ", number of replicates per fsequence (pees) "
            << p.Num_Of_Replicates;
        if (importanceSampling) {
            outfile << " with Importance Sampling, sample size "
            << importanceSampleSize << "\n";
        }
        else {
            outfile << " and sampling to make sum(P(fsequence)) = 1" << "\n";
        }
        outfile << "rgls seed is " << p.seed
            << "True rho per locus is " << true_rho_per_locus
            << ", true growth rate is " << true_growth_rate
            << ", true theta per locus is " << true_theta_per_locus
            << "\nGrid is " << (((p.rho_max - p.rho_min)/p.rho_incr) + 1)
                    << " x " << p.growth_number_of_points
                    << " x " << p.theta_number_of_points
                    << " = " << ((((p.rho_max - p.rho_min)/p.rho_incr) + 1)
                            * p.growth_number_of_points
                            * p.theta_number_of_points)
            << "\n";
        outfile << "DATA_id"
            << '\t' << "rep/t"
            << '\t' << "MaxI"
            << '\t' << "Param_rho"
            << '\t' << "Param_gr"
            << '\t' << "Param_th"
            << '\t' << "MaxX"
            << '\t' << "ShanInfoX"
            << '\t' << "Cov[gr]"
            << '\t' << "Cov[th]"
            << '\t' << "Cov[joint]"
            << std::fixed << std::showpoint << std::setprecision(2)
            << '\t' << "Cov[gr]<=" << coverageCriteria
            << '\t' << "Cov[th]<=" << coverageCriteria
            << '\t' << "Cov[joint]<=" << coverageCriteria
            << '\t' << "SfsSieveSize"
            << '\t' << "rho err"
            << '\t' << "gr err"
            << '\t' << "th err"
            << '\t' << "K-L div"        // Kullback-Leibler divergence
            << std::endl;
        outfile.close();
    }
    else {
        cout << "Error with likelihood output file " << filename << endl;
    }

    return filename;
}


/*! Send a map from xstars to ns to fsequence counts to an outfile.*/
void sendToOutfile(const PaddedXStarToN_FSequenceCountMap& paddedXStarMap,
                            const size_t min_n,
                            const size_t max_n,
                            const string& filename)
{

    ofstream outfile;

    outfile.open (filename.c_str(), ios::out);

    if (outfile.is_open()) {

        outfile << "\n\nSummarising fsequence counts for each padded xstar\n\n" << endl;

        outfile << "\tn is" << endl;
        outfile << "padded xstar" << endl;

        // header row
        for (size_t this_n = max_n; this_n >= min_n; this_n--) {
            outfile << "\t" << this_n;
        }
        outfile << endl;

        PaddedXStarToN_FSequenceCountMap::const_iterator pit;
        for (pit = paddedXStarMap.begin(); pit != paddedXStarMap.end(); pit++) {

            string xstar_label = stateToString(pit->first);

            outfile << xstar_label;

            N_FSequenceCountMap::const_reverse_iterator nit;
            for (nit = (pit->second).rbegin(); nit != (pit->second).rend(); nit++) {
                outfile << "\t" << nit->second;

            }
            outfile << endl;
        }

        outfile << "\n" << endl;
    }
    else {
        cout << "Error opening file " << filename << endl;
    }
}

/*! Output results for one SFS to likelihood output file.*/
void outputResultsForOneSFS(const string& filename,
                            const unsigned indexSFS,
                            const unsigned actual_Num_Of_Replicates_per_theta,
                            const vector<unsigned>& maxIndexX,
                            const vector<output_params>& output_GridParams,
                            const vector<double>& maxX,
                            const vector<double>& shanInfoX,
                            const vector<double>& coverageProb,
                            const double coverageCriteria,
                            const unsigned sfsSieveSize,
                            const double true_rho_per_locus,
                            const double true_growth_rate,
                            const double true_theta_per_locus,
                            const double divergenceKL)
{
    // send to the output file we started right at the top
    ofstream outfile;
    outfile.open (filename.c_str(), ios::app);
    if (outfile.is_open()) {
        outfile << std::scientific << std::showpoint << std::setprecision(6);

        outfile << indexSFS
            << '\t' << actual_Num_Of_Replicates_per_theta
            << '\t' << maxIndexX[2]
            << std::fixed << std::showpoint << std::setprecision(4)
            << '\t' << output_GridParams[maxIndexX[2]]
            << std::scientific << std::showpoint << std::setprecision(6)
            << '\t' << maxX[2]
            << '\t' << shanInfoX[2]
            << std::fixed << std::showpoint << std::setprecision(3)
            << '\t' << coverageProb[0]
            << '\t' << coverageProb[1]
            << '\t' << coverageProb[2]
            << '\t' << ( (coverageProb[0] <= coverageCriteria) ? "1" : "0")
            << '\t' << ( (coverageProb[1] <= coverageCriteria) ? "1" : "0")
            << '\t' << ( (coverageProb[2] <= coverageCriteria) ? "1" : "0")
            << '\t' << sfsSieveSize
            << '\t' << (output_GridParams[maxIndexX[2]].rho_per_locus
                            - true_rho_per_locus)
            << '\t' << (output_GridParams[maxIndexX[2]].growth_rate
                            - true_growth_rate)
            << '\t' << (output_GridParams[maxIndexX[2]].theta_per_locus
                            - true_theta_per_locus)
            << std::fixed << std::showpoint << std::setprecision(4)
            << '\t' << divergenceKL
            << std::endl;

        outfile.close();
    }
}

/*! Output a summary of all results to likelihood output file.*/
void outputSummaryResults(const string& filename,
                        const vector<output_params>& resultsMaxLikelihoods,
                        vector<double>& resultsKLDivergence, // gets sorted here
                        const size_t resultsGrowthCoverageAccumulation,
                        const size_t resultsThetaCoverageAccumulation,
                        const size_t resultsJointGrowthThetaCoverageAccumulation,
                        const double true_rho_per_locus,
                        const double true_growth_rate,
                        const double true_theta_per_locus,
                        const double coverageCriteria)
{
    // summarise results
    double sumRhoError = 0.0;
    double sumGrowthError = 0.0;
    double sumThetaError = 0.0;

    double sumSquaredRhoError = 0.0;
    double sumSquaredGrowthError = 0.0;
    double sumSquaredThetaError = 0.0;

    size_t numResults = resultsMaxLikelihoods.size();

    vector<output_params>::const_iterator pit;
    for (pit = resultsMaxLikelihoods.begin(); pit < resultsMaxLikelihoods.end(); pit++) {

        double thisRhoError = ((*pit).rho_per_locus) - true_rho_per_locus;
        double thisGrowthError = ((*pit).growth_rate) - true_growth_rate;
        double thisThetaError = ((*pit).theta_per_locus) - true_theta_per_locus;

        sumRhoError += thisRhoError;
        sumGrowthError += thisGrowthError;
        sumThetaError += thisThetaError;

        sumSquaredRhoError += thisRhoError * thisRhoError;
        sumSquaredGrowthError += thisGrowthError * thisGrowthError;
        sumSquaredThetaError += thisThetaError * thisThetaError;
    }

    if(numResults > 0 ) { // if we have any results

        // work out the bias and root mean squared errors
        double biasRho = sumRhoError / numResults;
        double biasGrowth = sumGrowthError / numResults;
        double biasTheta = sumThetaError / numResults;

        double rmseRho = sqrt((sumSquaredRhoError / numResults));
        double rmseGrowth = sqrt((sumSquaredGrowthError / numResults));
        double rmseTheta = sqrt((sumSquaredThetaError / numResults));

        // work out the Kullback-Leibler divergence stats
        // use percentiles method to interpolate between values either side
        // of desired percentile as Excel does (makes checking easy...)
        double resultsKLDivergence_LQ = 0;
        double resultsKLDivergence_Med = 0;
        double resultsKLDivergence_UQ = 0;

        // this hacked method for quartiles only works if number of results is a mult of 4
        if (resultsKLDivergence.size()%4 == 0) {
            // sort the KL divergences so that we can get order statistics
            sort(resultsKLDivergence.begin(), resultsKLDivergence.end());

            resultsKLDivergence_LQ =
                (resultsKLDivergence[((resultsKLDivergence.size()*25)/100) - 1] +
                resultsKLDivergence[((resultsKLDivergence.size()*25)/100)]*3.0) / 4.0;
            resultsKLDivergence_Med =
                (resultsKLDivergence[((resultsKLDivergence.size()*50)/100) - 1] +
                resultsKLDivergence[((resultsKLDivergence.size()*50)/100)]) / 2.0;
            resultsKLDivergence_UQ =
                (resultsKLDivergence[((resultsKLDivergence.size()*75)/100) - 1]*3.0 +
                resultsKLDivergence[((resultsKLDivergence.size()*75)/100)]) / 4.0;
        }
        // send summary to the output file
        ofstream outfile;
        outfile.open (filename.c_str(), ios::app);
        if (outfile.is_open()) {
            outfile << "\nSummary results:\n"
                << "rmseRho\tbiasRho\trmseTheta\tbiasTheta"
                << "\tcoverage" << coverageCriteria << "Theta"
                << "\trmseGrowth\tbiasGrowth"
                << "\tcoverage" << coverageCriteria << "Growth"
                << "\tcoverage" << coverageCriteria << "Joint"
                << "\tKL_LQ\tKL_Med\tKL_UQ" << endl;
            outfile << std::fixed << std::showpoint << std::setprecision(6);
            outfile << rmseRho << "\t" << biasRho
                << "\t" << rmseTheta << "\t" << biasTheta
                << "\t" << resultsThetaCoverageAccumulation
                << "\t" << rmseGrowth << "\t" <<biasGrowth
                << "\t" << resultsGrowthCoverageAccumulation
                << "\t" << resultsJointGrowthThetaCoverageAccumulation
                << "\t" << resultsKLDivergence_LQ
                << "\t" << resultsKLDivergence_Med
                << "\t" << resultsKLDivergence_UQ
                << endl;
            outfile.close();
        }
    } // end of summary output

}
