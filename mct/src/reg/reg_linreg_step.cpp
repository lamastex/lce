/* MCT - Markov Chains on Trees.

   Copyright (C) 2011 Jenny Harlow, Kevin Thornton
   Copyright (C) 2008-2010 Kevin Thornton	

   This program is closely based on the regression code in ABCreg,
   Thornton, K.R. (2009) Automating approximate Bayesian computation
   by local linear regression.  BMC Genetics 10: 35.	

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
\brief LinearRegressionStep definitions.
*/

//Headers for this project
#include "reg_linreg_step.hpp"

 //C++ headers
#include <cmath>
#include <sstream>
#include <cstdlib>
#include <algorithm>
#include <utility>
#include <functional>
#include <limits>
#include <stdexcept>
#include <iostream>
#include <cassert>

// gsl headers
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_blas.h>

// Boost headers
#include <boost/lambda/bind.hpp>
#include <boost/lambda/lambda.hpp>

#include <iterator>

//#define MYDEBUG
//#define DETAILDEBUG


using namespace reg;

LinearRegressionStep::LinearRegressionStep(const boost::shared_ptr < Transformer >& _t)
: trans(_t) {}

			
LinearRegressionStep::~LinearRegressionStep() {std::cout << "Linreg destroyed" << std::endl;}
		
//std::vector < boost::shared_ptr < const aabc::ParameterParticleSet > >
void			LinearRegressionStep::linreg(
				const boost::shared_ptr < const aabc::ParameterParticleSet>& pset,
				const boost::shared_ptr <mct::SummaryStatisticSet>& ss_set,
				const boost::shared_ptr <mct::SummaryStatisticSet>& ref_stat_set,
				const double tolerance,
				const std::string& basename)
{
	std::vector< std::vector <double> > prior = *(pset->getAllValues());
	std::vector< std::vector <double> > summaries = *(ss_set->getAllValues());
	
	size_t nparams = prior.size();
	size_t nsumm = summaries.size();
	
	std::string msg = std::string("LinearRegressionStep::linreg")
		+ std::string("(const boost::shared_ptr < const aabc::ParameterParticleSet>&, ")
		+ std::string("const boost::shared_ptr <mct::SummaryStatisticSet>&, ")
		+ std::string("const boost::shared_ptr <mct::SummaryStatisticSet>&, ")
		+ std::string("const double, const std::string&)");
	
	if (nparams < 1) {
		throw std::invalid_argument(msg + std::string (":\nnparams < 1") );
	}
			
	if (nsumm < 1) {
		throw std::invalid_argument(msg + std::string (":\nnsumm < 1") );
	}
	
	if (summaries.at(0).size() < 4) {
		throw std::invalid_argument(msg + std::string (
		":\nneed at least 4 simulations from the prior") );
	}
	
	#ifdef MYDEBUG
	{
		std::cout.precision(10);
			
		std::cout << "the parameters are" <<  std::endl;
		for (std::vector< std::vector < double > >::iterator it = prior.begin();
				it < prior.end();
				++it) {
					std::copy(it->begin(), it->end(), 
						std::ostream_iterator<double>(std::cout, "\t"));
					std::cout <<  std::endl;
				}
				std::cout <<  std::endl;
				
		std::cout << "the summary stats are" <<  std::endl;
		for (std::vector< std::vector < double > >::iterator it = summaries.begin();
				it < summaries.end();
				++it) {
					std::copy(it->begin(), it->end(), 
						std::ostream_iterator<double>(std::cout, "\t"));
					std::cout <<  std::endl;
				}
				std::cout <<  std::endl;
	}
	#endif
	
	//transform the prior
	prior = trans->data_transform(prior);
	
	#ifdef MYDEBUG
		std::cout << "After transformation, prior is" <<  std::endl;
		for (std::vector< std::vector < double > >::iterator it = prior.begin();
				it < prior.end();
				++it) {
					std::copy(it->begin(), it->end(), 
						std::ostream_iterator<double>(std::cout, "\t"));
					std::cout <<  std::endl;
			}
			std::cout <<  std::endl;
	#endif
	
	
	//normalize the simulated summary statistics
	std::vector< std::pair<double,double> > mean_and_sd 
				= scale_simulated_summaries(nsumm, &summaries);

	#ifdef MYDEBUG
	{
		std::cout << "After normalisation of summaries," <<  std::endl;
		for (std::vector< std::pair<double,double> >::iterator it = mean_and_sd.begin();
				it < mean_and_sd.end();
				++it) {
					std::cout << "Mean = " << it->first 
					<< " and sd = " << it->second 
					<<  std::endl;
				}
		std::cout << "and the normalised summary stats are," <<  std::endl;
		for (std::vector< std::vector < double > >::iterator it = summaries.begin();
				it < summaries.end();
				++it) {
					std::copy(it->begin(), it->end(), 
						std::ostream_iterator<double>(std::cout, "\t"));
					std::cout <<  std::endl;
				}
				std::cout <<  std::endl;
	}
	#endif
	

	//go through the observed summary stats,
	//normalizing each along the way and doing the regression,
	//and the output
	
	size_t run = 0;
	
	for (mct::SummaryStatisticSet::const_iterator it = ref_stat_set->begin();
			it < ref_stat_set->end();
			++it) {
		
		
		std::vector<double> observed = ( (*it)->getValues() );
		
		#ifdef MYDEBUG
			std::cout << "Run = " << run << " and observed is " << std::endl;
			std::copy(observed.begin(), observed.end(), 
						std::ostream_iterator<double>(std::cout, "\t"));
			std::cout <<  std::endl;
		#endif
	

		using namespace boost::lambda;
		transform(observed.begin(), observed.end(), mean_and_sd.begin(),
			observed.begin(), _1 -= bind(&std::pair<double, double>::first, _2 ) );
		transform(observed.begin(), observed.end(), mean_and_sd.begin(),
			observed.begin(), _1 /= bind(&std::pair<double, double>::second, _2 ) );
		
		#ifdef MYDEBUG
			std::cout << "After normalisation,  observed is " << std::endl;
			std::copy(observed.begin(), observed.end(), 
						std::ostream_iterator<double>(std::cout, "\t"));
			std::cout <<  std::endl;
		#endif
		
		std::vector< std::vector<double> > posterior = generate_posterior(
										nparams,
										nsumm,
										observed,
										prior,
										summaries,
										tolerance,
										run);
		
		
		#ifdef MYDEBUG
			{
			std::cout << "\nAfter generating posterior,  posterior size is " 
					<< posterior.size() << std::endl;
			std::cout << "and the (transformed) posterior elements are," <<  std::endl;
				for (std::vector< std::vector < double > >::iterator it = posterior.begin();
					it < posterior.end();
					++it) {
					std::copy(it->begin(), it->end(), 
						std::ostream_iterator<double>(std::cout, "\t"));
					std::cout <<  std::endl;
				}
				std::cout <<  std::endl;
			}
		#endif
	
		
		//reverse-transform the posterior
		posterior = trans->data_untransform(posterior);
		
		#ifdef MYDEBUG
			{
			std::cout << "Reverse-transformed posterior elements are," <<  std::endl;
				for (std::vector< std::vector < double > >::iterator it = posterior.begin();
					it < posterior.end();
					++it) {
					std::copy(it->begin(), it->end(), 
						std::ostream_iterator<double>(std::cout, "\t"));
					std::cout <<  std::endl;
				}
				std::cout <<  std::endl;
			}
		#endif

		outputPosterior(run, basename, nparams, posterior);
			
		run++;
	} // end loop through ref stats
	
	
}


// protected
std::vector< std::pair< double, double > > 
LinearRegressionStep::scale_simulated_summaries( const size_t nsummaries,
			   std::vector< std::vector< double > >* summaries )
{
	std::vector< std::pair< double,double> > mean_and_sd;
	const double n = double((*summaries)[0].size());
	
	for(unsigned stat = 0 ; stat < nsummaries ; ++stat) {
		double sum=0.,sumsq=0.;
		for( std::vector< double >::iterator itr1 = (*summaries)[stat].begin();
			itr1 != (*summaries)[stat].end() ;
			++itr1 ) {
	
			sum += *itr1;
			sumsq += (*itr1)*(*itr1);
		}
		
		double mean = sum/n;
		double sd = std::sqrt( (sumsq/(n-1)) - (sum*sum)/(n*(n-1.)) );
		
		mean_and_sd.push_back( std::make_pair(mean,sd) );
		
		for( std::vector< double >::iterator itr2 = (*summaries)[stat].begin() ;
		itr2 != (*summaries)[stat].end() ; ++itr2 ) {
		  *itr2 -= mean;
		  *itr2 /= sd;
		}
	}
	return mean_and_sd;
}

std::vector< double > LinearRegressionStep::euclidean_distances( 
					const size_t nsummaries,
				    const std::vector< double > & observed_summaries,
				    const std::vector< std::vector< double > >& simulated_summaries )
{
	std::vector< double > dist( simulated_summaries[0].size(),0. );
	
	for(unsigned i = 0 ; i < dist.size() ; ++i) {
		for(unsigned stat=0;stat<nsummaries;++stat) {
			
			#ifndef NDEBUG
				if ( ! std::isfinite(std::pow(observed_summaries[stat]-simulated_summaries[stat][i],2)) ) {
				  std::cerr << "error in distance calculation: " << observed_summaries[stat] << ' ' << simulated_summaries[stat][i] << '\n';		
				}
			#endif 
			dist[i] += std::pow(observed_summaries[stat]-simulated_summaries[stat][i],2);
		}
		dist[i] = sqrt(dist[i]);
	}
	return dist;
}

double LinearRegressionStep::get_distance_quantile(
				const std::vector< double > & dist, const double tolerance)
{
	std::vector< double > dist_copy(dist);
	
	sort(dist_copy.begin(),dist_copy.end());
	
	return ( gsl_stats_quantile_from_sorted_data(&dist_copy[0],1,
						   dist_copy.size(),
						   tolerance) );
}


/*! @todo split this up */
/*! @todo try catch blocks to free gsl stuff should be put in if we keep all this code */
std::vector<std::vector<double> > LinearRegressionStep::generate_posterior(
						const size_t nparams,
						const size_t nsumm,
						const std::vector< double > & scaled_observed,
						const std::vector<std::vector< double > > & prior,
						const std::vector<std::vector< double > > & scaled_summaries,
						const double tolerance,
						const size_t run)
{
	std::vector< std::vector<double> > posterior_params;
	
	//calculate the euclidean distances
	std::vector<double> dist = 
		euclidean_distances(nsumm, scaled_observed, scaled_summaries);
	
	#ifdef MYDEBUG
	{
		std::cout <<  "Distances are: " << std::endl;
		std::copy(dist.begin(), dist.end(), 
						std::ostream_iterator<double>(std::cout, "\t"));
					std::cout <<  std::endl;
  
	}
	#endif
	
	// and quantile
	double q = get_distance_quantile(dist, tolerance);
	
	#ifdef MYDEBUG
	{
		std::cout <<  "q = " << q << std::endl;
	}
	#endif

	//some storage space
	std::vector< std::vector<double> > accepted_params(nparams, std::vector<double>());
	std::vector< std::vector<double> > accepted_scaled_summs(nsumm, std::vector<double>());
	std::vector<double> regression_weights;

	std::vector<double>::iterator itr = dist.begin(),itr2;
	size_t naccepts = 0;

	if( std::fabs(q - 0.) > std::numeric_limits<double>::epsilon() )
	//then we do the linear regression bit
	{
		while( (itr2 = find_if(itr,
				 dist.end(),
				 std::bind2nd(std::less< double >(),q))) != dist.end() )
		{
			++naccepts;
			std::vector<double>::difference_type d = itr2-dist.begin();
			#ifdef DETAILDEBUG
				{
					std::cout <<  "found at pos = " << d << " and distance " << dist[d] 
					<<"\naccepted params\t";
				}
				#endif
			for(unsigned param=0; param<nparams; ++param)
			{
				accepted_params[param].push_back(prior[param][d]);
				#ifdef DETAILDEBUG
				{
					std::cout << prior[param][d] << "\t";
				}
				#endif
				
			}
			#ifdef DETAILDEBUG
				{
					std::cout <<  "\nwith stats\t";
				}
				#endif
			for(unsigned stat=0; stat<nsumm; ++stat)
			{
				accepted_scaled_summs[stat].push_back(scaled_summaries[stat][d]);
				#ifdef DETAILDEBUG
				{
					std::cout << scaled_summaries[stat][d] << "\t";
				}
				#endif
			}
			regression_weights.push_back( 1.0-(dist[d]*dist[d])/(q*q) ) ;
				#ifdef DETAILDEBUG
				{
					std::cout <<  "\nand weight = " << 1.0-(dist[d]*dist[d])/(q*q) << std::endl;
				}
				#endif
			itr = (itr2+1);
		}  // end finding acceptances
		
		if(naccepts==0)
		{
			return accepted_params;
		}
	
		#ifdef MYDEBUG
				{
					std::cout <<  "naccepts = " << naccepts << std::endl;
				}
			#endif
		
		// vector for posterior parameters
		
		posterior_params = std::vector< std::vector<double> > 
							(nparams, std::vector<double>(naccepts));
		
		// vector for 1,scaled reference summaries (the target)
		gsl_vector * target_plus_1 = gsl_vector_alloc(nsumm+1);
		gsl_vector_set(target_plus_1,0,1.);
		
		for(unsigned stat=0; stat<nsumm; ++stat)
		{
			gsl_vector_set(target_plus_1,stat+1,scaled_observed[stat]);
		}
      
		gsl_matrix * acc_summstat_matrix 
			= gsl_matrix_alloc(naccepts,nsumm+1);

		for(unsigned i = 0 ; i < naccepts ; ++i)
		{
			//add a dummy variable (column of 1s) so that intercept is 
			//estimated as the 1st parameter in vector b
			gsl_matrix_set(acc_summstat_matrix,i,0,1.);

			//WARNING: changed p.nparams to p.nsumm in next line...
			for(unsigned j = 0 ; j < nsumm ; ++j)
			{
				assert( std::isfinite( accepted_scaled_summs[j][i] ) );
				gsl_matrix_set(acc_summstat_matrix,
					 i,j+1,accepted_scaled_summs[j][i]);
			}
		}
		
		#ifdef MYDEBUG
		std::cout << "acc_summstat_matrix is " << std::endl;
		for(unsigned i = 0 ; i < acc_summstat_matrix->size1 ; ++i)
		{
			for(unsigned j = 0 ; j < acc_summstat_matrix->size2 ; ++j)
			{
				std::cout << gsl_matrix_get(acc_summstat_matrix,i,j) << "\t";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
		#endif
		      
		//free some memory - this reduces memory allocation to min poss
		// can't see we need this but a clear op was in Kevin's code
		// and at least this one really does free the memory allocation
		std::vector< std::vector<double> >().swap(accepted_scaled_summs); 
		
		//make the regression weights a gsl vector
		gsl_vector * gsl_weights_v = gsl_vector_alloc(naccepts);
		for(unsigned i=0; i<naccepts; ++i) {
			gsl_vector_set(gsl_weights_v,i,regression_weights[i]);
		}
		
		#ifdef MYDEBUG
		std::cout << "gsl_weights_v is " << std::endl;
		for(unsigned i = 0 ; i < gsl_weights_v->size ; ++i)
		{
			std::cout << gsl_vector_get(gsl_weights_v, i) << "\t";
		}
		std::cout << std::endl;
		#endif
		
		// can't see we need this but a clear op was in Kevin's code
		// and at least this one really does free the memory allocation
		std::vector<double>().swap(regression_weights);
		
		for(unsigned param=0; param<nparams; ++param)
		{
			//now, do the regression for each parameter
			gsl_vector * pvec = gsl_vector_alloc(naccepts);
			
			for(unsigned i=0; i<naccepts; ++i)
			{
				gsl_vector_set(pvec, i, accepted_params[param][i]);
			}
			
			#ifdef MYDEBUG
			std::cout << "For parameter index " << param << ", pvec is " << std::endl;
			for(unsigned i = 0 ; i < pvec->size ; ++i)
			{
				std::cout << gsl_vector_get(pvec,i) << "\t";
			}
			std::cout << std::endl;
			#endif

			//some gsl data types we will need
			gsl_matrix * cov_matrix 
				= gsl_matrix_alloc(nsumm+1, nsumm+1);

			gsl_vector * b = gsl_vector_alloc(nsumm+1);

			gsl_multifit_linear_workspace * work 
				= gsl_multifit_linear_alloc (naccepts, nsumm+1);

			double predicted_mean = 0;
			double chisq = 0;
		
			/* compute best fit parameters b for weighted model
			 * pvec = acc_summstat_matrix pvec
			 * using preallocated workspace work
			 * cov_matrix is computed with the given weights
			 * weighted sum of squares of residuals from best fit 
			 * is put into into chisq*/			
			gsl_multifit_wlinear(acc_summstat_matrix,
			       gsl_weights_v, pvec,
			       b, cov_matrix,
			       &chisq, work);
		    
		    #ifdef MYDEBUG
			std::cout << "Coefficients b are " << std::endl;
			for(unsigned i = 0 ; i < b->size ; ++i)
			{
				std::cout << gsl_vector_get(b,i) << "\t";
			}
			std::cout << std::endl;
			#endif
			
			//now, need to do predicted mean + residuals thing...
			//the b vector needs to be multiplied by the scaled target matrix..
			
			/* predicted mean is the predicted parameter value for 
			* for the actual observed */
			// note to me - what if observed is outside range of sims??
			gsl_blas_ddot(target_plus_1, b, &predicted_mean);
			
			#ifdef NDEBUG
			if (! isfinite(predicted_mean) )
			{
				cerr << "Warning: during analysis of data set " << (run+1) 
				<< " the predicted mean value of the regression is: " << predicted_mean 
				<< " for parameter " << (param+1) << '\n';
			}
			#endif
			
			#ifdef MYDEBUG
					std::cout << "predicted mean = " << predicted_mean << std::endl;
			#endif
			
			assert(std::isfinite(predicted_mean));
			
			//have to go through this and calculate residuals by hand
			//and generate the posterior
			// posterior is predicted mean + residual from regression
			for(unsigned i=0; i<naccepts; ++i)
			{
				double yi = gsl_vector_get(pvec,i);
				assert( std::isfinite(yi) );
				gsl_vector_const_view row = gsl_matrix_const_row(acc_summstat_matrix, i);
				double y_est = 0.0;
				gsl_blas_ddot(&row.vector, b ,&y_est);
				double ri = (yi-y_est);
				assert(std::isfinite(ri));
				assert(std::isfinite(y_est));
				posterior_params[param][i] = predicted_mean+ri;
				assert( std::isfinite(posterior_params[param][i]) );
				#ifdef DETAILDEBUG
					std::cout << "yi = " << yi << ", y_est = " << y_est << std::endl;
					std::cout << "ri = " << ri << std::endl;
				#endif	
			}
			#ifdef MYDEBUG
			std::cout << "posterior parameters for param index " << param << std::endl;
			for(unsigned i = 0 ; i < naccepts ; ++i)
			{
				std::cout << posterior_params[param][i] << "\t";
			}
			std::cout << std::endl;
			#endif
			
			gsl_vector_free(pvec);
			gsl_matrix_free(cov_matrix);
			gsl_vector_free(b);
			gsl_multifit_linear_free(work);
		
		}
		#ifdef MYDEBUG
			std::cout << "\nend of fitting parameters " << std::endl;
			#endif

		gsl_matrix_free(acc_summstat_matrix);
		gsl_vector_free(gsl_weights_v);
	}
	else //FIX FOR NANS/INFS PROBLEM IN THE R CODE
    /*
      We are in a situation where there are lots of exact matches to the summaries,
      or, we cannot tell the majority of euclidean distances from 0.  In this case, 
      the regression will not work stably, and so we just find the set of exact matches 
      and return it.
    */
	{
		while( (itr2 = find_if(itr,
			 dist.end(),
			 std::bind2nd(std::less_equal< double >(),q))) != dist.end() )
		{
			std::vector<double>::difference_type d = itr2-dist.begin();
			for(unsigned param=0;param<nparams;++param)
			{
			  accepted_params[param].push_back(prior[param][d]);
			}
			itr = (itr2+1);
		}
    }
	
	// just makes use of run if not debugging - prevent compiler warning
	assert(run+1 > 0 );
	
	return posterior_params;
}


void LinearRegressionStep::outputPosterior(
				size_t run, 
				const std::string& basename, 
				const size_t nparams,
				const std::vector< std::vector < double > >& posterior) const
{
	//creat suffix for output file names
	std::string suffix = "." + trans->getType() + ".post";
	
	//create output file name
	std::ostringstream ofname;
	ofname << basename << '.' << run  << suffix;
	FILE * ofp = fopen(ofname.str().c_str(), "w");
	
	//write posterior distribution out to file
	//as needed
	//!TODO clean this up
	for(unsigned i = 0 ; i < posterior[0].size() ; ++i)
	{
		unsigned j = 0;
		for(j = 0 ; j < nparams-1 ; ++j)
		{
			#ifndef NDEBUG
			
				//yell about any non-finite data
			
				if (!std::isfinite(posterior[j][i]))
					std::cerr << "Non-finite value for " << i 
					<< "th posterior parameter at element " << j << ": " 
					<<posterior[j][i] << std::endl;
				
				assert( std::isfinite(posterior[j][i]) );
			#endif
			fprintf(ofp,"%g\t", posterior[j][i] );
		}
		fprintf(ofp,"%g\n", posterior[j][i] );
	}
	#ifdef MYDEBUG
		{
		std::cout << "\nOutput file generated: " 
				<< ofname.str() << std::endl;
		std::cout <<  std::endl;
		}
	#endif

}
