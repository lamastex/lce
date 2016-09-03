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
\brief LinearRegressionStrp declarations.
*/

#ifndef __REG_LINREG_HPP__
#define __REG_LINREG_HPP__

//aabc headers
#include "parameter_particle_set.hpp"
#include "summary_statistic_set.hpp"

//Headers for this project
#include "reg_transformer.hpp" 


#include <boost/smart_ptr.hpp>

#include <string>
#include <vector>
#include <utility>

namespace aabc {
	class ParameterParticleSet;
	class SummaryStatisticSet;
	class Transformer;
	
}

namespace reg {
	
	class Transformer;
	
	class LinearRegressionStep {
		
		public:
		
			explicit LinearRegressionStep(const boost::shared_ptr < Transformer >& _t);
			
			virtual ~LinearRegressionStep();
		
			//virtual std::vector < boost::shared_ptr < const aabc::ParameterParticleSet > >
			/*! @todo change this to just take our whole validation results set (that has distances so no ref stat needed?) */
			virtual void 
				linreg(
				const boost::shared_ptr < const aabc::ParameterParticleSet>& pset,
				const boost::shared_ptr <mct::SummaryStatisticSet>& ss_set,
				const boost::shared_ptr <mct::SummaryStatisticSet>& ref_stat_set,
				const double tolerance,
				const std::string& basename);
		
		protected:
		
			explicit LinearRegressionStep ();
			
			explicit LinearRegressionStep (const LinearRegressionStep& other);
		
			LinearRegressionStep& operator=(LinearRegressionStep tmp);
			
			std::vector< std::pair < double, double > > 
				scale_simulated_summaries( 
						const size_t nsummaries,
						std::vector< std::vector< double > > * summaries );

			std::vector< double > euclidean_distances( const size_t nsummaries,
						const std::vector< double >& observed_summaries,
						const std::vector< std::vector< double > >& simulated_summaries );

			double get_distance_quantile(const std::vector< double>& dist,
											const double tolerance);
											
			std::vector<std::vector< double> > generate_posterior(
						const size_t nparams,
						const size_t nsumm,
						const std::vector< double> & scaled_observed,
						const std::vector< std::vector< double > > & prior,
						const std::vector< std::vector< double > > & scaled_summaries,
						const double tolerance,
						const size_t run);

			void outputPosterior(
				size_t run, 
				const std::string& basename, 
				const size_t nparams,
				const std::vector< std::vector < double > >& posterior) const;
		
			boost::shared_ptr < Transformer > trans;
			
		private:
		
	};
}

#endif

