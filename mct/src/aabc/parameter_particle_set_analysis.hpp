/* MCT - Markov Chains on Trees.

   Copyright (C) 2011 Jenny Harlow

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
\brief ParameterParticleSetAnalysis declarations.
*/

#ifndef _INC_PARAMETER_PARTICLE_SET_ANALYSIS_H
#define _INC_PARAMETER_PARTICLE_SET_ANALYSIS_H

#include "parameter_particle_set.hpp"
#include "parameter_particle.hpp"
#include "summary_statistic_set.hpp"
#include "summary_statistic.hpp"


#include <boost/function.hpp>
#include <boost/smart_ptr.hpp>


#include <vector>
#include <string>


namespace aabc {
	
	/** @brief A class for analysing a 
	ParameterParticle ParameterParticleSet.
	
	\note 
	This class keeps the analysis of data separate
	from the data itself, so that the %ParameterParticleSet
	is mostly just a way of being able to hold data. 
	This means that we can add to/change etc the analsysis,
	and in particular have augmented derived analysis classes using
	this as the base class, separate from the %ParameterParticleSet.
	However, you could argue that the %ParameterParticleSet
	should do its own analysis, in which case we could move these
	methods into that it.  At present I am in favour of the separate
	analysis class but could change my mind ...
	
	A %ParameterParticleSetAnalysis is constructed from a 
	%ParameterParticleSet and can produce various analyses or
	statistics for that %ParameterParticleSet.
	
	Changes the given %ParameterParticleSet after the
	%ParameterParticleSetAnalysis has been created will have
	no effect on the %ParameterParticleSetAnalysis. */
	 

	class ParameterParticleSetAnalysis {
		public:
		
			/*! \brief Constructor.
			
			\param pset A pointer to the ParameterParticleSet 
			to analyse.
			\pre \a pset must point to a ParameterParticleSet with
			at least one ParameterParticle in it, and each
			ParameterParticle in the set must have the same number of
			parameters.*/				
			explicit ParameterParticleSetAnalysis(
				const boost::shared_ptr< const ParameterParticleSet >& pset);
			
			virtual ~ParameterParticleSetAnalysis();
			
			/*! \brief Get the number of elements in the 
			\link aabc::ParameterParticle ParameterParticles\endlinkin the set being analysed.*/
			size_t size() const;
			
			/*! \brief Get whether the set being analysed is empty. */
			bool empty() const;
			
			/*! \brief Get the number of 
			\link aabc::ParameterParticle ParameterParticles\endlink in 
			the set being analysed.*/
			size_t nParticles() const;
			
			/*! \brief Make and return a SummaryStatisticSet
			holding the standardised values from the particles in 
			the set being analysed by this.
			* 
			Standardising means centering the values from each 
			ParameterParticle in the set being analysed by this on its
			mean (ie subtracting the mean for that element of the 
			particle over all particles in the set) and scaling each centred
			value by its sample standard deviation (ie for each
			centred value, dividing
			by the sample standard deviation of the values of 
			that element across all particles in the set).
			
			\return a pointer to a %SummaryStatisticSet holding	
			standardised values from the ParameterParticles in 
			the set being analysed by this.
			\pre If the set being analysed by this is not empty, it
			should contain at least two values.*/				
			virtual boost::shared_ptr < mct::SummaryStatisticSet >
				makeSummaryStatisticSetStandardisedParameterValues() const;

			
			/*! \brief Get all the parameter values as a container
			 of containers of doubles.
			 
			 If the given
			 ParameterParticleSet contains n
			 \link aabc::ParameterParticle ParameterParticles\endlink
			 with nparams elements in each, then the outer container of
			 the returned structure contains nparams inner containers,
			 and each inner container contains n values, and the
			 inner containers are in the same order as the 
			 elements are in the ParameterParticles.  ie, the 
			 i<sup>th</sup> value in the j<sup>th</sup> 
			 %ParameterParticle in the set will be the j<sup>th</sup>
			 element of the i<sup>th</sup> inner container.
			
			\return a pointer to a copy of 
			the container of containers of values from the 
			ParameterParticleSet used to construct this.*/
			virtual boost::shared_ptr < std::vector< std::vector< double > > > 
						getAllParameters() const;
			
			/*! \brief Get all the parameter values, standardised.
			
			Standardising means centering each value on its
			mean (ie subtracting the mean over all %ParameterParticles
			in the set being analysed from each
			%ParameterParticle in the set) and scaling each centred
			%ParameterParticle by its sample standard deviation (ie for each
			centred values, dividing
			the each element by the
			corresponding element of the sample standard deviations of
			the ParameterParticles in the set being analysed.
			
			 If this is analysing a set containing n
			 \link aabc::ParameterParticle ParameterParticles\endlink
			 with nparam elements in each, then the outer container in
			 the returned structure contains nparam inner containers,
			 and each inner container contains n values, and the
			 inner containers are in the same order as the 
			 elements are in the ParameterParticles.  ie, the 
			 i<sup>th</sup> value in the j<sup>th</sup> 
			 %ParameterParticle in the set will be the 
			 result of standardising the j<sup>th</sup>
			 element of the i<sup>th</sup> inner container
			 using mean and standard deviation of the j<sup>th</sup>
			 elements over all inner containers.
			\return a pointer to a container of containers of	
			standardised values from the set being analysed by this.	
			\pre If this is not empty, there should be at least two
			elements in the set being analysed..*/		
			virtual boost::shared_ptr < std::vector < std::vector < double > > >
						getAllParametersStandardised() const;
			
			/*! \brief Get a mct::SummaryStatistic containing 
			the means of the values in the ParameterParticleSet
			used to construct this analysis.  
						
			\return A pointer to a SummaryStatistic
			whose values are the means of the values in 
			the ParameterParticles in the ParameterParticleSet
			used to construct this analysis, i.e. the i<sup>th</sup>
			value in the returned SummaryStatistic is the mean
			of the i<sup>th</sup> values of all the 
			ParameterParticles in the given set. 
			\pre There should be at least one element in the set
			being analysed.*/				
			virtual boost::shared_ptr < const mct::SummaryStatistic > 
							getMeans() const;
			
			/*! \brief Get a mct::SummaryStatistic containing 
			the sample variances of the values in the ParameterParticleSet
			used to construct this analysis.  
						
			\return A pointer to a SummaryStatistic
			whose values are the sample variances of the values in 
			the ParameterParticles in the ParameterParticleSet
			used to construct this analysis, i.e. the i<sup>th</sup>
			value in the returned SummaryStatistic is the sample variance
			of the i<sup>th</sup> values of all the 
			ParameterParticles in the given set. 
			\pre There should be at least one element in the set
			being analysed.*/	
			virtual boost::shared_ptr < const mct::SummaryStatistic > 
							getSampleVars() const;
			
			/*! \brief Get a mct::SummaryStatistic containing 
			the sample standard deviations of the values in the ParameterParticleSet
			used to construct this analysis.  
						
			\return A pointer to a SummaryStatistic
			whose values are the sample standard deviations of the values in 
			the ParameterParticles in the ParameterParticleSet
			used to construct this analysis, i.e. the i<sup>th</sup>
			value in the returned SummaryStatistic is the sample 
			standard deviation
			of the i<sup>th</sup> values of all the 
			ParameterParticles in the given set. 
			\pre There should be at least one element in the set
			being analysed.*/	
			virtual boost::shared_ptr < const mct::SummaryStatistic > 
							getSampleSDs() const;
			
			/*! \brief Get a SummaryStatistic containing 
			the root mean squared errors of the values in the
			ParameterParticleSet
			used to construct this analysis against the
			values in the given ParameterParticle \a p.  
						
			\return A pointer to a SummaryStatistic
			whose values are the root mean squared errors
			of the values in 
			the ParameterParticles in the ParameterParticleSet
			used to construct this analysis against the
			values in the given %ParameterParticle \a p,
			i.e. the i<sup>th</sup>
			value in the returned SummaryStatistic is the RMSE 
			of the i<sup>th</sup> values of all the 
			ParameterParticles in the given set against the 
			i<sup>th</sup> value of \a p.
			
			For a single parameter the RMSE is calculated as 
			\f[ \frac{1}{n} \sum_{i=1}^2\left(\frac{\widehat{\theta}_i - \theta}{\theta} \right)^2  \f]
			where \f$ \widehat{\theta}_i \f$ is the value of the
			parameter for the i<sup>th</sup> estimate (i.e., 
			%ParameterParticle in the set) and
			\f$ \theta \f$ is the 'true' value, 
			(i.e. the value in \a p) 
			
			\pre There should be at least one element in the set
			being analysed.*/	
			virtual boost::shared_ptr < const mct::SummaryStatistic > 
					getRMSEs(
					const boost::shared_ptr < const ParameterParticle >& p) 
																const;
					
		protected:
			
			ParameterParticleSetAnalysis();
			
			explicit ParameterParticleSetAnalysis(
							const ParameterParticleSetAnalysis& other);
			
			ParameterParticleSetAnalysis& operator=(
									ParameterParticleSetAnalysis tmp);
			
			std::vector < double >& vec_sampleVars(
						std::vector < double >& vec) const;
					
			std::vector < double >& vec_sampleMeans(
						std::vector < double >& vec) const;
		
		
			/* \brief A pointer to a container of containers of doubles.
			 The i-th inner container is a container of the parameter values
			 at position i in the particles in the set.  There
			 are as many inner containers as there are parameters in
			 the particles. */
			boost::shared_ptr < std::vector< std::vector< double > > > parameters;
			
			std::vector< double > sums;
			
			std::vector< double > sumsqs;
			
			// number of particles in the parameterparticleset
			size_t n;
			
		private:
		
			static double standardise(double value, double mean, double sd);
			
	};
	
	//std::ostream& operator<< (std::ostream& out, 
							//const ParameterParticleSetAnalysis& p);

} // end namespace aabc

#endif 
