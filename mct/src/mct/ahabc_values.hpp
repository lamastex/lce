/* MCT - Markov Chains on Trees.

   Copyright (C) 2012 Jenny Harlow

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
\brief AHABCValues static class declarations.
*/

#ifndef _AHABC_VALUES_H
#define _AHABC_VALUES_H

#include "summary_statistic_set.hpp"
#include "parameter_particle_set.hpp"

#include <boost/smart_ptr.hpp>


#include <vector>
#include <string>

namespace mct {
	
	/*! \brief A structure for values for summary statistics and
	parameter particles for ahabc.
	
	Does not do a lot but enables these values to be organised and
	manipulated for passing on to the histograms.  
	
	 */
	class AHABCValues {
		public:
			/*! \brief Constructor.
			 * 
			 The constructor takes the sets of summary statistics and
			 parameter particles that this will manipulate.  The 
			 set of summary statistics and parameter particles will be 
			 treated as one set of ordered pairs, ie the i<sub>th</sub>
			 summary statistic in \a s_ptr will be paired with the
			 i<sub>th</sub> parameter particle in \a p_ptr.
			 
			 \param p_ptr Pointer to the set of parameter particles.
			 \param s_ptr Pointer to the set of summary statistics.
			 \pre Both \a p_ptr and \a s_ptr are not empty pointers and 
			 both point to sets of the same size (ie we have the same 
			 number of summary statistics as parameter particles).
			 \pre The size of the summary statistics within \a s_ptr are all
			 consistent and the size of the parameter particles within
			 \a p_ptr are all consistent. */
			AHABCValues(
				const boost::shared_ptr< const aabc::ParameterParticleSet > p_ptr,
				const boost::shared_ptr < const SummaryStatisticSet > s_ptr);
			
			/*! \brief Constructor.
			 * 
			 The constructor takes sets of summary statistics and
			 parameter particles and information about the required
			 parts of these that this will manipulate.  
			 
			 The required parts specify the positions of the 
			 required parameters and summary statistics to be included
			 in this. Positions are given as index numbers with the 
			 index of the first position being 1.

			 For example, each %ParameterParticle in 
			 \a p_ptr may have 3 elements, and \a reqParams may specify
			 the positions 1 and 3.  This will then only hold values 
			 for the first and third elements of all the ParameterParticles
			 in \a p_ptr.  Similarly, if each %SummaryStatistic in 
			 \a s_ptr has 2 elements, and \a reqStats specifies
			 just the position 2, then this will then only hold values 
			 for the second element of all the SummaryStatistics
			 in \a s_ptr. 
			 
			 The required summary statistics and required parameter 
			 particles will be 
			 treated as one set of ordered pairs, ie the required elements
			 for the i<sub>th</sub>
			 summary statistic in \a s_ptr will be paired with the
			 required elements of the i<sub>th</sub> parameter particle
			 in \a p_ptr.
			 
			\param p_ptr Pointer to the set of parameter particles.
			\param reqParams specifies the indices of the required 
			ParameterParticles for the new set.  The indices 
			can take values from 1 upwards.  The maximum index
			is the number of parameters in each %ParameterParticle
			in \a p_ptr.
			\param s_ptr Pointer to the set of summary statistics.
			\param reqStats specifies the indices of the required 
			SummaryStatistics for the new set.  The indices 
			can take values from 1 upwards.  The maximum index
			is the number of statistics in each %SummaryStatistic
			in \a s_ptr.
			 \pre Both \a p_ptr and \a s_ptr are not empty pointers and 
			 both point to sets of the same size (ie we have the same 
			 number of summary statistics as parameter particles).
			 \pre The size of the summary statistics within \a s_ptr are all
			 consistent and the size of the parameter particles within
			 \a p_ptr are all consistent. 
			 \pre \a reqParams and \a reqStats are non-empty and contain
			 only values that are valid indices into the sequences of 
			 parameters and statistics represented in \a p_ptr and \a s_ptr
			 respectively.*/
			AHABCValues(
			const boost::shared_ptr< const aabc::ParameterParticleSet > p_ptr,
			const std::vector <size_t>& reqParams,
			const boost::shared_ptr < const SummaryStatisticSet > s_ptr,
			const std::vector <size_t>& reqStats);
			
			/*! \brief Constructor.
			 * 
			 The constructor takes sets of summary statistics and
			 parameter particles and information about the required
			 summary statistics that this will manipulate.  
			 
			 The required summary statistics specify the positions of the 
			 summary statistics to be included
			 in this. Positions are given as index numbers with the 
			 index of the first position being 1.

			 For example, if each %SummaryStatistic in 
			 \a s_ptr has 2 elements, and \a reqStats specifies
			 just the position 2, then this will then only hold values 
			 for the second element of all the SummaryStatistics
			 in \a s_ptr. 
			 
			 The required summary statistics and parameter 
			 particles will be 
			 treated as one set of ordered pairs, ie the required elements
			 for the i<sub>th</sub>
			 summary statistic in \a s_ptr will be paired with the
			 elements of the i<sub>th</sub> parameter particle
			 in \a p_ptr.
			 
			\param p_ptr Pointer to the set of parameter particles.
			\param s_ptr Pointer to the set of summary statistics.
			\param reqStats specifies the indices of the required 
			SummaryStatistics for the new set.  The indices 
			can take values from 1 upwards.  The maximum index
			is the number of statistics in each %SummaryStatistic
			in \a s_ptr.
			 \pre Both \a p_ptr and \a s_ptr are not empty pointers and 
			 both point to sets of the same size (ie we have the same 
			 number of summary statistics as parameter particles).
			 \pre The size of the summary statistics within \a s_ptr are all
			 consistent and the size of the parameter particles within
			 \a p_ptr are all consistent. 
			 \pre \a reqStats is non-empty and contains
			 only values that are valid indices into the sequence of 
			 statistics represented in \a s_ptr.*/
			AHABCValues(
			const boost::shared_ptr< const aabc::ParameterParticleSet > p_ptr,
			const boost::shared_ptr < const SummaryStatisticSet > s_ptr,
			const std::vector <size_t>& reqStats);
			
			/*! \brief Constructor.
			 * 
			 The constructor takes sets of summary statistics and
			 parameter particles and information about the required
			 parts of the parameter particles that this will manipulate.  
			 
			 The required parts specify the positions of the 
			 required parameters to be included
			 in this. Positions are given as index numbers with the 
			 index of the first position being 1.

			 For example, each %ParameterParticle in 
			 \a p_ptr may have 3 elements, and \a reqParams may specify
			 the positions 1 and 3.  This will then only hold values 
			 for the first and third elements of all the ParameterParticles
			 in \a p_ptr.  
			 
			 The summary statistics and required parameter 
			 particles will be 
			 treated as one set of ordered pairs, ie the elements
			 for the i<sub>th</sub>
			 summary statistic in \a s_ptr will be paired with the
			 required elements of the i<sub>th</sub> parameter particle
			 in \a p_ptr.
			 
			\param p_ptr Pointer to the set of parameter particles.
			\param reqParams specifies the indices of the required 
			ParameterParticles for the new set.  The indices 
			can take values from 1 upwards.  The maximum index
			is the number of parameters in each %ParameterParticle
			in \a p_ptr.
			\param s_ptr Pointer to the set of summary statistics.
			\pre Both \a p_ptr and \a s_ptr are not empty pointers and 
			 both point to sets of the same size (ie we have the same 
			 number of summary statistics as parameter particles).
			 \pre The size of the summary statistics within \a s_ptr are all
			 consistent and the size of the parameter particles within
			 \a p_ptr are all consistent. 
			 \pre \a reqParams is non-empty and contains
			 only values that are valid indices into the sequence of 
			 parameters represented in \a p_ptr.*/
			AHABCValues(
			const boost::shared_ptr< const aabc::ParameterParticleSet > p_ptr,
			const std::vector <size_t>& reqParams,
			const boost::shared_ptr < const SummaryStatisticSet > s_ptr);
				
			~AHABCValues(){}
			
			/*! \brief Get size of the data collection in this.
			
			The size is the number of summary statistic - parameter
			value pairs held in this.*/	
			size_t size() const;
			
			
			/*! \brief Get the number of parameters, ie the
			number of dimensions we will need for parameters.*/	
			size_t getNumberParameters() const;
			
			/*! \brief Get the number of summary statistics, ie the
			number of dimensions we will need for summary statistics.*/	
			size_t getNumberSummaryStatistics() const;
			
			/*! \brief Get a collection of the dimensions relating
			to the parameters.
			
			Returns a collection of the positions of the parameter
			values for each row of the data comprising the
			parameter-summary statistic pairs.  
			
			\note The positions are labeled from 1 upwards, ie the 
			first dimension is labelled 1, etc.
			 
			\return A collection of the positions of the parameter
			values for each row of the data comprising the
			parameter-summary statistic pairs.*/	
			std::vector < int >  getParametersDims() const;
			
			/*! \brief Get a collection of the dimensions relating
			to the summary statistics.
			
			Returns a collection of the positions of the summary statistic
			values for each row of the data comprising the
			parameter-summary statistic pairs.  
			
			\note The positions are labeled from 1 upwards, ie the 
			first dimension is labelled 1, etc (this would normally
			be a parameter position so that the first position 
			we could expect for summary statistics would be 2, if there
			was only one parameter).
			 
			\return A collection of the positions of the summary statistic
			values for each row of the data comprising the
			parameter-summary statistic pairs.*/	
			std::vector < int >  getSummaryStatisticsDims() const;
			
			/*! \brief Fill a provided container with our data.
			 * 
			The i<sub>th</sub> inner container holds the 
			values for the i<sub>th</sub> parameter particle and the
			i<sub>th</sub> summary statistic.
			\pre \a container is empty.
			\post \a container holds size() inner vectors, each 
			containing first the parameter particle values and then
			the summary statistic values for a summary.*/	
			void fillDataContainer(
				std::vector < std::vector < double > >& container) const;
			
			/*! \brief Fill a provided container with transformed data.
			
			The data is transformed using the Householder Transformation.
			
			The i<sub>th</sub> inner container holds the transformed
			values for the i<sub>th</sub> parameter particle and the
			i<sub>th</sub> summary statistic.
			\pre \a container is empty.
			\post \a container holds size() inner vectors, each 
			containing first the transformed parameter particle values and then
			the transformed  summary statistic values for a summary.*/	
			void fillTransformedDataContainer(
				std::vector < std::vector < double > >& container) const;
				
			void fillStandardisedDataContainer(
				std::vector < std::vector < double > >& container) const;
			
			void fillTransformedStandardisedDataContainer(
				std::vector < std::vector < double > >& container) const;

			
			/*! \brief Get the Holdholder transformation (reflection)
			 * for the data in this.
			 * 
			 * The Householder matrix gives a transformation
			 * that aligns the dominant direction the data
			 * onto the first axis in the transformed 
			 * coordinates. */
			std::vector < std::vector < double > > 
								getHouseholderMatrix() const;
			
			std::vector < std::vector < double > > 
								getHouseholderMatrixStandardised() const;
			
			std::vector < double > getMeans() const;
			
			std::vector < double > getSDs() const;
			
			
			std::vector < double > standardise(
							const std::vector < double >& invec) const;
							
			/*! \brief Get the Holdholder transformation of some
			 * specific data.
			 * 
			 * \param invec The data to transform
			 * \return The transformed data.
			 * \pre \a invec must have the same dimensions as the data
			 * in this, ie invec.size() == getNumberParameters() +
			 * getNumberSummaryStatistics(). */			
			std::vector < double >householderTransform(
							const std::vector < double >& invec) const;
							
			std::vector < double > standardiseAndHouseholderTransform(
							const std::vector < double >& invec) const;
			
			void dataSVD() const;
			
		private:
		
			AHABCValues();
			
			AHABCValues(const AHABCValues& other);
			
			AHABCValues& operator=(AHABCValues rhs);
			
			boost::shared_ptr < std::vector < std::vector < double > > > 
							extractData() const;
			
			boost::shared_ptr < std::vector < std::vector < double > > > 
							extractStandardisedData() const;
			
			boost::shared_ptr < std::vector < std::vector < double > > > 
							extractTransformedData() const;
			
			boost::shared_ptr < std::vector < std::vector < double > > > 
							extractTransformedStandarisedData() const;
			
			void createHouseholderMatrix() const;
			
			void createHouseholderMatrix(
				boost::shared_ptr < 
					std::vector < std::vector < double > > 
										> dataPtr) const;
			
			void createHouseholderMatrixStandardised() const;
			
			void createHouseholderMatrixStandardised(
				boost::shared_ptr < 
						std::vector < std::vector < double > > 
											> dataPtr) const;
			
			std::vector < std::vector < double > > createHouseholderMatrix(
				boost::shared_ptr 
					< std::vector < std::vector < double > > > dataPtr,
					const std::vector < double >& meansForDataPtr) const;
			
			// transform with householder transformation			
			boost::shared_ptr < std::vector < std::vector < double > > >
				transformData(
				boost::shared_ptr < 
						std::vector < std::vector < double > > 
									> dataPtr) const;
			
			// transform with householder transformation on standardised data		
			boost::shared_ptr < std::vector < std::vector < double > > >
				transformStandardisedData(
				boost::shared_ptr <
						std::vector < std::vector < double > > 
									> dataPtr) const;
			
			//the summary statistics 
			boost::shared_ptr < const SummaryStatisticSet > ss_ptr;
			
			//the parameter values
			boost::shared_ptr < const SummaryStatisticSet > ps_ptr;
			
			//a container of the extracted values
			//boost::shared_ptr < std::vector < std::vector < double > > > dataPtr;
			
			//a container of the extracted values with Householder transformation
			//mutable boost::shared_ptr < std::vector < std::vector < double > > > transformedDataPtr;
			
			// inner vectors are the rows of the Householder matrix
			mutable std::vector < std::vector < double > > 
										householderMatrix;
			
			// Householder matrix for standardised data
			mutable std::vector < std::vector < double > > 
										householderMatrixStandardised;
			
			mutable std::vector < double > means;
			mutable std::vector < double > sds;
		
	};

} // end namespace mct

#endif 

