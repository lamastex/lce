/* mscplusplus - ms for c++.

   Copyright (C) 2011, 2012 Jenny Harlow
   
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
\brief SubPopulationStructure and PopulationStructure declarations.
*/

#ifndef _INC_POPULATION_STRUCTURE_H
#define _INC_POPULATION_STRUCTURE_H

#include "recomb_params.hpp"
#include "segsite_params.hpp"

#include "event.hpp"

#include <map>
#include <vector>
#include <string>
#include <set>
#include <sstream>

namespace hudson_ms {
	
	/** @brief A class for describing a population structure.
	 
	 Essentially, a manager for an ordered collection of sub-populations,
	 including their inter-subpopulation migration parameters,
	 and a population event list.
	 
	 Each subpopulation can have a text (string) label.
	  
	 The population event list describes demographic events
	 which happened to the population backwards in time from the
	 present (the population structure as it is now).  No
	 checks are carried out on the demographic events given
	 with respect to the current structure.  The demographic
	 events are mutable, ie can be changed or added to while this 
	 (ie the current structure) remains const.
	 
	 All indexes into the population structure should be
	 given, and will be interpreted, in terms of
	 array indexing {0, 1, ... [size of array - 1]}, ie 
	 the index for the first subpopulation in the structure
	 is 0.
	  
	*/

	class PopulationStructure {
				
		
		public:
		
		
			/*! \brief Constructor.
			
			Constructs an empty population structure.
			*/
			explicit PopulationStructure();
			
			/*! \brief Constructor.
			
			Constructs a population structure with just one subpopulation.
			
			\param _sn The number of samples in the subpopulation.
			\pre \a _sn >= 0.*/
			explicit PopulationStructure(int _sn);
			
			/*! \brief Copy constructor. */
			PopulationStructure(const PopulationStructure& other);
			
			virtual ~PopulationStructure();
			
			/*!Assignment operator.*/
			PopulationStructure& operator=(PopulationStructure tmp);
			
			void swap(PopulationStructure& other);
			
			/*! \brief Get number of samples for an indexed subpopulation.
			
			\param index The index of the subpopulation to get.
			\return Value of the sample size in the subpopulation
			at \a index. */
			size_t getSubNsamAt(int index) const;
			
			/*! \brief Get number of samples for a labeled subpopulation.
			
			\param label The label of the subpopulation.
			\return Value of the sample size for the 
			subpopulation labeled \a label. */
			size_t getSubNsamAt(const std::string& label) const;
			
			
			/*! \brief Get size ratio for an indexed subpopulation.
			
			\param index The index of the subpopulation to get.
			\return Value of the size ratio in the subpopulation
			at \a index. */
			double getSizeRatioAt(int index) const;
			
			/*! \brief Get size ratio for a labeled subpopulation.
			
			\param label The label of the subpopulation.
			\return Value of the size ratio for the 
			subpopulation labeled \a label. */
			double getSizeRatioAt(const std::string& label) const;
			
			/*! \brief Get growth rate for an indexed subpopulation.
			
			\param index The index of the subpopulation.
			\return Value of the growth rate for the 
			subpopulation at \a index. */
			double getGrowthRateAt(int index) const;
			
			/*! \brief Get growth rate for a labeled subpopulation.
			
			\param label The label of the subpopulation.
			\return Value of the growth rate ratio for the 
			subpopulation labeled \a label. */
			double getGrowthRateAt(const std::string& label) const;
			
			/*! \brief Get an ordered collection of labels for
			all the subpopulations.
			
			\return an ordered collection of labels for
			all the subpopulations. */
			std::vector < std::string > getLabels() const;
			
			/*! \brief Get label for an indexed subpopulation.
			
			\param index The index of the subpopulation.
			\return Label for the  
			subpopulation at index \a index. */
			std::string getLabelAt(int index) const;
			
			/*! \brief Get whether a label exists in the population.
			
			\param label The label to check.
			\return True if some subpopulation currently has this label,
			false otherwise. */
			bool hasLabel(const std::string& label) const;
			
			/*! \brief Add a subpopulation to the 
			end of the collection.
			
			Only the number of samples in the subpopulation is specified.
			* 
			All migration parameters into and out of the new 
			subpopulation are 0.0 by default; 
			
			\param _sn The number of samples in the subpopulation to add.
			\pre \a _sn >= 0.*/
			void add(int _sn);
			
			/*! \brief Add a subpopulation to the 
			end of the collection.
			
			All migration parameters into and out of the new 
			subpopulation are 0.0 by default; 
			
			\param _sn The number of samples in the subpopulation to add.
			\param _sr The size ratio of the subpopulation to add.
			\param _gr The growth rate of the subpopulation to add.
			\pre \a _sn >= 0, \a _sr > 0.0.*/
			void add(int _sn, double _sr, double _gr);
			
			/*! \brief Add a subpopulation to the 
			end of the collection.
			
			All migration parameters into and out of the new 
			subpopulation are 0.0 by default; 
			
			\param _sn The number of samples in the subpopulation to add.
			\param _gr The growth rate of the subpopulation to add.
			\pre \a _sn >0 0.*/
			void add(int _sn, double _gr);
			
			/*! \brief Erase a subpopulation from the collection.
			
			All migration parameters into other subpopulatinos
			out of the erased subpopulation are erased
			with the erasure of the subpopulation from the collection.  
			
			\param eraseLabel The label of the subpopulation to erase.
			\pre The list of demographic events must be empty
			for erasure of a subpopulation to be possible.*/
			void erase(const std::string& eraseLabel);
			
			/*! \brief Erase a subpopulation from the collection.
			
			All migration parameters into other subpopulatinos
			out of the erased subpopulation are erased
			with the erasure of the subpopulation from the collection.  
			
			\param orderErase The index of the subpopulation to erase.
			\pre The list of demographic events must be empty
			for erasure of a subpopulation to be possible.*/
			void erase(int orderErase);
			
			/*! \brief Reset the ordering for the subpopulations.
		
			Reorders by taking a parameter representing the required
			order, referencing each existing subpopulation by name, ie
			\a newLabelOrder should be a permuation of the set of 
			labels for the subpopulations in this.
			
			\param newLabelOrder An ordered collection of
			unique labels representing a permutation of the existing
			set of labels for the subpopulations.
			\pre The list of demographic events must be empty
			for reordering to be possible.
			\pre ie every existing
			label should be in \a newLabelOrder and the total 
			size of \a newLabelOrder should be the same as the number
			of subpopulations in this population.
			\post The subpopulations in this will be reordered so
			that the subpopulation at index i has label = \a newLabelOrder[i].
			*/
			void resetOrdering(
					const std::vector < std::string >& newLabelOrder);
			
			
			/*! \brief Relabel a subpopulation from the collection.
			
			The new label (\a newLabel) must be unique in the currently
			existing set of labels (but can be the same as \a oldLabel).
			
			All other subpopulations with migration parameters 
			referencing the label \a oldLabel will be updated
			with the new label information.
			
			\param oldLabel The label of the subpopulation to relabel.
			\param newLabel The new label of the subpopulation currently
			labeled \a oldLabel.*/
			void relabel(const std::string& oldLabel,
									const std::string& newLabel);
									
			/*! \brief Relabel a subpopulation from the collection.
			
			The new label (\a newLabel) must be unique in the currently
			existing set of labels (but can be the same as 
			the existing label for the subpopulation at \a index).
			
			All other subpopulations with migration parameters 
			referencing the label \a oldLabel will be updated
			with the new label information.
			
			\param order The index of the subpopulation to relabel.
			\param newLabel The new label of the subpopulation currently
			labeled \a oldLabel.*/
			void relabel(int order,
									const std::string& newLabel);
			
			/*! \brief Relabel all subpopulations in the collection.
			
			The new labels must each be unique within the collection
			\newLabels (but any subpopulation can be 'relabelled'
			with its existing name subject to this).
			
			All subpopulation migration parameters 
			referencing the old labels will be updated
			with the new label information.
			
			\param newLabels The ordered collection of new labels
			for the subpopulations.  ie the label at
			newLabels[i] will be applied to the subpopulation indexed i.
			\pre \a newLabels must contain the same number of 
			elements as there are subpopulations in this population.*/
			void relabel(
						const std::vector < std::string >& newLabels);
			
			/*! \brief Set the sample size for a subpopulation 
			in the population.
			
			\param index The index of the subpopulation.
			\param _sr The sample size to set for indexed 
			subpopulation.*/
			virtual void setSubNsamAt(int index, 
											int _sn);
			
			/*! \brief Set the sample size for a subpopulation 
			in the population.
			
			\param label The label for the subpopulation.
			\param _sr The sample size to set for labeled  
			subpopulation.*/
			virtual void setSubNsamAt(const std::string& label, 
											int _sn);
			
			/*! \brief Set the size ratio for a subpopulation 
			in the population.
			
			\param index The index of the subpopulation.
			\param _sr The size ratio to set for indexed 
			subpopulation.*/
			virtual void setSizeRatioAt(int index, 
											double _sr);
			
			/*! \brief Set the size ratio for a subpopulation 
			in the population.
			
			\param label The label for the subpopulation.
			\param _sr The size ratio to set for labeled  
			subpopulation.*/
			virtual void setSizeRatioAt(const std::string& label, 
											double _sr);
			
			/*! \brief Set the growth rate for a subpopulation 
			in the population.
			
			\param index The index of the subpopulation.
			\param _gr The growth rate to set for indexed 
			subpopulation.*/
			virtual void setGrowthRateAt(int index, 
											double _gr);
			
			/*! \brief Set the growth rate for a subpopulation 
			in the population.
			
			\param label The label for the subpopulation.
			\param _gr The growth rate to set for labeled  
			subpopulation.*/
			virtual void setGrowthRateAt(const std::string& label, 
											double _gr);
			
			/*! \brief Set the growth rate for all subpopulations
			 in the population.
			
			\param _gr The growth rate for all the
			subpopulations.*/
			void setAllGrowthRates(double _gr);
			
			/*! \brief Get the number of subpopulations
			 in the population.	*/
			size_t nSubPops() const;
			
			/*! \brief Return whether this contains no 
			subpopulations.	*/
			bool empty() const;
			
			/*! \brief Return whether the demographic event list
			 for this is empty.	*/
			bool eListEmpty() const;
			
			/*! \brief Get the total number of samples across
			all subpopulations.*/
			size_t totalNsam() const;
			
			
			/*! \brief Set the migration matrix between all subpopulations.
			 
			The migration matrix is a backwards migration matrix where
			\f$ M_{ij} = 4N_0m_{ij} \f$ and \f$ m_{ij} \f$ is the fraction
			of subpopulation \f$ i \f$ made up of migrants <strong>from</strong>
			subpopulation \f$ j \f$ each generation.
			
			\a _mm is specified such that the inner vectors are the 'rows', ie
			row \f$ i \f$ specifies \f$ M_{ij} = 4N_0m_{ij} \f$ for 
			\f$ j \f$ = each subpop in turn and \f$ m_{ij} \f$ is the fraction
			of subpopulation \f$ i \f$ made up of migrants <strong>from</strong>
			subpopulation \f$ j \f$ each generation (ie an inwards
			migration parameter).  All elements on the diagonal are
			ignored.
			 
			\param _mm The migration matrix to use.
			\pre _mm has appropriate dimensions for the population structure
			and contains no negative values.*/
			void setMigrationMatrix(const std::vector < std::vector <double> >& _mm);
		
			/*! \brief Set the migration parameter between all subpopulations.
			 
			Sets the migration parameter between all subpopulations to be
			\a _m.
			
			The migration matrix is a backwards migration matrix where
			\f$ M_{ij} = 4N_0m_{ij} \f$ and \f$ m_{ij} \f$ is the fraction
			of subpopulation \f$ i \f$ made up of migrants <strong>from</strong>
			subpopulation \f$ j \f$ each generation.
			
			\a _m will be interpreted as the fraction
			of each subpopulation \f$ i \f$ made up of migrants <strong>from</strong>
			every other subpopulation \f$ j \neq i \f$ each generation (ie an inwards
			migration parameter).  All elements on the diagonal are
			ignored.
			 
			\param _m The migration parameter to use.
			\pre _mm >= 0.0.*/
			void setMigrationMatrix(double _m);
		
			
			/*! \brief Set the migration parameter between two specified 
			subpopulations.
			 
			The migration matrix is a backwards migration matrix where
			\f$ M_{ij} = 4N_0m_{ij} \f$ and \f$ m_{ij} \f$ is the fraction
			of subpopulation \f$ i \f$ made up of migrants <strong>from</strong>
			subpopulation \f$ j \f$ each generation.
			 
			\param into The index of the subpopulation to migration
			 is into.
			\param from The index of the subpopulation to migration
			 is from. 
			\param _m The value to set this migration parameter to.
			\pre \a _m >= 0, \a into is a valid index, 
			\a from is a valid index.*/
			void setMigration(int into, 
								int from,
								double _m);
			
			/*! \brief Set the migration parameter between two specified 
			subpopulations.
			 
			The migration matrix is a backwards migration matrix where
			\f$ M_{ij} = 4N_0m_{ij} \f$ and \f$ m_{ij} \f$ is the fraction
			of subpopulation \f$ i \f$ made up of migrants <strong>from</strong>
			subpopulation \f$ j \f$ each generation.
			 
			\param toLabel The label of the subpopulation migration
			 is into.
			\param fromLabel The label of the subpopulation migration
			 is from. 
			\param _m The value to set this migration parameter to.
			\pre \a _m >= 0, \a toLabel is a valid label, 
			\a fromLabel is a valid label.*/
			void setMigration(const std::string& toLabel, 
								const std::string& fromLabel,
								double _m);
			
			
			/*! \brief Get the migration parameter between two specified 
			subpopulations.
			
			The migration matrix is a backwards migration matrix where
			\f$ M_{ij} = 4N_0m_{ij} \f$ and \f$ m_{ij} \f$ is the fraction
			of subpopulation \f$ i \f$ made up of migrants <strong>from</strong>
			subpopulation \f$ j \f$ each generation.
			
			\param into The index of the sub-population to migration
			 is into.
			\param from The index of the sub-population to migration
			 is from.
			\pre \a into is a valid index, \a from is a valid index.*/
			double getMigration(int into, 
									int from) const;
			
			
			/*! \brief Get the migration parameter between two specified 
			subpopulations.
			
			The migration matrix is a backwards migration matrix where
			\f$ M_{ij} = 4N_0m_{ij} \f$ and \f$ m_{ij} \f$ is the fraction
			of subpopulation \f$ i \f$ made up of migrants <strong>from</strong>
			subpopulation \f$ j \f$ each generation.
			
			\param toLabel The label of the subpopulation migration
			 is into.
			\param fromLabel The label of the subpopulation migration
			 is from. 
			\pre \a toLabel is a valid label, \a fromLabel is a valid label.*/
			double getMigration(const std::string& toLabel, 
									const std::string& fromLabel) const;
			
			/*! \brief Set an event list of demographic events
			(backwards in time) for this population.
			
			Elements in the demographic event list may reference 
			subpopulations according to the existing
			subpopulation order:  no reordering will 
			be possible once the event list is set.
			
			\param _eList a collection of demographic events to become
			the demographic event list for this. */
			void setEventList(
				const std::vector< hudson_ms::HudsonMSDemogEventPtr >& _eList) const;
			
			/*! \brief Add to the list of demographic events
			(backwards in time) for this population.
			
			Elements in the demographic event list may reference 
			subpopulations according to the existing
			subpopulation order:  no reordering will 
			be possible once the event list contains at least one element.
			
			\param e a demographic event to add to the event list for this. */
			void addToEventList(const hudson_ms::HudsonMSDemogEventPtr& e) const;
			
			/*! \brief Get a string representation of the migration
			 matrix.
			
			The migration matrix is a backwards migration matrix where
			\f$ M_{ij} = 4N_0m_{ij} \f$ and \f$ m_{ij} \f$ is the fraction
			of subpopulation \f$ i \f$ made up of migrants <strong>from</strong>
			subpopulation \f$ j \f$ each generation.
			*/
			std::string migrationMatrixToString() const;
			
			/*! @name Make a an MS command string using this 
			 population structure.
			 
			 \param recombParams is an object able to provide
			 chromosone recombination (crossover and gene conversion)
			 parameters.
			 \param ssParams is an object able to provide
			 parameters controlling mutation.
			 \param seed is the seed for the random number generator.
			 \param time is an indicator for whether the 
			 ms command line switch -L (output to include time to most
			 recent common ancestor (MRCA) and total tree length) is
			 included in the resulting string.  Defaults to false
			 (-L switch not included). 
			 \return a string command which can be fed to Hudsons' MS
			 to generate an arg (and segsites if suitable
			 mutation parameters are contained in \a ssParams).*/
			//@{
			std::string makeMScommand(	
					int seed, bool time = false ) const;
			
			std::string makeMScommand(	
					const hudson_ms::HudsonMSRecombParams& recombParams,
					int seed, bool time = false ) const;
			
			std::string makeMScommand(	
					const hudson_ms::HudsonMSRecombParams& recombParams,
					const hudson_ms::HudsonMSSegSiteParams& ssParams,
					int seed, bool time = false) const;
			//@}

			/*! \brief String summary of this population structure.
			 
			 \return A string representation of this.*/
			std::string toString() const;
			
			/*! \brief Return an ordered collection of the subpopulation
			labels for the individuals in the population structure.
			
			if subpop i has label \f$ label_i \f$ and \f$ n_i \f$ individuals, \f$ i = 0,\,\ldots,\,m-1 \f$
			for \f$ m \f$ subpopulations, then this method returns an
			ordered collection of labels
			\f$ label_0^{0},\,\ldots,\,label_0^{[n_0-1]},\,\, label_1^{0},\,\ldots,\,label_1^{[n_1-1]},\,\,\dots,\,\,label_{m-1}^{0},\,\ldots,\,label_{m-1}^{[n_{m-1}-1]} \f$
			i.e. \f$ n_0 \f$ label_0 then \f$ n_1 \f$ label_1, then \f$ n_2 \f$ label_2 , etc, and finally \f$ n_{m-1} \f$ label_{m-1}.
					
			\return a collection of population labels for the subpopulation
			labels for the individuals in the population structure.	*/
			std::vector < std::string > 
						makePopulationLabels() const;
						
			/*! \brief Make a collection of subpopulation sample 
			sizes for this.
			
			\return an ordered collection of the number of samples from
			each subpopulation.
			\note uses int rather than double for compatibility with\
			ms-hacked code.*/
			std::vector < int > makePopnsams() const;
			
			/*! \brief Make a collection of subpopulation size ratios
			for this.
			
			\return an ordered collection of the size ratios for
			each subpopulation.*/
			std::vector < double > makePopsizeratios() const;
			
			/*! \brief Make a collection of subpopulation growth
			rates for this.
			
			\return an ordered collection of the growth rates for
			each subpopulation.*/
			std::vector < double > makePopgrowthrates() const;
			
			/*! \brief Make a collection of subpopulation growth
			rates for this.
			
			\param growthRate the growth rate to use for all
			subpopulations.
			\return an ordered collection of growth rates, each set to 
			\a growthRate, for each subpopulation.*/
			std::vector < double > makePopgrowthrates(double growthRate) const;
			
			/*! \brief Make a representation of the migration matrix
			for this.
			
			The migration matrix is a backwards migration matrix where
			\f$ M_{ij} = 4N_0m_{ij} \f$ and \f$ m_{ij} \f$ is the fraction
			of subpopulation \f$ i \f$ made up of migrants <strong>from</strong>
			subpopulation \f$ j \f$ each generation.
			
			\return an ordered collection of the migration rates where
			inner vectors represent rows in the migration matrix.*/
			std::vector< std::vector < double > > makeMigration_matrix() const;
			
			std::vector< hudson_ms::HudsonMSDemogEventPtr > makeEventList() const;
			
			/*! \brief Check the migration matrix for this.
			
			Clients can check whether the return value is empty, indicating 
			that the check has been passed, or use the string
			as an error message if the check has not been passed.
			
			\return a string which is empty if the migration matrix
			passes the test, or which summarises the problem if the
			migration matrix would result in
			infinite coalesence because there are non-communicating
			groups of subpopulations.*/
			std::string doMigCheck() const;
							
			
			
		protected:
		
					
						
		private:
		
		
			typedef std::map < std::string, double > migMap;
			
			/* A structure for describing a sub-population structure.
			 * 
			 * A sub-population structure contains:
			 * <ul>
			 * <li>Number of samples in that sub-population</li>
			 * <li>Size of subpopulation as a ratio of assumed N<sub>0</sub> size.
			 * <li>Exponential growth rate for subpopulation.
			 * <li>A string label.
			 * </ul>
			  
			*/
			struct SubPopulationStructure {
				
				size_t subNsam;
				
				double sizeRatio;
				
				double growthRate;
				
				std::string label;
				
				bool isDefaultLabel;
				
				migMap inwardMigration;
				
				SubPopulationStructure(size_t _sn, 
										double _sr,
										double _gr,
										std::string _label,
										bool _defLab);
				
				SubPopulationStructure(size_t _sn,
										double _gr,
										std::string _label,
										bool _defLab);
				
				SubPopulationStructure(size_t _sn,
										std::string _label,
										bool _defLab);
				
				SubPopulationStructure(std::string _label, //sn defaults to 0
										bool _defLab);
				
				SubPopulationStructure& operator=(SubPopulationStructure tmp);
				
				void swap(SubPopulationStructure& other);
					
				
				size_t getSubNsam() const;
				
				void setSubNsam(size_t _sn);

				double getSizeRatio() const;

				void setSizeRatio(double _sr);

				double getGrowthRate() const;
				
				void setGrowthRate(double _gr);
				
				std::string getLabel() const;
				
				void setLabel(const std::string& _label);
				
				bool getIsDefaultLabel() const;
				
				void setIsDefaultLabel(bool newSetting);
				
				// return 0.0 if there is no matching fromLabel
				double getInwardMig(const std::string& fromLabel) const;
				
				// change if fromLabel is there already, otherwise add
				void setInwardMig(const std::string& fromLabel, double mig);
				
				// do nothing if there is no matching fromLabel
				void eraseInwardMig(const std::string& fromLabel);
				
				/* do nothing if there is no matching oldFromLabel
				 * else throw an exception if there is already a newFromLabel
				 * and replace oldFromLabel with newFromLabel if there is not 
				 * already a newFromLabel.*/
 				void relabelInwardMig(
						const std::string& oldFromLabel,
						const std::string& newFromLabel);
						
				std::string toString() const;
				
			}; 
	
	
			/* Iterator into the ordered collection of parameter values.*/
			typedef std::vector< SubPopulationStructure >::const_iterator
											const_iterator;
		
			/* Return a const iterator to the first 
			SubPopulationStructure in the population.*/
			const_iterator begin() const;
			
			/* Return a const iterator to the
			'just past the end' SubPopulationStructure in the population.
			
			\note The idiom for iterating through a 
			%PopulationStructure pop is 
			for (PopulationStructure::const_iterator it = pop.begin();
						it < pop.end();
						++it) {...}.
			*/
			const_iterator end() const;
			
			std::vector < std::vector < double > >
										makeMigMatrixFromSubPops() const;
			
			
			std::set < size_t > & traceMig(
				std::set < size_t > & been,
				const std::vector < std::vector < double > >& migration_matrix, 
				size_t at) const;
			
			std::string migrationMatrixStringForMS() const;
			
			std::string eventStringForMS() const;
			
			// estimating precision for printing ms string
			int precEst() const;
			
			int precEst(double val) const;
			
			//get a default label for a new subpop given its order
			static std::string getNewDefaultLabel(int newOrder);
						
			
			void completeErase(	size_t orderErase);
			
			void completeRelabel(size_t order,
									const std::string& newLabel);
				
			// throw "label" if label not found			
			size_t getOrder(const std::string& label) const;
			
			bool checkUniqueLabel(size_t order,
								const std::string& label) const;
			
			static std::string& addToErrorMsg(std::string& msg,
												const std::exception& e);

			
			// outputing the bits
			std::ostringstream& outputStart(
								std::ostringstream& stm,
								size_t reps)	const;

								
			std::ostringstream& outputSubpops(
								std::ostringstream& stm)	const;


			
			/* The collection of SubPopulationStructure held.*/
			std::vector< SubPopulationStructure > container;
			
			mutable std::vector< hudson_ms::HudsonMSDemogEventPtr > eList;
			
	};

} 

/*! A specialisation of std::swap for mct::PopulationStructure types.*/
	namespace std
	{
		template <>
		void swap (hudson_ms::PopulationStructure & p1, 
				hudson_ms::PopulationStructure & p2);
	}
#endif 

