/* MCT - Markov Chains on Trees.

   Copyright (C) 2009 Brendan Bycroft <brb44@student.canterbury.ac.nz>

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
\brief ModelSet declarations.
*/

#ifndef _INC_MODELSET_H
#define _INC_MODELSET_H

#include "model_data.hpp"
#include "model.hpp"
#include "mct_locus_map.hpp"

#include "prng.hpp"

#include <boost/smart_ptr.hpp>

#include <string>
#include <vector>


namespace mct {
	
	/** @brief A class managing an ordered collection of 
	\link mct::Model Models\endlink.
	  */
	 

	class ModelSet {
		public:
		
			/*! \brief Constructor.
			
			\param _mdata_list A collection of model configuration
			options specifying the models to use.
			\param _r_ptr A pointer to a random number generator.
			\note no checking of the configurations options in 
			\a _mdata_list is assumed to have taken place; any
			problems encountered in using these options for
			model creation will result in an exception being thrown
			from this constructor.	*/
			ModelSet(const std::vector<ModelData>& _mdata_list,
						const boost::shared_ptr < PRNGen >& _r_ptr);
			
			virtual ~ModelSet();
			
			/*! \brief Return a string representation of this. */
			virtual std::string toString() const;
			
			/*! \brief Return the number of Models in this set. */
			virtual size_t size() const;
			
			/*! \brief Return the number of sites covered by this. */
			virtual size_t getNsites() const;
			
			/*! \brief Returns the next state for the 
			\a index<sup>th</sup> Model in the collection given a 
			current state.
			
			\param index The index of the required Model in this. 
			\param current_state The current state to use to generate
			next state for the indexed Model.
			\param len The branch length between the parent state and 
			and the current state.
			\return The next state for the indexed Model given the
			\a current_state.*/
			virtual int step_at(const size_t index,
							const int current_state,
							const double len) const;

			/*! \brief Returns the ancestor state for the 
			\a index<sup>th</sup> Model in the collection.
			
			\param index The index of the required Model in this. 
			\return the state at time t=0 for the indexed Model.*/
			virtual int ancestor_at(const size_t index) const;
			
			/*! \brief Set the thetas for each  
			Model in the collection, separate theta for each model.
			
			\param new_thetas The new thetas to set the models to. 
			\pre thetas.size() == this.size().*/
			virtual void resetThetas(const std::vector < double >& new_thetas);

			/*! \brief Set the thetas for each  
			Model in the collection, same theta for each model.
			
			\param new_theta The new theta to set all the models to. */
			virtual void resetThetas(double new_theta);


			/*! \brief Return the map of the locus this models. */
			virtual const boost::shared_ptr <MCTLocusMap> 
												getLocusMap() const;
/*
			virtual std::string print_phylip_format(
						const std::vector<mct::site_seq> &sample_seqs) const;
*/
			
				

		protected:

		
		private:
			/** A container of Model pointers. */
			typedef std::vector< boost::shared_ptr < Model > > model_set;

			std::vector<ModelData> mdata_list; /** < Model data for chosen models.*/

			boost::shared_ptr < PRNGen > r_rand_ptr;
			
			model_set m_set; /** < one slot per model */
			model_set m_repeated_set; /** < one slot per site (for evolution of sequences) */
		
			boost::shared_ptr <MCTLocusMap> locusMap;
			
			ModelSet(const ModelSet& other);
			ModelSet& operator=(ModelSet tmp);
			
	};

} // end namespace mct

#endif /* _INC_MODELSET_H */

