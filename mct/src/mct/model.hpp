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
\brief Model and Model subclass declarations.
*/

#ifndef _INC_MODEL_H
#define _INC_MODEL_H

#include "model_data.hpp"
#include "mct_model_data_types.hpp"

#include <iostream>
#include <string>

//#include <boost/bind.hpp>
#include <boost/function.hpp>
//#include <boost/any.hpp>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

namespace mct {

	/*! \brief Base class for Markov chain models.

	A Model applies a time-continuous markov process that 
	defines how the site state in a child
	node in a coalescence tree is simulated given the site state in its 
	parent node and the time between the nodes (branch lengths in the 
	tree).

	The site state is described by an integer code.
	
	Each derived class will have its own specific set of parameters.

	 */
	class Model {
		
		/*! \brief Output the model to the given stream.*/
		friend std::ostream& operator<<(std::ostream& s, 
											const Model& m);
			
		public:

			/*! \brief Destructor. */
			virtual ~Model();

			/*! \brief Determines the state of a site at time t+len, 
			 given the state at time t using the continuous time 
			 Markov chain model.
			 
			 The probability of going from the ith state at time t to 
			 the jth state at time t+len is given by
			 
			 P(t+len) = exp(Q/c * len * theta)
			 
			 where Q is the continuous-time Markov chain;
			 c is the scaling factor for Q;
			 and theta is the model's per site coalescent scale 
			 mutation rate (\f$ 4 N_0 \mu_{site} \f$).
			 
			 This function uses previously determined eigen values and 
			 vectors of Q, and c, to calculate the matrix exponential.
			 
			 \param current_state The current state of the site.
			 \param len The length of time between states.
			 \return the state at time t+len
			 */
			int step(int current_state, double len) const;

			/*! \brief Selects the state at time 0 based on the long term
			state of the model.
			
			\return the state at time t=0.
			*/
			int ancestor() const;

			/*! \brief Get the name of the Model.
			
			\return The Model name.
			*/
			std::string getModelName() const;
			
			/*! \brief Get whether the Model is DNA model.
			
			\return True if the Model is a DNA model, false otherwise.
			*/
			DATATYPE modelType() const;
			
			/*! \brief Get the lower bound of the Model.
			
			\return The lower bound of the Model.
			*/
			int getLowerBound() const;
			
			/*! \brief Reset theta for the Model.
			
			\param new_theta The new theta for the Model.
			*/
			void resetTheta(double new_theta);
			
			/*! \brief Get a string representation of the Model.
			
			\return A string representation of the Model.
			*/
			std::string toString() const;
			
			/*! \brief Get theta for the Model.
			
			\return The theta for the Model.
			*/
			double getTheta() const;
											
											

		protected:
		
			/*! \brief Constructor.
			 *
			 * This constructor can only be used by subclasses.
			 * 
			 \param _rand_num_func The function with signature 
			 double func(void). Must return a random number in the 
			 range [0, 1) with a uniform distribution.
			 \param _ss The sample space of the Model.
			 \param _name The name of the Model. This should correspond 
			 to the name of the subclass.
			 \param _dna The type of data the model deals in.
			 \param _lb The lower bound for the Model.
			 \param _th The per-site mutation rate \f$ \theta \f$
			 for the Model (\f$ 4 N_0 \mu_{site} \f$).
			 \param _param The model-specific parameters for the Model.
			 */
			Model(boost::function<double ()> _rand_num_func,
						int _ss,
						std::string& _name,
						DATATYPE _dna,
						int _lb,
						double _th,
						ModelData::ParamD& _param);

			/*
			Not sure what this function does, or how it relates to the standard gamma
			function.
			 */
			static double gamma_function(int i,int j, int ss, double m);

			double get_param(const std::string& name) const;

			/* Calculates the eigen values, eigen vectors and steady state
			 *         of Q.
			 * Q and steady state are then used to calculate c
			 * This function reads Model::Q, and modifies Model::Q_evec, Model::Q_eval,
			 * and the inverses of evec and eval,
			 * and Model::pi, Model::c.
			 *
			 */
			void calc_Q_eigen_and_pi_and_c();


			/* Calculates the eigen values, eigen vectors of Q.
			 * This function reads Model::Q, and modifies Model::Q_evec, Model::Q_eval,
			 * and the inverses of evec and eval.
			 */
			void calc_Q_eigen();
			
			/* Calculates steady state of Q.
			 * Q and steady state are then used to calculate c
			 * This function reads Model::Q, and modifies
			 * Model::pi, Model::c.
			 *
			 */
			void calc_pi_and_c();
			
			
			/* Picks a randomly chosen integer in the range [0,n), where n
			 *  is the length of vec, and each integer is weighted by its value in vec.
			 *
			 * vec is a vector of weights. The sum of its components must be 1.0
			 */
			int pick_from_freq_dist(const gsl_vector *vec) const;

			/* Calculates the embedded markov chain from a continuous time
			 *         markov chain.
			 *
			 * \param S The resulting embedded markov chain.
			 * */
			void calc_EMC_from_Q(gsl_matrix *S) const;

			
			// catches and swallows exceptions on freeing dynamic memory
			void freeDynamicMemory();
			
			// set default values for pointers
			void setDefaultPointers();
			

			
			/*! \brief Name of the model. */
			const std::string model_name; 
			
			/*! \brief State space of the model */
			const int ss;  //Defines the size of all matrices & vectors. */
			
			/*! \brief The type of data the model deals in.*/
			
			const DATATYPE DNA_model; 
			
			/*! \brief Lower bound for integer-coded states. */
			const int lower_bound; 
			
			/*! \brief Function which generates a uniform distribution
			in the range (0, 1]. */
			boost::function<double ()> rand_num_func;
			
			/*! \brief The \f$ \theta \f$ or per site 
			 continuous time mutation rate.
			 
			 \f$ \theta  = 4N_0 \mu_{site} \f$ 
			 where \f$ N_0 \f$ is the effective 
			 population size and \f$ \mu_{site} \f$ is the neutral
			 mutation rate per site.
			 
			 \note This data member is not const.  This allows the same 
			 generic model to be run with different mutatin rates
			 */
			double theta; 
			
			/*! \brief The supplied model-specific parameters. */
			const ModelData::ParamD params; 
			
			/*! \brief Mutation rate scaling factor.
			 * 
			 * Calculated as
			 * \f[ c = -\sum_{i=0}^{ss} Q_{ii}pi_i \f]
			 */
			double c; 

			/*! \brief  Continuous-time markov chain. */
			gsl_matrix* Q; 

			gsl_matrix* Q_evec; /*< Eigen-vectors of Q. */
			gsl_matrix* Q_evec_inv; /*< Eigen-vectors of Q inverse. */
			gsl_vector* Q_eval; /*< Eigen-values of Q. */
			
			/*! \brief  Steady state of Q. */
			gsl_vector* pi;
			
			/* Intermediary vectors used in step().
			 * Allocate once at construction and deallocate once
			 * at destruction to save lots of individual 
			 * allocations/deallocations every time we call step().  */
			gsl_vector* prob_vector;
			gsl_vector* scaled;
		

			double avgmr; /* TODO: No idea what this does. */

			
			
		private:
			Model(); //no no-argument constructor 
	};

	/* Jukes and Cantor 1969 Model.

	\note State space is always 4 for this Model.

	 */
	class ModelJC69 : public Model {
	public:
		ModelJC69(boost::function<double ()> _func,
							int _ss, 
							std::string& _name,
							DATATYPE _dna,
							int _lb,
							double _th,
							ModelData::ParamD& _param);
	private:
		void create_Q();
	};

	/* The Stepwise Mutation model. */
	class ModelSMM : public Model {
	public:
		ModelSMM(boost::function<double ()> _func,
							int _ss, 
							std::string& _name,
							DATATYPE _dna,
							int _lb,
							double _th,
							ModelData::ParamD& _param);
	private:
		void create_Q();
	};

	/*General model for nucleotide substitution subject
	to the time-reversibility constraint. */
	class ModelREV : public Model {
	public:
		ModelREV(boost::function<double ()> _func,
							int _ss, 
							std::string& _name,
							DATATYPE _dna,
							int _lb,
							double _th,
							ModelData::ParamD& _param);
	private:
		void create_Q();

		gsl_vector *r_mat;
	};

	/* AQB mode. */
	/* A and B marker loci, Q a locus influencing the expression
	 *  of a quantitative trait; AQB a genotype?*/
	class ModelAQB : public Model {
	public:
		ModelAQB(boost::function<double ()> _func,
							int _ss, 
							std::string& _name,
							DATATYPE _dna,
							int _lb,
							double _th,
							ModelData::ParamD& _param);
	private:
		void create_Q();

		double threskappa;
		double u0, u1, u2;
		double d0, d1, d2;
	};

	/* EXP mode. */
	class ModelEXP : public Model {
	public:
		ModelEXP(boost::function<double ()> _func,
							int _ss, 
							std::string& _name,
							DATATYPE _dna,
							int _lb,
							double _th,
							ModelData::ParamD& _param);
	private:
		void create_Q();

		double gammaU, gammaD;
		double alphaU, alphaD;
		double lambdaU, lambdaD;
	};

	/* EXP mode. */
	class ModelPL : public Model {
	public:
		ModelPL(boost::function<double ()> _func,
							int _ss, 
							std::string& _name,
							DATATYPE _dna,
							int _lb,
							double _th,
							ModelData::ParamD& _param);
	private:
		void create_Q();

		double u;
		double v;
		double m;
		double s;
		double p;
		double bdryPLo, bdryPHi;
	};

	std::ostream& operator<<(std::ostream& s, const gsl_matrix* m);
	std::ostream& operator<<(std::ostream& s, const gsl_vector* m);
	
} // end namespace mct

#endif /* _INC_MODEL_H */

