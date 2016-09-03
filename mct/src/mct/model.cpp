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
\brief Model and Model subclass definitions.
*/

#include "model.hpp"


#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#include <ostream>
#include <sstream>
#include <cmath>
#include <stdexcept>
#include <cassert>

//#define MYDEBUG
//#define MYDEBUGDESTRUCTOR

#if defined (MYDEBUG) || defined (MYDEBUGDESTRUCTOR)
	#include <iostream>
	
#endif

using namespace mct;

//obselete
#if 0
std::string model_types[] = {
    "REV", // * Revised
    "SMM", // Stepwise Mutation Model
    "AQB", // *
    "EXP", // * Exponential model
    "PL1", //
    "PL2", //
    ""};
#endif


std::ostream& mct::operator<<(std::ostream& s, const gsl_matrix* m)
{
    unsigned i, j;
    for (i=0; i<m->size1; ++i) {
        for (j=0; j<m->size2; ++j) {
            s.width(10);
            s.precision(4);
            s << gsl_matrix_get(m, i, j) << " ";
        }
        s << "\n";
    }
    return s;
}

std::ostream& mct::operator<<(std::ostream& s, const gsl_vector* m)
{
    unsigned j;
    for (j=0; j<m->size; ++j) {
        s.width(10);
        s.precision(4);
        s << gsl_vector_get(m, j) << " ";
    }
    s << "\n";
    return s;
}


std::ostream& mct::operator<<(std::ostream& s, const Model& m) {

    s << m.toString();
    
    if (m.Q->size1 < 6)
	s << "\nContinuous-time markov chain:\n" << m.Q;

    gsl_matrix* S = gsl_matrix_alloc(m.ss, m.ss);

    m.calc_EMC_from_Q(S);

    if (S->size1 < 6)
	s << "\nDiscrete-time markov chain:\n" << S;
    //if (m.pi->size < 6)
	s << "\nStationary Probability Distribution:\n" << m.pi;

    gsl_matrix_free(S);

    return s;
}

Model::Model(boost::function<double ()> _rand_num_func,
				int _ss,
				std::string& _name,
				DATATYPE _dna,
				int _lb,
				double _th,
				ModelData::ParamD& _param)
 : model_name(_name), ss(_ss), 
 DNA_model(_dna), lower_bound(_lb), 
 rand_num_func(_rand_num_func),
 theta(_th), params(_param) 
{
	try {
		setDefaultPointers();
		Q =           gsl_matrix_calloc(ss, ss);
		Q_evec =      gsl_matrix_alloc(ss, ss);
		Q_evec_inv =  gsl_matrix_alloc(ss, ss);
		Q_eval =      gsl_vector_alloc(ss);
		pi =          gsl_vector_alloc(ss);
		scaled =          gsl_vector_alloc(ss);
		prob_vector =          gsl_vector_alloc(ss);
		
	}
	catch (std::exception& e) {
		freeDynamicMemory();
		throw;
	}
}

Model::~Model()
{
	#ifdef MYDEBUGDESTRUCTOR
		std::cout << "Destroying model " << getModelName() << std::endl;
	#endif
	
	freeDynamicMemory();
}


std::string Model::getModelName() const
{
	return model_name;
}


DATATYPE Model::modelType() const
{
	return DNA_model;
}


int Model::getLowerBound() const
{
	return lower_bound;
}


double Model::getTheta() const
{
	return theta;
}

void Model::resetTheta(double new_theta)
{
	theta = new_theta;
}

std::string Model::toString() const
{
	std::ostringstream stm;

    stm << "Model: " << model_name << "\n";
	
	stm << "ss: " << ss << "\t lower_bound: " << lower_bound
		<< "\ttheta: "<< theta;
	
	for (ModelData::ParamD::const_iterator pit = params.begin();
			pit != params.end();
			++pit) {
		if (params.begin()== pit) stm << "\n";
		stm << pit->first << ": " << pit->second << "\t";
		
	}
	stm << "\n";
	
	return stm.str();
}
    
void Model::calc_EMC_from_Q(gsl_matrix *S) const
{
    int n = Q->size1;
    int i, j;
    double total_inv, val;
    for (i=0; i<n; ++i) {
        total_inv = -1.0/gsl_matrix_get(Q, i, i);
        for (j=0; j<n; ++j) {
            if (i != j) {
                val = gsl_matrix_get(Q, i, j);
                gsl_matrix_set(S, i, j, val*total_inv);
            } else {
                gsl_matrix_set(S, i, j, 0);
            }
        }
    }

}


void Model::calc_Q_eigen_and_pi_and_c()
{
	Model::calc_Q_eigen();
	Model::calc_pi_and_c();
}

void Model::calc_Q_eigen()
{
	gsl_matrix* dummy = NULL;
	gsl_matrix_complex* evec_temp = NULL;
	gsl_vector_complex* eval_temp = NULL;
	gsl_eigen_nonsymmv_workspace* ws = NULL;
	gsl_permutation* perm = NULL;
		
	try {
		int signum = 0;
		double temp = 0.0;

		/* calculate the eigen values and eigen vectors of Q */
		
		 ws = gsl_eigen_nonsymmv_alloc(ss);
		dummy = gsl_matrix_alloc(ss, ss);
		evec_temp = gsl_matrix_complex_alloc(ss, ss);
		eval_temp = gsl_vector_complex_alloc(ss);
		perm = gsl_permutation_alloc(ss);


		gsl_matrix_memcpy(dummy, this->Q);
		gsl_eigen_nonsymmv(dummy, eval_temp, evec_temp, ws);
		
		// evals and evecs are all real. Convert them to real from complex
		for (int i = 0; i < ss; ++i) {
			for (int j = 0; j < ss; ++j) {
				temp = GSL_REAL(gsl_matrix_complex_get(evec_temp, i, j));
				gsl_matrix_set(this->Q_evec, i, j, temp);
			}
			temp = GSL_REAL(gsl_vector_complex_get(eval_temp, i));
			gsl_vector_set(this->Q_eval, i, temp);
		}
		
		/* calculating the inverse of Q_evec */

		gsl_matrix_memcpy(dummy, this->Q_evec);
		gsl_linalg_LU_decomp(dummy, perm, &signum);
		gsl_linalg_LU_invert(dummy, perm, this->Q_evec_inv);
		
		gsl_eigen_nonsymmv_free(ws);
		ws = NULL;
		gsl_matrix_free(dummy);
		dummy = NULL;
		gsl_matrix_complex_free(evec_temp);
		evec_temp = NULL;
		gsl_vector_complex_free(eval_temp);
		eval_temp = NULL;
		gsl_permutation_free(perm);
		perm = NULL;
		
	}
	// catch any exceptions, try to deal with all memory allocations,
	catch (...) {
		
		try { if (ws != NULL) gsl_eigen_nonsymmv_free(ws); }
		catch (...) {}
		try { if (dummy != NULL) gsl_matrix_free(dummy); }
		catch (...) {}
		try { if (evec_temp != NULL) gsl_matrix_complex_free(evec_temp); }
		catch (...) {}
		try { if (eval_temp != NULL) gsl_vector_complex_free(eval_temp); }
		catch (...) {}
		try { if (perm != NULL) gsl_permutation_free(perm); }
		catch (...) {}
		
		throw;
	}
}


void Model::calc_pi_and_c()
{
	gsl_matrix* U = NULL;
	gsl_matrix* V = NULL;
	gsl_vector* S = NULL;
	gsl_vector* work = NULL;
		
	try {
		double length = 0.0;

		U = gsl_matrix_alloc(ss, ss);
		V = gsl_matrix_alloc(ss, ss);
		S = gsl_vector_alloc(ss);
		work = gsl_vector_alloc(ss);

		/* Calculate the steady state, pi, from the Q matrix. */

		gsl_matrix_memcpy(U, this->Q);
		// U is input mtx Q, and also output U
		gsl_linalg_SV_decomp(U, V, S, work);
		// null space of Q is last col of U.
		gsl_matrix_get_col(this->pi, U, ss-1);

		length = 0.0;
		for (int i=0; i<ss; ++i)
			length += gsl_vector_get(this->pi, i);

		// U has norm(U, 2) = 1. We want norm(U, 1) = 1
		gsl_vector_scale(this->pi, 1.0 / length);

		gsl_matrix_free(U); 
		U = NULL;
		gsl_matrix_free(V);
		V = NULL;
		gsl_vector_free(S);
		S = NULL;
		gsl_vector_free(work);
		work = NULL;
		
		// now do c
		double temp = 0.0;
		//sum of Q[ii]*pi[i] over 0 <= i <= ss
		for (int i = 0; i < ss; ++i) {
			temp += gsl_matrix_get(this->Q, i, i) 
					* gsl_vector_get(this->pi, i);
		}
		this->c = -temp; // take the negative
	}

	// catch any exceptions, try to deal with all memory allocations,
	catch (...) {
		try { if (U != NULL) gsl_matrix_free(U); }
		catch (...) {}
		try { if (U != NULL) gsl_matrix_free(U); }
		catch (...) {}
		try { if (V != NULL) gsl_matrix_free(V); }
		catch (...) {}
		try { if (S != NULL) gsl_vector_free(S); }
		catch (...) {}
		try { if (work != NULL) gsl_vector_free(work); }
		catch (...) {}
		
		throw;
	}
}



int Model::pick_from_freq_dist(const gsl_vector *vec) const
{
    unsigned i;
    double accum = 0.0;
    double rand_num = rand_num_func();
   
    assert(
		abs(gsl_blas_dasum(vec) - 1.0) < 1e-14);
	
    accum = 0.0;
    for (i=0; i < vec->size; ++i) {
        accum += gsl_vector_get(vec, i);
        if (accum >= rand_num) break; // break out of loop
            
    }
	if (i >= vec->size) {
		throw std::logic_error(
		"Model::pick_from_freq_dist(const gsl_vector*)");
	}
	return i;
}

int
Model::step(int current_state_code, double len) const
{
    /* If j = current_state, then this function calculates
     *
     * prob_std::vector = exp(Q * len)[j, :]
     *
     * The matrix exponentiation is done using the eigen method, where the
     * eigen std::vectors and values have been precomputed. Only multiplication
     * is neccessary.
     *
     * exp(Q*t)  =  Q_evec * exp(Q_eval * len) * Q_evec_inv
     *
     * This is simplified further by only doing the calculations neccessary
     * to compute a single row of exp(Q * len).
     *
     * prob_std::vector = Q_evec[j, :] * exp(Q_eval * len) * Q_evec_inv
     *
     */
	
	#ifdef MYDEBUG
		std::cout << "In step" << std::endl;
		std::cout << "current state code is " << current_state_code << ", lower_bound is " << lower_bound << " current_state_code - lower_bound " << current_state_code - lower_bound << std::endl;
	#endif

	
	int current_state = current_state_code - lower_bound;
	assert (current_state >= 0);
	
	// scaled[i] = exp(Q_eval[i] * len * theta/c)
	for (int i=0; i<ss; ++i) {
		gsl_vector_set(scaled, i, 
			exp(gsl_vector_get(Q_eval, i) * len * theta/c));
	}

	// scaled[i] = scaled[i] * Q_evec[current_state, i]
	gsl_vector_mul(scaled, &gsl_matrix_const_row(Q_evec, current_state).vector);

	// prob_vector = scaled * Q_evec_inv  (vector/matrix mult)
	gsl_blas_dgemv(CblasTrans, 1.0, Q_evec_inv, scaled, 0.0, prob_vector);
	
	
	return (pick_from_freq_dist(prob_vector) + lower_bound);
}

int
Model::ancestor() const
{
		return (pick_from_freq_dist(this->pi) + lower_bound);
}

/*
Brendan's note:
Not sure what this function does, or how it relates to the standard gamma
function.
 */
double Model::gamma_function(int i,int j, int ss, double m)
{
    double gamma = 0.0;

    if(i < j && j <= ss)
        gamma = m * gsl_pow_int(1 - m, abs(i - j) - 1) / (1 - gsl_pow_int(1 - m, ss - i));
    else if(j < i && i <= ss)
        gamma = m * gsl_pow_int(1 - m, abs(i - j) - 1) / (1 - gsl_pow_int(1 - m, i));

    return gamma;
}


double Model::get_param(const std::string& name)  const
{
    ModelData::ParamD::const_iterator it = params.find(name);

    if (it == params.end()) {
		throw std::invalid_argument("Model::get_param(const std::string&)");
	}
        
    return it->second;
}

// set default values for pointers
void Model::setDefaultPointers()
{
	Q = NULL;
	Q_evec = NULL;
	Q_evec_inv = NULL;
	Q_eval = NULL;
	pi = NULL;
	prob_vector = NULL;
	scaled = NULL;
	
}


// catches and swallows exceptions on freeing dynamic memory
void Model::freeDynamicMemory()
{
	try {
		if (Q != NULL) gsl_matrix_free(Q);
	}
	catch (...) {}	
	try {
		if (Q_evec != NULL) gsl_matrix_free(Q_evec);
	}
	catch (...) {}	
	try {
		if (Q_evec_inv != NULL) gsl_matrix_free(Q_evec_inv);
	}
	catch (...) {}	
	try {
		if (Q_eval != NULL) gsl_vector_free(Q_eval);
	}
	catch (...) {}	
	try {
		if (pi != NULL) gsl_vector_free(pi);
	}
	catch (...) {}	
	try {
		if (prob_vector != NULL) gsl_vector_free(prob_vector);
	}
	catch (...) {}	
	try {
		if (scaled != NULL) gsl_vector_free(scaled);
	}
	catch (...) {}	
	
}



ModelJC69::ModelJC69(boost::function<double ()> _func,
						int _ss, 
						std::string& _name,
						DATATYPE _dna,
						int _lb,
						double _th,
						ModelData::ParamD& _param)
 : Model(_func, (_ss == 4 ? _ss : 4), _name, _dna, _lb, _th, _param ) 
{
    /* No parameters. */
	
    create_Q();
    Model::calc_Q_eigen_and_pi_and_c();
}

void ModelJC69::create_Q()
{
    int i, j;

    //set non-diagonals to 0.25
    for (i=0; i<ss; ++i) {
        for (j=0; j<ss; ++j) {
            if (i != j)
                gsl_matrix_set(Q, i, j, 0.25);
        }
    }

    /* q_ii = -sum(q_ij, j!=i) */
    for (i=0; i<ss; ++i) {
        double sum = 0.0;
        for (j=0; j<ss; ++j) {
            if (j != i)
                sum += gsl_matrix_get(Q, i, j);
        }
        gsl_matrix_set(Q, i, i, -sum);
    }
}


ModelSMM::ModelSMM(boost::function<double ()> _func,
						int _ss, 
						std::string& _name,
						DATATYPE _dna,
						int _lb,
						double _th,
						ModelData::ParamD& _param)
 : Model(_func, _ss, _name, _dna, _lb, _th, _param ) 
 {
    /* No parameters. */

    create_Q();
    Model::calc_Q_eigen_and_pi_and_c();
}

void ModelSMM::create_Q()
{
    int i, j;

    //set off-diagonals to 0.5
    for (i=0; i<ss; ++i) {
        for (j=0; j<ss; ++j) {
            if (abs(i - j) == 1)
                gsl_matrix_set(Q, i, j, 0.5); /* *freq[j] */
        }
    }

    /* q_ii = -sum(q_ij, j!=i) */
    for (i=0; i<ss; ++i) {
        double sum = 0.0;
        for (j=0; j<ss; ++j) {
            if (j != i)
                sum += gsl_matrix_get(Q, i, j);
        }
        gsl_matrix_set(Q, i, i, -sum);
    }
}

ModelREV::ModelREV(boost::function<double ()> _func,
						int _ss, 
						std::string& _name,
						DATATYPE _dna,
						int _lb,
						double _th,
						ModelData::ParamD& _param)
 : Model(_func, _ss, _name, _dna, _lb, _th, _param ) 
{
    /* No parameters */

    /* number of non-zero elements in an upper triangular matrix of size
       (ss x ss) (Wiki says 1+2+..+n = (ss(ss+1))/2 = ss*ss/2 + ss/2 )
       FIXME ? */
    int up_tri = ss*ss-ss/2;

    r_mat = gsl_vector_alloc(up_tri);

    if (1) { /* default parameters */
        int i;
        for (i=0; i<up_tri; ++i)
            gsl_vector_set(r_mat, i, 1.0);
        avgmr = 1.0;
    }

    create_Q();
    Model::calc_Q_eigen_and_pi_and_c();
}

void ModelREV::create_Q()
{
    int i, j, k=0;
    int last = (ss - 1) * ss;

    for (i=0; i<(ss-2); ++i) {
        for (j=i+1; j<ss; ++j) {
            if (i * ss + j != last) {
                gsl_matrix_set(Q, i, j, gsl_vector_get(r_mat, k));
                gsl_matrix_set(Q, j, i, gsl_vector_get(r_mat, k));
                ++k;
            }
        }
    }

    gsl_matrix_set(Q, ss-1, ss-2, 1.0);
    gsl_matrix_set(Q, ss-1, ss-1, 1.0);

    for (i=0; i<ss; ++i) {
        for (j=0; j<ss; ++j) {
            gsl_matrix_set(Q, i, j, 1.0); /* *freq[j] */
        }
    }

    /* q_ii = -sum(q_ij, j!=i) */
    for (i=0; i<ss; ++i) {
        double row_sum = 0.0;
        for (k=0; k<ss; ++k) {
            if (k != i) {
                row_sum += gsl_matrix_get(Q, i, k);
            }
        }
        gsl_matrix_set(Q, i, i, -row_sum);
    }
}


ModelAQB::ModelAQB(boost::function<double ()> _func,
						int _ss, 
						std::string& _name,
						DATATYPE _dna,
						int _lb,
						double _th,
						ModelData::ParamD& _param)
 : Model(_func, _ss, _name, _dna, _lb, _th, _param) 
{
    u0 = get_param("u0");
    u1 = get_param("u1");
    u2 = get_param("u2");
    d0 = get_param("d0");
    d1 = get_param("d1");
    d2 = get_param("d2");
    threskappa = get_param("threskappa");
    avgmr = get_param("avgmr");

    #if 0
    if (1) { /* default parameters */
        u0 = 0.5;
        u1 = 0.005;
        u2 = 0.03;
        d0 = 0.5;
        d1 = 0.005;
        d2 = 0.03;
        threskappa = 0.0;
        avgmr = 1.0;
    }
    #endif

    create_Q();
    Model::calc_Q_eigen_and_pi_and_c();
}

void ModelAQB::create_Q()
{
    int i, j, k;

    for (i=0; i<ss; ++i) {
        for (j=0; j<ss; ++j) {
            double temp = std::max(i - threskappa, 0.0);
            if (i < j)
                gsl_matrix_set(Q, i, j, u0 + u1 * temp + gsl_pow_2(u2 * temp));
            else if (i > j)
                gsl_matrix_set(Q, i, j, d0 + d1 * temp + gsl_pow_2(d2 * temp));
            else if (i == j)
                gsl_matrix_set(Q, i, j, 0.0); /* redundant */
        }
    }

    /* q_ii = -sum(q_ij, j!=i) */
    for (i=0; i<ss; ++i) {
        double row_sum = 0.0;
        for (k=0; k<ss; ++k) {
            if (k != i)
                row_sum += gsl_matrix_get(Q, i, k);
        }
        gsl_matrix_set(Q, i, i, -row_sum);
    }
}

ModelEXP::ModelEXP(boost::function<double ()> _func,
						int _ss, 
						std::string& _name,
						DATATYPE _dna,
						int _lb,
						double _th,
						ModelData::ParamD& _param)
 : Model(_func, _ss, _name, _dna, _lb, _th, _param ) 
{
    gammaU = get_param("gammaU");
    gammaD = get_param("gammaD");
    alphaU = get_param("alphaU");
    alphaD = get_param("alphaD");
    lambdaU = get_param("lambdaU");
    lambdaD = get_param("lambdaD");
    avgmr = get_param("avgmr");


    #if 0
    if (1) { /* default parameters */
        gammaU = 3.1e-6;
        gammaD = 4e-7;
        alphaU = 0.2;
        alphaD = 0.302;
        lambdaU = 1.06;
        lambdaD = 1.06;
        avgmr = 1.0;
    }
    #endif

    create_Q();
    Model::calc_Q_eigen_and_pi_and_c();
}

void ModelEXP::create_Q()
{
    int i, j, k;

    for (i=0; i<ss; ++i) {
        for (j=0; j<ss; ++j) {
            if (i < j)
                gsl_matrix_set(Q, i, j, gammaU * exp(alphaU * i) * exp(-lambdaU * (j - i)));
            else if (i > j)
                gsl_matrix_set(Q, i, j, gammaD * exp(alphaD * i) * exp(-lambdaD * (i - j)));
            else if (i == j)
                gsl_matrix_set(Q, i, j, 0.0); /* redundant */
        }
    }

    /* q_ii = -sum(q_ij, j!=i) */
    for (i=0; i<ss; ++i) {
        double row_sum = 0.0;
        for (k=0; k<ss; ++k) {
            if (k != i)
                row_sum += gsl_matrix_get(Q, i, k);
        }
        gsl_matrix_set(Q, i, i, -row_sum);
    }
}

ModelPL::ModelPL(boost::function<double ()> _func,
						int _ss, 
						std::string& _name,
						DATATYPE _dna,
						int _lb,
						double _th,
						ModelData::ParamD& _param)
 : Model(_func, _ss, _name, _dna, _lb, _th, _param )
{
    u = get_param("u");
    v = get_param("v");
    m = get_param("m");
    s = get_param("s");
    avgmr = get_param("avgmr");
    p = get_param("p");
    bdryPLo = get_param("bdryPLo");
    bdryPHi = get_param("bdryPHi");


    #if 0
    if (1) { /* default parameters (PL1) */
        u = 0.6246;
        v = 0.01542;
        m = 1.0;
        s = 0.8752;
        avgmr = 1.0;
    }
    else if (2) { /* default parameters (PL2) */
        u = 0.8158;
        v = 0.03947;
        m = 0.5475;
        s = 0.7638;
        avgmr = 0.761;
    }

    p = 0.0;
    bdryPLo = 0.5;
    bdryPHi = 0.5;
    #endif


    create_Q();
    Model::calc_Q_eigen_and_pi_and_c();
}

void ModelPL::create_Q()
{
    int i, j, k;

    double p2 = (1 - p);
    double beta, alpha, gamma;

    for (i=0; i<ss; ++i) {
        beta = avgmr * (1 + i * s);
        alpha = std::max(0.0, std::min(1.0, u - v * i));
        for (j=0; j<ss; ++j) {
            gamma = gamma_function(i, j, ss, m);
            if (j - i == 1)
                gsl_matrix_set(Q, i, j, beta * alpha * (p + p2 * gamma));
            else if (i < j - 1)
                gsl_matrix_set(Q, i, j, beta * alpha * p2 * gamma);
            else if (i == j + 1)
                gsl_matrix_set(Q, i, j, beta * (1.0 - alpha) * (p + p2 * gamma));
            else if (i > j + 1)
                gsl_matrix_set(Q, i, j, beta * (1.0 - alpha) * p2 * gamma);
            else if (i == j)
                gsl_matrix_set(Q, i, j, 0.0); /* redundant */
        }
    }

    double *d_ptr;
    d_ptr = gsl_matrix_ptr(Q, 0, 1);
    *d_ptr += bdryPLo;

    d_ptr = gsl_matrix_ptr(Q, ss-1, ss-2);
    *d_ptr += bdryPHi;

    /* q_ii = -sum(q_ij, j!=i) */
    for (i=0; i<ss; ++i) {
        double row_sum = 0.0;
        for (k=0; k<ss; ++k) {
            if (k != i)
                row_sum += gsl_matrix_get(Q, i, k);
        }
        gsl_matrix_set(Q, i, i, -row_sum);
    }
}
