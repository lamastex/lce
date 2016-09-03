/***************************************************************************************
 *  Multivariate Normal density function and random number generator
 *  Multivariate Student t density function and random number generator
 * 
 *  Using GSL -> www.gnu.org/software/gsl
 *
 *  Based on the original by Ralph dos Santos Silva
 *  see http://www.mail-archive.com/help-gsl@gnu.org/msg00631.html
 *  with minor adaptations by Jenny Harlow, 2011
 *  Copyright (C) 2006  Ralph dos Santos Silva
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 *  AUTHOR 
 *     Ralph dos Santos Silva,  [EMAIL PROTECTED]
 *     March, 2006       
***************************************************************************************/
#include "rmv.h"

#include <cmath>

#include <boost/function.hpp>

#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>



using namespace aabc;

/*****************************************************************************************************************/
/*****************************************************************************************************************/
void aabc::rmvnorm(boost::function<double ()> ugaussian_func,
							const gsl_vector *mean,
							const gsl_matrix *var,
							gsl_vector *result)
{
/* multivariate normal distribution random number generator */
/*
*	mean	vector of means of size n
*	var	variance matrix of dimension n x n
*	result	output variable with a sigle random vector normal distribution generation
*/
	gsl_matrix *work = NULL;
	try {
		int n = mean->size; 
		int k;
		work = gsl_matrix_alloc(n,n);

		gsl_matrix_memcpy(work,var);
		
		gsl_linalg_cholesky_decomp(work);

		for(k=0; k<n; ++k)
			gsl_vector_set( result, k, ugaussian_func() );

		gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, work, result);
		gsl_vector_add(result,mean);

		gsl_matrix_free(work);
		work = NULL;
	}
	catch(...) {
		try {
			if (work != NULL) gsl_matrix_free(work);
		}
		catch(...) {}
		try {
			int n = result->size;
			gsl_vector_free(result);
			result = gsl_vector_calloc(n);
		}
		catch(...) {}
		throw;
	}

}

/*****************************************************************************************************************/
/*****************************************************************************************************************/
double aabc::dmvnorm(const gsl_vector *x, const gsl_vector *mean, const gsl_matrix *var){
/* multivariate normal density function    */
/*
*	x	vector to calculate pdf of
*	mean	vector of means of size n
*	var	variance matrix of dimension n x n
*/
	gsl_vector *ym = NULL;
	gsl_vector *xm = NULL;
	gsl_matrix *work = NULL;
	gsl_matrix *winv = NULL;
	gsl_permutation *p = NULL;

	try {
		int n = mean->size; 
		int s;
		double ax,ay;
		work = gsl_matrix_alloc(n,n), 
		winv = gsl_matrix_alloc(n,n);
		p = gsl_permutation_alloc(n);

		gsl_matrix_memcpy( work, var );
		gsl_linalg_LU_decomp( work, p, &s );
		gsl_linalg_LU_invert( work, p, winv );
		ax = gsl_linalg_LU_det( work, s );
		gsl_matrix_free( work );
		work = NULL;
		gsl_permutation_free( p );
		p = NULL;

		xm = gsl_vector_alloc(n);
		gsl_vector_memcpy( xm, x);
		gsl_vector_sub( xm, mean );
		ym = gsl_vector_alloc(n);
		gsl_blas_dsymv(CblasUpper,1.0,winv,xm,0.0,ym);
		gsl_matrix_free( winv );
		winv = NULL;
		gsl_blas_ddot( xm, ym, &ay);
		gsl_vector_free(xm);
		xm = NULL;
		gsl_vector_free(ym);
		ym = NULL;
		
		ay = exp(-0.5*ay)/sqrt( pow((2*M_PI),n)*ax );

		return ay;
	}
	catch (...) {
		try { if (xm != NULL) gsl_vector_free(xm); }
		catch (...) {}
		try { if (ym != NULL) gsl_vector_free(ym); }
		catch (...) {}
		try { if (work != NULL) gsl_matrix_free(work); }
		catch (...) {}
		try { if (winv != NULL) gsl_matrix_free(winv); }
		catch (...) {}
		try { if (p != NULL) gsl_permutation_free(p); }
		catch (...) {}
		
		throw;
	}
}

/*****************************************************************************************************************/
/*****************************************************************************************************************/

#if(0)
void rmvt(const gsl_rng *r,
		const gsl_vector *location,
		const gsl_matrix *scale,
		const int dof,
		gsl_vector *result)
{
/* multivariate Student t distribution random number generator */
/*
*	location vector of locations of size n
*	scale	 scale matrix of dimension n x n
*	dof	 degrees of freedom
*	result	 output variable with a single random vector normal distribution generation
*/

	int n = location->size; 
	int k;
	gsl_matrix *work = gsl_matrix_alloc(n,n);
	double ax = 0.5*dof; 

	ax = gsl_ran_gamma(r,ax,(1/ax));     // gamma distribution

	gsl_matrix_memcpy(work,scale);
	gsl_matrix_scale(work,(1/ax));       // scaling the matrix
	gsl_linalg_cholesky_decomp(work);

	for(k=0; k<n; ++k)
		gsl_vector_set( result, k, gsl_ran_ugaussian(r) );

	gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, work, result);
	gsl_vector_add(result, location);

	gsl_matrix_free(work);

}

/*****************************************************************************************************************/
/*****************************************************************************************************************/

double dmvt(const gsl_vector *x, const gsl_vector *location, const gsl_matrix *scale, const int dof){
/* multivariate Student t density function */
/*
*	x vector to calculate pdf of
*	location vector of locations of size n
*	scale	 scale matrix of dimension n x n
*	dof	 degrees of freedom
*/
	
	int n = location->size; 
	int s;
	double ax,ay,az=0.5*(dof + n);
	gsl_vector *ym, *xm;
	gsl_matrix *work = gsl_matrix_alloc(n,n), 
			   *winv = gsl_matrix_alloc(n,n);
	gsl_permutation *p = gsl_permutation_alloc(n);

	gsl_matrix_memcpy( work, scale );
	gsl_linalg_LU_decomp( work, p, &s );
	gsl_linalg_LU_invert( work, p, winv );
	ax = gsl_linalg_LU_det( work, s );
	gsl_matrix_free( work );
	gsl_permutation_free( p );

	xm = gsl_vector_alloc(n);
	gsl_vector_memcpy( xm, x);
	gsl_vector_sub( xm, location );
	ym = gsl_vector_alloc(n);
	gsl_blas_dsymv(CblasUpper,1.0,winv,xm,0.0,ym);
	gsl_matrix_free( winv );
	gsl_blas_ddot( xm, ym, &ay);
	gsl_vector_free(xm);
	gsl_vector_free(ym);

	ay = std::pow((1+ay/dof),-az)*gsl_sf_gamma(az)/(gsl_sf_gamma(0.5*dof)*std::sqrt( std::pow((dof*M_PI),n)*ax ));

	return ay;
}
#endif
