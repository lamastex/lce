
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

#ifndef _INC_RMV_H
#define _INC_RMV_H

#include <boost/function.hpp>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

namespace aabc {

	void rmvnorm(boost::function<double ()> ugaussian_func, const gsl_vector *mean, const gsl_matrix *var, gsl_vector *result);

	double dmvnorm(const gsl_vector *x, const gsl_vector *mean, const gsl_matrix *var);

	//int rmvt(const gsl_rng *r, const gsl_vector *location, const gsl_matrix *scale, const int dof, gsl_vector *result);

	//double dmvt(const gsl_vector *x, const gsl_vector *locationn, const gsl_matrix *scale, const int dof);
} // end namespace aabc

#endif
