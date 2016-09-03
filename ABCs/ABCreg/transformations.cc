#include <transformations.hpp>
#include <cmath>
#include <gsl/gsl_math.h>

using namespace std;

inline double tan_transform( const double & x, const double & minval, const double & maxval)
//Apply the tangent transformation of Hamilton et al. 2005 PNAS 7476
{
  return -log(1./(tan( ((x-minval)/(maxval-minval))*(M_PI/2.) )));
}

inline double tan_untransform( const double & y, const double & minval, const double & maxval)
//undo the tangent transformation
{
  return minval + (2./M_PI)*(maxval-minval)*atan(exp(y));
}

double data_transform(const double & x, const params & p,const double  & minval, const double & maxval )
{
  if(p.transform_data)
    {
      if( p.transformation == params::LOG )
	{
	  return log(x);
	}
      else if (p.transformation == params::TAN)
	{
	  return tan_transform(x,minval,maxval);
	}
    }
  return x;
}

double data_untransform( const double & x, const params & p,  const double  & minval, const double & maxval )
{
  if(p.transform_data)
    {
      if( p.transformation == params::LOG )
	{
	  return exp(x);
	}
      else if (p.transformation == params::TAN)
	{
	  return tan_untransform(x,minval,maxval);
	}
    }
  return x;
}
