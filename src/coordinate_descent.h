#ifndef __COORDINATE_DESCENT__
#define __COORDINATE_DESCENT__

#include <RcppArmadillo.h>

using namespace arma;

void coordinate_descent(const mat& X, const vec& y, vec& beta, const double lambda, const double alpha, 
                        const mat& XtX, const vec& Xty, const double& tol /*= 1e-5*/);

void strong_coordinate_descent(const mat& X, const vec y, vec& beta, const double lambda, const double alpha, 
                               const mat& XtX, const vec& Xty, const double& tol /*= 1e-5*/);


#endif
