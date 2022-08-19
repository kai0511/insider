#ifndef __FEATURE_SIGN__
#define __FEATURE_SIGN__

#include <RcppArmadillo.h>

using namespace arma;

vec feature_sign_with_screening(const mat& X, const vec& y, const vec sugg_start, 
                                const double lambda, const double alpha, 
                                const mat& XtX, const vec& Xty, 
                                const int& max_iter /*= 1000*/);

#endif
