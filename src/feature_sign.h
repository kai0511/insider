#ifndef __FEATURE_SIGN__
#define __FEATURE_SIGN__

#include <RcppArmadillo.h>

using namespace arma;

vec strong_feature_sign(const mat& X, const vec& y, const vec& wstart, 
                        const double& lambda, const double& alpha, 
                        const mat& XtX, const vec& Xty, const unsigned int& max_iter /*= 1000*/);

#endif
