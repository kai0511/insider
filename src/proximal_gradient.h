#ifndef __PROXIMAL_GRADIENT__
#define __PROXIMAL_GRADIENT__

#include <RcppArmadillo.h>


using namespace arma;

void prox_operator(vec& v, const double& theta);

void proximal_gradient(const mat& X, const vec& y, vec& beta, 
                       const double& lambda, const double& alpha, 
                       mat XtX, const vec& Xty, const double tol /*= 1e-5*/, const int max_iter /*= 100*/);


#endif