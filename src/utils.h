#ifndef __UTILS__
#define __UTILS__

#include <RcppArmadillo.h>

using namespace arma;

double objective(const mat& X, const vec& y, const vec& beta,
                 const double& lambda, const double& alpha);

double compute_loss(const vec& residual, const vec& beta, const double& lambda, const double& alpha);

void predict(const mat& row_factor, const mat& column_factor, mat& predictions);

void evaluate(mat& residual, const uvec& train_idx, const uvec& test_idx, double& sum_residual, 
              double& train_rmse, double& test_rmse, const int& tuning, 
              const int& iter /* = 0*/, const int& verbose /*= 1*/);

double compute_loss(const field<mat>& cfd_factor, const mat& column_factor, 
                    const double& lambda, const double& alpha, 
                    double& sum_residual, const int& verbose /*= 1*/);

#endif
