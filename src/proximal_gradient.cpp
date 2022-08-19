// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins("cpp11")]]

#include "../inst/include/insider_types.h"
#include <iostream>
#include <RcppArmadillo.h>
#include <omp.h>
#include "proximal_gradient.h"
#include "utils.h"


void prox_operator(vec& v, const double& theta){
    // proximal operator for lasso

    // vec zero_vec = zeros<vec>(v.n_elem);
    // return (max(zero_vec, (v - theta)) - max(zero_vec, (- v - theta)));
    vec::row_iterator it = v.begin();
    for(; it != v.end(); ++it){
        *it = std::max(0.0, *it - theta) - std::max(0.0, - *it - theta);
    }
}

// [[Rcpp::export]]
vec proximal_gradient(const mat& X, const vec& y, const vec& wstart, const double& lambda, const double& alpha, 
                      const mat& XtX, const vec& Xty, const double& tol = 1e-5, const int& max_iter = 100){
    
    // reference see http://www.princeton.edu/~yc5/ele538b_sparsity/lectures/lasso_algorithm_extension.pdf
    // https://web.stanford.edu/~boyd/papers/prox_algs/lasso.html

    int col_num = X.n_cols, k = 1;
    double step_size, prev_loss, ridge_loss, loss, neg_eigen, pos_eigen;
    vec beta = wstart, pre_beta, old_beta, vals, gradient, pre_gradient;
    
    // suggested start
    prev_loss = objective(X, y, beta, lambda, alpha);
        
    vals = eig_sym(XtX);
    neg_eigen = vals(0);
    pos_eigen = vals(col_num-1);

    if(std::abs(neg_eigen) > std::abs(pos_eigen)){
        step_size = 1.0/neg_eigen;
    }else{
        step_size = 1.0/pos_eigen;
    }
    
    // XtX.diag() += (1 - alpha) * lambda;
    // pre_gradient = XtX * beta  - Xty;
    
    do { 
        pre_beta = beta; 

        if(k > 2){
            beta = pre_beta + (k - 2)/(k + 1)*(pre_beta - old_beta);
        }
        
        gradient = (XtX + (1 - alpha) * lambda * eye(size(XtX))) * beta  - Xty;
        beta -= step_size * gradient;
        prox_operator(beta, lambda * alpha * step_size);
        
        // loss = objective(X, y, beta, lambda, alpha);
        // cout << "Delta loss: " << prev_loss - loss << endl;
        // if(std::abs(prev_loss - loss) < tol){
        if(norm(pre_gradient, 1) * tol < norm(gradient, 1)){
            break;
        }
        old_beta = pre_beta;
        // prev_loss = loss;
        k++;
    }
    while(k <= max_iter);
    return beta;
}
