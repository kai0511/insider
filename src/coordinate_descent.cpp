// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins("cpp11")]]

#include <iostream>
#include <RcppArmadillo.h>
#include "coordinate_descent.h"
#include "utils.h"

using namespace arma;

void coordinate_descent(const mat& X, const vec& y, vec& beta, const double& lambda, const double& alpha, 
                        const mat& XtX, const vec& Xty, const double& tol = 1e-5){
    /* updates (20200829)
        Detail of optimization algorithm, refer to https://ieeexplore.ieee.org/document/6413853.
    Args:
        @beta warm start using the solution from the previous iteration
    */
    
    unsigned int i, k = 0;
    vec residual;
    double upper_scala, update, pre_loss, iter_loss;
    uvec ord(X.n_cols);

    // compute residual with warm start
    residual = y - X * beta;

    do {
        pre_loss = iter_loss;

        ord = randperm(X.n_cols);
        for(i = 0; i < beta.n_elem; i++){
            k = ord(i);

            upper_scala = dot(residual,  X.col(k)) + beta(k) * XtX(k, k);
            // right_scala = XtX(k, k) + lambda * (1 - alpha);

            // soft-thresholding
            if(std::abs(upper_scala) > lambda * alpha){
                update = sign(upper_scala) * std::max(std::abs(upper_scala) - lambda * alpha, 0.0) / (XtX(k, k) + lambda * (1 - alpha));
            } else{
                update = 0.0;
            }

            if(update != beta(k)){
                residual -= (update - beta(k)) * X.col(k);
                beta(k) = update;
		    }
        }
        iter_loss = compute_loss(residual, beta, lambda, alpha);
    } 
    while(std::abs(pre_loss - iter_loss) > tol);

    // return beta;
}

void strong_coordinate_descent(const mat& X, const vec y, vec& beta, const double& lambda, const double& alpha, 
                               const mat& XtX, const vec& Xty, const double& tol = 1e-5){
    /* update (20200829)
        1. Detail of optimization algorithm, refer to https://ieeexplore.ieee.org/document/6413853.
        2. screening rules to hasten computation(https://statweb.stanford.edu/~tibs/ftp/strong.pdf)
        3. add screening rules to filter out coefficients shrunk to zeros. 
    Args:
        @beta warm start using the solution from the previous iteration
    */
    
    unsigned int i, k;
    vec residual, grad, active_set = ones(beta.n_elem); 
    double upper_scala, update, pre_loss, iter_loss;
    uvec ord, inc_idx, violated_idx;

    // double c = alpha * (2 * lambda - max(abs(Xty)));
    // uvec ex_idx = find(abs(Xty) < c);
    uvec ex_idx = find(abs(Xty) < alpha * (2 * lambda - max(abs(Xty))));
    active_set(ex_idx).zeros();

    // suggested start
    beta.elem(ex_idx).zeros();
    residual = y - X * beta;
    iter_loss = compute_loss(residual, beta, lambda, alpha);

    while(true){
        inc_idx = find(active_set);
        ex_idx = find(active_set == 0);

        do {
            pre_loss = iter_loss;
            // pre_beta = beta;
            ord = randperm(inc_idx.n_elem);
            
            for(i = 0; i < inc_idx.n_elem; i++){
                k = inc_idx(ord(i));

                upper_scala = dot(residual,  X.col(k)) + beta(k) * XtX(k, k);
                // a trick to accelerate computation
                // right_scala = XtX(k, k) + lambda * (1 - alpha);

                // soft-thresholding
                if(std::abs(upper_scala) > lambda * alpha){
                    // update = sign(upper_scala) * std::max(std::abs(upper_scala) - lambda * alpha, 0.0) / right_scala;
                    update = sign(upper_scala) * std::max(std::abs(upper_scala) - lambda * alpha, 0.0) / (XtX(k, k) + lambda * (1 - alpha));
                } else{
                    update = 0.0;
                }
                
                if(update != beta(k)){
                    residual -= (update - beta(k)) * X.col(k);
                    beta(k) = update;
		        }
            }

            iter_loss = compute_loss(residual, beta, lambda, alpha);
        }
        while(std::abs(pre_loss - iter_loss) > tol);
        // while(accu(abs(pre_beta - beta)) > tol);

        // the key line to find the variable violating KKT condition
        grad = XtX(ex_idx, inc_idx) * beta(inc_idx) - Xty(ex_idx);
        violated_idx = find(abs(grad) > alpha * lambda);
        if(violated_idx.is_empty()){
            break;
        } else {
            active_set(ex_idx(violated_idx)).ones();
        }
    }
}
