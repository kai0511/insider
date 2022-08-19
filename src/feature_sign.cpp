// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins("cpp11")]]

#include <iostream>
#include <RcppArmadillo.h>
#include <omp.h>
#include "feature_sign.h"
#include "utils.h"

using namespace arma;

vec strong_feature_sign(const mat& X, const vec& y, const vec& wstart, 
                        const double& lambda, const double& alpha, 
                        const mat& XtX, const vec& Xty, const int& max_iter = 1000){
    /* update (20220610)
        adopt screening rules to filter out coefficients that shrink to zeros.
        for details, see https://statweb.stanford.edu/~tibs/ftp/strong.pdf. 
        use analytical solution to the elastic net.
    */
    
    unsigned int k = 0, iter = 0, col_num = X.n_cols, optimality1 = 0, optimality0 = 0;
    double pre_loss = 0.0, inner_loss = 0.0, line_search = 0.0;
    
    uvec exc_idx, inc_idx, search_idx, violate_idx, remove_idx;
    vec tmp_beta, new_beta, grad, progress;
    vec beta = wstart, active_set = ones(col_num), theta = sign(beta);
    
    mat X2, Gram;

    double cutoff = alpha * (2 * lambda - max(abs(Xty)));
    exc_idx = find(abs(Xty) < cutoff);
    active_set.elem(exc_idx).zeros();
    beta.elem(exc_idx).zeros();

    while(iter < max_iter){
        
        while(true) {
            line_search = 0.0;
            inc_idx = find(active_set);
            X2 = X.cols(inc_idx);

            Gram = XtX(inc_idx, inc_idx) + (1 - alpha) * lambda * eye(inc_idx.n_elem, inc_idx.n_elem);
            new_beta = solve(Gram, Xty(inc_idx) - lambda * alpha * theta(inc_idx), solve_opts::likely_sympd);

            if (!any(sign(new_beta) - theta(inc_idx))){
                beta(inc_idx) = new_beta;
                optimality1 = 1;
                line_search = 1.0;
                break; 
            }

            progress = -beta(inc_idx)/(new_beta - beta(inc_idx));
            progress.insert_rows(progress.n_elem, 1);
            

            // a = sum(square(X.cols(inc_idx) * (beta(inc_idx) - pre_beta(inc_idx))))/2;
            // b = trans(beta(inc_idx) - pre_beta(inc_idx)) * (XtX(inc_idx, inc_idx) * pre_beta(inc_idx) - Xty(inc_idx));


            // pre_loss = lambda * (1 - alpha) * sum(square(pre_beta(inc_idx)))/2 + lambda * alpha * sum(abs(pre_beta(inc_idx)));
            pre_loss = objective(X2, y, new_beta, lambda, alpha);

            // [sort_lsearch, ix_lsearch] = sort([progress',1]);
            search_idx = sort_index(progress);
            for(int i = 0; i < search_idx.n_elem; i++){
                k = search_idx(i);
                
                if (k <= 0) {
                    continue;
                }

                if(k >= 1){
                    break;
                }

                tmp_beta = beta(inc_idx) + (new_beta - beta(inc_idx)) * k;
                // inner_loss = a * std.square(k) + b * k + lambda * (1 - alpha) * sum(square(tmp_beta))/2 + lambda * alpha * sum(abs(tmp_beta));
                inner_loss = objective(X2, y, tmp_beta, lambda, alpha);

                if (inner_loss <= pre_loss) {
                    pre_loss = inner_loss;
                    new_beta = tmp_beta;
                    line_search = k;
                } else {
                    break;
                } 
            }
            
            beta(inc_idx) = new_beta;
            theta(inc_idx) = sign(new_beta);

            // if beta encounters zero along the line search, then remove it from active set
            
            remove_idx = find(abs(new_beta) < datum::eps);
            if(!remove_idx.is_empty()){
                beta.elem(inc_idx(remove_idx)).zeros();
                theta.elem(inc_idx(remove_idx)).zeros();
                active_set.elem(inc_idx(remove_idx)).zeros();
            }
        }

        exc_idx = find(active_set == 0);
        grad = XtX(exc_idx, inc_idx) * beta(inc_idx) - Xty(exc_idx);
        violate_idx = find(abs(grad) > lambda * alpha);

        if(violate_idx.is_empty()){
            optimality0 = 1;
            if(optimality1 == 1){
                break;
            }
        } else {
            active_set.elem(exc_idx(violate_idx)).ones();
            theta(exc_idx(violate_idx)) = - sign(grad(exc_idx(violate_idx)));
        }
        iter++;
    }
    return beta;
}
