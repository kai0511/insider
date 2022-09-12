// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins("cpp11")]]

#include "../inst/include/insider_types.h"
#include <iostream>
#include <RcppArmadillo.h>
#include "coordinate_descent.h"
#include "utils.h"

// [[Rcpp::export]]

void row_optimize(const mat& residual, const mat& indicator, const mat& c_factor, 
                  mat& updating_factor, const uvec& updating_confd, const mat& fixed_factor, const uvec& fixed_confd, 
                  const double& lambda, const int& tuning, const int& n_cores = 10) {
    /*
        fix column parameters and update row factors 
    args:
        @dimension must be in (1, 2), denoting which row_factor will be updated
        @n_cores number of cores for parallel computation
    */
    uvec levels = unique(updating_confd);
    int ncores = std::min(levels.n_elem, n_cores);

    if(tuning == 1){
        
        #pragma omp parallel for num_threads(ncores) schedule(dynamic, 1)
        for(unsigned int i = 0; i < levels.size(); i++) {

            uvec non_zeros;
            mat feaures, XtX = zeros(c_factor.n_rows, c_factor.n_rows);
            vec outcome, Xty = zeros<vec>(c_factor.n_rows);

            uvec ids = find(updating_confd == seq(i));
            uvec selected = fixed_confd(ids) - 1;

            for(unsigned int k = 0; k < ids.n_elem; k++){
                
                non_zeros = find(trans(indicator.row(ids(k))));
                feaures = c_factor.cols(non_zeros);
                feaures.each_col() %= trans(fixed_factor.row(selected(k)));
                
                outcome = trans(residual.row(ids(k)));
                outcome = outcome(non_zeros);

                XtX += feaures * trans(feaures);
                Xty += feaures * outcome;
            }
            
            XtX.diag() += lambda;
            updating_factor.row(levels(i) - 1) = trans(solve(XtX, Xty, solve_opts::likely_sympd));
        }

    } else if(tuning == 0){
        
        #pragma omp parallel for num_threads(n_cores) schedule(dynamic, 1)
        for(unsigned int i = 0; i < levels.size(); i++) {
            uvec non_zeros; 
            mat XtX = zeros(c_factor.n_rows, c_factor.n_rows);
            vec Xty = zeros<vec>(c_factor.n_rows);

            uvec ids = find(updating_confd == levels(i));
            uvec selected = fixed_confd(ids) - 1;
            
            for(unsigned int k = 0; k < ids.n_elem; k++){
                
                mat X = c_factor.each_col() % trans(fixed_factor.row(selected(k)));
                XtX += trans(X) * X;
                Xty += c_factor * trans(residual.row(ids(k)));
            }
            XtX.diag() += lambda;
            updating_factor.row(levels(i) - 1) = trans(solve(XtX, Xty, solve_opts::likely_sympd));
        }

    }else{
        cout << "Parameter tuning should be either 0 or 1!" << endl;
        exit(1);
    }
}

void column_optimize(const mat& residual, const mat& indicator, const mat& row_factor, mat& c_factor, 
                     const double& lambda, const double& alpha, const int tuning, const int n_cores = 10, const double tol = 1e-5){
    
    if(tuning == 1){

        cube feature_space(c_factor.n_rows, c_factor.n_rows, residual.n_rows);
        for(unsigned int i = 0; i < row_factor.n_rows; i++){
            feature_space.slice(i) = row_factor.row(i).t() * row_factor.row(i);
        }

        #if defined(_OPENMP)
            #pragma omp parallel for num_threads(n_cores) schedule(dynamic, 100)
        #endif
        for(unsigned int i = 0; i < residual.n_cols; i++) {
            uvec row_selected = find(indicator.col(i));
            mat feature = row_factor.rows(row_selected);
            mat XtX = sum(feature_space.slices(row_selected), 2);
            vec outcome = residual.col(i);
            outcome = outcome(row_selected);
            vec Xty = trans(feature) * outcome;

            if(alpha == 0){
                XtX.diag() += lambda; 
                c_factor.col(i) = solve(XtX, Xty, solve_opts::likely_sympd);
            }else{
                c_factor.col(i) = strong_coordinate_descent(feature, outcome, c_factor.col(i), lambda, alpha, XtX, Xty, tol);
            }
        }

    }else if(tuning == 0) {
        
        mat XtX = trans(row_factor) * row_factor;
        mat Xty = trans(row_factor) * residual;

        if(alpha == 0){
            XtX.diag() += lambda; 
        }

        #if defined(_OPENMP)
            #pragma omp parallel for num_threads(n_cores) schedule(dynamic, 100)
        #endif
        for(unsigned int i = 0; i < residual.n_cols; i++) {
            
            if(alpha == 0){
                c_factor.col(i) = solve(XtX, Xty, solve_opts::likely_sympd);
            }else{
                c_factor.col(i) = strong_coordinate_descent(row_factor, residual.col(i), c_factor.col(i), lambda, alpha, XtX, Xty.col(i), tol);
            }
        }

    } else {
        cout << "Parameter tuning should be either 0 or 1!" << endl;
        exit(1);
    }
}