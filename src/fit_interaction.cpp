// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins("cpp11")]]

#include "../inst/include/insider_types.h"
#include <iostream>
#include <omp.h>
#include "coordinate_descent.h"
#include "utils.h"

using Rcpp::List;
using Rcpp::Named;


void fit_intraction(const mat& residual, const mat& train_indicator, const umat& inc_cfd_indicators,, const mat& columm_factor, 
                    const double lambda, const double alpha, const int tuning, const int n_cores = 10){
    /*
        fix column parameters and update row factors 
    args:
        @dimension must be in (1, 2), denoting which row_factor will be updated
        @n_cores number of cores for parallel computation
    */
    uvec seq = unique(updating_confd);
    // mat updates = zeros(size(updating_factor));

    if(tuning == 1){
        
        #pragma omp parallel for num_threads(n_cores) schedule(dynamic, 1)
        for(unsigned int i = 0; i < seq.size(); i++) {

            uvec non_zeros; 
            mat feaures, XtX = zeros(columm_factor.n_rows, columm_factor.n_rows);
            vec outcome, Xty = zeros<vec>(columm_factor.n_rows);

            uvec ids = find(updating_confd == seq(i));

            for(unsigned int k = 0; k < ids.n_elem; k++){
                non_zeros = find(trans(indicator.row(ids(k))));
                feaures = columm_factor.cols(non_zeros); 

                outcome = trans(residual.row(ids(k)));
                outcome = outcome(non_zeros);
                
                XtX += feaures * trans(feaures);
                Xty += feaures * outcome;
            }
            
            XtX.diag() += lambda;
            updating_factor.row(seq(i) - 1) = trans(solve(XtX, Xty, solve_opts::likely_sympd));
        }

    }else if(tuning == 0){

        mat gram = columm_factor * trans(columm_factor);
        mat Xtys = columm_factor * trans(residual);

        #pragma omp parallel for num_threads(n_cores) schedule(dynamic, 1)
        for(unsigned int i = 0; i < seq.size(); i++) {

            uvec non_zeros; 
            mat XtX = zeros(columm_factor.n_rows, columm_factor.n_rows);
            vec Xty = zeros<vec>(columm_factor.n_rows);

            uvec ids = find(updating_confd == seq(i));
            for(unsigned int k = 0; k < ids.n_elem; k++){
                XtX += gram;
                Xty += Xtys.col(ids(k));
            }
            
            XtX.diag() += lambda;
            updating_factor.row(seq(i) - 1) = trans(solve(XtX, Xty, solve_opts::likely_sympd));
        }

    }else{
        cout << "Parameter tuning should be either 0 or 1!" << endl;
        exit(1);
    }
    
}