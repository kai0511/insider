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


mat fit_intraction(const mat& residual, const mat& train_indicator, const umat& cfd_indicators, const mat& columm_factor, const umat& unique_cfd,
                   const double lambda, const double alpha, const int tuning, const double& tol = 1e-10, const int n_cores = 10){
    /*
        fix column parameters and update row factors 
    args:
        @dimension must be in (1, 2), denoting which row_factor will be updated
        @n_cores number of cores for parallel computation
    */
    mat interactions = zeros(unique_cfd.n_rows, columm_factor.n_rows);

    if(tuning == 1){
        
        #pragma omp parallel for num_threads(n_cores) schedule(dynamic, 1)
        for(unsigned int i = 0; i < unique_cfd.size(); i++) {

            uvec non_zeros; 
            uvec cfd = unique_cfd.row(i);
            uvec ids = find((inc_cfd_indicators.col(1) == cfd(1)) && (inc_cfd_indicators.col(2) == cfd(2)));
            //uvec ids = inc_cfd_indicators.each_row( [](const rowvec& a){ return approx_equal(a, cfd, "absdiff", 1e-8)? 1 : 0; }); 

            vec wstart = zeros<vec>(columm_factor.n_rows);
            int st_idx = 0, ed_idx, nonzero_num = accu(train_indicator.rows(ids));
            vec sub_outcome, outcomes = zeros(nonzero_num), Xty = zeros<vec>(columm_factor.n_rows);
            mat sub_feature, features = zeros(nonzero_num, columm_factor.n_rows), XtX = zeros(columm_factor.n_rows, columm_factor.n_rows);

            for(unsigned int k = 0; k < ids.n_elem; k++){
                non_zeros = find(trans(train_indicator.row(ids(k))));
                ed_idx = st_idx + non_zeros.n_elem - 1;

                sub_feature = columm_factor.cols(non_zeros); 
                features.rows(st_idx, end_idx) = trans(sub_feature);

                sub_outcome = trans(residual.row(ids(k)));
                sub_outcome = sub_outcome(non_zeros);
                outcomes.subvec(st_idx, end_idx) = sub_outcome;

                XtX += sub_feature * features.rows(st_idx, end_idx);
                Xty += sub_feature * sub_outcome;

                st_idx = end_idx + 1;
            }
            interactions.row(i) = strong_coordinate_descent(features, outcomes, wstart, lambda, alpha, XtX, Xty, tol);
        }

    }else if(tuning == 0){

        mat gram = columm_factor * trans(columm_factor);
        mat Xtys = columm_factor * trans(residual);

        #pragma omp parallel for num_threads(n_cores) schedule(dynamic, 1)
        for(unsigned int i = 0; i < unique_cfd.size(); i++) {

            uvec non_zeros; 
            uvec cfd = unique_cfd.row(i);
            uvec ids = find((inc_cfd_indicators.col(1) == cfd(1)) && (inc_cfd_indicators.col(2) == cfd(2)));
            //uvec ids = inc_cfd_indicators.each_row( [](const rowvec& a){ return approx_equal(a, cfd, "absdiff", 1e-8)? 1 : 0; }); 

            vec wstart = zeros<vec>(columm_factor.n_rows);
            int st_idx = 0, ed_idx, nonzero_num = ids.n_elem * columm_factor.n_cols;
            vec sub_outcome, outcomes = zeros(nonzero_num), Xty = zeros<vec>(columm_factor.n_rows);
            mat sub_feature, features = zeros(nonzero_num, columm_factor.n_rows), XtX = zeros(columm_factor.n_rows, columm_factor.n_rows);

            for(unsigned int k = 0; k < ids.n_elem; k++){
                ed_idx = st_idx + columm_factor.n_cols - 1;
                features.rows(st_idx, ed_idx) = trans(columm_factor);
                outcomes.subvec(st_idx, ed_idx) = trans(residual.row(ids(k)));
                st_idx = ed_idx + 1;

                XtX += gram;
                Xty += Xtys.col(ids(k));
            }
            interactions.row(i) = strong_coordinate_descent(features, outcomes, wstart, lambda, alpha, XtX, Xty, tol);
        }

    }else{
        cout << "Parameter tuning should be either 0 or 1!" << endl;
        exit(1);
    }
    return interactions;
    
}
