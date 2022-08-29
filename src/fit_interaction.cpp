// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins("cpp11")]]

#include "../inst/include/insider_types.h"
#include "fit_interaction.h"
#include <iostream>
#include <omp.h>

void fit_interaction(const mat& residual, const mat& train_indicator, mat& interactions, const uvec& interaction_indicator, const mat& column_factor, 
                     const double lambda, const int tuning, const double& tol = 1e-10, const int n_cores = 10){
    /*
        fix column parameters and update row factors 
    args:
        @dimension must be in (1, 2), denoting which row_factor will be updated
        @n_cores number of cores for parallel computation
    */
    uvec unique_ita = unique(interaction_indicator);
    // mat interactions = zeros(unique_ita.n_elem, column_factor.n_rows);

    if(tuning == 1){
        
        #pragma omp parallel for num_threads(n_cores) schedule(dynamic, 1)
        for(unsigned int i = 0; i < unique_ita.n_elem; i++) {

            uvec non_zeros; 
            uvec ids = find(interaction_indicator == unique_ita(i));

            // vec wstart = zeros<vec>(column_factor.n_rows);
            // int st_idx = 0, ed_idx, nonzero_num = accu(train_indicator.rows(ids));
            vec sub_outcome, Xty = zeros<vec>(column_factor.n_rows);
            mat sub_feature, XtX = zeros(column_factor.n_rows, column_factor.n_rows);

            // vec outcomes = zeros(nonzero_num);
            // mat features = zeros(nonzero_num, column_factor.n_rows);

            for(unsigned int k = 0; k < ids.n_elem; k++){
                non_zeros = find(trans(train_indicator.row(ids(k))));
                // ed_idx = st_idx + non_zeros.n_elem - 1;

                sub_feature = column_factor.cols(non_zeros); 
                // features.rows(st_idx, ed_idx) = trans(sub_feature);

                sub_outcome = trans(residual.row(ids(k)));
                sub_outcome = sub_outcome(non_zeros);
                // outcomes.subvec(st_idx, ed_idx) = sub_outcome;

                XtX += sub_feature * trans(sub_feature);
                Xty += sub_feature * sub_outcome;

                // st_idx = ed_idx + 1;
            }
            // interactions.row(i) = trans(strong_coordinate_descent(features, outcomes, wstart, lambda, alpha, XtX, Xty, tol));
            interactions.row(i) = solve(XtX, Xty, solve_opts::likely_sympd);
        }

    }else if(tuning == 0){

        mat gram = column_factor * trans(column_factor);
        mat Xtys = column_factor * trans(residual);

        #pragma omp parallel for num_threads(n_cores) schedule(dynamic, 1)
        for(unsigned int i = 0; i < unique_ita.n_elem; i++) {

            uvec ids = find(interaction_indicator == unique_ita(i));
            // vec wstart = zeros<vec>(column_factor.n_rows);
            // int st_idx = 0, ed_idx, nonzero_num = ids.n_elem * column_factor.n_cols;
            // mat features = zeros(nonzero_num, column_factor.n_rows); 
            // vec outcomes = zeros(nonzero_num);
            mat XtX = zeros(column_factor.n_rows, column_factor.n_rows);
            vec Xty = zeros<vec>(column_factor.n_rows);

            for(unsigned int k = 0; k < ids.n_elem; k++){
                // ed_idx = st_idx + column_factor.n_cols - 1;
                // features.rows(st_idx, ed_idx) = trans(column_factor);
                // outcomes.subvec(st_idx, ed_idx) = trans(residual.row(ids(k)));
                // st_idx = ed_idx + 1;
                XtX += gram;
                Xty += Xtys.col(ids(k));
            }
            // interactions.row(i) = trans(strong_coordinate_descent(features, outcomes, wstart, lambda, alpha, XtX, Xty, tol));
            interactions.row(i) = solve(XtX, Xty, solve_opts::likely_sympd);
        }

    }else{
        cout << "Parameter tuning should be either 0 or 1!" << endl;
        exit(1);
    }
    // return interactions;
}
