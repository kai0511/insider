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

void optimize_row(const mat& residual, const mat& indicator, mat& updating_factor, const mat& c_factor, 
                  const uvec& updating_confd, const double lambda, const int tuning, const int n_cores = 10){
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
            mat feaures, XtX = zeros(c_factor.n_rows, c_factor.n_rows);
            vec outcome, Xty = zeros<vec>(c_factor.n_rows);

            uvec ids = find(updating_confd == seq(i));

            for(unsigned int k = 0; k < ids.n_elem; k++){
                non_zeros = find(trans(indicator.row(ids(k))));
                feaures = c_factor.cols(non_zeros); 

                outcome = trans(residual.row(ids(k)));
                outcome = outcome(non_zeros);
                
                XtX += feaures * trans(feaures);
                Xty += feaures * outcome;
            }
            
            XtX.diag() += lambda;
            updating_factor.row(seq(i) - 1) = trans(solve(XtX, Xty, solve_opts::likely_sympd));
        }

    }else if(tuning == 0){

        mat gram = c_factor * trans(c_factor);
        mat Xtys = c_factor * trans(residual);

        #pragma omp parallel for num_threads(n_cores) schedule(dynamic, 1)
        for(unsigned int i = 0; i < seq.size(); i++) {

            uvec non_zeros; 
            mat XtX = zeros(c_factor.n_rows, c_factor.n_rows);
            vec Xty = zeros<vec>(c_factor.n_rows);

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

void optimize_col(const mat& residual, const mat& indicator, const mat& row_factor, mat& c_factor, 
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
                // mat left_matrix = XtX + lambda * eye(c_factor.n_rows, c_factor.n_rows);
                XtX.diag() += lambda; 
                c_factor.col(i) = solve(XtX, Xty, solve_opts::likely_sympd);
            }else{
                c_factor.col(i) = strong_coordinate_descent(feature, outcome, c_factor.col(i), lambda, alpha, XtX, Xty, tol);
            }
        }

    }else if(tuning == 0){
        
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

    } else{
        cout << "Parameter tuning should be either 0 or 1!" << endl;
        exit(1);
    }
}

// [[Rcpp::export]]
List optimize(const mat& data, List cfd_factors, mat& column_factor, const umat& cfd_indicators, const mat& train_indicator, 
              const int latent_dim, const double lambda = 1.0, const double alpha = 0.1, const int tuning = 1, const double global_tol=1e-10, const double sub_tol = 1e-6, const unsigned int max_iter = 10000){
    
    cout.precision(12);
    double delta_loss;
    unsigned int i, iter = 0, cfd_num = cfd_factors.size();
    uvec train_idx, test_idx;
    double loss, pre_loss, sum_residual, train_rmse, test_rmse, optimal_rmse, decay = 1.0; 
    mat residual, sub_matrix; 
    mat row_factor = zeros(data.n_rows, latent_dim) , predictions = zeros(size(data));
    List row_matrices;

    // check whether the number of the confounding matrices is equal to the number of confounding indicators.
    if(cfd_num != cfd_indicators.n_cols){
        cout << "The number of confounding matrices should be the same the number of indicators." << endl;
        exit(1);
    } 

    // put confounding matrices from Rcpp::List into arma::field
    field<mat> cfd_matrices(cfd_num);
    for(i = 0; i < cfd_num; i ++) {
        Rcpp::NumericMatrix temp = cfd_factors[i];
        cfd_matrices(i) = mat(temp.begin(), temp.nrow(), temp.ncol(), false);
        row_factor += cfd_matrices(i).rows(cfd_indicators.col(i) - 1);
    }

    // find the indices for training and testing elements 
    train_idx = find(train_indicator);
    test_idx = find(train_indicator == 0);

    // check the fitting with initial values
    predict(row_factor, column_factor, predictions);
    residual = data - predictions;
    evaluate(residual, train_idx, test_idx, sum_residual, train_rmse, test_rmse, tuning, iter, 1);
    loss = compute_loss(cfd_matrices, column_factor, lambda, alpha, sum_residual, 1);

    while(iter < max_iter) {

        if(iter % 10 == 0){
            cout << "Iteration " << iter << " ---------------------------------" << endl;
        }

        // update all confonding matrices
        for(i = 0; i < cfd_num; i++){
            sub_matrix = cfd_matrices(i) * column_factor;
            residual += sub_matrix.rows(cfd_indicators.col(i) - 1);

            optimize_row(residual, train_indicator, cfd_matrices(i), column_factor, cfd_indicators.col(i), lambda, tuning);

            sub_matrix = cfd_matrices(i) * column_factor;
            residual -= sub_matrix.rows(cfd_indicators.col(i) - 1);
        }

        // compute row_factor
        row_factor.zeros();
        for(i = 0; i < cfd_num; i ++) {
            row_factor += cfd_matrices(i).rows(cfd_indicators.col(i) - 1);
        }

        // update columm_factor
        optimize_col(data, train_indicator, row_factor, column_factor, lambda, alpha, tuning, 30, sub_tol * decay);
        predict(row_factor, column_factor, predictions);
        residual = data - predictions;

        // check the fitting every 10 steps
        if(iter % 10 == 0){
            pre_loss = loss;
            evaluate(residual, train_idx, test_idx, sum_residual, train_rmse, test_rmse, tuning, iter, 1);
            loss = compute_loss(cfd_matrices, column_factor, lambda, alpha, sum_residual, 1);

            cout << "Delta loss for iter " << iter << ":" << pre_loss - loss << endl;

            if((pre_loss - loss)/1000 <= 1e-6){
                decay = 1e-6;
            }else if((pre_loss - loss)/1000 <= 1e-5){
                decay = 1e-5;
            }else if((pre_loss - loss)/1000 <= 1e-4){
                decay = 1e-4;
            }else if((pre_loss - loss)/1000 <= 1e-3){
                decay = 1e-3;
            }else if((pre_loss - loss)/1000 <= 1e-2){
                decay = 1e-2;
            }else if((pre_loss - loss)/1000 <= 1e-1){
                decay = 1e-1;
            }else{
                decay = 1.0;
            }

            if((tuning == 1) & (optimal_rmse > test_rmse)){
                optimal_rmse = test_rmse;
            }

            if((pre_loss - loss)/pre_loss < global_tol){
                break;
            }
        }
        iter++;
    }
    
    // put the updated confounding matrices into a R list
    for(i = 0; i < cfd_num; i ++) {
        row_matrices["factor" + std::to_string(i)] = cfd_matrices(i);
    }

    return List::create(Named("row_matrices") = row_matrices,
                        Named("column_factor") = column_factor, 
                        Named("optimal_rmse") = optimal_rmse, 
                        Named("loss") = loss);
}
