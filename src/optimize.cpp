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

// [[Rcpp::export]]
void optimize_continuous(const mat& data, const mat& indicator, rowvec& updating_factor, const mat& c_factor,
                         const vec& updating_confd, const mat& gram, const double lambda, const int tuning){

    if(tuning == 1){

        uvec non_zeros;
        vec factor, outcome;

        mat resid = data - updating_confd * updating_factor * c_factor;
        // cout << size(resid) << endl;
        double XtX, Xty, pre_loss = 0.0, delta_loss = 0.0;
        double loss = 0.5 * (sum(square(resid.elem(find(indicator)))) + lambda * pow(norm(updating_factor, "F"), 2));

        while(1){

            for(int i = 0; i < c_factor.n_rows; i++) {

                XtX = 0.0, Xty = 0.0;
                resid += updating_factor(i) * updating_confd * c_factor.row(i);
                factor = trans(c_factor.row(i));
                // cout << 1 << endl;
                for(int k = 0; k < indicator.n_rows; k++){
                    non_zeros = find(indicator.row(k));

                    outcome = trans(resid.row(k));
                    outcome = outcome(non_zeros);

                    XtX += std::pow(updating_confd(k), 2) * dot(factor(non_zeros), factor(non_zeros));
                    Xty += updating_confd(k) * dot(factor(non_zeros), outcome);
                }
                updating_factor(i) = Xty/(XtX + lambda);
                resid -= updating_factor(i) * updating_confd * c_factor.row(i);
            }

            pre_loss = loss;
            loss = 0.5 * (sum(square(resid.elem(find(indicator)))) + lambda * pow(norm(updating_factor, "F"), 2));
            delta_loss = pre_loss - loss;

            if(delta_loss < 1e-5){
                // cout << "Delta loss: " << delta_loss << endl;
                break;
            }
        }
    } else if(tuning == 0){
        vec Xty = c_factor * trans(data) * updating_confd ;
        mat XtX = dot(updating_confd, updating_confd) * gram;
        XtX.diag() += lambda;
        updating_factor = trans(solve(XtX, Xty, solve_opts::likely_sympd));

    } else {
        cout << "Parameter tuning should be either 0 or 1!" << endl;
        exit(1);
    }
}

void optimize_row(const mat& residual, const mat& indicator, mat& updating_factor, const mat& c_factor, 
                  const uvec& updating_confd, const mat& gram, const double lambda, const int tuning, const int n_cores = 10){
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

            uvec non_zeros, zero_idx; 
            mat feaures, XtX = zeros(c_factor.n_rows, c_factor.n_rows);
            vec outcome, Xty = zeros<vec>(c_factor.n_rows);

            uvec ids = find(updating_confd == seq(i));

            for(unsigned int k = 0; k < ids.n_elem; k++){
                non_zeros = find(trans(indicator.row(ids(k))));
                zero_idx = find(trans(indicator.row(ids(k))) == 0);  
                // feaures = c_factor.cols(non_zeros); 

                outcome = trans(residual.row(ids(k)));
                outcome = outcome(non_zeros);
                
                // XtX += feaures * trans(feaures);
                XtX += gram - c_factor.cols(zero_idx) * trans(c_factor.cols(zero_idx));
                Xty += c_factor.cols(non_zeros) * outcome;
            }
            
            XtX.diag() += lambda;
            updating_factor.row(seq(i) - 1) = trans(solve(XtX, Xty, solve_opts::likely_sympd));
        }

    }else if(tuning == 0){

        mat Xtys = c_factor * trans(residual);

        #pragma omp parallel for num_threads(n_cores) schedule(dynamic, 1)
        for(unsigned int i = 0; i < seq.size(); i++) {

            uvec ids = find(updating_confd == seq(i));
            mat XtX = ids.n_elem * gram;
            XtX.diag() += lambda;
            vec Xty = sum(Xtys.cols(ids), 1);
            
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
        
        mat gram = trans(row_factor) * row_factor;

        cube feature_space(c_factor.n_rows, c_factor.n_rows, residual.n_rows);
        for(unsigned int i = 0; i < row_factor.n_rows; i++){
            feature_space.slice(i) = row_factor.row(i).t() * row_factor.row(i);
        }

        #if defined(_OPENMP)
            #pragma omp parallel for num_threads(n_cores) schedule(dynamic, 100)
        #endif
        for(unsigned int i = 0; i < residual.n_cols; i++) {
            uvec selected = find(indicator.col(i));
            mat feature = row_factor.rows(selected);
            mat XtX = sum(feature_space.slices(find(indicator.col(i) == 0.0)), 2);
            XtX = gram - XtX;
            vec outcome = residual.col(i);
            outcome = outcome(selected);
            vec Xty = trans(feature) * outcome;

            if(alpha == 0.0){
                XtX.diag() += lambda; 
                c_factor.col(i) = solve(XtX, Xty, solve_opts::likely_sympd);
            }else{
                c_factor.col(i) = strong_coordinate_descent(feature, outcome, c_factor.col(i), lambda, alpha, XtX, Xty, tol);
            }
        }

    }else if(tuning == 0){
        
        mat XtX = trans(row_factor) * row_factor;
        mat Xty = trans(row_factor) * residual;

        if (alpha == 0.0) {
            XtX.diag() += lambda; 
            // c_factor.col(i) = solve(XtX, Xty.col(i), solve_opts::likely_sympd);
            c_factor = solve(XtX, Xty, solve_opts::likely_sympd);
        } else {
            #if defined(_OPENMP)
                #pragma omp parallel for num_threads(n_cores) schedule(dynamic, 100)
            #endif
            for(unsigned int i = 0; i < residual.n_cols; i++) {
                c_factor.col(i) = strong_coordinate_descent(row_factor, residual.col(i), c_factor.col(i), lambda, alpha, XtX, Xty.col(i), tol);
            }
        }
    } else{
        cout << "Parameter tuning should be either 0 or 1!" << endl;
        exit(1);
    }
}

// [[Rcpp::export]]
List optimize(const mat& data, List cfd_factors, mat& column_factor, const umat& cfd_indicators, const mat&ctns_confounder, const mat& train_indicator, const mat& test_indicator, const int& inc_continuous,
              const int latent_dim, const double lambda1 = 1.0, const double lambda2 = 1.0, const double alpha = 0.1, const int tuning = 1, const double global_tol=1e-10, const double sub_tol = 1e-5, const unsigned int max_iter = 10000) {
    
    cout.precision(12);
    double delta_loss;
    unsigned int i, j, iter = 0, cfd_num = cfd_indicators.n_cols;
    uvec train_idx, test_idx, ord;
    double loss, pre_loss, sum_residual, train_rmse, test_rmse, decay = 1.0; 
    mat gram, residual, sub_matrix; 
    mat row_factor = zeros(data.n_rows, latent_dim) , predictions = zeros(size(data));
    List row_matrices;

    // to support continuous covariates
    if(inc_continuous != 0 || inc_continuous != 1){
        cout << "The value of prarameter inc_continuous can only be 0 or 1." << endl;
        exit(1);
    }

    // to support continuous covariates
    if(inc_continuous == 1){
        cfd_num += 1;
    }

    // move confounding matrices from Rcpp::List into arma::field
    field<mat> cfd_matrices(cfd_num);
    for(i = 0; i < cfd_num; i ++) {
        Rcpp::NumericMatrix temp = cfd_factors[i];
        cfd_matrices(i) = mat(temp.begin(), temp.nrow(), temp.ncol(), false);
        if(i < cfd_indicators.n_cols){
            row_factor += cfd_matrices(i).rows(cfd_indicators.col(i) - 1);
        }else{
            // to support continuous covariates
            row_factor += ctns_confounder * cfd_matrices(i);
        }
    }

    // place indices of confounders into arma::field for computational consideration
    field<mat> index_matrices(cfd_num);
    // field<vec> confd_counts(cfd_num);
    for(i = 0; i < cfd_num; i ++) {
        if(i < cfd_indicators.n_cols){
            uvec levels = unique(cfd_indicators.col(i));
            mat cfd_idx = zeros(cfd_indicators.n_rows, levels.n_elem);

            for(unsigned int k = 0; k < levels.n_elem; k++) {
                vec tmp = cfd_idx.col(k);
                tmp.elem(find(cfd_indicators.col(i) == levels(k))).ones();
                cfd_idx.col(k) = tmp;
            }
            index_matrices(i) = cfd_idx;
            // confd_counts(i) = trans(sum(cfd_idx));
        } else {
            // to support continuous covariates
            index_matrices(i) = ctns_confounder;
            // confd_counts(i) = ones(cfd_indicators.n_rows);
        }
    }

    // find the indices for training and testing elements 
    train_idx = find(train_indicator);
    test_idx = find(test_indicator);

    // check the fitting with initial values
    predict(row_factor, column_factor, predictions);
    residual = data - predictions;
    evaluate(residual, train_idx, test_idx, sum_residual, train_rmse, test_rmse, tuning, iter, 1);
    loss = compute_loss(cfd_matrices, column_factor, lambda1, lambda2, alpha, sum_residual, 1);

    while(iter <= max_iter) {

        if(iter % 10 == 0){
            cout << "Iteration " << iter << " ---------------------------------" << endl;
        }

        // update all confonding matrices
        gram = column_factor * column_factor.t();
        // ord = randperm(cfd_num);

        for(i = 0; i < cfd_num; i++){
            // i = ord(j);
            if(i < cfd_indicators.n_cols){
                residual += index_matrices(i) * cfd_matrices(i) * column_factor;
                optimize_row(residual, train_indicator, cfd_matrices(i), column_factor, cfd_indicators.col(i), gram, lambda1, tuning);
            } else {
                // to support continuous covariates
                for (j = 0; j < ctns_confounder.n_cols; j++){
                    residual += ctns_confounder.col(j) * cfd_matrices(i).row(j) * column_factor;
                    optimize_continuous(residual, train_indicator, cfd_matrices(i).row(j), column_factor, ctns_confounder.col(j), gram, lambda1, tuning);
                    if(j != ctns_confounder.n_cols - 1){
                        residual -= ctns_confounder.col(j) * cfd_matrices(i).row(j) * column_factor;
                    }
                }
            }
            
            if(i != cfd_num - 1){
                residual -= index_matrices(i) * cfd_matrices(i) * column_factor;
            }
	        // pre_loss = loss;
            // evaluate(residual, train_idx, test_idx, sum_residual, train_rmse, test_rmse, tuning, iter, 1);
            // loss = compute_loss(cfd_matrices, column_factor, lambda, alpha, sum_residual, 1);

            // delta_loss = pre_loss - loss;
            // cout << "Delta loss for updating cfd factor in iter " << iter << ":" << delta_loss << endl;
        }
        
        // compute row_factor
        row_factor.zeros();
        for(i = 0; i < cfd_indicators.n_cols; i ++) {
            // row_factor += cfd_matrices(i).rows(cfd_indicators.col(i) - 1);
            row_factor += index_matrices(i) * cfd_matrices(i);
        }

        if(inc_continuous == 1){
            row_factor += ctns_confounder * cfd_matrices(cfd_num-1);
        }

        // update columm_factor
        optimize_col(data, train_indicator, row_factor, column_factor, lambda2, alpha, tuning, 30, sub_tol * decay);
        predict(row_factor, column_factor, predictions);
        residual = data - predictions;
        
        // check the fitting every 10 steps
        if(iter % 10 == 0){
            pre_loss = loss;
            evaluate(residual, train_idx, test_idx, sum_residual, train_rmse, test_rmse, tuning, iter, 1);
            loss = compute_loss(cfd_matrices, column_factor, lambda1, lambda2, alpha, sum_residual, 1);
            
            delta_loss = pre_loss - loss;
            cout << "Delta loss for iter " << iter << ":" << delta_loss << endl;

            if(delta_loss/1000 <= 1e-6){
                decay = 1e-6;
            }else if(delta_loss/1000 <= 1e-5){
                decay = 1e-5;
            }else if(delta_loss/1000 <= 1e-4){
                decay = 1e-4;
            }else if(delta_loss/1000 <= 1e-3){
                decay = 1e-3;
            }else if(delta_loss/1000 <= 1e-2){
                decay = 1e-2;
            }else if(delta_loss/1000 <= 1e-1){
                decay = 1e-1;
            }else{
                decay = 1.0;
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
                        Named("train_rmse") = train_rmse, 
                        Named("test_rmse") = test_rmse, 
                        Named("loss") = loss);
}
